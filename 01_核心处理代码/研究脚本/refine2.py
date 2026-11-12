#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Bellegarde 精配 v2：把尺度当未知量的联合搜索。

v1（refine_pair.py）失败的教训：拿全分辨率 ACI 特征模板去匹配重采样过的
WATSON，频段对不上，匹配全盘失败。这说明两机采样率差（几何说 ~2.75 倍，
纹理测量说 1.01）必须在搜索里当成未知量，不能假设几何模型给的对。

v2 的做法：
  1. 从 CAHVOR 的 |H| 直接算两台相机的角分辨率（µrad/px），
     连同预测尺度比曲线一起打印——这是"模型说什么"；
  2. 锚点不再靠目测报坐标，程序化检测 ACI 里的亮斑簇；
  3. 固定几何解 standoff=129 mm，联合搜索 4 种翻转 × 40 档缩放
     × 3 档模糊，每档对全部锚点做模板匹配；
  4. 选择判据：≥3 个锚点得分过关的前提下，得分中位数高、且各锚点
     残余平移一致（离散小）——缩放对了平移才会一致，这是自带的方向感；
  5. 赢家做零假设检验（±60 px 平移后 NCC 应显著掉），
     并反解"等效 standoff"：如果落在物理合理区间外，
     结论就是相机模型的尺度可疑，而不是距离猜错。

用法：python refine2.py
"""

from __future__ import annotations

import os
import sys

import numpy as np
import cv2

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
sys.path.insert(0, HERE)

from colorize_pair import (WAT_XML, ACI_XML, WAT_IMG, ACI_IMG,   # noqa: E402
                           load_cam, read_bsq, to8, feature,
                           homography, sample_watson, ncc)

D0 = 129                 # 几何解 standoff（视轴交点），搜索围绕它做
SEARCH_PX = 300          # 模板匹配邻域
TPL = 280                # 锚点模板边长
MUS = np.geomspace(0.30, 3.3, 40)   # 相对几何解的缩放搜索
# WATSON 角分辨率粗 5.4 倍，投影到 ACI 网格要上采样 ~2.8x，
# ACI 一侧的模板必须低通到 WATSON 的有效带宽，σ 要给到 5-7 才够
SIGMAS = (0.0, 1.5, 3.0, 5.0, 7.0)
FLIPS = {"none": None, "h": "h", "v": "v", "hv": "hv"}


def flip_matrix(kind, ns, nl):
    """坐标空间的翻转矩阵（作用在采样坐标上，不动图像数组）。"""
    if kind == "h":
        return np.float64([[-1, 0, ns - 1], [0, 1, 0], [0, 0, 1]])
    if kind == "v":
        return np.float64([[1, 0, 0], [0, -1, nl - 1], [0, 0, 1]])
    if kind == "hv":
        return np.float64([[-1, 0, ns - 1], [0, -1, nl - 1], [0, 0, 1]])
    return np.eye(3)


def zoom_matrix(mu, ns, nl):
    """绕画面中心的缩放（坐标形式：内容放大 mu 倍）。"""
    c = np.array([ns / 2.0, nl / 2.0])
    return np.float64([[mu, 0, (1 - mu) * c[0]],
                       [0, mu, (1 - mu) * c[1]],
                       [0, 0, 1]])


def ang_per_px(cam):
    """CAHV 的 |H|、|V| 模长 -> 每像素张角（µrad）。"""
    return (1.0 / np.linalg.norm(cam["H"]) * 1e6,
            1.0 / np.linalg.norm(cam["V"]) * 1e6)


def detect_anchors(aci8, want=6):
    """程序化找亮斑簇：局部极大点 -> 网格聚类 -> 取对比度×数量最高的簇。

    返回 [(cx, cy, 分数), ...]，按分排序，且互相至少隔 260 px。
    """
    g = cv2.GaussianBlur(aci8.astype(np.float32), (0, 0), 1.2)
    dil = cv2.dilate(g, np.ones((25, 25), np.float32))
    peaks = (g == dil) & (g > g.mean() + 2.2 * g.std())
    ys, xs = np.nonzero(peaks)
    if len(xs) == 0:
        return []
    # 100px 网格聚类
    bins = {}
    for x, y in zip(xs, ys):
        bins.setdefault((x // 100, y // 100), []).append((x, y))
    clusters = []
    for (_bx, _by), pts in bins.items():
        xs_c = np.array([p[0] for p in pts])
        ys_c = np.array([p[1] for p in pts])
        bright = float(np.mean(g[ys_c, xs_c]))
        clusters.append((float(xs_c.mean()), float(ys_c.mean()),
                         bright * len(pts)))
    clusters.sort(key=lambda c: -c[2])
    picked = []
    for cx, cy, sc in clusters:
        if all(np.hypot(cx - px, cy - py) > 260 for px, py, _s in picked):
            picked.append((cx, cy, sc))
        if len(picked) >= want:
            break
    return picked


def edge_ratio(Hx, ns, nl):
    corners = np.float32([[0, 0], [ns - 1, 0], [ns - 1, nl - 1], [0, nl - 1]])
    proj = cv2.perspectiveTransform(corners.reshape(-1, 1, 2),
                                    np.float64(Hx)).reshape(-1, 2)
    d_w = sum(np.hypot(*(proj[(i + 1) % 4] - proj[i])) for i in range(4))
    d_a = sum(np.hypot(*(corners[(i + 1) % 4] - corners[i])) for i in range(4))
    return d_w / d_a


def main():
    wat, aci = load_cam(WAT_XML), load_cam(ACI_XML)
    wr, ar = read_bsq(WAT_IMG, wat), read_bsq(ACI_IMG, aci)
    idx = [i - 1 for i in wat["rgb"]]
    wcol8 = np.stack([to8(wr[i]) for i in idx], -1)
    aci8 = to8(ar[0])
    wgray = cv2.cvtColor(wcol8, cv2.COLOR_RGB2GRAY)
    ns, nl = aci["ns"], aci["nl"]

    # ---- 0) 模型说了什么：角分辨率 + 预测尺度比 ----
    (ah, av), (wh, wv) = ang_per_px(aci), ang_per_px(wat)
    print("角分辨率（µrad/px）：ACI H=%.1f V=%.1f  WATSON H=%.1f V=%.1f"
          % (ah, av, wh, wv))
    sep = float(np.linalg.norm(wat["C"] - aci["C"])) * 1000
    print("光心间距 %.1f mm（WATSON 在 ACI 后方）" % sep)
    H0, _ = homography(aci, wat, D0)
    r0 = edge_ratio(H0, ns, nl)
    print("几何解 %d mm 的尺度比（WATSON px / ACI px，四边平均）= %.3f" % (D0, r0))
    print("  · 若 ACI=10.1 µm/px → WATSON ≈ %.1f µm/px" % (r0 * 10.1))
    print("  · 若 ACI= 6.66 µm/px → WATSON ≈ %.1f µm/px" % (r0 * 6.66))
    for d in (60, 90, 129, 190, 260):
        Hx, _ = homography(aci, wat, d)
        print("    预测比 @%3d mm = %.3f" % (d, edge_ratio(Hx, ns, nl)))

    # ---- 1) 锚点 ----
    anchors = detect_anchors(aci8)
    print("\n程序化锚点（亮斑簇）：%d 个" % len(anchors))
    for cx, cy, sc in anchors:
        print("   (%4.0f, %4.0f)  分 %.0f  [分数坐标 (%.2f, %.2f)]"
              % (cx, cy, sc, cx / ns, cy / nl))
    if len(anchors) < 3:
        print("锚点太少，没法继续")
        return 1

    af = feature(aci8)
    wf = feature(wgray)

    def make_tpl(cx, cy):
        x0 = int(np.clip(cx - TPL / 2, 0, ns - TPL))
        y0 = int(np.clip(cy - TPL / 2, 0, nl - TPL))
        return x0, y0

    boxes = [make_tpl(*a[:2]) for a in anchors]

    # ---- 2) 联合搜索 ----
    Zs = {mu: zoom_matrix(mu, ns, nl) for mu in MUS}
    # 分辨率匹配是主角：模板必须从模糊后的 ACI 特征图里取
    afs = {s: (af if s == 0 else cv2.GaussianBlur(af, (0, 0), s))
           for s in SIGMAS}
    tpls = {s: [afs[s][y0:y0 + TPL, x0:x0 + TPL] for x0, y0 in boxes]
            for s in SIGMAS}
    combos = []
    for fname, Fk in FLIPS.items():
        Fm = flip_matrix(Fk, wat["ns"], wat["nl"]) if Fk else np.eye(3)
        for mu in MUS:
            Mz = H0 @ Zs[mu] if Fk is None else Fm @ H0 @ Zs[mu]
            wfl = sample_watson(wf, Mz, (ns, nl))
            for s in SIGMAS:
                row = []
                for (x0, y0), tpl in zip(boxes, tpls[s]):
                    sx0, sy0 = max(0, x0 - SEARCH_PX), max(0, y0 - SEARCH_PX)
                    sx1 = min(ns, x0 + TPL + SEARCH_PX)
                    sy1 = min(nl, y0 + TPL + SEARCH_PX)
                    if sx1 - sx0 < TPL or sy1 - sy0 < TPL:
                        row.append((0.0, 0, 0))
                        continue
                    res = cv2.matchTemplate(wfl[sy0:sy1, sx0:sx1], tpl,
                                            cv2.TM_CCOEFF_NORMED)
                    _, mx, _, loc = cv2.minMaxLoc(res)
                    row.append((float(mx), (loc[0] + sx0) - x0,
                                (loc[1] + sy0) - y0))
                scores = np.array([r[0] for r in row])
                offs = np.array([[r[1], r[2]] for r in row])
                nok = int((scores > 0.18).sum())
                spread = (float(np.hypot(offs[:, 0] - offs[:, 0].mean(),
                                         offs[:, 1] - offs[:, 1].mean()).mean())
                          if nok else 9999.0)
                med = float(np.median(scores)) if nok else 0.0
                crit = med - 0.004 * spread if nok >= 3 else -9
                combos.append((crit, med, spread, nok, fname, mu, s, row))

    combos.sort(key=lambda c: -c[0])
    print("\nTop 8 组合（判据 = 得分中位数 − 0.004×平移离散；需 ≥3 锚点过 0.25）：")
    print("   flip    mu    sigma  n_ok  median  spread(px)")
    for c in combos[:8]:
        print("   %-4s  %5.3f  %4.1f   %d    %.3f   %.0f"
              % (c[4], c[5], c[6], c[3], c[1], c[2]))

    crit, med, spread, nok, fname, mu, sig, row = combos[0]
    if crit <= 0:
        print("\n没有任何组合过线——这对图在当前方法下配不起来。")
        print("结论倾向：换 Montpezat / Alfalfa 配对试，或核查 WATSON 模型。")
        return 1

    print("\n==> 赢家：flip=%s  mu=%.3f  sigma=%.1f  (%d 锚点过线，"
          "中位得分 %.3f，平移离散 %.0f px)" % (fname, mu, sig, nok, med, spread))
    print("    各锚点：")
    for a, r in zip(anchors, row):
        print("      (%4.0f,%4.0f)  得分 %.3f  平移 (%+d,%+d)"
              % (a[0], a[1], r[0], r[1], r[2]))

    # ---- 3) 组装校正后的完整映射，零假设检验 ----
    Fm = flip_matrix(FLIPS[fname], wat["ns"], wat["nl"]) if FLIPS[fname] \
        else np.eye(3)
    good = [(a, r) for a, r in zip(anchors, row) if r[0] > 0.18]
    offs = np.array([[r[1], r[2]] for _a, r in good], np.float64)
    t = offs.mean(axis=0)
    T = np.float64([[1, 0, t[0]], [0, 1, t[1]], [0, 0, 1]])
    M = Fm @ H0 @ Zs[mu] @ T        # corrected(p) = wat(M·p)

    wfl = sample_watson(wf, M, (ns, nl))
    afs_best = afs[sig]
    # 全局 NCC（阴影掩模）
    lit_raw = (wgray > np.percentile(wgray, 25)).astype(np.uint8) * 255
    lit_raw = cv2.erode(lit_raw, np.ones((31, 31), np.uint8))
    cover = sample_watson(np.full_like(wgray, 255), M, (ns, nl)) > 127
    mask = cover & (sample_watson(lit_raw, M, (ns, nl)) > 127)
    mb = np.zeros_like(mask)
    mb[60:-60, 60:-60] = True
    mask &= mb
    n_after = ncc(afs_best, wfl, mask)
    nulls = []
    for dx, dy in ((60, 0), (-60, 0), (0, 60), (0, -60)):
        Tn = np.float64([[1, 0, dx], [0, 1, dy], [0, 0, 1]])
        w2 = sample_watson(wf, M @ np.linalg.inv(Tn), (ns, nl))
        nulls.append(ncc(afs_best, w2, mask))
    print("\n阴影掩模全局 NCC：校正后 %.4f" % n_after)
    print("零假设（±60 px）：%s  中位 %.4f  → %s"
          % (" ".join("%.3f" % v for v in nulls), float(np.median(nulls)),
             "有分辨力" if n_after - float(np.median(nulls)) > 0.05 else "仍无分辨力"))

    # ---- 4) 等效 standoff 与尺度判决 ----
    r_emp = edge_ratio(M, ns, nl)
    print("\n经验尺度比（含校正）= %.3f" % r_emp)
    best_d, best_err = None, 9e9
    for d in range(40, 601, 5):
        Hx, _ = homography(aci, wat, d)
        e = abs(edge_ratio(Hx, ns, nl) - r_emp)
        if e < best_err:
            best_d, best_err = d, e
    print("反解等效 standoff ≈ %d mm（几何解是 %d mm）" % (best_d, D0))
    if best_d <= 45 or best_d >= 590:
        print("  → 落在物理合理区间外：不是距离猜错，是相机模型尺度可疑"
              "（对照 WATSON ECM/FDR 两版模型的 H/V）")

    # ---- 5) 渲染 + 上色 ----
    warped = sample_watson(wcol8, M, (ns, nl))
    cover_c = sample_watson(np.full_like(wgray, 255), M, (ns, nl)) > 127
    lit_c = sample_watson(lit_raw, M, (ns, nl)) > 127
    colmask = cover_c & lit_c

    lab_w = cv2.cvtColor(warped, cv2.COLOR_RGB2LAB).astype(np.float32)
    lab_a = cv2.cvtColor(cv2.cvtColor(aci8, cv2.COLOR_GRAY2RGB),
                         cv2.COLOR_RGB2LAB).astype(np.float32)
    k = max(int(round(r_emp * 3)) | 1, 3)
    ab = cv2.GaussianBlur(lab_w[..., 1:], (k, k), 0)
    out = lab_a.copy()
    out[..., 1:] = ab
    out[~colmask, 1:] = 128.0
    color = cv2.cvtColor(np.clip(out, 0, 255).astype(np.uint8), cv2.COLOR_LAB2RGB)

    def draw(im):
        im2 = im.copy()
        for (cx, cy, _s), r in zip(anchors, row):
            cv2.circle(im2, (int(cx), int(cy)), 26, (255, 90, 0), 3)
            if r[0] > 0.25:
                cv2.circle(im2, (int(cx + r[1]), int(cy + r[2])), 26,
                           (0, 220, 0), 3)
                cv2.line(im2, (int(cx), int(cy)),
                         (int(cx + r[1]), int(cy + r[2])), (0, 220, 0), 2)
        return im2

    h = 760
    def rs(im):
        s = h / im.shape[0]
        return cv2.resize(im, (int(im.shape[1] * s), h))

    def panel(img, text):
        bar = np.full((44, img.shape[1], 3), 255, np.uint8)
        cv2.putText(bar, text, (10, 30), cv2.FONT_HERSHEY_SIMPLEX, 0.65,
                    (0, 0, 0), 2, cv2.LINE_AA)
        return np.vstack([bar, img])

    comp = np.hstack([
        panel(rs(draw(cv2.cvtColor(aci8, cv2.COLOR_GRAY2RGB))),
              "ACI + detected anchors"),
        panel(rs(draw(warped)),
              "WATSON @%dmm flip=%s mu=%.2f (corrected)" % (D0, fname, mu)),
        panel(rs(color), "colorized (shadow neutral)"),
    ])
    cv2.imwrite(os.path.join(OUT, "bellegarde_v2_compare.png"),
                cv2.cvtColor(comp, cv2.COLOR_RGB2BGR))
    cv2.imwrite(os.path.join(OUT, "bellegarde_v2_color.png"),
                cv2.cvtColor(color, cv2.COLOR_RGB2BGR))
    print("\n写出 out/bellegarde_v2_compare.png / bellegarde_v2_color.png")
    print("锚点圈：橙=检测位置，绿=匹配位置（重合=配准对）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
