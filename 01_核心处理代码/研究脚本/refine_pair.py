#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Bellegarde 配对精配：掩掉阴影的多锚点模板匹配。

上一轮留下的三个问题，这个脚本一次回答：

  1. standoff 哪档对——几何解 129 mm 和联络表目测（129 对齐）、
     三联图目测（偏移像 90 mm）互相矛盾。两次目测至少一次是错的。
  2. 残余偏移多大——三联图目测说整体统一偏 (−16%, −17%)。
  3. 经验尺度比是多少——几何预测 0.35 px/px，纹理测量 1.01，差 2.85 倍。

方法：不再做全局 NCC（WATSON 约三分之一画幅是阴影，全局相关性被稀释成
噪声，上一轮已证实没有分辨力）。改用 ACI 上 5 个结构丰富的锚点区域，
每个 standoff 下把 WATSON 特征图投过来，在锚点邻域做模板匹配。
峰值和偏移是局部的、阴影进不了锚点就伤不了它。

用法：python refine_pair.py
输出：控制台报告 + out/bellegarde_refined_compare.png 等图
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

# ACI 锚点（分数坐标 + 模板边长 px）。位置来自目检：亮斑簇、亮斑对、
# 大亮块、上中亮区、画面中心。分布在四角与中心，便于分辨平移/尺度/旋转。
ANCHORS = [
    ("spot-cluster", 0.38, 0.31, 420),
    ("upper-mid",    0.55, 0.18, 360),
    ("spot-pair",    0.39, 0.71, 380),
    ("big-patch",    0.82, 0.77, 400),
    ("center",       0.55, 0.50, 420),
]
SEARCH_PX = 280                    # 模板匹配的邻域半径
SWEEP = list(range(70, 181, 2))    # standoff 扫描（mm，2 mm 步进）


def anchor_box(a, ns, nl):
    _, fx, fy, side = a
    cx, cy = int(fx * ns), int(fy * nl)
    x0 = min(max(cx - side // 2, 0), ns - side)
    y0 = min(max(cy - side // 2, 0), nl - side)
    return x0, y0, x0 + side, y0 + side


def match_anchor(wfl, tpl, box, ns, nl):
    """在 wfl 里 box 邻域搜 tpl，返回 (得分, dx, dy)。"""
    x0, y0, x1, y1 = box
    sx0, sy0 = max(0, x0 - SEARCH_PX), max(0, y0 - SEARCH_PX)
    sx1, sy1 = min(ns, x1 + SEARCH_PX), min(nl, y1 + SEARCH_PX)
    if sx1 - sx0 < x1 - x0 or sy1 - sy0 < y1 - y0:
        return 0.0, 0, 0
    res = cv2.matchTemplate(wfl[sy0:sy1, sx0:sx1], tpl, cv2.TM_CCOEFF_NORMED)
    _, mx, _, loc = cv2.minMaxLoc(res)
    # 响应面的峰噪比：峰再高、如果周围一样高，也只是"通用相似"
    med = float(np.median(res))
    return float(mx), (loc[0] + sx0) - x0, (loc[1] + sy0) - y0, mx - med


def main():
    wat, aci = load_cam(WAT_XML), load_cam(ACI_XML)
    wr, ar = read_bsq(WAT_IMG, wat), read_bsq(ACI_IMG, aci)
    idx = [i - 1 for i in wat["rgb"]]
    wcol8 = np.stack([to8(wr[i]) for i in idx], -1)
    aci8 = to8(ar[0])
    wgray = cv2.cvtColor(wcol8, cv2.COLOR_RGB2GRAY)

    af = feature(aci8)
    wf = feature(wgray)

    boxes = [anchor_box(a, aci["ns"], aci["nl"]) for a in ANCHORS]
    tpls = [af[y0:y1, x0:x1] for x0, y0, x1, y1 in boxes]

    # ---- 第 1 步：standoff 扫描，锚点模板匹配 ----
    print("standoff 扫描（每档 5 个锚点的匹配得分；判据取第 2 好的锚点）")
    print("  mm   " + "".join("%-14s" % a[0][:12] for a in ANCHORS) + "  2nd-best")
    table = {}
    for d in SWEEP:
        Hm, _ = homography(aci, wat, d)
        wfl = sample_watson(wf, Hm, (aci["ns"], aci["nl"]))
        row = []
        for box, tpl in zip(boxes, tpls):
            r = match_anchor(wfl, tpl, box, aci["ns"], aci["nl"])
            row.append(r)
        table[d] = row
        scores = sorted((r[0] for r in row), reverse=True)
        cells = "".join("  %5.3f(%+4d,%+4d)" % (r[0], r[1], r[2]) for r in row)
        print("  %3d %s   %5.3f" % (d, cells, scores[1]))

    # 稳健判据：第二名锚点的得分（单锚点高分可能是巧合）
    best_d = max(SWEEP, key=lambda d: sorted((r[0] for r in table[d]),
                                             reverse=True)[1])
    print("\n==> standoff 最优：%d mm" % best_d)
    print("    各锚点残余偏移（px）：")
    good = []
    for a, r in zip(ANCHORS, table[best_d]):
        print("      %-14s 得分 %.3f  峰噪 %.3f  偏移 (%+d, %+d)"
              % (a[0], r[0], r[3], r[1], r[2]))
        if r[0] > 0.20:
            good.append((a, r))

    # ---- 第 2 步：从锚点对估计位置校正 S ----
    src_pts, dst_pts = [], []
    for a, r in good:
        fx, fy = a[1], a[2]
        x0, y0 = int(fx * aci["ns"]), int(fy * aci["nl"])
        src_pts.append([x0, y0])               # ACI 里应在的位置 r
        dst_pts.append([x0 + r[1], y0 + r[2]])  # 投影图里实际在的位置 q
    src_pts = np.float32(src_pts)
    dst_pts = np.float32(dst_pts)

    offs = dst_pts - src_pts
    spread = float(np.hypot(offs[:, 0] - offs[:, 0].mean(),
                            offs[:, 1] - offs[:, 1].mean()).mean()) \
        if len(offs) > 1 else 0.0
    print("    偏移一致性：mean=(%+.0f,%+.0f) px  离散=%.0f px  (n=%d)"
          % (offs[:, 0].mean(), offs[:, 1].mean(), spread, len(offs)))

    Hm, _ = homography(aci, wat, best_d)
    if len(good) >= 2 and spread > 25:
        S, _ = cv2.estimateAffinePartial2D(src_pts, dst_pts,
                                           method=cv2.RANSAC)
        S = np.vstack([S, [0, 0, 1]])
        kind = "相似变换"
    else:
        t = offs.mean(axis=0)
        S = np.float32([[1, 0, t[0]], [0, 1, t[1]], [0, 0, 1]])
        kind = "纯平移"
    M = Hm @ S
    print("    位置校正：%s\n    校正矩阵 M =\n%s" % (kind, np.round(M, 3)))

    # ---- 第 3 步：阴影掩模下的全局 NCC（校正前 vs 后 vs 零假设）----
    lit_raw = (wgray > np.percentile(wgray, 25)).astype(np.uint8) * 255
    lit_raw = cv2.erode(lit_raw, np.ones((31, 31), np.uint8))
    cover = sample_watson(np.full_like(wgray, 255), M, (aci["ns"], aci["nl"])) > 127
    mask = cover & (sample_watson(lit_raw, M, (aci["ns"], aci["nl"])) > 127)
    m = np.zeros_like(mask)
    m[60:-60, 60:-60] = True
    mask &= m
    wfl = sample_watson(wf, M, (aci["ns"], aci["nl"]))
    ncc_after = ncc(af, wfl, mask)
    wfl0 = sample_watson(wf, Hm, (aci["ns"], aci["nl"]))
    mask0 = mask  # 同一阴影掩模下比较才有意义
    ncc_before = ncc(af, wfl0, mask0)
    nulls = []
    for dx, dy in ((60, 0), (-60, 0), (0, 60), (0, -60)):
        T = np.float32([[1, 0, dx], [0, 1, dy], [0, 0, 1]])
        w2 = sample_watson(wf, M @ np.linalg.inv(T), (aci["ns"], aci["nl"]))
        nulls.append(ncc(af, w2, mask))
    print("\n阴影掩模下全局 NCC：校正前 %.4f → 校正后 %.4f" % (ncc_before, ncc_after))
    print("零假设（±60 px 平移）：%s  中位 %.4f"
          % (" ".join("%.3f" % v for v in nulls), float(np.median(nulls))))

    # ---- 第 4 步：经验尺度比（了结 2.85 倍悬案） ----
    corners = np.float32([[0, 0], [aci["ns"] - 1, 0],
                          [aci["ns"] - 1, aci["nl"] - 1], [0, aci["nl"] - 1]])
    def edge_ratio(Hx):
        proj = cv2.perspectiveTransform(
            corners.reshape(-1, 1, 2), np.float64(Hx)).reshape(-1, 2)
        d_w = sum(np.hypot(*(proj[(i + 1) % 4] - proj[i])) for i in range(4))
        d_a = sum(np.hypot(*(corners[(i + 1) % 4] - corners[i])) for i in range(4))
        return d_w / d_a
    r_corr = edge_ratio(M)
    r_geom = edge_ratio(Hm)
    print("\n经验尺度比（WATSON px / ACI px，四边平均）：")
    print("    纯几何 %d mm：%.3f" % (best_d, r_geom))
    print("    几何 + 锚点校正：%.3f" % r_corr)
    print("    对照：纹理测量 1.01；上一轮几何预测 0.35")
    print("    若 ACI=6.66 µm/px → WATSON ≈ %.1f µm/px；"
          "若 ACI=10.1 → WATSON ≈ %.1f µm/px"
          % (r_corr * 6.66, r_corr * 10.1))

    # ---- 第 5 步：用 M 重新上色（阴影区色度置中性，不再涂黑） ----
    warped = sample_watson(wcol8, M, (aci["ns"], aci["nl"]))
    cover_c = sample_watson(np.full_like(wgray, 255), M,
                            (aci["ns"], aci["nl"])) > 127
    lit_c = sample_watson(lit_raw, M, (aci["ns"], aci["nl"])) > 127
    colmask = cover_c & lit_c

    lab_w = cv2.cvtColor(warped, cv2.COLOR_RGB2LAB).astype(np.float32)
    lab_a = cv2.cvtColor(cv2.cvtColor(aci8, cv2.COLOR_GRAY2RGB),
                         cv2.COLOR_RGB2LAB).astype(np.float32)
    # WATSON 像素粗于 ACI 时色度要模糊，避免色块锯齿；这里按尺度比定核
    k = max(int(round(r_corr * 3)) | 1, 3)
    ab = cv2.GaussianBlur(lab_w[..., 1:], (k, k), 0)
    out = lab_a.copy()
    out[..., 1:] = ab
    out[~colmask, 1:] = 128.0            # 阴影/画幅外：中性色，保留 ACI 亮度
    color = cv2.cvtColor(np.clip(out, 0, 255).astype(np.uint8), cv2.COLOR_LAB2RGB)

    def panel(img, text):
        bar = np.full((44, img.shape[1], 3), 255, np.uint8)
        cv2.putText(bar, text, (10, 30), cv2.FONT_HERSHEY_SIMPLEX, 0.7,
                    (0, 0, 0), 2, cv2.LINE_AA)
        return np.vstack([bar, img])

    def draw_anchors(img):
        im = img.copy()
        for (x0, y0, _x1, _y1), (a, _r) in zip(boxes, [(a, r) for a, r in
                                                       zip(ANCHORS, table[best_d])]):
            cv2.rectangle(im, (x0, y0), (x0 + 60, y0 + 60), (255, 80, 0), 3)
        return im

    h = 760
    def rs(im):
        s = h / im.shape[0]
        return cv2.resize(im, (int(im.shape[1] * s), h))

    comp = np.hstack([
        panel(draw_anchors(rs(cv2.cvtColor(aci8, cv2.COLOR_GRAY2RGB))),
              "ACI gray + anchors"),
        panel(draw_anchors(rs(warped)), "WATSON@%dmm corrected (+%s)" % (best_d, kind)),
        panel(rs(color), "colorized: L=ACI, ab=WATSON (shadow neutral)"),
    ])
    cv2.imwrite(os.path.join(OUT, "bellegarde_refined_compare.png"),
                cv2.cvtColor(comp, cv2.COLOR_RGB2BGR))
    cv2.imwrite(os.path.join(OUT, "bellegarde_refined_color.png"),
                cv2.cvtColor(color, cv2.COLOR_RGB2BGR))
    print("\n写出 out/bellegarde_refined_compare.png / bellegarde_refined_color.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
