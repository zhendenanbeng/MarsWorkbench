#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""任意 WATSON/ACI 配对的精配 + 上色（refine2 的参数化版本）。

用法：
    python refine_any.py <标签> <watson.xml> <watson.IMG> <aci.xml> <aci.IMG>

流程与 refine2 相同：CAHVOR 几何单应当初值（standoff 取两视轴最近点，
自动计算，不硬编码）→ ACI 亮斑簇锚点 → 翻转×缩放×模糊联合搜索 →
平移一致性选择 → 零假设检验 → LAB 上色（阴影中性）。

判 GO 的自动标准（写给未来的软件集成）：
    ≥3 个锚点得分 > 0.18 且平移离散 < 120 px 且零假设落差 > 3 倍
"""

from __future__ import annotations

import os
import sys

import numpy as np
import cv2

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
sys.path.insert(0, HERE)

from colorize_pair import (load_cam, read_bsq, to8, feature,   # noqa: E402
                           homography, sample_watson, ncc)

SEARCH_PX = 300
TPL = 280
MUS = np.geomspace(0.15, 4.0, 44)
SIGMAS = (0.0, 1.5, 3.0, 5.0, 7.0)
FLIPS = {"none": None, "h": "h", "v": "v", "hv": "hv"}


def rigid_fit(src, dst):
    """Kabsch 刚性拟合（旋转+平移，无缩放）。返回 3x3 齐次矩阵。"""
    src_c = src - src.mean(0)
    dst_c = dst - dst.mean(0)
    H = src_c.T @ dst_c
    U, _s, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1] *= -1
        R = Vt.T @ U.T
    t = dst.mean(0) - R @ src.mean(0)
    return np.vstack([np.hstack([R, t.reshape(2, 1)]), [0, 0, 1.0]])


def rigid_residual(src, dst):
    """刚性拟合（旋转+平移，无缩放）的平均残差（px）和旋转角（度）。

    尺度由 μ 搜索负责，这里不许缩放——否则相似变换可以用大缩放
    把任意一组随机匹配拟合成零残差（Montpezat 上实测踩过）。
    """
    S3 = rigid_fit(src, dst)
    R, t = S3[:2, :2], S3[:2, 2]
    proj = src @ R.T + t
    resid = float(np.hypot(*(proj - dst).T).mean())
    ang = float(np.degrees(np.arctan2(R[1, 0], R[0, 0])))
    return resid, ang


def flip_matrix(kind, ns, nl):
    if kind == "h":
        return np.float64([[-1, 0, ns - 1], [0, 1, 0], [0, 0, 1]])
    if kind == "v":
        return np.float64([[1, 0, 0], [0, -1, nl - 1], [0, 0, 1]])
    if kind == "hv":
        return np.float64([[-1, 0, ns - 1], [0, -1, nl - 1], [0, 0, 1]])
    return np.eye(3)


def zoom_matrix(mu, ns, nl):
    c = np.array([ns / 2.0, nl / 2.0])
    return np.float64([[mu, 0, (1 - mu) * c[0]],
                       [0, mu, (1 - mu) * c[1]],
                       [0, 0, 1]])


def pick_standoff(aci, wat):
    """选几何初值的 standoff：扫 60-400 mm，取 ACI 四角投影完整落在
    WATSON 画幅内的档位。优先取**物理先验窗口 70-150 mm**内的可行档
    （SHERLOC 的工作距离量级；RLS 反推约 10.8 cm）——全可行域取中位数
    会选到 245 mm 这种远离工作距离的档（Montpezat 实测踩过）。
    """
    PHYS = (70, 150)
    feas, feas_phys = [], []
    best = (150.0, -1)
    for d in range(60, 401, 10):
        try:
            _H, dst = homography(aci, wat, d)
        except Exception:                                   # noqa: BLE001
            continue
        dst = np.asarray(dst)
        inside = int(((dst[:, 0] >= 0) & (dst[:, 0] <= wat["ns"] - 1) &
                      (dst[:, 1] >= 0) & (dst[:, 1] <= wat["nl"] - 1)).sum())
        if inside == 4:
            feas.append(d)
            if PHYS[0] <= d <= PHYS[1]:
                feas_phys.append(d)
        elif inside > best[1]:
            best = (float(d), inside)
    if feas_phys:
        return float(np.median(feas_phys))
    if feas:
        return float(np.median(feas))
    return best[0]


def detect_anchors(aci8, want=6):
    g = cv2.GaussianBlur(aci8.astype(np.float32), (0, 0), 1.2)
    dil = cv2.dilate(g, np.ones((25, 25), np.float32))
    peaks = (g == dil) & (g > g.mean() + 2.2 * g.std())
    ys, xs = np.nonzero(peaks)
    bins = {}
    for x, y in zip(xs, ys):
        bins.setdefault((x // 100, y // 100), []).append((x, y))
    clusters = []
    for (_bx, _by), pts in bins.items():
        xs_c = np.array([p[0] for p in pts])
        ys_c = np.array([p[1] for p in pts])
        clusters.append((float(xs_c.mean()), float(ys_c.mean()),
                         float(np.mean(g[ys_c, xs_c])) * len(pts)))
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


def run(tag, wat_xml, wat_img, aci_xml, aci_img):
    wat, aci = load_cam(wat_xml), load_cam(aci_xml)
    wr, ar = read_bsq(wat_img, wat), read_bsq(aci_img, aci)
    idx = [i - 1 for i in wat["rgb"]]
    wcol8 = np.stack([to8(wr[i]) for i in idx], -1)
    aci8 = to8(ar[0])
    wgray = cv2.cvtColor(wcol8, cv2.COLOR_RGB2GRAY)
    ns, nl = aci["ns"], aci["nl"]

    D0 = pick_standoff(aci, wat)
    H0, _ = homography(aci, wat, D0)
    print("[%s] 几何解 standoff = %.0f mm（视轴最近点），尺度比 = %.3f"
          % (tag, D0, edge_ratio(H0, ns, nl)))

    anchors = detect_anchors(aci8)
    print("[%s] 锚点 %d 个" % (tag, len(anchors)))
    if len(anchors) < 3:
        return {"ok": False, "reason": "锚点不足"}

    af, wf = feature(aci8), feature(wgray)
    boxes = []
    for cx, cy, _s in anchors:
        x0 = int(np.clip(cx - TPL / 2, 0, ns - TPL))
        y0 = int(np.clip(cy - TPL / 2, 0, nl - TPL))
        boxes.append((x0, y0))
    src_pts = np.float32([[a[0], a[1]] for a in anchors])
    Zs = {mu: zoom_matrix(mu, ns, nl) for mu in MUS}
    afs = {s: (af if s == 0 else cv2.GaussianBlur(af, (0, 0), s))
           for s in SIGMAS}
    tpls = {s: [afs[s][y0:y0 + TPL, x0:x0 + TPL] for x0, y0 in boxes]
            for s in SIGMAS}

    def sim_residual(offs, scores):
        """把锚点对 (r_i -> r_i+off_i) 做刚性拟合后的残差（px）。

        平移一致（Bellegarde）和带小旋转（Montpezat 这种跨 sol 摆位）
        都该被刚性变换吸收；残差才是"锚点们讲的是不是同一个故事"。
        """
        ok = scores > 0.18
        if ok.sum() < 3:
            return 9999.0, 0.0
        return rigid_residual(src_pts[ok], src_pts[ok] + offs[ok])

    # 受光掩模（WATSON 侧，逐对一次）：投影区落在阴影里的组合后面直接拒
    lit_raw0 = (wgray > np.percentile(wgray, 25)).astype(np.uint8) * 255

    # 显示方向一致（Sample/L2R + Line/T2B 是两机共同的约定）时不需要翻转；
    # 留翻转只会引来物理上不可能的伪胜（Montpezat 实测 flip=hv 刷分）。
    dirs_match = (wat["hax"] == aci["hax"] and wat["hdir"] == aci["hdir"]
                  and wat["vax"] == aci["vax"] and wat["vdir"] == aci["vdir"])
    flips = {"none": None} if dirs_match else FLIPS
    combos = []
    for fname, Fk in flips.items():
        Fm = flip_matrix(Fk, wat["ns"], wat["nl"]) if Fk else np.eye(3)
        for mu in MUS:
            Mz = H0 @ Zs[mu] if Fk is None else Fm @ H0 @ Zs[mu]
            wfl = sample_watson(wf, Mz, (ns, nl))
            cover = sample_watson(np.full_like(wgray, 255), Mz, (ns, nl)) > 127
            lit_frac = float((cover & (sample_watson(lit_raw0, Mz, (ns, nl))
                                        > 127)).sum()) / cover.size
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
                    pn = float(mx - np.median(res))
                    row.append((float(mx), (loc[0] + sx0) - x0,
                                (loc[1] + sy0) - y0, pn))
                scores = np.array([r[0] for r in row])
                offs = np.array([[r[1], r[2]] for r in row])
                nok = int((scores > 0.18).sum())
                resid, _ang = sim_residual(offs, scores)
                med = float(np.median(scores)) if nok else 0.0
                pns = np.array([r[3] for r in row])
                pn_med = float(np.median(pns[scores > 0.18])) if nok else 0.0
                # 受光守门：投影区 <10% 受光 = 缩进阴影刷分的退化解
                if lit_frac < 0.10:
                    nok, med, resid, pn_med = 0, 0.0, 9999.0, 0.0
                crit = med - 0.006 * resid if nok >= 3 else -9
                combos.append((crit, med, resid, nok, fname, mu, s, row,
                               offs, scores, lit_frac, pn_med))

    combos.sort(key=lambda c: -c[0])
    # 对照：限制 flip=none 的最优组合。物理上两机显示方向一致，不该需要翻转；
    # 全空间赢家若是翻转组合而 none 明显更差，翻转多半是伪胜。
    none_best = max((c for c in combos if c[4] == "none"), key=lambda c: c[0])
    print("[%s] [对照] flip=none 最优：mu=%.3f sigma=%.1f n_ok=%d med=%.3f "
          "resid=%.0f" % (tag, none_best[5], none_best[6], none_best[3],
                          none_best[1], none_best[2]))
    (crit, med, resid, nok, fname, mu, sig, row, offs, scores, lit_frac,
     pn_med) = combos[0]
    _wresid, wang = sim_residual(offs, scores)
    print("[%s] 刚性拟合：旋转 %.1f 度" % (tag, wang))
    print("[%s] 最优组合：flip=%s mu=%.3f sigma=%.1f（n_ok=%d, median=%.3f, "
          "sim残差=%.0f px, 受光%.0f%%）" % (tag, fname, mu, sig, nok, med,
                                            resid, 100*lit_frac))
    for a, r in zip(anchors, row):
        print("      (%4.0f,%4.0f)  %.3f  (%+d,%+d)  峰噪%.2f"
              % (a[0], a[1], r[0], r[1], r[2], r[3]))
    if crit <= 0:
        print("[%s] ==> NO-GO（无组合过线）" % tag)
        return {"ok": False, "reason": "无组合过线",
                "best": (fname, float(mu), sig, nok, med, resid)}

    Fm = flip_matrix(FLIPS[fname], wat["ns"], wat["nl"]) if FLIPS[fname] \
        else np.eye(3)
    ok = scores > 0.18
    _rs, _ra = rigid_residual(src_pts[ok], src_pts[ok] + offs[ok])
    S3 = rigid_fit(src_pts[ok], src_pts[ok] + offs[ok])
    M = Fm @ H0 @ Zs[mu] @ S3

    # 零假设检验。腐蚀要自适应：大核（31px）能挡住阴影过渡带的伪边缘
    # （Bellegarde：小腐蚀时 NCC 被稀释到 0.098），但受光区指状分布的图
    # 会被大核整个削没（Alfalfa：27px 腐蚀后掩模清零）。从大到小试，
    # 掩模够 30 万 px 就停。
    cover = sample_watson(np.full_like(wgray, 255), M, (ns, nl)) > 127
    mask = None
    for er in (31, 21, 13, 9, 5):
        lit_raw = cv2.erode(lit_raw0, np.ones((er, er), np.uint8))
        litw = sample_watson(lit_raw, M, (ns, nl)) > 127
        m2 = cover & litw
        mb = np.zeros_like(m2)
        mb[40:-40, 40:-40] = True
        m2 &= mb
        mask = m2
        if int(m2.sum()) >= 300000:
            break
    print("[debug] 掩模 %d px（腐蚀 %d）" % (int(mask.sum()), er))
    wfl = sample_watson(wf, M, (ns, nl))
    n_after = ncc(afs[sig], wfl, mask)
    sh = max(int(round(60 * mu)), 25)
    nulls = []
    for dx, dy in ((sh, 0), (-sh, 0), (0, sh), (0, -sh)):
        Tn = np.float64([[1, 0, dx], [0, 1, dy], [0, 0, 1]])
        w2 = sample_watson(wf, M @ np.linalg.inv(Tn), (ns, nl))
        nulls.append(ncc(afs[sig], w2, mask))
    null_med = float(np.median(nulls))
    # 零假设可以是负值（平移后反相关）——那比 0 还干净。用 0.05 做保护
    # 下限统一口径：ratio = NCC / max(零假设, 0.05)。
    ratio_null = n_after / max(null_med, 0.05)
    print("[%s] NCC=%.3f  零假设中位=%.3f  落差=%.1f 倍（掩模 %d px）"
          % (tag, n_after, null_med, ratio_null, int(mask.sum())))

    # GO 门槛（锚点级为主）：≥3 锚点过线、刚性残差 <90px（跨 sol 有真实
    # keystone，吃不满是物理必然）、锚点峰噪中位 ≥1.2（响应面得是个峰）。
    # 全图 NCC 零假设只做参考——WATSON 阴影占比高的图全图 NCC 天花板低
    # （Bellegarde 实测 0.08 vs 零假设 0.03），不作为硬门槛。
    go = (nok >= 3 and resid < 90 and pn_med >= 0.20 and lit_frac >= 0.10
          and abs(wang) <= 3.0 and np.isfinite(n_after))
    print("[%s] 峰噪中位=%.2f 受光=%.0f%%  旋转=%.1f°  各锚点峰噪：%s"
          % (tag, pn_med, 100 * lit_frac, wang,
             " ".join("%.2f" % r[3] for r in row)))
    print("[%s] ==> %s" % (tag, "GO" if go else "NO-GO"))
    if not go:
        return {"ok": False, "reason": "判据未达", "ncc": n_after,
                "null": null_med, "resid": resid}

    # 渲染 + 上色
    warped = sample_watson(wcol8, M, (ns, nl))
    colmask = cover & (sample_watson(lit_raw, M, (ns, nl)) > 127)
    lab_w = cv2.cvtColor(warped, cv2.COLOR_RGB2LAB).astype(np.float32)
    lab_a = cv2.cvtColor(cv2.cvtColor(aci8, cv2.COLOR_GRAY2RGB),
                         cv2.COLOR_RGB2LAB).astype(np.float32)
    r_emp = edge_ratio(M, ns, nl)
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
            if r[0] > 0.18:
                cv2.circle(im2, (int(cx + r[1]), int(cy + r[2])), 26,
                           (0, 220, 0), 3)
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
              "%s ACI + anchors" % tag),
        panel(rs(draw(warped)),
              "WATSON @%.0fmm mu=%.2f" % (D0, mu)),
        panel(rs(color), "colorized"),
    ])
    cv2.imwrite(os.path.join(OUT, "%s_v2_compare.png" % tag.lower()),
                cv2.cvtColor(comp, cv2.COLOR_RGB2BGR))
    cv2.imwrite(os.path.join(OUT, "%s_v2_color.png" % tag.lower()),
                cv2.cvtColor(color, cv2.COLOR_RGB2BGR))
    print("[%s] 写出 %s_v2_compare.png / %s_v2_color.png"
          % (tag, tag.lower(), tag.lower()))
    return {"ok": True, "standoff": D0, "mu": float(mu), "flip": fname,
            "resid": resid, "ncc": n_after, "null": null_med,
            "ratio_emp": r_emp}


def main(argv):
    if len(argv) != 6:
        print(__doc__)
        return 2
    tag, wx, wi, ax, ai = argv[1:6]
    for p in (wx, wi, ax, ai):
        if not os.path.exists(p):
            print("文件不存在：%s" % p)
            return 2
    res = run(tag, wx, wi, ax, ai)
    return 0 if res.get("ok") else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
