# -*- coding: utf-8 -*-
"""把"能不能对上"这件事放到**物理尺度对等**的条件下来看。

为什么不沿用上一轮的做法：
  * 上一轮把 WATSON 放大 2.83 倍填到 ACI 网格上（几何要求如此），然后
    在 ACI 自己的像素尺度上比 NCC（blur σ ≤ 3 ACI px = 20 µm）。
    可是两台相机的**光学分辨率**并不相等：ACI 在岩面上约 25 µm，
    WATSON 约 60 µm。在 20 µm 的尺度上比，等于要求 WATSON 给出它
    根本没有的信息——NCC 自然平坦（0.32 vs 0.34 无分辨力）。
  * 这一轮改成：两边先带通到**同一个物理波段**（约 0.06–0.5 mm），
    再比。这是这两台相机**理论上还共享**的那部分信息。

三样东西：
  1. 几何落点标注图：ACI 的视场到底落在 WATSON 画面的哪一块
  2. 同物理尺度的比对图（WATSON 原生网格，两边都 21.4 µm/px）
  3. standoff 扫描 + 平移零假设（都在同一物理波段上）

用法：python visual_check.py [standoff_mm ...]
"""
from __future__ import annotations

import os
import sys

import numpy as np
import cv2

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from colorize_pair import (ACI_IMG, ACI_XML, WAT_IMG, WAT_XML, OUT,      # noqa: E402
                           homography, load_cam, read_bsq, to8,
                           sample_watson)

# 物理波段（岩面上的 mm）。两台相机共享的波段取 0.06–0.5 mm：
# 下限 = WATSON 的光学分辨率（约 60 µm），上限 = 观测量级。
BAND_LO_MM = 0.06
BAND_HI_MM = 0.50


def zscore(a, m=None):
    a = a.astype(np.float32)
    v = a[m] if m is not None else a
    s = float(v.std())
    return (a - float(v.mean())) / (s if s > 1e-9 else 1.0)


def bandpass(img, lo_px, hi_px):
    """带通：Gauss(lo) - Gauss(hi)，两个 σ 都以像素为单位。"""
    g = img.astype(np.float32)
    return cv2.GaussianBlur(g, (0, 0), lo_px) - cv2.GaussianBlur(g, (0, 0), hi_px)


def ncc(a, b, m):
    if m.sum() < 500:
        return float("nan")
    av = a[m] - a[m].mean()
    bv = b[m] - b[m].mean()
    den = np.sqrt((av ** 2).sum() * (bv ** 2).sum())
    return float((av * bv).sum() / den) if den > 0 else float("nan")


def mutual_info(a, b, m, bins=32):
    x, y = a[m], b[m]
    if x.size < 500:
        return float("nan")
    x = (x - x.min()) / max(float(np.ptp(x)), 1e-9)
    y = (y - y.min()) / max(float(np.ptp(y)), 1e-9)
    h, _, _ = np.histogram2d(x, y, bins=bins, range=((0, 1), (0, 1)))
    p = h / h.sum()
    px = p.sum(axis=1)[:, None]
    py = p.sum(axis=0)[None, :]
    nz = p > 0
    return float((p[nz] * np.log2(p[nz] / (px @ py)[nz])).sum())


def label(img, text, h=36):
    bar = np.full((h, img.shape[1], 3), 255, np.uint8)
    cv2.putText(bar, text, (8, h - 10), cv2.FONT_HERSHEY_SIMPLEX, 0.62,
                (0, 0, 0), 2, cv2.LINE_AA)
    return np.vstack([bar, img])


def gray_rgb(g):
    return cv2.cvtColor(g, cv2.COLOR_GRAY2RGB)


def main():
    stands = [float(v) for v in sys.argv[1:]] or [100.0, 129.0, 160.0]

    wat = load_cam(WAT_XML)
    aci = load_cam(ACI_XML)
    wr = read_bsq(WAT_IMG, wat)
    ar = read_bsq(ACI_IMG, aci)
    b = [i - 1 for i in wat["rgb"]]
    wcol = np.stack([to8(wr[b[i]]) for i in range(3)], axis=-1)
    wgray = cv2.cvtColor(wcol, cv2.COLOR_RGB2GRAY)
    aci8 = to8(ar[0])

    # --- 岩面上的像素尺度（几何给出的是视场比；绝对尺度用 RLS 反推的
    #     ACI 6.66 µm/px 当锚，这样"物理波段"才落得下来）---
    aci_um = 6.66                      # RLS: HDR_500 的 500 µm 步长 / 75.1 px
    print("ACI 锚定 %.2f µm/px（RLS HDR_500 反推）；WATSON 按几何比推算" % aci_um)

    # 两个网格上的带通 σ（像素）
    def sigmas(px_per_mm):
        lo = max(1.0, BAND_LO_MM * px_per_mm / 2.355)
        hi = max(lo * 1.5, BAND_HI_MM * px_per_mm / 2.355)
        return lo, hi

    aci_ppm = 1000.0 / aci_um                    # px per mm on the rock
    sla, sha = sigmas(aci_ppm)
    afull = zscore(bandpass(aci8, sla, sha))
    print("ACI  %.3f px/mm → 带通 σ = %.1f–%.1f px" % (aci_ppm, sla, sha))

    # ---------- 1. 几何落点标注 ----------
    annot = cv2.cvtColor(wgray, cv2.COLOR_GRAY2RGB)
    colors = [(0, 255, 255), (0, 255, 0), (255, 128, 0)]
    for d, col in zip(stands, colors):
        Hm, dst = homography(aci, wat, d)
        pts = dst.astype(np.int32).reshape(-1, 1, 2)
        cv2.polylines(annot, [pts], True, col, 3)
        cv2.putText(annot, "%d mm" % d, (int(dst[0, 0]) + 6, int(dst[0, 1]) + 34),
                    cv2.FONT_HERSHEY_SIMPLEX, 1.1, col, 3, cv2.LINE_AA)
    cv2.imwrite(os.path.join(OUT, "vc_annot.png"),
                cv2.cvtColor(annot, cv2.COLOR_RGB2BGR))
    print("写到 out/vc_annot.png（ACI 视场在 WATSON 里的落点）")

    # WATSON 像素尺度：由 ACI 单应推出（每个 standoff 不同），用 129 mm 那档
    Hm129, dst129 = homography(aci, wat, 129.0)
    wat_um = aci_um / 0.353                      # ≈ 18.9 µm/px（几何比）
    wat_ppm = 1000.0 / wat_um
    slw, shw = sigmas(wat_ppm)
    print("WATSON %.3f px/mm（几何比 0.353）→ 带通 σ = %.1f–%.1f px"
          % (wat_ppm, slw, shw))

    # ---------- 2. standoff 扫描（同一物理波段）----------
    wfull = zscore(bandpass(wgray, slw, shw))
    # WATSON 的"阴影/暗区"掩膜：局部均值明显偏低的地方
    loc = cv2.blur(wgray.astype(np.float32), (151, 151))
    shadow = loc < 0.55 * float(np.median(loc))
    print("WATSON 暗区占画面 %.0f%%" % (100.0 * shadow.mean()))

    print("\nstandoff 扫描（两边都带通到 0.06–0.5 mm；掩膜=有效覆盖∩非暗区）")
    print("  stand(mm)  覆盖率%%   NCC      NCC(含暗区)   MI(bit)   MI(含暗区)")
    best = None
    for d in [80, 90, 100, 110, 120, 129, 140, 160, 180, 210, 250, 300]:
        Hm, _ = homography(aci, wat, d)
        warpf = sample_watson(wfull, Hm, (aci["ns"], aci["nl"]))
        valid = sample_watson(np.full((wat["nl"], wat["ns"]), 255, np.uint8), Hm,
                              (aci["ns"], aci["nl"]),
                              ) > 0
        # 阴影也投到 ACI 网格上
        sh = sample_watson(shadow.astype(np.uint8) * 255, Hm,
                           (aci["ns"], aci["nl"])) > 127
        m_all = valid.copy()
        m_all[:60, :] = m_all[-60:, :] = False
        m_all[:, :60] = m_all[:, -60:] = False
        m_lit = m_all & ~sh
        cov = 100.0 * valid.sum() / valid.size
        n1 = ncc(afull, warpf, m_lit)
        n2 = ncc(afull, warpf, m_all)
        i1 = mutual_info(afull, warpf, m_lit)
        i2 = mutual_info(afull, warpf, m_all)
        print("  %7.0f %9.0f %8.3f %12.3f %9.4f %11.4f"
              % (d, cov, n1, n2, i1, i2))
        if np.isfinite(n1) and (best is None or n1 > best[0]):
            best = (n1, d)
    if best:
        print("  → 最优 standoff %.0f mm（NCC %.3f）" % (best[1], best[0]))

    # ---------- 3. 平移零假设 + 出图 ----------
    for d in stands:
        Hm, dst = homography(aci, wat, d)
        warpf = sample_watson(wfull, Hm, (aci["ns"], aci["nl"]))
        sh = sample_watson(shadow.astype(np.uint8) * 255, Hm,
                           (aci["ns"], aci["nl"])) > 127
        valid = sample_watson(np.full((wat["nl"], wat["ns"]), 255, np.uint8), Hm,
                              (aci["ns"], aci["nl"])) > 0
        m = valid & ~sh
        m[:60, :] = m[-60:, :] = False
        m[:, :60] = m[:, -60:] = False

        # 平移扫描（在 ACI 网格上，±300 px，步长 15）
        offs = list(range(-300, 301, 15))
        grid = np.full((len(offs), len(offs)), np.nan)
        for i, dy in enumerate(offs):
            for j, dx in enumerate(offs):
                T = np.array([[1, 0, dx], [0, 1, dy], [0, 0, 1.0]])
                w2 = cv2.warpPerspective(wfull, T @ Hm, (aci["ns"], aci["nl"]),
                                         flags=cv2.INTER_LINEAR)
                grid[i, j] = ncc(afull, w2, m)
        k = np.unravel_index(np.nanargmax(grid), grid.shape)
        print("\nstandoff %.0f mm 平移扫描 ±300 px：峰 %.3f 在 (dx=%+d, dy=%+d)；"
              "中位 %.3f；(0,0) %.3f"
              % (d, grid[k], offs[k[1]], offs[k[0]], float(np.nanmedian(grid)),
                 grid[offs.index(0), offs.index(0)]))

        # 图 1：棋盘格（结构图，80 px 块）
        blk = 80
        checker = afull.copy()
        for i in range(0, aci["nl"], blk):
            for j in range(0, aci["ns"], blk):
                if ((i // blk) + (j // blk)) % 2:
                    checker[i:i + blk, j:j + blk] = \
                        warpf[i:i + blk, j:j + blk]
        ch = np.clip(128 + 45 * checker, 0, 255).astype(np.uint8)
        ch[~valid] = 0
        cv2.imwrite(os.path.join(OUT, "vc_checker_%d.png" % d),
                    cv2.cvtColor(ch, cv2.COLOR_GRAY2BGR))

        # 图 2：红=ACI，绿=WATSON（同物理波段、各自零均值单位方差）
        rg = np.zeros((aci["nl"], aci["ns"], 3), np.uint8)
        rg[..., 0] = np.clip(128 + 60 * afull, 0, 255)
        rg[..., 1] = np.clip(128 + 60 * warpf, 0, 255)
        rg[..., 2] = 128
        rg[~valid] = 0
        cv2.imwrite(os.path.join(OUT, "vc_anaglyph_%d.png" % d),
                    cv2.cvtColor(rg, cv2.COLOR_RGB2BGR))
        print("  写 out/vc_checker_%d.png, out/vc_anaglyph_%d.png" % (d, d))

    # ---------- 4. WATSON 原生网格上的比对（两边都 21 µm/px，物理上最公平）----------
    aci_in_w = cv2.warpPerspective(afull.astype(np.float32), Hm129,
                                   (wat["ns"], wat["nl"]),
                                   flags=cv2.INTER_LINEAR)
    # ACI 到 WATSON 网格后，ACI 自己的 PSF 只剩 1.2 px，WATSON 是 3 px：把 ACI 补糊
    aci_match = cv2.GaussianBlur(aci_in_w, (0, 0), 2.7)
    cov = aci_in_w != 0
    aci_match[~cov] = 0
    x0, x1 = int(dst129[:, 0].min()), int(dst129[:, 0].max())
    y0, y1 = int(dst129[:, 1].min()), int(dst129[:, 1].max())
    pad = 12
    x0, y0 = max(0, x0 - pad), max(0, y0 - pad)
    x1, y1 = min(wat["ns"], x1 + pad), min(wat["nl"], y1 + pad)
    wc = 255 - wgray[y0:y1, x0:x1]                    # 让结构在黑底上看得清
    ac = 255 - np.clip(128 + 90 * aci_match[y0:y1, x0:x1], 0, 255).astype(np.uint8)
    out = np.hstack([
        label(cv2.cvtColor(wc, cv2.COLOR_GRAY2RGB),
              "WATSON native %.1f um/px  [%d..%d, %d..%d]" % (wat_um, x0, x1, y0, y1)),
        label(cv2.cvtColor(ac, cv2.COLOR_GRAY2RGB),
              "ACI mapped to WATSON grid, blurred to match")])
    cv2.imwrite(os.path.join(OUT, "vc_native_129.png"),
                cv2.cvtColor(out, cv2.COLOR_RGB2BGR))
    print("写到 out/vc_native_129.png（%.0fx%.0f）" % (out.shape[1], out.shape[0]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
