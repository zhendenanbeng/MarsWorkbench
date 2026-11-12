# -*- coding: utf-8 -*-
"""在"几何预测的位置"附近把平移空间搜一遍，看有没有真正的相关峰。

为什么需要这一步：`colorize_pair.py` 里直接用几何预测的单应算 NCC，
得到 0.34；但把配准**平移 60–200 像素**后 NCC 还是 0.32——说明那个数
根本没有分辨力（两张岩石纹理图，模糊之后随便怎么对都差不多）。

所以改成：固定几何给出的尺度/旋转，只扫平移 (dx, dy)，画出一张 NCC 曲面。
判读：
  * 曲面在 (0,0) 附近有尖峰  → 几何对，配准成立
  * 峰在别处且很尖          → 几何有系统偏差，用那个偏移当修正量
  * 到处都差不多（无峰）      → 这两张图没有对应关系，别再往下做

用法：python ncc_search.py [standoff_mm]
输出：out/ncc_surface.txt
"""

from __future__ import annotations

import os
import sys

import numpy as np
import cv2

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from colorize_pair import (ACI_IMG, ACI_XML, WAT_IMG, WAT_XML, OUT,   # noqa: E402
                           feature, homography, load_cam, read_bsq, to8)

STEP = 8            # 平移搜索步长（全分辨率像素）
SPAN = 200          # 搜索范围 ±SPAN 像素
LEVEL = 4           # 降采样倍数


def mutual_info(a, b, mask, bins=48):
    """互信息（bit）。跨传感器配准的标准判据：只看"两边灰度是否共同变化"，
    不在乎一边亮一边暗（ACI 是环形照明，WATSON 是蓝光 LED，明暗关系本来就不同）。
    对得准时联合直方图集中 → MI 高；对不准时散开 → MI 低。"""
    x = a[mask]
    y = b[mask]
    if x.size < 2000:
        return float("nan")
    # 各自按分位数归一化到 0-1，避免两边的动态范围差异影响分箱
    # numpy 2.x 去掉了 ndarray.ptp()，用 np.ptp()
    x = (x - x.min()) / max(float(np.ptp(x)), 1e-9)
    y = (y - y.min()) / max(float(np.ptp(y)), 1e-9)
    h, _, _ = np.histogram2d(x, y, bins=bins, range=((0, 1), (0, 1)))
    p = h / h.sum()
    px = p.sum(axis=1)[:, None]
    py = p.sum(axis=0)[None, :]
    nz = p > 0
    return float((p[nz] * np.log2(p[nz] / (px @ py)[nz])).sum())


def main():
    stand = float(sys.argv[1]) if len(sys.argv) > 1 else 100.0
    wat = load_cam(WAT_XML)
    aci = load_cam(ACI_XML)
    wr = read_bsq(WAT_IMG, wat)
    ar = read_bsq(ACI_IMG, aci)
    b = [i - 1 for i in wat["rgb"]]
    wcol = np.stack([to8(wr[b[i]]) for i in range(3)], axis=-1)
    aci8 = to8(ar[0])

    wf = feature(cv2.cvtColor(wcol, cv2.COLOR_RGB2GRAY))
    af = feature(aci8)
    # 分辨率对齐：两侧都降到 1/LEVEL，再把细的糊一下，避免"模糊换分"
    wf = cv2.GaussianBlur(cv2.resize(wf, None, fx=1 / LEVEL, fy=1 / LEVEL,
                                     interpolation=cv2.INTER_AREA), (0, 0), 1.2)
    af = cv2.GaussianBlur(cv2.resize(af, None, fx=1 / LEVEL, fy=1 / LEVEL,
                                     interpolation=cv2.INTER_AREA), (0, 0), 1.2)
    h, w = af.shape
    m = np.zeros_like(af, bool)
    m[20:-20, 20:-20] = True

    H0, dst = homography(aci, wat, stand)
    H0 = H0 / LEVEL
    H0[2, 2] = 1.0
    if dst[:, 0].min() < 0 or dst[:, 0].max() >= wat["ns"] or \
       dst[:, 1].min() < 0 or dst[:, 1].max() >= wat["nl"]:
        print("standoff %.0f mm 时 ACI 投出 WATSON 画面外" % stand)
        return 1

    offs = list(range(-SPAN, SPAN + 1, STEP))
    grid = np.full((len(offs), len(offs)), np.nan)
    mi = np.full((len(offs), len(offs)), np.nan)
    for i, dy in enumerate(offs):
        for j, dx in enumerate(offs):
            T = np.array([[1, 0, dx / LEVEL], [0, 1, dy / LEVEL], [0, 0, 1.0]])
            warp = cv2.warpPerspective(wf, T @ H0, (w, h), flags=cv2.INTER_LINEAR)
            av = af[m] - af[m].mean()
            bv = warp[m] - warp[m].mean()
            den = np.sqrt((av ** 2).sum() * (bv ** 2).sum())
            grid[i, j] = float((av * bv).sum() / den) if den > 0 else np.nan
            mi[i, j] = mutual_info(af, warp, m)

    km = np.unravel_index(np.nanargmax(mi), mi.shape)
    print("互信息：峰值 %.4f bit 在 (dx=%+d, dy=%+d)；中位 %.4f；(0,0) %.4f；峰-中位 %.4f"
          % (mi[km], offs[km[1]], offs[km[0]], np.nanmedian(mi),
             mi[offs.index(0), offs.index(0)], mi[km] - np.nanmedian(mi)))
    print("互信息曲面（行=dy，列=dx，每格 %d px；数字是 (MI - 中位)×1000）：" % STEP)
    print("      " + " ".join("%4d" % o for o in offs))
    for i, dy in enumerate(offs):
        row = " ".join("%4.0f" % ((mi[i, j] - np.nanmedian(mi)) * 1000)
                       for j in range(len(offs)))
        print("%4d  %s" % (dy, row))
    print()

    k = np.unravel_index(np.nanargmax(grid), grid.shape)
    peak = grid[k]
    med = float(np.nanmedian(grid))
    print("standoff %.0f mm，搜索 ±%d px（步长 %d，%d 倍降采样）"
          % (stand, SPAN, STEP, LEVEL))
    print("  NCC 峰值 %.4f 出现在 (dx=%+d, dy=%+d)" % (peak, offs[k[1]], offs[k[0]]))
    print("  中位数 %.4f  最小值 %.4f  峰值/中位差 %.4f"
          % (med, float(np.nanmin(grid)), peak - med))
    print("  (0,0) 处 %.4f" % grid[offs.index(0), offs.index(0)])
    print("\n曲面（行=dy，列=dx，每格 %d px；数字是 (NCC-中位)×1000）：" % STEP)
    header = "      " + " ".join("%4d" % o for o in offs)
    print(header)
    for i, dy in enumerate(offs):
        row = " ".join("%4.0f" % ((grid[i, j] - med) * 1000) for j in range(len(offs)))
        print("%4d  %s" % (dy, row))

    with open(os.path.join(OUT, "ncc_surface.txt"), "w", encoding="utf-8") as fh:
        fh.write("standoff=%g peak=%.4f at dx=%d dy=%d median=%.4f at00=%.4f\n"
                 % (stand, peak, offs[k[1]], offs[k[0]], med,
                    grid[offs.index(0), offs.index(0)]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
