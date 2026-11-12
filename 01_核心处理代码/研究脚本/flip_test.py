# -*- coding: utf-8 -*-
"""试像素坐标的朝向约定：CAHVOR 给的 (x,y) 是不是和数组的 (列,行) 同向？

动机：几何明明自洽（两条视轴在目标处相交），但按模型投出来两图对不上。
最可能的系统性错误就是坐标朝向——只要有一边是上下或左右反的，
投影出来的"形状"仍然落在画面内（所以之前的覆盖测试会通过），但内容完全错位。

做法：把 ACI 的 x/y 和 WATSON 的 x/y 各自允许翻转，4×4 = 16 种组合，
每种算一遍 NCC，看有没有哪一组显著高出一截。

用法：python flip_test.py [standoff_mm]
"""

from __future__ import annotations

import os
import sys

import numpy as np
import cv2

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from colorize_pair import (ACI_IMG, ACI_XML, WAT_IMG, WAT_XML, OUT,   # noqa: E402
                           feature, load_cam, read_bsq, to8)


def homog_flip(aci, wat, stand_mm, fa_x, fa_y, fw_x, fw_y):
    """带朝向开关的单应：ACI 像素 → WATSON 像素。"""
    d = stand_mm / 1000.0
    P0 = aci["C"] + d * aci["A"]
    n = aci["A"]
    src = np.array([[0, 0], [aci["ns"] - 1, 0],
                    [aci["ns"] - 1, aci["nl"] - 1], [0, aci["nl"] - 1]], np.float32)
    dst = []
    for x, y in src:
        xa = (aci["ns"] - 1 - x) if fa_x else x
        ya = (aci["nl"] - 1 - y) if fa_y else y
        dirv = np.cross(aci["H"] - xa * aci["A"], aci["V"] - ya * aci["A"])
        t = float((P0 - aci["C"]) @ n) / float(dirv @ n)
        v = (aci["C"] + t * dirv) - wat["C"]
        den = float(v @ wat["A"])
        wx = float(v @ wat["H"]) / den
        wy = float(v @ wat["V"]) / den
        if fw_x:
            wx = wat["ns"] - 1 - wx
        if fw_y:
            wy = wat["nl"] - 1 - wy
        dst.append((wx, wy))
    return cv2.getPerspectiveTransform(src, np.array(dst, np.float32)), np.array(dst)


def main():
    stand = float(sys.argv[1]) if len(sys.argv) > 1 else 129.0
    wat = load_cam(WAT_XML)
    aci = load_cam(ACI_XML)
    wr = read_bsq(WAT_IMG, wat)
    ar = read_bsq(ACI_IMG, aci)
    b = [i - 1 for i in wat["rgb"]]
    wcol = np.stack([to8(wr[b[i]]) for i in range(3)], axis=-1)
    aci8 = to8(ar[0])
    wf = feature(cv2.cvtColor(wcol, cv2.COLOR_RGB2GRAY))
    af = feature(aci8)
    # 分辨率对齐后再比
    wf = cv2.GaussianBlur(wf, (0, 0), 1.2)
    af = cv2.GaussianBlur(af, (0, 0), 1.2)
    m = np.zeros_like(aci8, bool)
    m[80:-80, 80:-80] = True
    print("standoff %.0f mm，16 种朝向组合的 NCC（越接近 1 越像）：" % stand)
    best = None
    for fx in (0, 1):
        for fy in (0, 1):
            row = []
            for wx in (0, 1):
                for wy in (0, 1):
                    Hm, dst = homog_flip(aci, wat, stand, fx, fy, wx, wy)
                    if dst[:, 0].min() < -50 or dst[:, 0].max() > wat["ns"] + 50 or \
                       dst[:, 1].min() < -50 or dst[:, 1].max() > wat["nl"] + 50:
                        row.append("  --  ")
                        continue
                    warp = cv2.warpPerspective(
                        wf, Hm, (aci["ns"], aci["nl"]),
                        flags=cv2.INTER_LINEAR | cv2.WARP_INVERSE_MAP)
                    av = af[m] - af[m].mean()
                    bv = warp[m] - warp[m].mean()
                    den = np.sqrt((av ** 2).sum() * (bv ** 2).sum())
                    s = float((av * bv).sum() / den) if den > 0 else float("nan")
                    row.append("%6.3f" % s)
                    if np.isfinite(s) and (best is None or s > best[0]):
                        best = (s, fx, fy, wx, wy)
            print("   ACI翻转(x=%d,y=%d)  %s" % (fx, fy, " ".join(row)))
    print("   （列依次是 WATSON翻转 x=0,y=0 / 0,1 / 1,0 / 1,1）")
    if best:
        print("=> 最好 %.4f，用 ACI(%d,%d) + WATSON(%d,%d)"
              % (best[0], best[1], best[2], best[3], best[4]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
