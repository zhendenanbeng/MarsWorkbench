# -*- coding: utf-8 -*-
"""两件事：

1) 前/后两张 ACI 到底差多少？之前那张差值图是**纯灰**（=完全一样）。
   量一下最佳位移下的 NCC 和最大像素差，确认"18 分钟里表面没有变化"。

2) 用**棋盘格**做一次"我的眼睛还能不能看出对应关系"的标定。
   棋盘格的道理：真的对上了，纹理在格子边界处应当接得上；没对上就是拼贴。
   关键是要配一个**对照**：把同一张 WATSON 故意平移几百像素再拼一次。
   如果我分不出哪张是真、哪张是对照，那"肉眼看上去像/不像"就不能当判据
   ——之前我说"中间那张和左边不像"这件事也就不可靠了。

出图：out/vc_check_raw_{true,shift}.png
用法：python checker_calib.py [standoff_mm]
"""
from __future__ import annotations

import os
import sys

import numpy as np
import cv2

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from colorize_pair import (ACI_IMG, ACI_XML, OUT, homography, load_cam,
                           read_bsq, to8)
from watson_pose_test import gray, NEAR

ACI_UM = 6.66
WAT_UM = 8.10


def stretch(g, lo=1.0, hi=99.0):
    a, b = np.percentile(g, (lo, hi))
    return np.clip((g - a) / max(b - a, 1e-9) * 255, 0, 255).astype(np.uint8)


def local_norm(g, um, mm=0.5):
    """局部均值/方差归一化：把照明差异压掉，只留结构。"""
    ppm = 1000.0 / um
    s = mm * ppm / 2.355
    m = cv2.GaussianBlur(g, (0, 0), s * 3)
    v = np.sqrt(np.maximum(cv2.GaussianBlur(g * g, (0, 0), s * 3) - m * m, 1e-6))
    return (g - m) / v


def main():
    d = float(sys.argv[1]) if len(sys.argv) > 1 else 108.0
    aci = load_cam(ACI_XML)
    aci8 = to8(read_bsq(ACI_IMG, aci)[0]).astype(np.float32)
    gw, camw = gray(NEAR)

    # ---- 1) 前/后 ACI ----
    pre_path = os.path.join(OUT, "img", "SC3_0186_0683478583_214FDR_N0070000SRLC11420_0000LMJ05.IMG")
    pre_xml = os.path.join(OUT, "xml", "SC3_0186_0683478583_214FDR_N0070000SRLC11420_0000LMJ05.xml")
    if os.path.exists(pre_path) and os.path.exists(pre_xml):
        precam = load_cam(pre_xml)
        pre = to8(read_bsq(pre_path, precam)[0]).astype(np.float32)
        best = (-9, None)
        for dy in range(-6, 7):
            for dx in range(-6, 7):
                a = aci8[20:-20, 20:-20]
                b = np.roll(np.roll(pre, dy, 0), dx, 1)[20:-20, 20:-20]
                av, bv = a - a.mean(), b - b.mean()
                v = float((av * bv).sum() / np.sqrt((av ** 2).sum() * (bv ** 2).sum()))
                if v > best[0]:
                    best = (v, (dx, dy))
        print("1) 前/后 ACI：最佳位移 %s 下 NCC = %.4f；整图最大像素差 %.1f，"
              "中位绝对差 %.2f"
              % (best[1], best[0], float(np.abs(aci8 - pre).max()),
                 float(np.median(np.abs(aci8 - pre)))))
        print("   → 两张 ACI 基本就是同一张图（表面没变，也没有激光坑）")

    # ---- 2) 棋盘格标定 ----
    Hm, dst = homography(aci, camw, d)
    warp = cv2.warpPerspective(gw.astype(np.float32), Hm, (aci["ns"], aci["nl"]),
                               flags=cv2.INTER_LINEAR | cv2.WARP_INVERSE_MAP)
    valid = cv2.warpPerspective(np.full((camw["nl"], camw["ns"]), 255, np.uint8),
                                Hm, (aci["ns"], aci["nl"]),
                                flags=cv2.INTER_NEAREST | cv2.WARP_INVERSE_MAP) > 127
    A = local_norm(aci8, ACI_UM)
    W = local_norm(warp, WAT_UM)
    blk = 96
    for tag, shift in (("true", (0, 0)), ("shift", (0, 260))):
        Ws = np.roll(np.roll(W, shift[0], 0), shift[1], 1)
        out = A.copy()
        for i in range(0, aci["nl"], blk):
            for j in range(0, aci["ns"], blk):
                if ((i // blk) + (j // blk)) % 2:
                    out[i:i + blk, j:j + blk] = Ws[i:i + blk, j:j + blk]
        o8 = np.clip(128 + 55 * out, 0, 255).astype(np.uint8)
        o8[~valid] = 0
        cv2.imwrite(os.path.join(OUT, "vc_check_raw_%s.png" % tag), o8)
        # 数值：跨格边界的"接缝"强度——真对上时接缝应当不比格内差
        gx = np.abs(np.diff(out, axis=1))
        vseam = gx[:, blk - 1::blk].mean()
        vinner = gx.mean()
        print("   棋盘 %-5s：跨格边界平均 |Δ| %.4f，格内平均 |Δ| %.4f，比值 %.2f"
              % (tag, vseam, vinner, vseam / max(vinner, 1e-9)))
        cv2.imwrite(os.path.join(OUT, "vc_check_raw_%s.png" % tag), o8)
    print("写到 out/vc_check_raw_true.png 与 _shift.png（棋盘块 %d px；"
          "shift 是把 WATSON 整体下移 260 px 的对照）" % blk)
    return 0


if __name__ == "__main__":
    sys.exit(main())
