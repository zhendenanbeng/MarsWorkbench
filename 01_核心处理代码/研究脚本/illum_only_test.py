# -*- coding: utf-8 -*-
"""决定性对照：**同一位姿、只换照明/曝光**的两帧 WATSON，还能不能对上？

A = SIF_0186_0683466024_335FDR_..._0000LMJ (1648×1200 全画幅)
B = SIF_0186_0683466037_832FDR_..._0000LUJ (子帧 first_line=17 first_sample=33, 1584×1184)

两者同属 srlc00621 落位（C 应当相同）。后缀 LMJ / LUJ 是照明/曝光配置，
所以这一对比**只差照明**。判读：
  * 若 NCC 仍接近 1  → 换照明不破坏相关，那 ACI↔WATSON 的零就不能用照明解释；
  * 若 NCC 掉到 0.1–0.3 → **照明差异本身就足以毁掉相关**，
    于是"ACI 和 WATSON 对不上"完全不能说明配对有问题。

（离线运行：两张图与标签都已在 out/ 里。）
"""
from __future__ import annotations

import os
import re
import sys

import numpy as np
import cv2

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from colorize_pair import OUT, read_bsq, to8

IMG = os.path.join(OUT, "img")
XML = os.path.join(OUT, "xml")
A_XML = os.path.join(XML, "SIF_0186_0683466024_335FDR_N0070000SRLC00621_0000LMJ04.xml")
A_IMG = os.path.join(IMG, "SIF_0186_0683466024_335FDR_N0070000SRLC00621_0000LMJ04.IMG")
B_XML = os.path.join(XML, "SIF_0186_0683466037_832FDR_N0070000SRLC00621_0000LUJ05.xml")
B_IMG = os.path.join(IMG, "SIF_0186_0683466037_832FDR_N0070000SRLC00621_0000LUJ05.IMG")


def parse(path):
    t = open(path, encoding="utf-8", errors="replace").read()

    def vec(tag):
        i = t.find("<geom:%s>" % tag)
        seg = t[i:t.index("</geom:%s>" % tag, i)]
        return np.array([float(x) for x in
                         re.findall(r">\s*(-?[\d.eE+-]+)\s*<", seg)][:3])
    c = dict(C=vec("Vector_Center"), A=vec("Vector_Axis"),
             H=vec("Vector_Horizontal"), V=vec("Vector_Vertical"))
    c["ns"] = int(re.search(r"<img:samples>\s*(\d+)", t).group(1))
    c["nl"] = int(re.search(r"<img:lines>\s*(\d+)", t).group(1))
    c["fl"] = int(re.search(r"<img:first_line>\s*(\d+)", t).group(1))
    c["fs"] = int(re.search(r"<img:first_sample>\s*(\d+)", t).group(1))
    m = re.search(r"<img:exposure_duration unit=\"ms\">([\d.]+)", t)
    c["exp"] = float(m.group(1)) if m else None
    m = re.search(r"<img:exposure_count>\s*(\d+)", t)
    c["nexp"] = int(m.group(1)) if m else 1
    states = dict(re.findall(r"<img:device_name>([^<]+)</img:device_name>\s*"
                             r"<img:device_state>([^<]+)</img:device_state>", t))
    c["led"] = {k: v for k, v in states.items() if "LED" in k}
    m = re.search(r"<img:bands>\s*(\d+)", t)
    c["bands"] = int(m.group(1)) if m else 3
    rgb = []
    for ch in ("red", "green", "blue"):
        m = re.search(r"<disp:%s_channel_band>\s*(\d+)" % ch, t)
        rgb.append(int(m.group(1)) if m else 1)
    c["rgb"] = rgb
    return c


def gray(path, cam):
    r = read_bsq(path, cam)
    b = [k - 1 for k in cam["rgb"]]
    g = cv2.cvtColor(np.stack([to8(r[b[k]]) for k in range(3)], axis=-1),
                     cv2.COLOR_RGB2GRAY)
    return g[:cam["nl"], :cam["ns"]].astype(np.float32)


def main():
    ca, cb = parse(A_XML), parse(B_XML)
    for nm, c in (("LMJ", ca), ("LUJ", cb)):
        print("%s：%dx%d  子帧起点 (%d,%d)  曝光 %.1f ms × %d  LED %s"
              % (nm, c["ns"], c["nl"], c["fs"], c["fl"], c["exp"] or -1,
                 c["nexp"], c["led"]))
        print("     C=%s  A=%s" % (np.round(c["C"], 5), np.round(c["A"], 5)))
    print("两帧光心相距 %.4f mm；视轴夹角 %.5f°"
          % (np.linalg.norm(ca["C"] - cb["C"]) * 1000,
             np.degrees(np.arccos(np.clip(ca["A"] @ cb["A"], -1, 1)))))

    ga = gray(A_IMG, ca)
    gb = gray(B_IMG, cb)
    fl, fs = cb["fl"] - 1, cb["fs"] - 1
    crop = ga[fl:fl + cb["nl"], fs:fs + cb["ns"]]
    print("A 的对应块 %s  vs  B %s" % (crop.shape, gb.shape))
    if crop.shape != gb.shape:
        print("形状不一致，取公共部分")
        h = min(crop.shape[0], gb.shape[0])
        w = min(crop.shape[1], gb.shape[1])
        crop, gb = crop[:h, :w], gb[:h, :w]

    def bp(g, lo=1.5, hi=15.0):
        return cv2.GaussianBlur(g, (0, 0), lo) - cv2.GaussianBlur(g, (0, 0), hi)

    def ncc(a, b):
        a = a - a.mean()
        b = b - b.mean()
        return float((a * b).sum() / np.sqrt((a ** 2).sum() * (b ** 2).sum()))

    xa, xb = bp(crop), bp(gb)
    m = 80
    print("\n只差照明/曝光的两帧：零位移 NCC = %.3f"
          % ncc(xa[m:-m, m:-m], xb[m:-m, m:-m]))
    best = (-9, None)
    for dy in range(-30, 31, 2):
        for dx in range(-30, 31, 2):
            v = ncc(xa[m:-m, m:-m], np.roll(np.roll(xb, dy, 0), dx, 1)[m:-m, m:-m])
            if v > best[0]:
                best = (v, (dx, dy))
    print("±30 px 平移搜索最好 NCC = %.3f 在 %s" % best)

    # 直方图/亮度关系：同一位姿下两帧的像素是否线性相关
    a, b = xa[m:-m, m:-m].ravel(), xb[m:-m, m:-m].ravel()
    print("（参考：同落位同照明 = 0.999；ACI ↔ WATSON = 0.05–0.17）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
