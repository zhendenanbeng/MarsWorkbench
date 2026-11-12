#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SHERLOC 练习脚本 03：把三个 Region 画出来，看清「拉曼尖峰」和「荧光宽包」
====================================================================
用法：
    python 03_plot_regions.py <RRS文件> [输出png]

图中上排是 Region 1（拉曼区，250-277 nm），下排是 Region 2/3（荧光区，275-357 nm）。
两排的纵轴刻度不同 —— 这正是关键：荧光比拉曼强好几个数量级。

依赖：pip install pillow numpy pandas
"""

import os
import sys

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont

LASER_NM = 248.6
W, H = 1400, 920


def load_rrs(path):
    with open(path) as f:
        contents = f.readlines()
    n = int((len(contents) - 9) / 3)
    wl = np.array([float(v) for v in contents[2].split(",")])
    r1 = pd.read_csv(path, nrows=n, skiprows=4).to_numpy()
    r2 = pd.read_csv(path, nrows=n, skiprows=4 + n + 2).to_numpy()
    r3 = pd.read_csv(path, nrows=n, skiprows=4 + 2 * n + 4).to_numpy()
    return wl, r1, r2, r3


def get_font(size):
    for name in (r"C:\Windows\Fonts\msyh.ttc", r"C:\Windows\Fonts\simhei.ttf",
                 r"C:\Windows\Fonts\arial.ttf"):
        try:
            return ImageFont.truetype(name, size)
        except Exception:
            continue
    return ImageFont.load_default()


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    rrs = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.splitext(rrs)[0] + "_regions.png"

    wl, r1, r2, r3 = load_rrs(rrs)
    m1 = np.median(r1, axis=0)
    m2 = np.median(r2, axis=0)
    m3 = np.median(r3, axis=0)

    # 按波长排序（RRS 里 region 表与波长表未必严格同序）
    o = np.argsort(wl)
    wl, m1, m2, m3 = wl[o], m1[o], m2[o], m3[o]

    img = Image.new("RGB", (W, H), "white")
    d = ImageDraw.Draw(img)
    f_title = get_font(30)
    f_label = get_font(20)
    f_small = get_font(17)

    d.text((40, 20), "SHERLOC 三个光谱区：拉曼（尖峰）vs 荧光（宽包）", fill="black", font=f_title)
    d.text((40, 58), os.path.basename(rrs), fill=(110, 110, 110), font=f_small)

    X0, X1 = 130, 1350

    def xmap(nm):
        return X0 + (nm - 248.0) / (358.0 - 248.0) * (X1 - X0)

    def panel(y0, y1, title, ymax):
        d.rectangle([X0, y0, X1, y1], outline=(60, 60, 60), width=2)
        d.text((X0, y0 - 30), title, fill="black", font=f_label)
        for frac in (0.25, 0.5, 0.75, 1.0):
            yy = y1 - frac * (y1 - y0)
            d.line([X0, yy, X1, yy], fill=(225, 225, 225), width=1)
            d.text((X0 - 78, yy - 11), "%.3g" % (frac * ymax), fill=(90, 90, 90), font=f_small)

    def curve(nm_arr, y_arr, y0, y1, ymax, color, width=3):
        pts = []
        for nm, val in zip(nm_arr, y_arr):
            px = xmap(nm)
            if px < X0 or px > X1:
                continue
            v = min(max(val, 0.0), ymax)
            py = y1 - (v / ymax) * (y1 - y0)
            pts.append((px, py))
        if len(pts) > 1:
            d.line(pts, fill=color, width=width)

    sel1 = (wl >= 250.0) & (wl <= 276.0)
    A0, A1 = 110, 400
    ymaxA = float(np.percentile(m1[sel1], 99.0)) * 1.15
    panel(A0, A1, "上排：Region 1（拉曼区）250-277 nm  —— 又尖又窄", ymaxA)
    curve(wl[sel1], m1[sel1], A0, A1, ymaxA, (200, 30, 30))

    selB = (wl >= 275.0) & (wl <= 357.0)
    B0, B1 = 500, 790
    ymaxB = float(max(m2[selB].max(), m3[selB].max())) * 1.15
    panel(B0, B1, "下排：Region 2 / 3（荧光区）275-357 nm  —— 又宽又缓", ymaxB)
    curve(wl[selB], m2[selB], B0, B1, ymaxB, (30, 90, 210))
    curve(wl[selB], m3[selB], B0, B1, ymaxB, (20, 150, 90))

    # 252.93 nm 激光线
    xl = xmap(252.93)
    if X0 < xl < X1:
        d.line([xl, A0, xl, A1], fill=(255, 140, 0), width=2)
        d.text((xl + 6, A0 + 8), "252.93 nm\n激光线", fill=(200, 100, 0), font=f_small)

    # 荧光峰位
    k2 = int(np.argmax(np.where(selB, m2, -1)))
    xf = xmap(wl[k2])
    d.ellipse([xf - 6, B0 + (B1 - B0) * (1 - m2[k2] / ymaxB) - 6,
               xf + 6, B0 + (B1 - B0) * (1 - m2[k2] / ymaxB) + 6], fill=(30, 90, 210))
    d.text((xf - 60, B0 + 6), "荧光峰 %.1f nm" % wl[k2], fill=(20, 70, 180), font=f_small)

    # 坐标轴刻度
    for nm in range(250, 359, 10):
        xx = xmap(nm)
        d.line([xx, A1, xx, A1 + 6], fill=(60, 60, 60), width=2)
        d.line([xx, B1, xx, B1 + 6], fill=(60, 60, 60), width=2)
        d.text((xx - 16, B1 + 10), str(nm), fill=(60, 60, 60), font=f_small)
    d.text(((X0 + X1) / 2 - 40, B1 + 42), "波长 (nm)", fill="black", font=f_label)

    d.text((40, 850), "蓝=Region2  绿=Region3  红=Region1      两排纵轴刻度不同：拉曼比荧光弱约 3 个数量级",
           fill=(80, 80, 80), font=f_small)

    img.save(out)

    # 控制台里也给出可比较的数字
    print("荧光峰位: Region2 = %.1f nm, Region3 = %.1f nm" % (wl[np.argmax(np.where(selB, m2, -1))],
                                                              wl[np.argmax(np.where(selB, m3, -1))]))
    print("拉曼区（Region1, 253.5-276 nm）中位强度范围: %.1f - %.1f"
          % (np.median(m1[(wl >= 253.5) & (wl <= 276.0)]),
             np.percentile(m1[(wl >= 253.5) & (wl <= 276.0)], 99)))
    print("荧光区（Region2）中位强度最大值: %.0f" % m2[selB].max())
    print("已保存:", out)


if __name__ == "__main__":
    main()
