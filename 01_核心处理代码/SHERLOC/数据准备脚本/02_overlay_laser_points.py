#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SHERLOC 练习脚本 02：把激光点叠到 ACI 图上
==============================================
对应操作手册「第 7 步 图像数据处理」里最核心的那一步：
用 RLS（激光点坐标表）把 100 个光谱点的位置画到 ACI 灰度图上。

如果同时提供候选峰表（脚本 01 的输出），还会把有高置信度矿物峰的点
高亮出来 —— 这就是 SHERLOC 论文里的那张核心成果图。

用法：
    python 02_overlay_laser_points.py <ACI图像> <RLS文件> [候选峰表] [输出文件]

依赖：
    pip install pillow numpy pandas
"""

import os
import sys

import pandas as pd
from PIL import Image, ImageDraw, ImageFont

# 画点用的样式
SPOT_RADIUS = 9          # ACI 是 10.1 µm/像素，激光斑约 100-110 µm ≈ 10 像素直径
HIT_RADIUS = 22          # 有矿物峰的点，画大一圈

COLOR_BASE = (0, 200, 255, 220)      # 普通点：青色
COLOR_HIT = (255, 60, 60, 240)       # 有矿物峰的点：红色
COLOR_RING = (255, 255, 0, 230)      # 高亮点外面的黄圈
COLOR_TEXT = (255, 255, 0, 255)


def load_rls(path, image_keyword=None):
    """读 RLS，返回 (点号, samp, line)。同名文件里常有两个图像块，按关键字挑一块。"""
    df = pd.read_csv(path, skiprows=1)
    if image_keyword:
        df = df[df["image_name"].str.contains(image_keyword, case=False, na=False)]
    # 同一个点号只保留一行
    df = df.drop_duplicates(subset=["number"]).sort_values("number")
    return df


def load_peak_points(path):
    """读候选峰表，返回有高置信度矿物峰的点号集合。"""
    if not path or not os.path.exists(path):
        return set()
    df = pd.read_csv(path, encoding="utf-8-sig")
    df = df[df["区间"] == "矿物峰"]
    return set(int(v) for v in df["点号"].unique())


def get_font(size=18):
    for name in ("msyh.ttc", "simhei.ttf", "arial.ttf"):
        try:
            return ImageFont.truetype(name, size)
        except Exception:
            continue
    return ImageFont.load_default()


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    img_path, rls_path = sys.argv[1], sys.argv[2]
    peak_path = sys.argv[3] if len(sys.argv) > 3 and sys.argv[3] else None
    out_path = sys.argv[4] if len(sys.argv) > 4 else os.path.splitext(img_path)[0] + "_overlay.png"

    # 用图像文件名里的 SCLK 段（形如 SC3_0181_0683033276_164ECM 的第 3 段）
    # 去 RLS 的 image_name 列里匹配，挑出这张图对应的一块激光点坐标
    img_stem = os.path.splitext(os.path.basename(img_path))[0]
    key = img_stem.split("_")[2] if len(img_stem.split("_")) > 2 else None
    pts = load_rls(rls_path, image_keyword=key)
    if pts.empty:
        pts = load_rls(rls_path)

    hits = load_peak_points(peak_path)

    img = Image.open(img_path).convert("RGB")
    W, H = img.size
    overlay = Image.new("RGBA", img.size, (0, 0, 0, 0))
    draw = ImageDraw.Draw(overlay)
    font = get_font(18)

    print("图像尺寸 : %d x %d" % (W, H))
    print("点数     : %d" % len(pts))
    print("有矿物峰的点: %s" % (sorted(hits) if hits else "（未提供峰表）"))

    out_of_range = 0
    for _, r in pts.iterrows():
        x, y = float(r["samp"]), float(r["line"])
        n = int(r["number"])
        if not (0 <= x < W and 0 <= y < H):
            out_of_range += 1
            continue

        if n in hits:
            draw.ellipse([x - HIT_RADIUS, y - HIT_RADIUS, x + HIT_RADIUS, y + HIT_RADIUS],
                         outline=COLOR_RING, width=4)
            draw.ellipse([x - 5, y - 5, x + 5, y + 5], fill=COLOR_HIT)
            draw.text((x + HIT_RADIUS + 3, y - 12), str(n), fill=COLOR_TEXT, font=font)
        else:
            draw.ellipse([x - SPOT_RADIUS / 2, y - SPOT_RADIUS / 2,
                          x + SPOT_RADIUS / 2, y + SPOT_RADIUS / 2],
                         outline=COLOR_BASE, width=2)

    print("越界点数 : %d" % out_of_range)

    Image.alpha_composite(img.convert("RGBA"), overlay).convert("RGB").save(out_path)
    print("已保存   :", out_path)


if __name__ == "__main__":
    main()
