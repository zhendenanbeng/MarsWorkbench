# -*- coding: utf-8 -*-
"""给"上色"挑来源：哪张 WATSON 的视场里，ACI 那块地方是**有光的**？

既然：
  * 两张图的明暗在物理上无关（ACI 环光 vs WATSON 阳光），配准无法用图像验证；
  * 颜色只在"有光"的像素上有意义——阴影里的 a/b 是照明条件，不是岩石颜色；
  * 而 WATSON 的阴影位置随 standoff 变（太阳不动、相机在动），
    所以不同 standoff 的候选里，有的可能整块都在阳光下。

对每张候选算三件事：预测采样尺度、有效覆盖比例、覆盖区里落在阴影内的比例。
用法：python color_source_rank.py
"""
from __future__ import annotations

import glob
import os
import re
import sys

import numpy as np
import cv2

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from colorize_pair import ACI_IMG, ACI_XML, OUT, homography, load_cam, read_bsq, to8

ACI_UM = 6.66
D_ACI_MM = 108.0

CAND = [
    ("WATSON 论文那张 sol185 (43.8 mm)",
     "sif_0185_0683368184_652fdr_n0070000srlc00720_0000lmj"),
    ("srlc00621 sol186 (23.4 mm)",
     "sif_0186_0683466024_335fdr_n0070000srlc00621_0000lmj"),
    ("srlc01021 sol186 (63.4 mm)",
     "sif_0186_0683466603_441fdr_n0070000srlc01021_0000lmj"),
]


def um_of(stem):
    for line in open(os.path.join(OUT, "screen_watson.txt"), encoding="utf-8"):
        p = line.split()
        if len(p) >= 7 and p[6] == stem:
            return float(p[2])
    return None


def load_cached(stem):
    """从 out/xml 与 out/img 按文件名前缀取（离线）。"""
    up = stem.upper()
    xs = [p for p in glob.glob(os.path.join(OUT, "xml", "*.xml"))
          if os.path.basename(p).upper().startswith(up)
          and "_T0" not in os.path.basename(p).upper()]
    if not xs:
        raise FileNotFoundError("缓存里没有标签")
    t = open(sorted(xs)[0], encoding="utf-8", errors="replace").read()

    def vec(tag):
        i = t.find("<geom:%s>" % tag)
        seg = t[i:t.index("</geom:%s>" % tag, i)]
        return np.array([float(v) for v in
                         re.findall(r">\s*(-?[\d.eE+-]+)\s*<", seg)][:3])
    cam = dict(C=vec("Vector_Center"), A=vec("Vector_Axis"),
               H=vec("Vector_Horizontal"), V=vec("Vector_Vertical"))
    cam["ns"] = int(re.search(r"<img:samples>\s*(\d+)", t).group(1))
    cam["nl"] = int(re.search(r"<img:lines>\s*(\d+)", t).group(1))
    rgb = []
    for ch in ("red", "green", "blue"):
        m = re.search(r"<disp:%s_channel_band>\s*(\d+)" % ch, t)
        rgb.append(int(m.group(1)) if m else 1)
    cam["rgb"] = rgb
    mb = re.search(r"<img:bands>\s*(\d+)", t)
    cam["bands"] = int(mb.group(1)) if mb else 3
    ips = [p for p in glob.glob(os.path.join(OUT, "img", "*.IMG"))
           if os.path.basename(p).upper().startswith(up)
           and "_T0" not in os.path.basename(p).upper()]
    if not ips:
        raise FileNotFoundError("缓存里没有图")
    ip = sorted(ips)[0]
    r = read_bsq(ip, cam)
    b = [k - 1 for k in cam["rgb"]]
    g = cv2.cvtColor(np.stack([to8(r[b[k]]) for k in range(3)], axis=-1),
                     cv2.COLOR_RGB2GRAY)
    return g[:cam["nl"], :cam["ns"]], cam, os.path.basename(ip)


def main():
    aci = load_cam(ACI_XML)
    rows = []
    for tag, stem in CAND:
        um = um_of(stem)
        try:
            g, cam, fn = load_cached(stem)
        except Exception as exc:                                   # noqa: BLE001
            print("%-34s 取不到（%s）" % (tag, exc))
            continue
        Hm, dst = homography(aci, cam, D_ACI_MM)
        valid = cv2.warpPerspective(np.full((cam["nl"], cam["ns"]), 255, np.uint8),
                                    Hm, (aci["ns"], aci["nl"]),
                                    flags=cv2.INTER_NEAREST | cv2.WARP_INVERSE_MAP) > 127
        loc = cv2.blur(g.astype(np.float32), (151, 151))
        shadow = loc < 0.55 * float(np.median(loc))
        sh_w = cv2.warpPerspective((shadow * 255).astype(np.uint8), Hm,
                                   (aci["ns"], aci["nl"]),
                                   flags=cv2.INTER_NEAREST | cv2.WARP_INVERSE_MAP) > 127
        nv = int(valid.sum())
        dark = int((sh_w & valid).sum())
        usable = 100.0 * (nv - dark) / max(nv, 1)
        scale = ACI_UM / um if um else float("nan")
        rows.append((usable, um, tag))
        print("%-34s %.2f µm/px（ACI 1 px → 它 %.2f px） 覆盖 %.0f%%  阴影 %.0f%% "
              "→ **可用 %.0f%%**   [%s]"
              % (tag, um, scale, 100.0 * nv / valid.size,
                 100.0 * dark / max(nv, 1), usable, fn[:26]))
    print("\n按『可用比例优先、同样可用时取采样更细』排序：")
    for usable, um, tag in sorted(rows, key=lambda r: (-round(r[0] / 10), r[1])):
        print("  可用 %3.0f%%   %5.2f µm/px   %s" % (usable, um, tag))
    print("\n注：可用比例低不是『不能上色』，而是『那部分没有颜色信息』——"
          "诚实做法是留灰，不要用阴影里的 a/b 顶替（试过：会糊出一块假橙色）。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
