# -*- coding: utf-8 -*-
"""三张候选各自上一次色，看**色相**结论是否与"用哪张"无关。

关键点：既然逐像素配准无法验证，那就退一步只问"**颜色结论稳不稳**"：
如果三张不同 standoff/不同 sol 的 WATSON（论文那张 43.8 mm、srlc00621 23.4 mm、
srlc01021 63.4 mm）投到 ACI 网格后，**有光区域的平均/中位色相基本一致**，
那么"这块岩石是橙褐色、某些颗粒偏白"这个结论就不依赖于挑哪张、
也不依赖于那个无法验证的逐像素对应；反之就要把不确定性写清楚。

做法：ab 取白平衡后的 WATSON，阴影区留灰（不臆造），平滑到 0.4 mm，L 用 ACI。
用法：python color_variants.py
输出：out/v3_color_{paper,near,far}.png 与 out/v3_compare.png
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

D_ACI_MM = 108.0
ACI_UM = 6.66
CAND = [
    ("paper", "论文那张 sol185 (43.8 mm)",
     "sif_0185_0683368184_652fdr_n0070000srlc00720_0000lmj"),
    ("near", "srlc00621 sol186 (23.4 mm)",
     "sif_0186_0683466024_335fdr_n0070000srlc00621_0000lmj"),
    ("far", "srlc01021 sol186 (63.4 mm)",
     "sif_0186_0683466603_441fdr_n0070000srlc01021_0000lmj"),
]


def load_cached(stem):
    up = stem.upper()
    xs = [p for p in glob.glob(os.path.join(OUT, "xml", "*.xml"))
          if os.path.basename(p).upper().startswith(up)
          and "_T0" not in os.path.basename(p).upper()]
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
    ip = sorted(ips)[0]
    r = read_bsq(ip, cam)
    b = [k - 1 for k in cam["rgb"]]
    col = np.stack([to8(r[b[k]]) for k in range(3)], axis=-1)          # RGB
    return col[:cam["nl"], :cam["ns"]], cam, os.path.basename(ip)


def white_balance(rgb, mask):
    """用有光区的中位色做一次简单白平衡（灰世界做不到，就按中位色归一到中性）。"""
    med = np.median(rgb[mask], axis=0)
    g = float(med.mean())
    scale = np.clip(g / np.maximum(med, 1.0), 0.5, 2.0)
    return np.clip(rgb * scale[None, None, :], 0, 255).astype(np.uint8), scale


def main():
    aci = load_cam(ACI_XML)
    aci8 = to8(read_bsq(ACI_IMG, aci)[0])
    panels = []
    stats = []
    for key, tag, stem in CAND:
        col, cam, fn = load_cached(stem)
        Hm, _ = homography(aci, cam, D_ACI_MM)
        wcol = cv2.warpPerspective(col.astype(np.float32), Hm,
                                   (aci["ns"], aci["nl"]),
                                   flags=cv2.INTER_LINEAR | cv2.WARP_INVERSE_MAP)
        valid = cv2.warpPerspective(np.full((cam["nl"], cam["ns"]), 255, np.uint8),
                                    Hm, (aci["ns"], aci["nl"]),
                                    flags=cv2.INTER_NEAREST | cv2.WARP_INVERSE_MAP) > 127
        gw = cv2.cvtColor(wcol.astype(np.uint8), cv2.COLOR_RGB2GRAY).astype(np.float32)
        loc = cv2.blur(gw, (151, 151))
        dark = loc < 0.55 * float(np.median(loc))
        lit = valid & ~dark
        col8, scale = white_balance(wcol.astype(np.uint8), lit)
        lab = cv2.cvtColor(col8, cv2.COLOR_RGB2LAB).astype(np.float32)
        ab = lab[..., 1:3].copy()
        s = 0.4 * (1000.0 / ACI_UM) / 2.355
        ab = cv2.GaussianBlur(ab, (0, 0), s)
        out = cv2.cvtColor(aci8, cv2.COLOR_GRAY2RGB)
        lo = cv2.cvtColor(out, cv2.COLOR_RGB2LAB).astype(np.float32)
        lo[..., 1] = np.where(lit, ab[..., 0], 128.0)
        lo[..., 2] = np.where(lit, ab[..., 1], 128.0)
        rgb = cv2.cvtColor(np.clip(lo, 0, 255).astype(np.uint8), cv2.COLOR_LAB2RGB)
        rgb[~valid] = 0
        cv2.imwrite(os.path.join(OUT, "v3_color_%s.png" % key),
                    cv2.cvtColor(rgb, cv2.COLOR_RGB2BGR))
        a_med = float(np.median(ab[..., 0][lit]))
        b_med = float(np.median(ab[..., 1][lit]))
        stats.append((tag, fn, a_med, b_med, 100.0 * lit[valid].mean(),
                      float(np.std(ab[..., 0][lit])), scale))
        bar = np.full((30, aci["ns"], 3), 255, np.uint8)
        cv2.putText(bar, "%s  ab(median)=%.0f,%.0f" % (tag, a_med, b_med),
                    (8, 21), cv2.FONT_HERSHEY_SIMPLEX, 0.55, (0, 0, 0), 1, cv2.LINE_AA)
        panels.append(np.vstack([bar, cv2.cvtColor(rgb, cv2.COLOR_RGB2BGR)]))

    print("%-30s %-28s  ab 中位      有光占比  a 的标准差  白平衡系数" % ("候选", "文件"))
    for tag, fn, a_med, b_med, litp, sa, scale in stats:
        print("%-30s %-28s %5.0f,%-5.0f %8.0f%% %10.1f   %s"
              % (tag, fn[:28], a_med, b_med, litp, sa, np.round(scale, 2)))
    aa = [s[2] for s in stats]
    bb = [s[3] for s in stats]
    print("\n三张的 a 中位差 %.0f、b 中位差 %.0f（LAB 里 128 = 中性；"
          "差几以内说明色相结论与挑哪张无关）" % (max(aa) - min(aa), max(bb) - min(bb)))
    h = 700
    def rs(im):
        f = h / im.shape[0]
        return cv2.resize(im, (int(im.shape[1] * f), h))
    cv2.imwrite(os.path.join(OUT, "v3_compare.png"),
                np.vstack([rs(p) for p in panels]))
    print("写到 out/v3_color_{paper,near,far}.png 与 out/v3_compare.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
