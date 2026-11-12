# -*- coding: utf-8 -*-
"""WATSON 自己验自己：两张不同 standoff、瞄同一个靶点的 WATSON 该能对上。

逻辑：ACI 侧已经被两个独立事实钉住了——
  * ACI 的相机模型与标签自洽（FOV 差 0.4%），像素原点约定查过；
  * 打激光前 / 后两张 ACI（落位完全相同，C 差 0.00 mm）互相比 NCC = 0.672，
    说明"同一台相机、同一个地点"的比对管线是work的。
可是这两张 ACI 谁都配不上任何一张 WATSON（五张候选全 ≤0.17）。
那就只剩一种可能：**WATSON 标签里的位姿和它的图对不上**。

怎么验：sol 186 上有两张瞄同一个靶点的 WATSON
    srlc00621   standoff 23.4 mm  8.1 µm/px  视场 3.9 mm
    srlc01021   standoff 63.4 mm  21.2 µm/px 视场 10.5 mm
近的那张的 3.9 mm 视场，应当整体落在远的那张的 10.5 mm 视场里，
且采样相差 21.2/8.1 = 2.6 倍。如果模型对，把近图缩到远图的采样后在远图里搜，
应当在**预测位置（两张都瞄 T，所以是远图中心附近）**出现高相关。

用法：python watson_watson_test.py
"""
from __future__ import annotations

import os
import sys

import numpy as np
import cv2

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from colorize_pair import OUT, read_bsq, to8
from screen_watson import urls_from_lidvid, read_cam_from, collect
from sign_test import bandpass, zs, search

NEAR = "sif_0186_0683466001_667fdr_n0070000srlc00621_0000lmj"
FAR = "sif_0186_0683466593_199fdr_n0070000srlc01021_0000lmj"
CACHE = {}


def load_stem(stem):
    if stem in CACHE:
        return CACHE[stem]
    lv = None
    for s, variants in collect([stem[4:8]]):
        if s == stem:
            lv = variants
    cam = read_cam_from(lv)
    url = urls_from_lidvid(cam["lidvid"])[1]
    img = os.path.join(OUT, "img", os.path.basename(url))
    r = read_bsq(img, cam)
    if cam["bands"] == 1 or not cam["rgb"] or cam["rgb"][0] is None:
        g = to8(r[0])
    else:
        b = [i - 1 for i in cam["rgb"]]
        g = cv2.cvtColor(np.stack([to8(r[b[i]]) for i in range(3)], axis=-1),
                         cv2.COLOR_RGB2GRAY)
    g = g[:cam["nl"], :cam["ns"]]
    CACHE[stem] = (g, cam)
    return CACHE[stem]


def um_per_px(stem):
    for line in open(os.path.join(OUT, "screen_watson.txt"), encoding="utf-8"):
        p = line.split()
        if len(p) >= 7 and p[6] == stem:
            return float(p[2]), float(p[1]), float(p[3])   # um/px, standoff, offset
    raise KeyError(stem)


def main():
    gn, camn = load_stem(NEAR)
    gf, camf = load_stem(FAR)
    umn, dn, offn = um_per_px(NEAR)
    umf, df, offf = um_per_px(FAR)
    print("近：%s  %dx%d  %.2f µm/px  视场 %.1f mm  standoff %.1f mm  靶点偏移 %.1f mm"
          % (NEAR[-24:], camn["ns"], camn["nl"], umn, camn["ns"] * umn / 1000, dn, offn))
    print("远：%s  %dx%d  %.2f µm/px  视场 %.1f mm  standoff %.1f mm  靶点偏移 %.1f mm"
          % (FAR[-24:], camf["ns"], camf["nl"], umf, camf["ns"] * umf / 1000, df, offf))
    print("两张的光心相距 %.1f mm，视轴夹角 %.2f°"
          % (np.linalg.norm(camn["C"] - camf["C"]) * 1000,
             np.degrees(np.arccos(np.clip(camn["A"] @ camf["A"], -1, 1)))))

    # 两侧带通到同一物理波段（各自按自己的 µm/px 折算），再用同一个物理尺度比
    def band(g, um):
        ppm = 1000.0 / um
        hi = min(0.5, min(g.shape) * um / 1000.0 / 8.0)
        return np.ascontiguousarray(zs(bandpass(g, 0.06 * ppm / 2.355,
                                                hi * ppm / 2.355)).astype(np.float32))

    N = band(gn, umn)
    F = band(gf, umf)
    # 去掉各自的暗区
    for G, g in ((N, gn), (F, gf)):
        k = max(3, int(round(4.0 * 1000.0 / umn)) | 1)
        loc = cv2.blur(g.astype(np.float32), (k, k))
        G[loc < 0.55 * float(np.median(loc))] = 0.0
    ratio = umn / umf
    print("\n近图缩到远图采样：倍率 %.3f（%d px → %d px）"
          % (ratio, camn["ns"], int(camn["ns"] * ratio)))
    small = cv2.resize(N, None, fx=ratio, fy=ratio, interpolation=cv2.INTER_AREA)
    small = np.ascontiguousarray(zs(small).astype(np.float32))
    bmax, bmin, amax, amin = -9.0, 9.0, None, None
    for s in (ratio * 0.8, ratio, ratio * 1.25):
        for ang in range(0, 360, 10):
            mx, lm = search(F, small, ang, s, "max")
            mn, ln = search(F, small, ang, s, "min")
            if np.isfinite(mx) and mx > bmax:
                bmax, amax = mx, (s, ang, lm)
            if np.isfinite(mn) and mn < bmin:
                bmin, amin = mn, (s, ang, ln)
    print("全搜结果：max %+.3f%s   min %+.3f%s" %
          (bmax, (" @尺度%.3f/%d° 位置%s" % amax) if amax else "",
           bmin, (" @尺度%.3f/%d° 位置%s" % amin) if amin else ""))
    if amax:
        cx = camf["ns"] / 2 - small.shape[1] / 2
        cy = camf["nl"] / 2 - small.shape[0] / 2
        print("模型预测的位置（两张都瞄 T，远处的靶点偏移 %.1f mm → 中心附近）：(%.0f, %.0f)"
              % (offf, cx, cy))

    # ---- 直接把"模型预测的那一组"算出来：按单应把近图投到远图网格 ----
    from colorize_pair import homography
    Hm, dst = homography(camn, camf, dn / 1000.0)      # 近相机 → 远相机，靶面在近处 dn
    warp = cv2.warpPerspective(N, Hm, (camf["ns"], camf["nl"]),
                               flags=cv2.INTER_LINEAR)
    cov = np.abs(warp) > 0
    a = F[cov]
    b = warp[cov]
    if a.size > 500:
        av, bv = a - a.mean(), b - b.mean()
        nac = float((av * bv).sum() / np.sqrt((av ** 2).sum() * (bv ** 2).sum()))
        print("按模型直接配准（不搜）：有效重叠 %d px，NCC = %+.3f" % (a.size, nac))
    else:
        print("按模型直接配准：重叠只有 %d px（说明模型预测它俩几乎不重叠）" % a.size)

    # ---- 零假设对照：把远图按 16×16 块打乱，破坏一切对应但保留纹理统计 ----
    rng = np.random.default_rng(0)
    blk = 16
    Fs = F.copy()
    hh, ww = F.shape
    tiles = [(i, j) for i in range(0, hh - blk + 1, blk)
             for j in range(0, ww - blk + 1, blk)]
    idx = rng.permutation(len(tiles))
    tmp = F.copy()
    for (i, j), t in zip(tiles, idx):
        Fs[i:i + blk, j:j + blk] = tmp[tiles[t][0]:tiles[t][0] + blk,
                                        tiles[t][1]:tiles[t][1] + blk]
    nmax, nmin = -9.0, 9.0
    for s in (ratio * 0.8, ratio, ratio * 1.25):
        for ang in range(0, 360, 10):
            mx, _ = search(Fs, small, ang, s, "max")
            mn, _ = search(Fs, small, ang, s, "min")
            if np.isfinite(mx):
                nmax = max(nmax, mx)
            if np.isfinite(mn):
                nmin = min(nmin, mn)
    print("零假设对照（远图按 16×16 块打乱后同样的全搜）：max %+.3f  min %+.3f"
          % (nmax, nmin))
    print("\n判读：阳性对照（把一块内容贴进去再搜）是 0.95；"
          "模型预测的那一组若也只有 0.1 上下，说明**WATSON 的位姿与它的图对不上**。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
