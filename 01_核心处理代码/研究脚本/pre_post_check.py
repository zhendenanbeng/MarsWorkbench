# -*- coding: utf-8 -*-
"""这张 ACI 是**打激光之前**还是之后？取同一活动的另一张 ACI 来对。

sol 186 / srlc11420（= Bellegarde survey）的 ACI 产品时间线：
    0683478309 – 0683478397   44 张，间隔 1.75 s   ← 像自动对焦扫描
    0683478571.425            1 张
    0683478583.214            1 张   ← RLS 产品时刻 0683478599.650 之前 16 s
    0683479661.191            1 张
    0683479674.054            1 张   ← **现在上色用的就是这张**（RLS 之后 18 分钟）

如果"打完激光的表面"和"打之前的表面"是两块不同的东西，那么：
  * 两张 ACI 互相比应当**很像**（同一台相机、差不多同一个落位，只多了点坑）
  * 但只有**打之前那一张**才应该和 WATSON 对得上

用法：python pre_post_check.py
"""
from __future__ import annotations

import os
import sys
import urllib.request

import numpy as np
import cv2

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from colorize_pair import (ACI_IMG, ACI_XML, OUT, homography, load_cam,
                           read_bsq, to8)
from sign_test import bandpass, zs, search
from screen_watson import urls_from_lidvid, read_cam_from, collect

PRE = ("urn:nasa:pds:mars2020_imgops:data_aci_imgops:"
       "sc3_0186_0683478583_214fdr_n0070000srlc11420_0000lmj::5.0")
UA = {"User-Agent": "Mozilla/5.0"}
ACI_UM = 6.66
CAND = ["sif_0186_0683466001_667fdr_n0070000srlc00621_0000lmj",   # 8.1 µm/px
        "sif_0185_0683368174_617fdr_n0070000srlc00720_0000lmj",   # 14.7 µm/px
        "sif_0186_0683466593_199fdr_n0070000srlc01021_0000lmj"]   # 21.2 µm/px


def aci_urls(lidvid):
    lid, _, ver = lidvid.partition("::")
    stem = lid.rsplit(":", 1)[-1]
    sol = stem.split("_")[1].zfill(5)
    name = stem.upper() + ver.split(".")[0].zfill(2)
    base = ("https://pds-imaging.jpl.nasa.gov/archive/m20/cumulative/"
            "mars2020_imgops/data_aci_imgops/sol/%s/ids/fdr/shrlc/%s")
    return base % (sol, name + ".xml"), base % (sol, name + ".IMG")


def fetch(url, dst):
    if os.path.exists(dst) and os.path.getsize(dst) > 5000:
        return dst
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA),
                                timeout=600) as fh:
        open(dst, "wb").write(fh.read())
    return dst


def polar_vs(wat_stem, aci_img, aci_um, label):
    """把 aci_img 的结构图拿到某张 WATSON 上全搜。返回最好 |NCC|。"""
    lv = None
    for s, variants in collect([wat_stem[4:8]]):
        if s == wat_stem:
            lv = variants
            break
    cam = read_cam_from(lv)
    img = os.path.join(OUT, "img", os.path.basename(urls_from_lidvid(cam["lidvid"])[1]))
    fetch(urls_from_lidvid(cam["lidvid"])[1], img)
    r = read_bsq(img, cam)
    if cam["bands"] == 1 or not cam["rgb"] or cam["rgb"][0] is None:
        g = to8(r[0])
    else:
        b = [i - 1 for i in cam["rgb"]]
        g = cv2.cvtColor(np.stack([to8(r[b[i]]) for i in range(3)], axis=-1),
                         cv2.COLOR_RGB2GRAY)
    g = g[:cam["nl"], :cam["ns"]]
    um = aci_um / 0.353 * 0.353          # 占位，下面用 screen 表里的值覆盖
    for line in open(os.path.join(OUT, "screen_watson.txt"), encoding="utf-8"):
        p = line.split()
        if len(p) >= 7 and p[6] == wat_stem:
            um = float(p[2])
            break
    pred = aci_um / um
    ppm = 1000.0 / um
    hi_mm = min(0.5, min(cam["ns"], cam["nl"]) * um / 1000.0 / 8.0)
    W = np.ascontiguousarray(zs(bandpass(g, 0.06 * ppm / 2.355,
                                         hi_mm * ppm / 2.355)).astype(np.float32))
    k = max(3, int(round(min(cam["ns"] * um / 1000.0, 8.0) * ppm)) | 1)
    loc = cv2.blur(g.astype(np.float32), (k, k))
    W[loc < 0.55 * float(np.median(loc))] = 0.0
    max_t = int(min(cam["ns"], cam["nl"]) * 0.95 / pred)
    cy, cx = aci_img.shape[0] // 2, aci_img.shape[1] // 2
    h, w = min(max_t, aci_img.shape[0]), min(max_t, aci_img.shape[1])
    sub = aci_img[cy - h // 2:cy + h // 2, cx - w // 2:cx + w // 2]
    A = np.ascontiguousarray(zs(bandpass(
        sub, 0.06 * (1000.0 / aci_um) / 2.355,
        0.5 * (1000.0 / aci_um) / 2.355)).astype(np.float32))
    bmax, bmin = -9.0, 9.0
    for s in (pred * 0.75, pred * 0.9, pred, pred * 1.24, pred * 1.5):
        if s <= 0.05 or s > 10 or max(A.shape) * s > min(cam["ns"], cam["nl"]):
            continue
        for ang in range(0, 360, 10):
            mx, _ = search(W, A, ang, s, "max")
            mn, _ = search(W, A, ang, s, "min")
            if np.isfinite(mx):
                bmax = max(bmax, mx)
            if np.isfinite(mn):
                bmin = min(bmin, mn)
    print("    %-12s vs %-52s max %+.3f  min %+.3f"
          % (label, wat_stem[:52], bmax, bmin), flush=True)
    return max(bmax, -bmin)


def main():
    post = load_cam(ACI_XML)
    xu, iu = aci_urls(PRE)
    xp = fetch(xu, os.path.join(OUT, "xml", os.path.basename(xu)))
    ip = fetch(iu, os.path.join(OUT, "img", os.path.basename(iu)))
    pre = load_cam(xp)
    print("打之前那张：%s" % os.path.basename(xp))
    print("  C=(%.4f,%.4f,%.4f) A=(%.4f,%.4f,%.4f) |H|=%.1f"
          % (*pre["C"], *pre["A"], np.linalg.norm(pre["H"])))
    print("  （现在用的那张 C=(%.4f,%.4f,%.4f)）" % tuple(post["C"]))
    print("  光心相距 %.2f mm；视轴夹角 %.3f°"
          % (np.linalg.norm(pre["C"] - post["C"]) * 1000,
             np.degrees(np.arccos(np.clip(pre["A"] @ post["A"], -1, 1)))))

    gpre = to8(read_bsq(ip, pre)[0])
    gpost = to8(read_bsq(ACI_IMG, post)[0])

    # 两张 ACI 互比：如果落位不同，用模型把"打之前"投到"现在"的网格上
    H, _ = homography(pre, post, 0.108)
    wp = cv2.warpPerspective(gpre.astype(np.float32), H,
                             (post["ns"], post["nl"]), flags=cv2.INTER_LINEAR)
    m = wp > 0
    a = zs(bandpass(gpost, 4, 32))
    b = zs(bandpass(wp, 4, 32))
    n = len(a[m])
    r = float(((a[m] - a[m].mean()) * (b[m] - b[m].mean())).sum()
              / np.sqrt(((a[m] - a[m].mean()) ** 2).sum()
                        * ((b[m] - b[m].mean()) ** 2).sum()))
    print("\n两张 ACI 互比（按模型对齐后，同一物理尺度）：NCC %.3f" % r)
    diff = np.clip(128 + 0.6 * (gpre.astype(np.float32) - wp), 0, 255).astype(np.uint8)
    diff[~m] = 0
    cv2.imwrite(os.path.join(OUT, "vc_aci_prepost_diff.png"), diff)
    print("写到 out/vc_aci_prepost_diff.png（前 − 后）")

    print("\n拿两张 ACI 分别去 WATSON 里搜同一个判据（阳性对照 0.95）：")
    for stem in CAND:
        polar_vs(stem, gpre, ACI_UM, "打之前")
        polar_vs(stem, gpost, ACI_UM, "现在(后)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
