# -*- coding: utf-8 -*-
"""WATSON 自己的两个判决测试（都用**全画幅**产品，之前用的是 480×480 子帧）。

1) 阳性对照：同一落位两张全画幅（0683466024 与 0683466011，间隔 13 s）
   位姿应当几乎相同 → 模型预测"零位移"，NCC 该在 0.9 上下。
   若这一步就低，说明产品与标签层面就乱了。

2) 位姿测试：同一 sol、瞄同一个靶点、standoff 23.4 mm 与 63.4 mm 两张全画幅
   （视场 13.4 mm 与 36 mm，小的应当整体落在大的里面，采样比 0.369）。
   * 先按模型把近图投到远图网格，看预测重叠有多大、NCC 多少；
   * 再让判据自由搜（位置 × 尺度 × 转角 × 双极性），看有没有比零假设更高的峰；
   * 零假设用相位随机化（保功率谱、毁位置）。

用法：python watson_pose_test.py
"""
from __future__ import annotations

import os
import sys
import urllib.request

import numpy as np
import cv2

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from colorize_pair import OUT, homography, read_bsq, to8
from sign_test import bandpass, zs, search

UA = {"User-Agent": "Mozilla/5.0"}
NEAR = "urn:nasa:pds:mars2020_imgops:data_watson_imgops:sif_0186_0683466024_335fdr_n0070000srlc00621_0000lmj::4.0"
NEAR2 = "urn:nasa:pds:mars2020_imgops:data_watson_imgops:sif_0186_0683466011_449fdr_n0070000srlc00621_0000lmj::5.0"
FAR = "urn:nasa:pds:mars2020_imgops:data_watson_imgops:sif_0186_0683466603_441fdr_n0070000srlc01021_0000lmj::5.0"


def load_cam_xml(lidvid):
    """只用 XML 里的 CAHVOR + 像幅（不依赖 colorize_pair 的字段名）。"""
    import re
    lid, _, ver = lidvid.partition("::")
    stem = lid.rsplit(":", 1)[-1]
    sol = stem.split("_")[1].zfill(5)
    name = stem.upper() + ver.split(".")[0].zfill(2)
    base = ("https://pds-imaging.jpl.nasa.gov/archive/m20/cumulative/"
            "mars2020_imgops/data_watson_imgops/sol/%s/ids/fdr/shrlc/%s")
    xu, iu = base % (sol, name + ".xml"), base % (sol, name + ".IMG")
    dst = os.path.join(OUT, "xml", os.path.basename(xu))
    if not (os.path.exists(dst) and os.path.getsize(dst) > 5000):
        with urllib.request.urlopen(urllib.request.Request(xu, headers=UA),
                                    timeout=60) as fh:
            open(dst, "wb").write(fh.read())
    t = open(dst, encoding="utf-8", errors="replace").read()

    def vec(tag):
        i = t.find("<geom:%s>" % tag)
        seg = t[i:t.index("</geom:%s>" % tag, i)]
        return np.array([float(x) for x in
                         re.findall(r">\s*(-?[\d.eE+-]+)\s*<", seg)][:3])
    cam = dict(C=vec("Vector_Center"), A=vec("Vector_Axis"),
               H=vec("Vector_Horizontal"), V=vec("Vector_Vertical"))
    cam["ns"] = int(re.search(r"<img:samples>\s*(\d+)", t).group(1))
    cam["nl"] = int(re.search(r"<img:lines>\s*(\d+)", t).group(1))
    cam["bands"] = int(re.search(r"<img:bands>\s*(\d+)", t).group(1)) \
        if re.search(r"<img:bands>\s*(\d+)", t) else 3
    rgb = []
    for c in ("red", "green", "blue"):
        m = re.search(r"<disp:%s_channel_band>\s*(\d+)" % c, t)
        rgb.append(int(m.group(1)) if m else 1)
    cam["rgb"] = rgb
    return cam, iu, os.path.join(OUT, "img", os.path.basename(iu))


def gray(lidvid):
    cam, iu, p = load_cam_xml(lidvid)
    if not os.path.exists(p):
        print("    下载 %s …" % os.path.basename(iu), flush=True)
        with urllib.request.urlopen(urllib.request.Request(iu, headers=UA),
                                    timeout=900) as fh:
            open(p, "wb").write(fh.read())
    r = read_bsq(p, cam)
    b = [k - 1 for k in cam["rgb"]]
    if cam["bands"] == 1:
        g = to8(r[0])
    else:
        g = cv2.cvtColor(np.stack([to8(r[b[k]]) for k in range(3)], axis=-1),
                         cv2.COLOR_RGB2GRAY)
    return g[:cam["nl"], :cam["ns"]], cam


def band_at(g, um, lo=0.06, hi=0.5):
    ppm = 1000.0 / um
    hi = min(hi, min(g.shape) * um / 1000.0 / 8.0)
    return np.ascontiguousarray(zs(bandpass(g, lo * ppm / 2.355,
                                            hi * ppm / 2.355)).astype(np.float32))


def maskdark(g, um, k_mm=6.0):
    k = max(3, int(round(k_mm * 1000.0 / um)) | 1)
    loc = cv2.blur(g.astype(np.float32), (k, k))
    return loc < 0.55 * float(np.median(loc))


def main():
    # ---------- 1) 同落位阳性对照 ----------
    ga, ca = gray(NEAR)
    gb, cb = gray(NEAR2)
    print("1) 同落位（间隔 13 s）")
    print("   A %dx%d   B %dx%d   C 相差 %.4f mm   视轴夹角 %.5f°"
          % (ca["ns"], ca["nl"], cb["ns"], cb["nl"],
             np.linalg.norm(ca["C"] - cb["C"]) * 1000,
             np.degrees(np.arccos(np.clip(ca["A"] @ cb["A"], -1, 1)))))
    for lo, hi in ((0.06, 0.5), (0.3, 2.0)):
        xa = band_at(ga, 8.10, lo, hi)
        xb = band_at(gb, 8.10, lo, hi)
        m = maskdark(ga, 8.10)
        xa[m] = 0.0
        xb[m] = 0.0
        a = xa[100:-100, 100:-100]
        b = xb[100:-100, 100:-100]
        av, bv = a - a.mean(), b - b.mean()
        ncc = float((av * bv).sum() / np.sqrt((av ** 2).sum() * (bv ** 2).sum()))
        print("   波段 %.2f–%.1f mm：零位移 NCC = %.3f" % (lo, hi, ncc))

    # ---------- 2) 位姿测试：23.4 mm × 63.4 mm，全画幅 ----------
    gf, cf = gray(FAR)
    print("\n2) 不同落位：23.4 mm (8.10 µm/px) × 63.4 mm (21.96 µm/px)")
    print("   近 %dx%d  远 %dx%d  C 相距 %.2f mm  视轴夹角 %.3f°"
          % (ca["ns"], ca["nl"], cf["ns"], cf["nl"],
             np.linalg.norm(ca["C"] - cf["C"]) * 1000,
             np.degrees(np.arccos(np.clip(ca["A"] @ cf["A"], -1, 1)))))
    Hm, dst = homography(ca, cf, 23.4)      # 单位是 mm（homography 内部 /1000）
    print("   模型预测：近图四角投到远图 → %s"
          % np.array2string(dst, precision=0, suppress_small=True).replace("\n", " "))
    N = band_at(ga, 8.10)
    F = band_at(gf, 21.96)
    warp = cv2.warpPerspective(N, Hm, (cf["ns"], cf["nl"]), flags=cv2.INTER_LINEAR)
    cov = np.abs(warp) > 0
    a, b = F[cov], warp[cov]
    if a.size > 500:
        av, bv = a - a.mean(), b - b.mean()
        print("   按模型直接配准：重叠 %d px，NCC = %+.3f"
              % (a.size, float((av * bv).sum()
                               / np.sqrt((av ** 2).sum() * (bv ** 2).sum()))))
    else:
        print("   按模型直接配准：重叠只有 %d px" % a.size)

    # 自由搜：把近图整幅缩到远图采样（8.10/21.96 = 0.369）后到远图里搜
    ratio = 8.10 / 21.96
    small = np.ascontiguousarray(zs(cv2.resize(N, None, fx=ratio, fy=ratio,
                                               interpolation=cv2.INTER_AREA)
                                    ).astype(np.float32))
    Fm = F.copy()
    Fm[maskdark(gf, 21.96)] = 0.0
    rng = np.random.default_rng(11)
    ph = rng.uniform(0, 2 * np.pi, Fm.shape)
    ph[0, 0] = 0
    ph[0, :] = 0
    ph[:, 0] = 0
    Fr = np.real(np.fft.ifft2(np.abs(np.fft.fft2(Fm)) * np.exp(1j * ph)))
    Fr = np.ascontiguousarray(zs(Fr).astype(np.float32))
    Fr[maskdark(gf, 21.96)] = 0.0
    print("   诊断：Fm finite=%s std=%.3f  Fr finite=%s std=%.3f  small %s std=%.3f"
          % (np.isfinite(Fm).all(), Fm.std(), np.isfinite(Fr).all(), Fr.std(),
             small.shape, small.std()))
    for tag, WW in (("真图", Fm), ("相位随机化", Fr)):
        bmax, bmin, amax = -9.0, 9.0, None
        first = True
        for s in (ratio * 0.8, ratio, ratio * 1.25):
            for ang in range(0, 360, 10):
                mx, lm = search(WW, small, ang, s, "max")
                mn, ln = search(WW, small, ang, s, "min")
                if first:
                    print("   诊断 %s：第一组 s=%.3f ang=%d → mx=%s mn=%s"
                          % (tag, s, ang, mx, mn))
                    first = False
                if np.isfinite(mx) and mx > bmax:
                    bmax, amax = mx, (s, ang, lm)
                if np.isfinite(mn) and mn < bmin:
                    bmin, amin = mn, (s, ang, ln)
        print("   自由搜 %-12s max %+.3f  min %+.3f   (max@%s)"
              % (tag, bmax, bmin, amax))
    print("\n   模型预测的位置：近图中心应落在远图中心附近 —— 远图中心 (%.0f, %.0f)"
          % (cf["ns"] / 2, cf["nl"] / 2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
