# -*- coding: utf-8 -*-
"""把"尺度"这个自由度彻底放开，只在**物理波段**上对齐，再全搜一次。

为什么必须重做：前面几轮搜索的尺度范围都是拿"ACI standoff=108 mm"这个锚算出来的
（候选 ×0.75…×1.5）。可是 108 与几何解 129 之间差 20%，对应的尺度能差 1.25 倍；
对"视场最接近"的那一对（ACI 11.0 mm × WATSON 13.4 mm），129 档预测的尺度 0.53
**落在我搜过的 0.62–1.23 之外**——所以那次的"零"不算数。

这次的做法：
  * 模板 = ACI 中心一块（自带物理尺度 6.66 µm/px，带通 0.06–0.5 mm）
  * 对每个候选尺度 s（0.30–1.45，对数取点）：推出 WATSON 的 µm/px = 6.66/s，
    按它把 WATSON 带通到**同一个物理波段**，再让模板在整幅 WATSON 上搜
  * 36 个转角 × 双极性；零假设 = 把 WATSON 相位随机化后同样全搜

判读：真对上时应当在某个 s 上出现 >0.4 的尖峰（阳性对照 0.95）；
      若峰值 ≈ 零假设 → 这张 WATSON 里没有 ACI 的内容。

用法：python wide_scale_search.py [--stem sif_0186_0683466024_335...]
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import cv2

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from colorize_pair import ACI_IMG, ACI_XML, OUT, load_cam, read_bsq, to8
from sign_test import bandpass, zs, search
from watson_pose_test import gray, NEAR

ACI_UM = 6.66
BAND = (0.06, 0.5)
CROP = 1100          # ACI 中心取块（px）


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stem", default=NEAR)
    args = ap.parse_args()

    aci = load_cam(ACI_XML)
    aci8 = to8(read_bsq(ACI_IMG, aci)[0]).astype(np.float32)
    cy, cx = aci["nl"] // 2, aci["ns"] // 2
    sub = aci8[cy - CROP // 2:cy + CROP // 2, cx - CROP // 2:cx + CROP // 2]
    ppm_a = 1000.0 / ACI_UM
    A = np.ascontiguousarray(zs(bandpass(sub, BAND[0] * ppm_a / 2.355,
                                         BAND[1] * ppm_a / 2.355)).astype(np.float32))
    print("模板：ACI 中心 %dx%d px（%.2f mm 见方），带通 %.2f–%.2f mm"
          % (CROP, CROP, CROP * ACI_UM / 1000, BAND[0], BAND[1]))

    gw, camw = gray(args.stem)
    print("目标：%s  %dx%d px" % (args.stem[:52], camw["ns"], camw["nl"]))

    vals = []
    for s in np.geomspace(0.30, 1.45, 9):
        um = ACI_UM / s                      # 该尺度下 WATSON 的像素尺度
        ppm = 1000.0 / um
        hi = min(BAND[1], min(gw.shape) * um / 1000.0 / 8.0)
        W = np.ascontiguousarray(zs(bandpass(gw, BAND[0] * ppm / 2.355,
                                             hi * ppm / 2.355)).astype(np.float32))
        k = max(3, int(round(4.0 * ppm)) | 1)
        loc = cv2.blur(gw, (k, k))
        W[loc < 0.55 * float(np.median(loc))] = 0.0
        bmax, bmin, amax = -9.0, 9.0, None
        for ang in range(0, 360, 15):
            mx, lm = search(W, A, ang, s, "max")
            mn, ln = search(W, A, ang, s, "min")
            if np.isfinite(mx) and mx > bmax:
                bmax, amax = mx, (ang, lm)
            if np.isfinite(mn) and mn < bmin:
                bmin, amin = mn, (ang, ln)
        vals.append((s, um, bmax, bmin, amax, amin))
        print("  s=%.3f（WATSON %.2f µm/px，视场 %.1f mm）：max %+.3f%s  min %+.3f"
              % (s, um, camw["ns"] * um / 1000, bmax,
                 (" @%d°%s" % amax) if amax else "", bmin), flush=True)

    best = max(vals, key=lambda v: max(v[2], -v[3]))
    print("\n最好的一组：s=%.3f（%.2f µm/px）max %+.3f / min %+.3f"
          % (best[0], best[1], best[2], best[3]))

    # 零假设：相位随机化的 WATSON，同一套搜索
    rng = np.random.default_rng(5)
    ph = rng.uniform(0, 2 * np.pi, gw.shape)
    ph[0, 0] = 0
    ph[0, :] = 0
    ph[:, 0] = 0
    gr = np.real(np.fft.ifft2(np.abs(np.fft.fft2(gw)) * np.exp(1j * ph)))
    s = best[0]
    um = ACI_UM / s
    ppm = 1000.0 / um
    Wr = np.ascontiguousarray(zs(bandpass(gr, BAND[0] * ppm / 2.355,
                                          min(BAND[1], min(gw.shape) * um / 1000.0 / 8.0)
                                          * ppm / 2.355)).astype(np.float32))
    nb = -9.0
    for ang in range(0, 360, 15):
        mx, _ = search(Wr, A, ang, s, "max")
        if np.isfinite(mx):
            nb = max(nb, mx)
    print("零假设（相位随机化，同一 s）：max %+.3f" % nb)
    return 0


if __name__ == "__main__":
    sys.exit(main())
