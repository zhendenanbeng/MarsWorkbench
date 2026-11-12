# -*- coding: utf-8 -*-
"""线索 1 的判决实验：grad_ratio 这个尺度探针**能不能**测出采样率差异？

做法：把 ACI 降采样 2.83 倍（等于把"ACI 的采样率"改成 WATSON 的），
再用 scale_probe.py 里的同一个度量量它。如果降采样后数值几乎不变，
说明这个度量测的是"相对像素网格的相关长度"，对采样率是盲的——
那么 4.88 vs 4.82（比值 1.01）就没有任何关于两台相机采样率的信息。

同时给出正面对照：把 WATSON 升采样/降采样，看度量怎么动。

用法：python metric_blindness.py
"""
from __future__ import annotations

import os
import sys

import numpy as np
import cv2

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from colorize_pair import (ACI_IMG, ACI_XML, WAT_IMG, WAT_XML,             # noqa: E402
                           load_cam, read_bsq, to8)
from scale_probe import grad_ratio, autocorr_len                          # noqa: E402


def main():
    wat = load_cam(WAT_XML)
    aci = load_cam(ACI_XML)
    wr = read_bsq(WAT_IMG, wat)
    ar = read_bsq(ACI_IMG, aci)
    b = [i - 1 for i in wat["rgb"]]
    wgray = cv2.cvtColor(np.stack([to8(wr[b[i]]) for i in range(3)], axis=-1),
                         cv2.COLOR_RGB2GRAY)
    aci8 = to8(ar[0])

    # 用两块"纹理最干净"的区域：ACI 中心，WATSON 岩石上的亮区
    apatch = aci8[300:820, 540:1060]
    wpatch = wgray[200:720, 200:720]

    def probe(name, img):
        g = grad_ratio(img, 3.0)
        l, _ = autocorr_len(img, 3.0)
        print("  %-34s 梯度/标准差 %.3f   相关长度 %s px"
              % (name, g, ("%.0f" % l) if np.isfinite(l) else ">80"))
        return g

    print("原始：")
    ga = probe("ACI 原图", apatch)
    gw = probe("WATSON 原图", wpatch)

    print("\n把 ACI 按几何预测的倍率重采样（0.353 倍 = 变成 WATSON 的采样率）：")
    for f in (0.353, 0.5, 1.0, 2.0, 2.83):
        # 先按比例抗混叠模糊，再重采样——这才是"改变采样率"的正确做法
        if f < 1:
            sig = 0.5 * (1.0 / f - 1.0) ** 0.5 * 1.0
            small = cv2.resize(cv2.GaussianBlur(apatch.astype(np.float32), (0, 0), sig),
                               None, fx=f, fy=f, interpolation=cv2.INTER_AREA)
        else:
            small = cv2.resize(apatch, None, fx=f, fy=f, interpolation=cv2.INTER_CUBIC)
        probe("ACI ×%.3f（%d×%d）" % (f, small.shape[1], small.shape[0]), small)

    print("\n结论判读：如果这些数值几乎不动 → 度量对采样率是盲的。")
    print("（参考：%s）" % ("ACI %.2f vs WATSON %.2f，比值 %.2f"
                          % (ga, gw, ga / gw)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
