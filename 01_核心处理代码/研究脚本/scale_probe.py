# -*- coding: utf-8 -*-
"""不用配准，只比"纹理粗细"——判断两图之间的尺度差到底是多少。

原理：同一块岩石表面，两台相机采样率不同，纹理的特征尺度（自相关长度）
就会按采样率之比缩放。所以在两图上分别量出纹理相关长度 L，再跟几何预测的
倍率比一比，就能判断尺度对不对——**完全不需要先对齐**。

量法两条，互相印证：
  1. 高通图的自相关降到 0.5 时的半径（L_ac）
  2. 高通图的平均梯度 / 标准差（∝ 1/L，对噪声更稳）

用法：python scale_probe.py
"""

from __future__ import annotations

import os
import sys

import numpy as np
import cv2

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from colorize_pair import (ACI_IMG, ACI_XML, WAT_IMG, WAT_XML,          # noqa: E402
                           homography, load_cam, read_bsq, to8)


def highpass(img, sigma):
    g = img.astype(np.float32)
    return g - cv2.GaussianBlur(g, (0, 0), sigma)


def autocorr_len(patch, sigma=3.0):
    """高通后做二维自相关，取径向平均降到 0.5 的半径（像素）。"""
    g = highpass(patch, sigma)
    g = g - g.mean()
    win = np.outer(np.hanning(g.shape[0]), np.hanning(g.shape[1]))
    g = g * win
    F = np.fft.rfft2(g)
    ac = np.fft.irfft2(F * np.conj(F), s=g.shape)
    ac = np.fft.fftshift(ac, axes=0)
    ac = ac / ac[0, 0]
    h, w = ac.shape
    yy = np.arange(h) - h // 2
    xx = np.arange(w)
    r = np.sqrt(yy[:, None] ** 2 + xx[None, :] ** 2)
    prof = []
    for rad in np.arange(0, 80, 1.0):
        sel = (r >= rad) & (r < rad + 1)
        prof.append(ac[sel].mean() if sel.any() else np.nan)
    prof = np.array(prof)
    below = np.where(prof < 0.5)[0]
    return (float(below[0]) if len(below) else float("nan")), prof


def grad_ratio(patch, sigma=3.0):
    """平均梯度 / 标准差。对随机纹理，这个量正比于 1/相关长度。"""
    g = highpass(patch, sigma)
    gx = cv2.Sobel(g, cv2.CV_32F, 1, 0, ksize=3)
    gy = cv2.Sobel(g, cv2.CV_32F, 0, 1, ksize=3)
    mag = np.sqrt(gx ** 2 + gy ** 2)
    sd = float(g.std())
    return float(mag.mean() / sd) if sd > 0 else float("nan")


def main():
    wat = load_cam(WAT_XML)
    aci = load_cam(ACI_XML)
    wr = read_bsq(WAT_IMG, wat)
    ar = read_bsq(ACI_IMG, aci)
    b = [i - 1 for i in wat["rgb"]]
    wcol = np.stack([to8(wr[b[i]]) for i in range(3)], axis=-1)
    wgray = cv2.cvtColor(wcol, cv2.COLOR_RGB2GRAY)
    aci8 = to8(ar[0])

    H, dst = homography(aci, wat, 129.0)
    cx, cy = dst[:, 0].mean(), dst[:, 1].mean()
    print("几何预测：ACI 视场投到 WATSON 的 (%.0f, %.0f)，占 %.0f×%.0f px"
          % (cx, cy, dst[:, 0].max() - dst[:, 0].min(),
             dst[:, 1].max() - dst[:, 1].min()))
    half = 260
    x0, x1 = int(max(0, cx - half)), int(min(wat["ns"], cx + half))
    y0, y1 = int(max(0, cy - half)), int(min(wat["nl"], cy + half))
    wpatch = wgray[y0:y1, x0:x1]
    apatch = aci8[340:860, 560:1080]
    print("取块：WATSON %s，ACI %s" % (wpatch.shape, apatch.shape))

    for sigma in (2.0, 3.0, 4.0, 6.0):
        la, _ = autocorr_len(apatch, sigma)
        lw, _ = autocorr_len(wpatch, sigma)
        ga = grad_ratio(apatch, sigma)
        gw = grad_ratio(wpatch, sigma)
        ratio = (lw / la) if la else float("nan")
        print("  σ=%.0f  ACI 相关长度 %5.1f px  WATSON %5.1f px  → L_W/L_A = %.2f"
              % (sigma, la, lw, ratio))
        print("        梯度/标准差  ACI %.3f  WATSON %.3f  → 比值 %.2f（约等于 L_A/L_W）"
              % (ga, gw, (ga / gw) if gw else float("nan")))

    mag = (dst[:, 0].max() - dst[:, 0].min()) / aci["ns"]
    print("\n几何预测的采样倍率：ACI 1 px → WATSON %.3f px" % mag)
    print("  若尺度一致，量到的 L_W/L_A 应接近 %.2f，梯度比应接近 %.2f"
          % (mag, 1 / mag))
    return 0


if __name__ == "__main__":
    sys.exit(main())
