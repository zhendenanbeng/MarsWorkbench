# -*- coding: utf-8 -*-
"""粗筛：28 张 WATSON ECM 里，哪一张真正盖住了 ACI 的视野。

思路（都是为了让"跨相机、跨尺度、跨照明"的匹配能活下来）：

1. **用梯度幅值，不用亮度**。WATSON 是蓝色 LED 照明、ACI 是环形照明，
   两边的亮度直方图差得远，但岩石的**边缘纹理**应该一致。
2. **多尺度**。ACI 是 10.07 µm/px、WATSON 约 31.85 µm/px，线性尺度差约 3.16 倍。
   真实值可能有偏差，所以搜 2.2–5.0。
3. **多旋转**。两台相机都装在转塔上、都垂直朝下，但转塔姿态不同，
   图像之间可能有几度到十几度的相对旋转。
4. **先粗后细**。先在 1/8 分辨率上把（尺度 × 旋转 × 平移）全扫一遍，
   再对赢家在 1/2 和全分辨率上局部细化。

输出：一张排名表（每张帧的最好得分），以及
`screen_best.txt`（赢家的变换参数）。**不做颜色迁移**——那是下一步的事。

用法：python screen_frames.py
"""

from __future__ import annotations

import glob
import os
import sys

import numpy as np
import cv2

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "01_软件工作台"))
from sherloc.vicar import read_vicar          # noqa: E402  项目里已验证过的 VICAR 读取

ACI_DIR = os.path.join(os.path.dirname(HERE), "02_练习数据与验证基准", "images", "aci")
ACI_NAME = "SC3_0141_0679457520_796ECM_N0051812SRLC10203_0000LUJ04.IMG"

# 尺度搜索范围：WATSON 像素 / ACI 像素 的比值 k 的候选。
# k = 10.07 / 31.85 = 0.316 是估算值，两边各放宽一倍多。
K_CANDIDATES = [0.20, 0.22, 0.24, 0.26, 0.28, 0.30, 0.316, 0.33, 0.35,
                0.38, 0.41, 0.45]
ANGLES = [-12, -8, -4, 0, 4, 8, 12]        # ACI 相对 WATSON 的旋转（度）


def read_img(path):
    """读 ACI / WATSON 的 8 位单波段 IMG，返回 (float32 数组, 标签)。"""
    arr, lab = read_vicar(path)
    return arr.astype(np.float32), lab


def features(img, blur=0):
    """把图像变成"只留结构"的特征图：去低频 → CLAHE → 梯度幅值。"""
    g = img.astype(np.float32)
    if blur:
        g = cv2.GaussianBlur(g, (blur, blur), 0)
    # 去低频：大切片的局部均值，压掉大尺度的照明梯度
    bg = cv2.GaussianBlur(g, (0, 0), sigmaX=max(g.shape) / 12.0)
    g = g - bg
    g = cv2.normalize(g, None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
    g = cv2.createCLAHE(2.0, (8, 8)).apply(g)
    gx = cv2.Sobel(g, cv2.CV_32F, 1, 0, ksize=3)
    gy = cv2.Sobel(g, cv2.CV_32F, 0, 1, ksize=3)
    mag = cv2.magnitude(gx, gy)
    return cv2.GaussianBlur(mag, (0, 0), 0.8)


def coarse_scan(aci_f, frames_f, level):
    """在 1/level 分辨率上扫（尺度 × 旋转 × 平移），返回每张帧的最好结果。"""
    a_small = cv2.resize(aci_f, None, fx=1 / level, fy=1 / level,
                         interpolation=cv2.INTER_AREA)
    ah, aw = a_small.shape
    out = []
    for name, fr in frames_f:
        f_small = cv2.resize(fr, None, fx=1 / level, fy=1 / level,
                             interpolation=cv2.INTER_AREA)
        fh, fw = f_small.shape
        best = None
        for k in K_CANDIDATES:
            tw, th = int(round(aw * k)), int(round(ah * k))
            if tw < 12 or th < 12 or tw >= fw or th >= fh:
                continue
            tpl0 = cv2.resize(a_small, (tw, th), interpolation=cv2.INTER_AREA)
            for ang in ANGLES:
                if ang:
                    M = cv2.getRotationMatrix2D((tw / 2, th / 2), ang, 1.0)
                    tpl = cv2.warpAffine(tpl0, M, (tw, th),
                                         flags=cv2.INTER_LINEAR,
                                         borderMode=cv2.BORDER_REPLICATE)
                else:
                    tpl = tpl0
                res = cv2.matchTemplate(f_small, tpl, cv2.TM_CCOEFF_NORMED)
                _, mx, _, loc = cv2.minMaxLoc(res)
                if best is None or mx > best[0]:
                    best = (float(mx), k, ang, int(loc[0]), int(loc[1]))
        out.append((name, best))
    return out


def main():
    aci, alab = read_img(os.path.join(ACI_DIR, ACI_NAME))
    print("ACI: %s  %dx%d" % (ACI_NAME[:28], aci.shape[1], aci.shape[0]))
    aci_f = features(aci)

    paths = sorted(glob.glob(os.path.join(HERE, "ecm", "*.IMG")))
    frames = []
    for p in paths:
        arr, lab = read_img(p)
        frames.append((os.path.basename(p), features(arr)))
    print("WATSON ECM 帧: %d 张（%dx%d）" % (len(frames), frames[0][1].shape[1],
                                            frames[0][1].shape[0]))

    print("\n--- 1/8 分辨率粗筛：%d 尺度 × %d 旋转 ---"
          % (len(K_CANDIDATES), len(ANGLES)))
    coarse = coarse_scan(aci_f, frames, 8)
    coarse.sort(key=lambda t: -(t[1][0] if t[1] else -1))
    for name, b in coarse:
        if b is None:
            print("  %-58s （无有效尺度）" % name[:58])
            continue
        print("  %-58s NCC=%.3f  k=%.3f  rot=%+d  last=(%d,%d)"
              % (name[:58], b[0], b[1], b[2], b[3], b[4]))

    winner, wb = coarse[0]
    print("\n赢家：%s  NCC=%.3f" % (winner, wb[0]))
    with open(os.path.join(HERE, "screen_best.txt"), "w", encoding="utf-8") as f:
        f.write("frame=%s\nncc=%.4f\nk=%s\nangle=%s\nloc=%s,%s\nlevel=8\n"
                % (winner, wb[0], wb[1], wb[2], wb[3], wb[4]))
        for name, b in coarse[:5]:
            f.write("runner-up %s %.4f\n" % (name, b[0] if b else -1))
    print("参数已写入 screen_best.txt")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
