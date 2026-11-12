# -*- coding: utf-8 -*-
"""最后一个"什么都不信"的判据：位置 × 旋转 × 尺度 全搜。

为什么必须试旋转：盲搜只覆盖了二面体群 8 种朝向（含转置），**不含任意角度旋转**。
而 WATSON / ACI 都挂在转塔上，转塔角度（JOINT 5）的建模误差表现出来就是
一张相对另一张转了若干度——这恰恰是上一轮 16 种翻转组合永远测不到的。

搜法：把 ACI 的粗尺度结构图（0.06–0.5 mm）在 WATSON 同一波段的图上
做掩膜版归一化互相关，扫 72 个旋转角 × 4 个尺度 × 全部位置。
判读：
  * 某角度上出现 >0.3 的尖峰 → 两张图**是**同一块地方，是相机模型错了
    （转角/尺度就是误差量，可以直接拿来当修正）
  * 全部 <0.15 → 这两张图内容上没有共同结构

用法：python similarity_search.py
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np
import cv2

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from colorize_pair import (ACI_IMG, ACI_XML, WAT_IMG, WAT_XML, OUT,      # noqa: E402
                           homography, load_cam, read_bsq, to8)

ACI_UM = 6.66


def bandpass(img, lo, hi):
    g = img.astype(np.float32)
    return cv2.GaussianBlur(g, (0, 0), lo) - cv2.GaussianBlur(g, (0, 0), hi)


def zs(a):
    s = float(a.std())
    return (a - float(a.mean())) / (s if s > 1e-9 else 1.0)


def main():
    wat = load_cam(WAT_XML)
    aci = load_cam(ACI_XML)
    wr = read_bsq(WAT_IMG, wat)
    ar = read_bsq(ACI_IMG, aci)
    b = [i - 1 for i in wat["rgb"]]
    wgray = cv2.cvtColor(np.stack([to8(wr[b[i]]) for i in range(3)], axis=-1),
                         cv2.COLOR_RGB2GRAY)
    aci8 = to8(ar[0])

    aci_ppm = 1000.0 / ACI_UM
    wat_ppm = aci_ppm * 0.353
    A = zs(bandpass(aci8, 0.06 * aci_ppm / 2.355, 0.5 * aci_ppm / 2.355))
    W = np.ascontiguousarray(zs(bandpass(wgray, 0.06 * wat_ppm / 2.355,
                                         0.5 * wat_ppm / 2.355)).astype(np.float32))
    # 去掉阴影区，避免"两张都平的地方"刷出假峰
    loc = cv2.blur(wgray.astype(np.float32), (151, 151))
    W[(loc < 0.55 * float(np.median(loc)))] = 0.0

    H, dst = homography(aci, wat, 129.0)
    print("开始全搜：72 个角度 × 4 个尺度 × 全部位置", flush=True)
    t0 = time.time()
    rows = []
    for scale in (0.30, 0.353, 0.44, 0.55):
        for ang in range(0, 360, 5):
            M = cv2.getRotationMatrix2D((aci["ns"] / 2.0, aci["nl"] / 2.0),
                                        ang, scale)
            t = cv2.warpAffine(A, M, (aci["ns"], aci["nl"]),
                               flags=cv2.INTER_LINEAR)
            t = np.ascontiguousarray(t.astype(np.float32))
            if t.shape[0] > W.shape[0] or t.shape[1] > W.shape[1]:
                continue
            # 掩膜：只信旋转后仍在画面内的那块（外圈是 warpAffine 填的 0）
            m = cv2.warpAffine(np.full(A.shape, 255, np.uint8), M,
                               (aci["ns"], aci["nl"]),
                               flags=cv2.INTER_NEAREST)
            m = cv2.erode(m, np.ones((41, 41), np.uint8))
            if m.sum() == 0:
                continue
            res = cv2.matchTemplate(W, t, cv2.TM_CCORR_NORMED, mask=m.astype(np.float32))
            _, mx, _, ml = cv2.minMaxLoc(res)
            # 几何预测位置（只在恒等/几何尺度下有意义）
            gx = float(res[min(max(int(dst[0, 1]), 0), res.shape[0] - 1),
                           min(max(int(dst[0, 0]), 0), res.shape[1] - 1)]) \
                if abs(scale - 0.353) < 1e-6 and ang == 0 else float("nan")
            rows.append((mx, scale, ang, ml, gx))
        best = max(r for r in rows if abs(r[1] - scale) < 1e-6)
        print("  尺度 %.3f 完成，当前最好：峰 %.3f @ %d° 位置 %s（%.0f s）"
              % (scale, best[0], best[2], best[3], time.time() - t0), flush=True)

    rows.sort(reverse=True)
    print("\n最强的 10 个（全搜 288 组）：")
    print("     峰   尺度   转角   位置(x,y)            几何预测处")
    for mx, s, ang, ml, gx in rows[:10]:
        print("  %6.3f  %.3f  %4d°  (%5d,%5d)  %s"
              % (mx, s, ang, ml[0], ml[1], ("%.3f" % gx) if np.isfinite(gx) else "—"))
    geom = [r for r in rows if abs(r[1] - 0.353) < 1e-6 and r[2] == 0]
    if geom:
        print("\n恒等/几何尺度那一组（相机模型若正确，峰就该出现在这一组）：峰 %.3f 在 (%d,%d)"
              % (geom[0][0], geom[0][3][0], geom[0][3][1]))
    vals = np.array([r[0] for r in rows])
    print("\n全部 288 组峰值的分布：中位 %.3f，90%% 分位 %.3f，最大 %.3f"
          % (float(np.median(vals)), float(np.percentile(vals, 90)), float(vals.max())))
    print("（两张无关的岩石纹理、同一波段，峰值通常落在 0.05–0.15；"
          "真正对上时应有 0.4 以上的尖峰）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
