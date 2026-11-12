# -*- coding: utf-8 -*-
"""不看几何，直接把 ACI 的**粗尺度结构图**在整张 WATSON 上搜一遍。

上一轮的所有判据都建立在"几何把 ACI 投到了 WATSON 的正确位置"之上；
如果这个前提本身有问题（原点约定、模型版本、配对），在预测位置附近
怎么搜都搜不出来。所以这一步把前提撤掉：

  * 尺度：0.20–0.80（几何预测 0.353；1.0 = 两台相机采样率相同）
  * 位置：整张 WATSON 画幅（1648×1200）任意位置
  * 朝向：二面体群 8 种（含转置与 90°/180°/270° 旋转——翻转测试没覆盖转置）
  * 波段：只用两台相机**都**有的 0.06–0.5 mm

输出：每个尺度/朝向的最佳相关峰，以及几何预测位置处（不动）的响应值。
判读：
  * 某尺度上峰显著且落在岩石上 → 尺度/朝向对，几何位置错
  * 所有尺度/朝向都平    → 两张图在共享波段上没有共同结构

用法：python blind_search.py
"""
from __future__ import annotations

import os
import sys

import numpy as np
import cv2

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from colorize_pair import (ACI_IMG, ACI_XML, WAT_IMG, WAT_XML, OUT,      # noqa: E402
                           homography, load_cam, read_bsq, to8)

ACI_UM = 6.66          # RLS HDR_500 反推的 ACI 像素尺度
GEOM_SCALE = 0.353     # 129 mm 档几何预测的 WATSON px / ACI px


def bandpass(img, lo, hi):
    g = img.astype(np.float32)
    return cv2.GaussianBlur(g, (0, 0), lo) - cv2.GaussianBlur(g, (0, 0), hi)


def zs(a):
    s = float(a.std())
    return (a - float(a.mean())) / (s if s > 1e-9 else 1.0)


def variants(a):
    """8 种：恒等 / 左右翻 / 上下翻 / 180° / 转置 / 转置+翻 / 90° / 270°"""
    return [("identity", a), ("fliplr", np.fliplr(a)), ("flipud", np.flipud(a)),
            ("rot180", np.rot90(a, 2)), ("transpose", a.T),
            ("transpose+fliplr", np.fliplr(a.T)),
            ("rot90", np.rot90(a, 1)), ("rot270", np.rot90(a, 3))]


def main():
    wat = load_cam(WAT_XML)
    aci = load_cam(ACI_XML)
    wr = read_bsq(WAT_IMG, wat)
    ar = read_bsq(ACI_IMG, aci)
    b = [i - 1 for i in wat["rgb"]]
    wgray = cv2.cvtColor(np.stack([to8(wr[b[i]]) for i in range(3)], axis=-1),
                         cv2.COLOR_RGB2GRAY)
    aci8 = to8(ar[0])

    # 两边带通到同一物理波段（0.06–0.5 mm，σ 按 µm/px 折算）
    aci_ppm = 1000.0 / ACI_UM
    wat_ppm = aci_ppm * GEOM_SCALE
    a = zs(bandpass(aci8, 0.06 * aci_ppm / 2.355, 0.5 * aci_ppm / 2.355))
    w = zs(bandpass(wgray, 0.06 * wat_ppm / 2.355, 0.5 * wat_ppm / 2.355))

    # 阴影掩膜（WATSON 局部均值过低处）——只用于报告最佳位置的"是否在岩石上"
    loc = cv2.blur(wgray.astype(np.float32), (151, 151))
    dark = loc < 0.55 * float(np.median(loc))
    print("WATSON 基准亮度中位 %.0f；暗区 %.0f%%"
          % (np.median(loc), 100 * dark.mean()))

    _, dst = homography(aci, wat, 129.0)
    pred_x, pred_y = float(dst[0, 0]), float(dst[0, 1])
    print("几何预测：ACI 左上角落在 WATSON (%.0f, %.0f)，尺度 %.3f"
          % (pred_x, pred_y, GEOM_SCALE))

    print("\n  尺度  朝向              峰    位置(x,y)      预测处   峰/中位   峰处亮度")
    rows = []
    for s in (0.20, 0.25, 0.30, 0.353, 0.40, 0.50, 0.60, 0.80):
        tw, th = int(round(aci["ns"] * s)), int(round(aci["nl"] * s))
        if tw < 40 or th < 40 or tw > wat["ns"] or th > wat["nl"]:
            continue
        small = cv2.resize(a, (tw, th), interpolation=cv2.INTER_AREA)
        for name, v in variants(small):
            v = np.ascontiguousarray(zs(v).astype(np.float32))
            resp = cv2.matchTemplate(w.astype(np.float32), v, cv2.TM_CCOEFF_NORMED)
            _, mx, _, ml = cv2.minMaxLoc(resp)
            med = float(np.median(resp))
            px, py = int(round(pred_x)), int(round(pred_y))
            at_pred = float(resp[min(max(py, 0), resp.shape[0] - 1),
                                  min(max(px, 0), resp.shape[1] - 1)]) \
                if 0 <= px < resp.shape[1] and 0 <= py < resp.shape[0] else float("nan")
            lum = float(loc[ml[1] + th // 2, ml[0] + tw // 2])
            rows.append((mx, s, name, ml, med, at_pred, lum, dark[ml[1] + th // 2,
                                                                  ml[0] + tw // 2]))
            print("  %.3f  %-16s %6.3f  (%5d,%5d) %8.3f %8.3f  %8.0f%s"
                  % (s, name, mx, ml[0], ml[1], at_pred, mx - med, lum,
                     "  ←暗区" if dark[ml[1] + th // 2, ml[0] + tw // 2] else ""))

    rows.sort(reverse=True)
    print("\n最强的 6 个（不论尺度/朝向）：")
    for mx, s, name, ml, med, at_pred, lum, dk in rows[:6]:
        print("  峰 %.3f（中位 %.3f，超出 %+.3f） 尺度 %.3f 朝向 %-16s 位置 (%d,%d)%s"
              % (mx, med, mx - med, s, name, ml[0], ml[1], "  ←暗区" if dk else ""))

    # 同尺度、同朝向（恒等）下的响应面，直接出图看有没有"峰"
    s = GEOM_SCALE
    tw, th = int(round(aci["ns"] * s)), int(round(aci["nl"] * s))
    small = zs(cv2.resize(a, (tw, th), interpolation=cv2.INTER_AREA)).astype(np.float32)
    resp = cv2.matchTemplate(w.astype(np.float32), small, cv2.TM_CCOEFF_NORMED)
    vis = np.clip(128 + 900 * (resp - np.median(resp)), 0, 255).astype(np.uint8)
    vis = cv2.cvtColor(vis, cv2.COLOR_GRAY2BGR)
    cv2.rectangle(vis, (int(pred_x), int(pred_y)), (int(pred_x) + tw, int(pred_y) + th),
                  (0, 0, 255), 2)
    cv2.imwrite(os.path.join(OUT, "bs_response.png"), vis)
    print("\n写到 out/bs_response.png（尺度 %.3f 恒等朝向的响应面，红框=几何预测位置）"
          % s)
    return 0


if __name__ == "__main__":
    sys.exit(main())
