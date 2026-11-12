# -*- coding: utf-8 -*-
"""把 ACI 上最显眼的那块特征（右下那块 1.3×1.6 mm 亮斑）拿去 WATSON 里找。

这是给肉眼用的定点检验：如果一个 1.3×1.6 mm 的亮斑真在同一块地方，
它投到 WATSON 网格上应该有 ~90 px，绝对看得见。看不到就是"不是同一块地方"。

出图 out/vc_patch.png：每档 standoff 一行，
  左 = ACI 里那块（按该档尺度重采样到 WATSON 网格）
  右 = WATSON 同一位置同样大小的块
两边同一物理尺度、同样放大倍数。

用法：python patch_check.py
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

ACI_UM = 6.66
HALF_MM = 2.2          # 取 ±2.2 mm 的方块（岩面）


def stretch(g, lo=0.5, hi=99.5):
    a, b = np.percentile(g, (lo, hi))
    return np.clip((g - a) / max(b - a, 1e-9) * 255, 0, 255).astype(np.uint8)


def main():
    wat = load_cam(WAT_XML)
    aci = load_cam(ACI_XML)
    wr = read_bsq(WAT_IMG, wat)
    ar = read_bsq(ACI_IMG, aci)
    b = [i - 1 for i in wat["rgb"]]
    wgray = cv2.cvtColor(np.stack([to8(wr[b[i]]) for i in range(3)], axis=-1),
                         cv2.COLOR_RGB2GRAY)
    aci8 = to8(ar[0])

    # 找 ACI 里最亮的那块连通区（=那块亮斑）
    thr = np.percentile(aci8, 98.0)
    mask = (aci8 >= thr).astype(np.uint8)
    mask = cv2.morphologyEx(mask, cv2.MORPH_CLOSE, np.ones((9, 9), np.uint8))
    n, lab, stats, cent = cv2.connectedComponentsWithStats(mask)
    k = 1 + int(np.argmax(stats[1:, 4]))
    cx, cy = cent[k]
    print("ACI 最亮连通区：中心 (%.0f, %.0f)，面积 %d px，外接框 %dx%d px"
          % (cx, cy, stats[k, 4], stats[k, 2], stats[k, 3]))
    print("  → 岩面上的尺寸约 %.2f × %.2f mm（按 %.2f µm/px）"
          % (stats[k, 2] * ACI_UM / 1000, stats[k, 3] * ACI_UM / 1000, ACI_UM))

    half_aci = int(HALF_MM * 1000 / ACI_UM)
    x0a, x1a = int(cx - half_aci), int(cx + half_aci)
    y0a, y1a = int(cy - half_aci), int(cy + half_aci)
    aci_patch = aci8[max(0, y0a):y1a, max(0, x0a):x1a]
    print("ACI 取块 %dx%d px" % (aci_patch.shape[1], aci_patch.shape[0]))

    rows = []
    for d in (108.0, 129.0, 160.0):
        Hm, _ = homography(aci, wat, d)
        p = Hm @ np.array([cx, cy, 1.0])
        px, py = p[0] / p[2], p[1] / p[2]
        s = 0.353 * 129.0 / d          # 采样率之比随 standoff 变
        n_px = int(2 * HALF_MM * 1000 / ACI_UM * s)
        # ACI 那块重采样到 WATSON 网格尺度（先按比例糊一下再抽）
        sig = 0.55 * max(0.0, (1.0 / s) - 1.0) ** 0.5
        small = cv2.resize(cv2.GaussianBlur(aci_patch.astype(np.float32), (0, 0), sig),
                           (n_px, n_px), interpolation=cv2.INTER_AREA)
        # WATSON 同位置
        wx0, wy0 = int(px - n_px / 2), int(py - n_px / 2)
        wx1, wy1 = wx0 + n_px, wy0 + n_px
        if wx0 < 0 or wy0 < 0 or wx1 > wat["ns"] or wy1 > wat["nl"]:
            print("  standoff %.0f mm：预测位置 (%.0f,%.0f) 出界" % (d, px, py))
            continue
        wpatch = wgray[wy0:wy1, wx0:wx1]
        z = 4
        a_big = cv2.resize(stretch(small), (n_px * z, n_px * z),
                           interpolation=cv2.INTER_NEAREST)
        w_big = cv2.resize(stretch(wpatch), (n_px * z, n_px * z),
                           interpolation=cv2.INTER_NEAREST)
        lab_a = np.full((34, a_big.shape[1], 3), 255, np.uint8)
        cv2.putText(lab_a, "ACI bright patch -> WATSON grid (%.0f mm)" % d,
                    (6, 24), cv2.FONT_HERSHEY_SIMPLEX, 0.55, (0, 0, 0), 2, cv2.LINE_AA)
        lab_w = np.full((34, w_big.shape[1], 3), 255, np.uint8)
        cv2.putText(lab_w, "WATSON same place, same size (%dx%d px)"
                    % (n_px, n_px), (6, 24), cv2.FONT_HERSHEY_SIMPLEX, 0.55,
                    (0, 0, 0), 2, cv2.LINE_AA)
        row = np.vstack([np.hstack([lab_a, lab_w]),
                         np.hstack([cv2.cvtColor(a_big, cv2.COLOR_GRAY2RGB),
                                    cv2.cvtColor(w_big, cv2.COLOR_GRAY2RGB)])])
        rows.append(row)
        print("  standoff %.0f mm：预测位置 (%.0f, %.0f)，块 %d px，"
              "WATSON 该处亮度中位 %.0f（全图中位 %.0f）"
              % (d, px, py, n_px, np.median(wpatch), np.median(wgray)))

    if rows:
        h = max(r.shape[0] for r in rows)
        w = max(r.shape[1] for r in rows)
        rows = [cv2.copyMakeBorder(r, 0, h - r.shape[0], 0, w - r.shape[1],
                                   cv2.BORDER_CONSTANT, value=(255, 255, 255))
                for r in rows]
        out = np.vstack(rows)
        cv2.imwrite(os.path.join(OUT, "vc_patch.png"),
                    cv2.cvtColor(out, cv2.COLOR_RGB2BGR))
        print("写到 out/vc_patch.png (%dx%d)" % (out.shape[1], out.shape[0]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
