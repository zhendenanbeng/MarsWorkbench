# -*- coding: utf-8 -*-
"""直接看：ACI 与**视场几乎相同**的那张 WATSON（23.4 mm 档，13.4 mm vs 11 mm）。

之前的比对都在 2.83 倍采样差上做（论文那张），这次换一对**尺度最接近**的：
    ACI         6.66 µm/px，视场 11.0 mm，standoff ≈108 mm（RLS 锚定）
    WATSON srlc00621（0683466024 全画幅） 8.10 µm/px，视场 13.4 mm，standoff 23.4 mm
两者都瞄 T（偏移 0.8–0.9 mm），采样只差 1.2 倍，视场只差 1.2 倍。
如果配对是对的，**大尺度的形状**（暗带、亮斑轮廓、裂缝走向）应当肉眼可认。

出图：
  out/vc_pair_srlc00621.png  上排 原图（各自拉对比度）| 下排 结构图 | 末排 红绿叠加
  另存 out/vc_aci_in_watson.png（把 ACI 投到 WATSON 网格，与 WATSON 并排）

用法：python pair_srlc00621.py
"""
from __future__ import annotations

import os
import sys

import numpy as np
import cv2

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from colorize_pair import (ACI_IMG, ACI_XML, OUT, homography, load_cam,
                           read_bsq, to8, WARP_FLAGS)
from watson_pose_test import gray, NEAR

ACI_UM = 6.66
WAT_UM = 8.10
D_ACI_MM = 108.0        # RLS 锚定


def stretch(g, lo=1.0, hi=99.0):
    a, b = np.percentile(g, (lo, hi))
    return np.clip((g - a) / max(b - a, 1e-9) * 255, 0, 255).astype(np.uint8)


def lp(g, um, lo_mm, hi_mm):
    """按物理尺度带通（σ 单位 px）。"""
    ppm = 1000.0 / um
    return cv2.GaussianBlur(g, (0, 0), lo_mm * ppm / 2.355) - \
        cv2.GaussianBlur(g, (0, 0), hi_mm * ppm / 2.355)


def zs(a):
    s = float(a.std())
    return (a - float(a.mean())) / (s if s > 1e-9 else 1.0)


def bar(w, text):
    b = np.full((34, w, 3), 255, np.uint8)
    cv2.putText(b, text, (6, 24), cv2.FONT_HERSHEY_SIMPLEX, 0.55, (0, 0, 0), 2,
                cv2.LINE_AA)
    return b


def main():
    aci = load_cam(ACI_XML)
    aci8 = to8(read_bsq(ACI_IMG, aci)[0]).astype(np.float32)
    gw, camw = gray(NEAR)                       # WATSON 全画幅 (1648×1200, 8.10 µm/px)
    print("ACI %dx%d   WATSON %dx%d" % (aci["ns"], aci["nl"], camw["ns"], camw["nl"]))

    # ACI → WATSON 的单应（岩面在 ACI 前 108 mm），把 WATSON 采到 ACI 网格上
    Hm, dst = homography(aci, camw, D_ACI_MM)
    warp = cv2.warpPerspective(gw.astype(np.float32), Hm, (aci["ns"], aci["nl"]),
                               flags=cv2.INTER_LINEAR | cv2.WARP_INVERSE_MAP)
    valid = cv2.warpPerspective(np.full((camw["nl"], camw["ns"]), 255, np.uint8),
                                Hm, (aci["ns"], aci["nl"]),
                                flags=cv2.INTER_NEAREST | cv2.WARP_INVERSE_MAP) > 127
    print("几何落点：ACI 四角 → WATSON %s；有效覆盖 %.0f%%"
          % (np.array2string(dst, precision=0).replace("\n", " "),
             100.0 * valid.mean()))

    # 大尺度形状（0.3–2 mm）——照明差最可能保留下来的那一档
    Ab = zs(lp(aci8, ACI_UM, 0.3, 2.0))
    Wb = zs(lp(warp, WAT_UM, 0.3, 2.0))
    Wb[~valid] = 0
    rg = np.full((aci["nl"], aci["ns"], 3), 128, np.uint8)
    rg[..., 0] = np.clip(128 + 70 * Ab, 0, 255)
    rg[..., 1] = np.clip(128 + 70 * Wb, 0, 255)
    rg[~valid] = 0

    def g3(a):
        return cv2.cvtColor(np.clip(128 + 60 * a, 0, 255).astype(np.uint8),
                            cv2.COLOR_GRAY2RGB)

    raw_a = cv2.cvtColor(stretch(aci8), cv2.COLOR_GRAY2RGB)
    raw_w = cv2.cvtColor(stretch(np.where(valid, warp, 0)), cv2.COLOR_GRAY2RGB)
    rows = [np.hstack([bar(aci["ns"], "ACI 6.66 um/px, 11.0 mm field"),
                       bar(aci["ns"], "WATSON srlc00621 warped to ACI grid (108 mm)")]),
            np.hstack([raw_a, raw_w]),
            np.hstack([bar(aci["ns"], "structure 0.3-2 mm: ACI"),
                       bar(aci["ns"], "structure 0.3-2 mm: WATSON")]),
            np.hstack([g3(Ab), g3(Wb)]),
            np.hstack([bar(aci["ns"], "red = ACI, green = WATSON (yellow = agree)"),
                       bar(aci["ns"], "same physical band 0.3-2 mm")]),
            np.hstack([rg, rg * 0]),
            ]
    # 对半裁开，避免整图缩得太小
    out = np.vstack(rows)
    cv2.imwrite(os.path.join(OUT, "vc_pair_srlc00621.png"),
                cv2.cvtColor(out, cv2.COLOR_RGB2BGR))
    half = aci["ns"] // 2
    for i, (x0, x1, nm) in enumerate(((0, half, "L"), (half, aci["ns"], "R"))):
        sub = np.vstack([r[:, x0:x1] for r in rows])
        cv2.imwrite(os.path.join(OUT, "vc_pair_srlc00621_%s.png" % nm),
                    cv2.cvtColor(sub, cv2.COLOR_RGB2BGR))
    print("写到 out/vc_pair_srlc00621.png（%dx%d）及 _L/_R 半幅"
          % (out.shape[1], out.shape[0]))
    m = valid
    a, b = Ab[m], Wb[m]
    av, bv = a - a.mean(), b - b.mean()
    print("0.3–2 mm 波段、有效区内的 NCC = %+.3f"
          % float((av * bv).sum() / np.sqrt((av ** 2).sum() * (bv ** 2).sum())))
    return 0


if __name__ == "__main__":
    sys.exit(main())
