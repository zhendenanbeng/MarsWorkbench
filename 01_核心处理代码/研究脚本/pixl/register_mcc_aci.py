# -*- coding: utf-8 -*-
"""MCC ↔ ACI 图像配准（PIXL spike 的核心验证，可复跑）。

链路：SHERLOC RLS(ACI 像素) → [ACI↔MCC 单应] → MCC 像素 → rxl(PMC) → rqc 化学。

用法：python register_mcc_aci.py
依赖：本地已有 06_算法实验与依据/pixl/mcc_bellegarde.IMG 和练习数据里的 ACI 图。
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import cv2

# 旧开发树位置：重做实验时改成你的本地路径
DEV_ROOT = Path.home() / "Desktop" / "Sherlock"
sys.path.insert(0, str(DEV_ROOT / "01_软件工作台"))
from sherloc.vicar import load_aci, read_vicar          # noqa: E402

HERE = str(DEV_ROOT / "06_算法实验与依据" / "pixl")
ACI = str(DEV_ROOT / "02_练习数据与验证基准" / "images" / "aci"
          / "SC3_0186_0683476704_113ECM_N0070000SRLC10500_0000LMJ05.IMG")


def stretch(g, lo_pct=0.5, hi_pct=99.5):
    g = g.astype(np.float32)
    lo, hi = np.percentile(g, (lo_pct, hi_pct))
    return np.clip((g - lo) / max(hi - lo, 1e-9) * 255, 0, 255).astype(np.uint8)


def prep(g):
    g = g.astype(np.float32)
    bg = cv2.GaussianBlur(g, (0, 0), sigmaX=max(g.shape) / 12.0)
    g = cv2.normalize(g - bg, None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
    return cv2.createCLAHE(2.0, (8, 8)).apply(g)


def load_mcc_gray(path):
    img, _ = read_vicar(path)
    return stretch(img.astype(np.float64))


def main():
    mcc = load_mcc_gray(HERE + r"\mcc_bellegarde.IMG")
    aci = np.asarray(load_aci(ACI).convert("L"))
    sift = cv2.SIFT_create(20000, contrastThreshold=0.01)
    flann = cv2.FlannBasedMatcher(dict(algorithm=1, trees=5), dict(checks=64))
    k2, d2 = sift.detectAndCompute(prep(mcc), None)
    k1, d1 = sift.detectAndCompute(prep(aci), None)
    ms = flann.knnMatch(d2, d1, k=2)
    bt = {}
    for m, n in ms:
        if m.distance < 0.8 * n.distance:
            t = m.trainIdx
            if t not in bt or m.distance < bt[t].distance:
                bt[t] = m
    good = list(bt.values())
    src = np.float64([k2[m.queryIdx].pt for m in good]).reshape(-1, 1, 2)
    dst = np.float64([k1[m.trainIdx].pt for m in good]).reshape(-1, 1, 2)
    M, inl = cv2.estimateAffinePartial2D(src, dst, method=cv2.RANSAC,
                                         ransacReprojThreshold=4.0,
                                         maxIters=20000, confidence=0.999)
    mask = inl.ravel().astype(bool)
    proj = cv2.transform(src, M).reshape(-1, 2)
    err = np.hypot(*(proj[mask] - dst.reshape(-1, 2)[mask]).T)
    print("内点 %d/%d  尺度(MCC->ACI)=%.3f  重投影中位 %.2f px"
          % (mask.sum(), len(good), np.hypot(M[0, 0], M[1, 0]), np.median(err)))

    # MCC 上画 ACI 框（橙）和 PMC 扫描网格（绿，来自 rxl 的位置在别处填充）
    import csv
    vis = cv2.cvtColor(mcc, cv2.COLOR_GRAY2BGR)
    h, w = aci.shape
    Minv = cv2.invertAffineTransform(M)
    c = np.float64([[0, 0], [w - 1, 0], [w - 1, h - 1], [0, h - 1]])
    prj = cv2.transform(c.reshape(-1, 1, 2), Minv).reshape(-1, 2)
    cv2.polylines(vis, [prj.astype(int)], True, (0, 80, 255), 2)
    rows = list(csv.reader(open(HERE + r"\pe__0186_0683483744_000rxl__"
                                r"00700000699274310005___j07.csv",
                                encoding="utf-8")))
    hdr = rows[0]
    ci, cj = hdr.index("PMC_0005_MCC_i"), hdr.index("PMC_0005_MCC_j")
    for r in rows[1:]:
        try:
            x, y = float(r[ci]), float(r[cj])
        except ValueError:
            continue
        cv2.circle(vis, (int(x), int(y)), 2, (0, 220, 80), -1)
    cv2.imwrite(HERE + r"\mcc_with_aci_frame_pmcs.png", vis)
    print("saved mcc_with_aci_frame_pmcs.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
