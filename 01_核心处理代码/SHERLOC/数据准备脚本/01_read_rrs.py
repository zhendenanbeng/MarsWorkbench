#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SHERLOC 练习脚本 01：读取 RRS 光谱 + 按指南的双阈值筛候选矿物峰
==================================================================
对应操作手册「路线 A」的步骤 A-2 ~ A-5。

做什么：
  1. 把 RRS（扣暗背景的科学光谱）读进来，拆成 波长表 + Region1/2/3 三张强度表
  2. 把发射波长(nm)换算成拉曼位移(cm-1)
  3. 裁掉首尾各 50 个不感光的像素
  4. 迭代多项式减基线
  5. 按指南阈值筛峰：800-1200 cm-1 用 FWHM 25-200，3000-4000 cm-1 用 FWHM 40-300
    且 SNR > 3（噪声标准差取 2000-2100 cm-1 的 RMS）
  6. 导出 中位光谱 CSV 和 候选峰 CSV

重要提醒：
  SHERLOC 一次扫描里每个点的矿物信号可能差别很大，中位光谱经常把峰"抹平"。
  所以本脚本对 100 个点**逐个**筛峰，而不是只看中位光谱。

用法：
    python 01_read_rrs.py <RRS文件路径> [输出目录]

依赖：
    pip install numpy pandas

注意：这是帮你理解流程的入门脚本，不是 SHERLOC 团队的官方代码。
      基线和阈值都应该按你自己的数据继续调优。
"""

import os
import sys

import numpy as np
import pandas as pd

LASER_NM = 248.6                  # SHERLOC 的 DUV 激光波长
PIXEL_LO, PIXEL_HI = 50, 2098     # SCCD 首尾各 50 个像素不感光，作图/分析要去掉

# 指南里的阈值
FWHM_MIN_MINERAL, FWHM_MAX_MINERAL = 25.0, 200.0     # 800-1200 cm-1 矿物峰
FWHM_MIN_WATER, FWHM_MAX_WATER = 40.0, 300.0         # 3000-4000 cm-1 水合带
SNR_MIN = 3.0                                        # 信噪比下限
NOISE_WINDOW = (2000.0, 2100.0)                      # 用来算噪声 σ 的安静区间


def nm_to_raman_shift(nm):
    """发射波长(nm) -> 拉曼位移(cm-1)： 1/λ_exc − 1/λ_em"""
    nm = np.asarray(nm, dtype=float)
    return 1.0e7 / LASER_NM - 1.0e7 / nm


def load_rrs(path):
    """读 RRS。结构固定为：1 张波长表 + 3 张强度表（对应 3 个 region）。"""
    with open(path, "r") as f:
        contents = f.readlines()

    # 行数 = 光谱数×3(region) + 4(表头) + 4(列名) + 1(波长表)  ->  反推光谱数
    n_spectra = int((len(contents) - 9) / 3)

    wavelength = np.array([float(v) for v in contents[2].split(",")])
    r1 = pd.read_csv(path, nrows=n_spectra, skiprows=4).to_numpy()
    r2 = pd.read_csv(path, nrows=n_spectra, skiprows=4 + n_spectra + 2).to_numpy()
    r3 = pd.read_csv(path, nrows=n_spectra, skiprows=4 + 2 * n_spectra + 4).to_numpy()
    return wavelength, r1, r2, r3, n_spectra


def poly_baseline(x, y, degree=5, n_iter=10, clip=1.5):
    """迭代 sigma-clipping 多项式基线：反复拟合并丢掉明显高于基线的峰。"""
    keep = np.ones_like(y, dtype=bool)
    coef = np.polyfit(x[keep], y[keep], degree)
    for _ in range(n_iter):
        resid = y - np.polyval(coef, x)
        sigma = resid[keep].std()
        keep = resid < clip * sigma
        if keep.sum() < degree + 2:
            break
        coef = np.polyfit(x[keep], y[keep], degree)
    return np.polyval(coef, x)


def find_peaks(x, y, min_fwhm, max_fwhm, min_height):
    """找局部极大值并量半峰宽(FWHM)。返回 (位置, 高度, FWHM) 列表。"""
    peaks = []
    for i in range(1, len(y) - 1):
        if not (y[i] > y[i - 1] and y[i] >= y[i + 1]):
            continue
        if y[i] < min_height:
            continue

        half = y[i] / 2.0
        left = i
        while left > 0 and y[left] > half:
            left -= 1
        right = i
        while right < len(y) - 1 and y[right] > half:
            right += 1

        fwhm = x[right] - x[left]
        if min_fwhm <= fwhm <= max_fwhm:
            peaks.append((float(x[i]), float(y[i]), float(fwhm)))
    return peaks


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    rrs_path = sys.argv[1]
    out_dir = sys.argv[2] if len(sys.argv) > 2 else os.path.dirname(os.path.abspath(rrs_path))
    os.makedirs(out_dir, exist_ok=True)
    stem = os.path.splitext(os.path.basename(rrs_path))[0]

    # ---------- 1. 读文件 ----------
    wavelength_nm, r1, r2, r3, n_spectra = load_rrs(rrs_path)
    print("=" * 70)
    print("文件        :", os.path.basename(rrs_path))
    print("光谱条数    :", n_spectra)
    print("每条光谱通道:", wavelength_nm.size)
    print("波长范围    : %.2f - %.2f nm" % (wavelength_nm.min(), wavelength_nm.max()))

    # ---------- 2. 换算拉曼位移 + 裁掉不感光像素 ----------
    shift = nm_to_raman_shift(wavelength_nm)[PIXEL_LO:PIXEL_HI]
    r1 = r1[:, PIXEL_LO:PIXEL_HI]
    r2 = r2[:, PIXEL_LO:PIXEL_HI]
    r3 = r3[:, PIXEL_LO:PIXEL_HI]
    print("拉曼位移范围: %.0f - %.0f cm-1" % (shift.min(), shift.max()))

    order = np.argsort(shift)
    shift = shift[order]
    r1 = r1[:, order]
    wl_sorted = wavelength_nm[PIXEL_LO:PIXEL_HI][order]

    # ---------- 3. 取中位光谱（比单条光谱稳，适合先看整体）----------
    med1 = np.median(r1, axis=0)
    med2 = np.median(r2, axis=0)[order]
    med3 = np.median(r3, axis=0)[order]

    # 顺带看一眼 Region 2/3（荧光区）大概落在什么波长上
    for tag, med in (("Region2", med2), ("Region3", med3)):
        hot = med > 5.0 * np.median(med[med > 0]) if (med > 0).any() else None
        if hot is not None and hot.any():
            print("%s 强信号波长: %.1f - %.1f nm" % (tag, wl_sorted[hot].min(), wl_sorted[hot].max()))

    # ---------- 4. 只看 800-4000 cm-1 ----------
    band = (shift >= 800.0) & (shift <= 4000.0)
    x = shift[band]
    y = med1[band]

    # ---------- 5. 中位光谱减基线（先给个整体印象）----------
    baseline = poly_baseline(x, y)
    y_sub = y - baseline

    # ---------- 6. 导出中位光谱 ----------
    spec_df = pd.DataFrame({
        "raman_shift_cm-1": x,
        "R1_intensity": y,
        "R1_baseline": baseline,
        "R1_baseline_subtracted": y_sub,
    })
    spec_path = os.path.join(out_dir, stem + "_median.csv")
    spec_df.to_csv(spec_path, index=False)

    # ---------- 7. 逐点减基线 + 双阈值筛峰 ----------
    # 指南要求两个区间分别处理：
    #   矿物峰 800-1200  -> 用 800-4000 整体做多项式基线
    #   水合带 3200-3700 -> 单独在 3000-4000 里做低次基线，否则宽的水带会被基线吃掉
    mineral_band = (shift >= 800.0) & (shift <= 4000.0)
    water_band = (shift >= 3000.0) & (shift <= 4000.0)
    xm = shift[mineral_band]
    quiet = (xm >= NOISE_WINDOW[0]) & (xm <= NOISE_WINDOW[1])
    mineral_sel = (xm >= 800.0) & (xm <= 1200.0)

    xw = shift[water_band]
    water_sel = (xw >= 3200.0) & (xw <= 3700.0)

    rows = []
    for i in range(n_spectra):
        # --- 矿物峰 ---
        ym_sub = r1[i][mineral_band] - poly_baseline(xm, r1[i][mineral_band], degree=5)
        sigma = float(np.sqrt(np.mean(ym_sub[quiet] ** 2)))
        if not np.isfinite(sigma) or sigma <= 0:
            continue
        for pos, height, fwhm in find_peaks(
            xm[mineral_sel], ym_sub[mineral_sel],
            FWHM_MIN_MINERAL, FWHM_MAX_MINERAL, SNR_MIN * sigma,
        ):
            rows.append({
                "点号": i, "区间": "矿物峰",
                "峰位_cm-1": round(pos, 1), "峰高": round(height, 1),
                "FWHM_cm-1": round(fwhm, 1), "SNR": round(height / sigma, 2),
                "噪声σ": round(sigma, 2), "判定": "高置信度",
            })

        # --- 水合带（单独拟合）---
        yw_sub = r1[i][water_band] - poly_baseline(xw, r1[i][water_band], degree=1)
        for pos, height, fwhm in find_peaks(
            xw[water_sel], yw_sub[water_sel],
            FWHM_MIN_WATER, FWHM_MAX_WATER, SNR_MIN * sigma,
        ):
            rows.append({
                "点号": i, "区间": "水合带",
                "峰位_cm-1": round(pos, 1), "峰高": round(height, 1),
                "FWHM_cm-1": round(fwhm, 1), "SNR": round(height / sigma, 2),
                "噪声σ": round(sigma, 2), "判定": "高置信度",
            })

    peaks_df = pd.DataFrame(rows)
    if not peaks_df.empty:
        peaks_df = peaks_df.sort_values("SNR", ascending=False).reset_index(drop=True)

    peak_path = os.path.join(out_dir, stem + "_peaks.csv")
    peaks_df.to_csv(peak_path, index=False, encoding="utf-8-sig")

    print("\n逐点筛查结果（已同时满足 FWHM 和 SNR>3）：")
    if peaks_df.empty:
        print("  （这一次扫描没筛出峰 —— 可以放宽阈值，或换一次扫描再看）")
    else:
        n_pts = peaks_df["点号"].nunique()
        print("  共 %d 个峰，分布在 %d 个点（该扫描共 %d 个点）"
              % (len(peaks_df), n_pts, n_spectra))
        print()
        print(peaks_df.head(20).to_string(index=False))
        if len(peaks_df) > 20:
            print("  ...（完整结果见导出的 CSV）")

    print("\n已导出：")
    print("  中位光谱 :", spec_path)
    print("  候选峰表 :", peak_path)
    print("=" * 70)


if __name__ == "__main__":
    main()
