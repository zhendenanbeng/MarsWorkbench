# -*- coding: utf-8 -*-
"""检查基线拟合的数值条件：5 次多项式直接作用在 800-4000 cm-1 上是否病态。

背景：同一份 sol 293 数据，用 Python 3.14 跑出「碳酸盐 48 点」，
用 Python 3.12 跑出「碳酸盐 3 点」。如果差异只是数值条件问题，
那么把横轴先中心化 + 归一化再拟合，两个解释器就应该给出同样的结果。

本脚本用当前解释器做两件事：
  1. 报告范德蒙矩阵的条件数（原始 x vs 缩放后的 x）
  2. 用两种拟合方式各跑一遍逐点检测，比较「有合格波段的点数」
"""
from __future__ import annotations

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout

import os
import sys

import numpy as np

RRS = (_BootPath(__file__).resolve().parents[3]
       / "03_下载数据/SHERLOC/练习数据/processed/sol_00293"
       / "ss__0293_0692979640_115rrs__0090000srlc11360w108cgzj03.csv")

LASER_NM = 248.6
PIXEL_LO, PIXEL_HI = 50, 2098
NOISE = (2500.0, 3000.0)
FINGER = (800.0, 1500.0)
SNR_MIN = 3.0
MIN_FWHM = 30.0


def read_rrs(path):
    """按区块标题解析 RRS，不假设每个 Region 有多少个点。

    sol 181 的 RRS 每个 Region 100 行，sol 293 每段只有 50 行——
    写死行号会读错，所以这里按标题切块。
    """
    lines = open(path, encoding="utf-8", errors="replace").read().splitlines()

    def nums(line):
        return [float(v) for v in line.split(",") if v.strip()]

    wl = None
    for ln in lines[:6]:
        try:
            vals = nums(ln)
        except ValueError:
            continue
        if len(vals) > 100:
            wl = np.array(vals)
            break

    regions = []
    i = 0
    while i < len(lines):
        if "SPECTRA:_REGION" in lines[i] or "SPECTRA: _REGION" in lines[i]:
            i += 2                                     # 跳过通道表头
            rows = []
            while i < len(lines):
                try:
                    row = nums(lines[i])
                except ValueError:
                    break
                if not row:
                    break
                rows.append(row)
                i += 1
            regions.append(np.array(rows))
            continue
        i += 1
    return wl, regions


def baseline_raw(x, y, degree=5, n_iter=10, clip=1.5):
    """现在工具里的做法：直接对原始波数拟合。"""
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


def baseline_scaled(x, y, degree=5, n_iter=10, clip=1.5):
    """先把 x 中心化、归一化到 [-1,1]，再拟合同阶多项式。"""
    mid = 0.5 * (x[0] + x[-1])
    half = 0.5 * (x[-1] - x[0]) or 1.0
    xs = (x - mid) / half
    keep = np.ones_like(y, dtype=bool)
    coef = np.polyfit(xs[keep], y[keep], degree)
    for _ in range(n_iter):
        resid = y - np.polyval(coef, xs)
        sigma = resid[keep].std()
        keep = resid < clip * sigma
        if keep.sum() < degree + 2:
            break
        coef = np.polyfit(xs[keep], y[keep], degree)
    return np.polyval(coef, xs)


def find_bands(x, y, sigma):
    out = []
    n = len(y)
    for j in range(1, n - 1):
        if not (y[j] > y[j - 1] and y[j] >= y[j + 1]):
            continue
        if sigma <= 0 or y[j] < SNR_MIN * sigma:
            continue
        half = y[j] / 2.0
        k = m = j
        while k > 0 and y[k] > half:
            k -= 1
        while m < n - 1 and y[m] > half:
            m += 1
        fwhm = x[m] - x[k]
        if MIN_FWHM <= fwhm <= 200.0:
            out.append({"pos": float(x[j]), "fwhm": float(fwhm),
                        "snr": float(y[j] / sigma)})
    return out


def run(x_full, y_full, baseline):
    """返回 (有波段的点数, 波段总数, 噪声中位数)。"""
    wide = (x_full >= FINGER[0]) & (x_full <= 4000.0)
    xw = x_full[wide]
    quiet = (xw >= NOISE[0]) & (xw <= NOISE[1])
    in_win = (xw >= FINGER[0]) & (xw <= FINGER[1])
    n = y_full.shape[0]
    sig = np.empty(n)
    n_pts, n_bands = 0, 0
    for i in range(n):
        ys = y_full[i, wide] - baseline(xw, y_full[i, wide])
        s = float(np.sqrt(np.mean(ys[quiet] ** 2)))
        sig[i] = s
        bands = [b for b in find_bands(xw, ys, s) if in_win[np.searchsorted(xw, b["pos"])]]
        if bands:
            n_pts += 1
            n_bands += len(bands)
    return n_pts, n_bands, float(np.median(sig))


def main():
    wl, regs = read_rrs(RRS)
    if wl is None or len(regs) < 3:
        print("解析失败：波长 %s，Region 数 %d" % (wl is not None, len(regs)))
        return 1
    r1, r2, r3 = regs[0], regs[1], regs[2]
    shift = 1.0e7 / LASER_NM - 1.0e7 / wl
    order = np.argsort(shift)
    shift = shift[order]
    total = (r1 + r2 + r3)[:, order]
    sl = slice(PIXEL_LO, PIXEL_HI)
    x_full, y_full = shift[sl], total[:, sl]

    print("解释器 numpy:", np.__version__)
    print("点数 %d，波数范围 %.0f - %.0f" % (y_full.shape[0], x_full[0], x_full[-1]))

    # 1. 条件数
    for tag, xs in (("原始 x", x_full),
                    ("缩放 x", (x_full - x_full.mean()) / np.ptp(x_full))):
        v = np.vander(xs, 6)
        v = v / np.linalg.norm(v, axis=0)
        print("  范德蒙条件数（%s）: %.3e" % (tag, np.linalg.cond(v)))

    # 2. 两种拟合方式
    for tag, fn in (("原始 x 拟合（现状）", baseline_raw),
                    ("缩放 x 拟合", baseline_scaled)):
        n_pts, n_bands, sig = run(x_full, y_full, fn)
        print("  %-18s -> 有波段 %3d 点 / %3d 个波段，噪声中位 %.0f 计数"
              % (tag, n_pts, n_bands, sig))

    # 3. 同一根谱线上，两种基线差多少
    wide = (x_full >= FINGER[0]) & (x_full <= 4000.0)
    xw = x_full[wide]
    i = 20
    b1, b2 = baseline_raw(xw, y_full[i, wide]), baseline_scaled(xw, y_full[i, wide])
    d = np.abs(b1 - b2)
    print("  第 %d 点：两种基线最大差 %.1f 计数（该段信号幅度约 %.1f 计数）"
          % (i, d.max(), np.ptp(y_full[i, wide])))
    return 0


if __name__ == "__main__":
    sys.exit(main())
