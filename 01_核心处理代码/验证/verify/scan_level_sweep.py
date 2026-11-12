# -*- coding: utf-8 -*-
"""给"扫描级兜底"定阈值：把每个扫描的两个量都量出来，看真实的检出与已知的零能不能分开。

扫描级兜底要回答的问题是：**单个点的信噪比撑不起逐点判据，但整片扫描合起来
是不是有带？** 两个量：

  A      —— 全扫描平均谱上最显著的带，相对平均谱自身噪声的倍数（σ_mean）
  比值 r —— 同一个带的幅度 ÷ 单点噪声 σ_pt

r 是关键的护栏：平均谱的噪声按 1/√N 缩小，1296 点的扫描里哪怕只有 7 个计数
的起伏也够 4σ，光看 A 会把"论文说是零"的扫描也报成有矿物。r 要求这条带
至少达到单点噪声的某个比例，才有资格叫"弱信号检出"。

用法：python verify\\scan_level_sweep.py
"""

from __future__ import annotations

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = str(_layout.project_layout.APP)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from sherloc import minerals                       # noqa: E402
from sherloc.organize import scan_folder           # noqa: E402

DATA = os.path.join(os.path.dirname(ROOT), "03_下载数据/SHERLOC/练习数据")

# 已知的"论文说没有矿物"的对照，兜底绝不能在它们身上报出新类别
KNOWN_ZERO = {("141", "srlc10203"), ("181", "srlc15050"), ("349", "srlc11420")}


def scan_metrics(x, spectra):
    """返回 (A, sigma_mean, r, 峰位, FWHM, n_points, sigma_pt)。"""
    sel = (x >= 2500) & (x <= 3000)
    # 窗口跟着内核走（2026-09-15 扩到 800–1800），不要再写死 1500——
    # 阈值就是照这个窗口标定的，两处不一致的话标定等于白做。
    band = (x >= minerals.FINGERPRINT[0]) & (x <= minerals.FINGERPRINT[1])
    xs = x[band]
    n = spectra.shape[0]
    acc = np.zeros_like(x)
    sig = np.empty(n)
    for i in range(n):
        y = spectra[i] - minerals.baseline(x, spectra[i])
        acc += y
        sig[i] = np.sqrt(np.mean(y[sel] ** 2))
    sig = np.nan_to_num(sig)
    sig_med = float(np.median(sig))

    mean_all = acc / max(n, 1)           # 全谱长度
    # 与内核保持一致：掩峰迭代基线（普通多项式在 800–1800 上会放大两端的伪影）
    resid = (mean_all - minerals.baseline_masked(x, mean_all))[band]
    # 平均谱自身的噪声：用稳健尺度（落在带上的点只占窗口一小部分，不影响 MAD）
    sig_mean = 1.4826 * float(np.median(np.abs(resid - np.median(resid))))
    bands = minerals.find_bands(xs, resid, sig_mean, min_fwhm=30.0,
                                max_fwhm=200.0, snr_min=3.0)
    if not bands:
        return 0.0, sig_mean, 0.0, float("nan"), float("nan"), sig_med
    b = max(bands, key=lambda d: d["snr"])
    return (float(b["height"]), sig_mean, float(b["height"] / max(sig_med, 1e-9)),
            float(b["pos"]), float(b["fwhm"]), sig_med)


def main():
    _p, observations = scan_folder(DATA)
    rows = []
    for ob in observations:
        if ob.spectra is None:
            continue
        wl, r1, _r2, _r3 = ob.read_regions()
        if wl is None or r1 is None:
            continue
        x = 1.0e7 / minerals.__dict__.get("LASER_NM", 248.6) - 1.0e7 / wl
        # 与工具保持一致：只用 Region 1、按波数排序
        o = np.argsort(x)
        x, r1 = x[o], r1[:, o]
        A, sm, r, pos, fw, spt = scan_metrics(x, r1)
        rows.append((r, A, sm, pos, fw, spt, ob))
    rows.sort(key=lambda t: -t[0])
    print("%-6s %-9s %-16s %6s %7s %8s %9s %8s  %s"
          % ("sol", "seq", "sclk", "r=A/σpt", "A", "σ_mean", "A/σ_mean", "峰位",
             "FWHM"))
    for r, A, sm, pos, fw, spt, ob in rows:
        tag = "  ← 论文说没有矿物" if (ob.sol, ob.sequence_id) in KNOWN_ZERO else ""
        print("%-6s %-9s %-16s %6.2f %7.1f %8.2f %9.1f %8.0f  %5.0f%s"
              % (ob.sol, ob.sequence_id, ob.sclk[:16], r, A, sm,
                 A / sm if sm else float("nan"), pos, fw, tag))
    return 0


if __name__ == "__main__":
    sys.exit(main())
