# -*- coding: utf-8 -*-
"""PIXL SFS 原始 XRF 谱的自研读取与峰面积层（Phase 3；全自有实现）。

自算/对账边界声明：
  管线输入 = PDS data_raw_spectroscopy 的 SFS csv（逐 PMC 4096 通道 XRF 计数）；
  自研科学变换 = CSV 解析 → PMC 谱矩阵 → 能量轴自证（两点定标：钙富扫描的
               最强峰=Ca Kα 3.691 keV、次强=Fe Kα 6.404 keV）→ 元素 Kα 窗
               净峰面积（线性背景）+误差+检出判定；
  对账基准 = PIQUANT（官方开源，可选交叉验证）与官方 RQC wt%——只在台账。
  能量轴自证结果写入 energy_axis.json 供全部扫描复用；钙/铁身份假设仅在
  CaO 显著高的扫描上启用，且轴参数跨扫描一致性要校验。

不做：基本参数浓度反演（PIQUANT 覆盖）；谱去卷积（重叠峰解析留待
能量轴锚定后另立项）。
"""
from __future__ import annotations

import csv
import json
import os
import re

import numpy as np

# 元素 Kα 线（keV）——PIXL RQC 元素集（Na-Br）+ 主量
# P 2.0144（2026-11-14 增补，文献多点标定需要 P2O5；纯增量，不影响既有窗口）
KA_KEV = {"Na": 1.041, "Mg": 1.254, "Al": 1.487, "Si": 1.740, "P": 2.014,
          "S": 2.308,
          "Cl": 2.622, "K": 3.314, "Ca": 3.691, "Ti": 4.511, "Cr": 5.415,
          "Mn": 5.899, "Fe": 6.404, "Ni": 7.478, "Zn": 8.639, "Sr": 14.165,
          "Br": 13.291}
CAL_PAIR = (("Ca", 3.691), ("Fe", 6.404))   # 两点定标元素（按序）
N_CHANNELS = 4096


def read_sfs(path):
    """读一份 SFS csv：返回 (谱矩阵 [行 x 4096], 行元数据列表)。

    行结构（官方 SIS）：前若干列为遥测（SCLK/PMC/live_time...），其后
    4096 个通道计数。列位置由表头定位，不猜固定列号。
    """
    lines = open(path, "r", encoding="utf-8", errors="replace").read().splitlines()
    if not lines:
        raise ValueError("空 SFS 文件：%s" % path)
    header = [h.strip() for h in lines[0].split(",")]
    ch0 = None
    for i, h in enumerate(header):
        if re.fullmatch(r"c_?0+|ch0|channel[_ ]?0+", h, re.I) or (
                i > 0 and h == "" and len(header) - i == N_CHANNELS + 1):
            ch0 = i
            break
    if ch0 is None:
        # 兜底：首个"一行 4096 通道"的行结构由数据行长度推断
        data_row = lines[1].split(",")
        ch0 = len(data_row) - N_CHANNELS
        if ch0 <= 0:
            return None, None          # 非谱表（如遥测/配置行），调用方跳过
    spectra, meta = [], []
    for line in lines[1:]:
        parts = [p.strip() for p in line.split(",")]
        if len(parts) < ch0 + N_CHANNELS:
            continue
        spec = np.asarray(parts[ch0:ch0 + N_CHANNELS], dtype=float)
        if not np.isfinite(spec).all():
            continue
        spectra.append(spec)
        pm = {}
        for key in ("SCLK", "PMC", "live_time", "real_time"):
            if key in header:
                j = header.index(key)
                pm[key] = parts[j] if j < len(parts) else ""
        meta.append(pm)
    if not spectra:
        raise ValueError("SFS 无有效谱行")
    return np.vstack(spectra), meta


def find_peak_channels(spec, n=12, min_gap=20):
    """谱的局部极大通道清单（按计数降序，最小间隔 min_gap 通道）。"""
    s = np.convolve(spec, np.ones(5) / 5, mode="same")
    order = np.argsort(s)[::-1]
    picked = []
    for ch in order:
        if s[ch] <= 0:
            break
        if all(abs(ch - p) > min_gap for p in picked):
            picked.append(int(ch))
        if len(picked) >= n:
            break
    return sorted(picked)


def calibrate_energy(mean_spec, gain_lo=0.0030, gain_hi=0.0100):
    """能量轴自证 v4（多点线系两层；v3 单点启发保留为平票兜底）。

    第一层（v4 多点线系）：枚举强峰对 (p1<p2) 解两点 (g,b)，物理约束
    g∈[0.0030,0.0100] keV/ch（v4：真轴 ≈0.007884——2026-11-14 发现旧
    [0.0035,0.0070] 域漏掉真解；全量程 ≈32 keV/4096 ch）、b 近 0
    （|b|≤0.15）。每个候选按"多点线系支持"打分：KA_KEV 主线（Si/S/Cl/K/
    Ca/Ti/Cr/Mn/Fe）期望通道 ±5 ch 内存在实测峰的条数；取支持数最高者，
    平票时回退 v3 的 Ca/Fe 计数和启发。第二层（Fe 峰弱时兜底）：对每个
    强峰 p1 假设 p1=Ca Kα（g=3.691/p1，b=0），验证 Fe Kα 通道
    (6.404/g)±15 处存在显著峰（计数 > 全谱中位数的 1.5 倍）。
    无解抛 ValueError。

    修订记录（冻结纪律）：
      v4 2026-11-14 —— Máaz 四研磨斑 sum 谱 10 线最小二乘发现真轴
      g≈0.007884/b≈-0.028（线残差 ≤0.015 keV）；旧两点法把最强峰（富硅
      斑的 Si）误当 Ca，且旧增益域漏掉真解，冻结轴窗实际落在 Kβ/肩部
      （仍随元素浓度涨落但相关性退化：Ca 0.9712→0.9902、Fe 0.9624→0.9776）。
    """
    peaks = find_peak_channels(mean_spec, n=30)
    if len(peaks) < 2:
        raise ValueError("可辨识峰不足 2 个，无法能量自证")
    s = np.convolve(mean_spec, np.ones(5) / 5, mode="same")
    e1, e2 = CAL_PAIR[0][1], CAL_PAIR[1][1]
    # 线匹配只认强峰（≥3×中位数）：噪声假匹配会抬高错误候选的线系支持数
    med = float(np.median(s[s > 0])) if (s > 0).any() else 0.0
    strong = [p for p in peaks if s[p] >= 3 * med] or list(peaks)
    cands = []
    for i, p1 in enumerate(peaks):
        for p2 in peaks[i + 1:]:
            if p2 <= p1:
                continue
            g = (e2 - e1) / (p2 - p1)
            b = e1 - g * p1
            if gain_lo <= g <= gain_hi and abs(b) <= 0.15:
                # 多点线系支持：主线期望通道 ±5 ch 内有实测强峰的条数
                support = 0
                for el, kev in KA_KEV.items():
                    if el in ("Na", "Mg", "Sr", "Br"):
                        continue        # 低信噪/近满量程线不参与评分
                    ch_exp = (kev - b) / g
                    if any(abs(ch_exp - p) <= 5 for p in strong):
                        support += 1
                cands.append((support, s[p1] + s[p2], p1, p2, g, b))
    if cands:
        cands.sort(reverse=True)
        support, _, p1, p2, g, b = cands[0]
        return {"keV_per_ch": float(g), "zero_keV": float(b),
                "peaks": [[int(p1), e1], [int(p2), e2]],
                "mode": "两点" if support <= 2 else "两点+%d线支持" % support}
    # 第二层：单点 Ca 锚 + Fe 窗验证
    med = float(np.median(s[s > 0])) if (s > 0).any() else 0.0
    for p1 in peaks:
        g = e1 / p1
        if not (gain_lo <= g <= gain_hi):
            continue
        fe_ch = int(round(e2 / g))
        win = s[max(0, fe_ch - 15):fe_ch + 16]
        if len(win) and win.max() > 1.5 * med:
            return {"keV_per_ch": float(g), "zero_keV": 0.0,
                    "peaks": [[int(p1), e1], [fe_ch, e2]],
                    "mode": "单点Ca锚+Fe窗验证"}
    raise ValueError("无可行的能量自证解（两层约束内）")


def ch_to_kev(ch, axis):
    return axis["keV_per_ch"] * np.asarray(ch, dtype=float) + axis["zero_keV"]


def peak_areas(spec, axis, elements=("Ca", "Fe", "Si", "Mg", "Al", "S", "Ti", "K")):
    """元素 Kα 窗净面积（半宽 ±0.06 keV，线性背景取窗两端）。

    返回 {元素: {area, error, channel, kev, detected}}；窗越界→detected=False
    且 area=None（不猜值）。error=面积的开方近似（泊松）。
    """
    x = np.asarray(spec, dtype=float)
    out = {}
    half = 0.06
    for el in elements:
        kev = KA_KEV[el]
        ch = (kev - axis["zero_keV"]) / axis["keV_per_ch"]
        lo, hi = int(ch - half / axis["keV_per_ch"]), int(ch + half / axis["keV_per_ch"])
        if lo < 1 or hi >= len(x) - 1:
            out[el] = {"area": None, "error": None, "channel": int(ch),
                       "kev": kev, "detected": False,
                       "note": "窗越界"}
            continue
        win = x[lo:hi + 1]
        bg = float(np.median(np.concatenate([x[max(0, lo - 8):lo],
                                             x[hi + 1:hi + 9]]))) if lo > 8 and hi < len(x) - 9 \
            else float(np.median(win[:2].tolist() + win[-2:].tolist()))
        net = win - bg
        net[net < 0] = 0.0
        area = float(net.sum())
        out[el] = {"area": area, "error": float(np.sqrt(max(area, 0.0))),
                   "channel": int(ch), "kev": kev,
                   "detected": bool(area > 5 * np.sqrt(max(area, 1.0)))}
    return out
