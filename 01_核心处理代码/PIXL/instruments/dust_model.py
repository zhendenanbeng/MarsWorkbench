# -*- coding: utf-8 -*-
"""PIXL 火星尘厚度与成分估计（VanBommel 2026 方法的 Python 移植）。

来源：WUSTL Library DOI 10.7936/6rxs-108366 的 PIXL_CT_dust_thickness.m
（VanBommel, S. J., 用于 PNAS 投稿 "Fine and Ultrafine Martian Dust
Present Additional Risks for Future Crewed Missions to Mars"，
MATLAB R2024a，2026-06-10）。本模块逐式移植其向量化的混合模型与
蒙特卡洛成分搜索，不冒充官方代码——口径与公式按 .m 原文对照。

物理模型（两层混合，逐元素）：
    测量谱 = LF(T) × 尘成分信号 + (1 − LF(T)) × 底材信号
其中 LF(T) 是预计算的采样分数表（PIXL_Sampling_Fractions.csv）：
质量厚度 T 时信号中来自尘层的比例（逐元素不同——能量越低尘层
衰减/贡献越强）。表值由官方蒙特卡洛计算给出，本模块只读取。

约束途径（与 .m 一致）：
  * PTFE 定标靶（C/F 组成，与 8 元素尘无重叠）→ 测到的 8 元素全来自
    尘层：model = LF × 尘信号（归一到 100），配 T 先验（3e-4 g/cm² ± 2e-4）
  * BHVO-2 定标靶（已知成分玄武岩玻璃）→ 只用 P、S 两元素（岩石几乎
    不含、尘富含——最灵敏的尘信号）：model = LF×尘 + (1−LF)×BHVO
  * 8 元素尘成分由蒙特卡洛逐步增减扰动搜索（保留 χ² 改善的步）

支持数据（仪器灵敏度、PTFE/BHVO 飞行计数、初始成分猜测）按 .m 的
loadData() 原样转录——官方将这些值硬编码在脚本内以方便使用。
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

ELEMENTS_DUST = ("Mg", "Al", "Si", "P", "S", "Ca", "Ti", "Fe")
PRIOR_T = 3e-4          # g/cm²，.m: prior_T = 3e-4
PRIOR_SIGMA = 2e-4      # .m: sigma_T = 2e-4
RHO = 2.0               # g/cm³，.m: rho = 2.0
T_TO_UM = lambda t: (t / RHO) * 1e4          # noqa: E731

# loadData() 的硬编码支持数据（.m 原文转录）
SENS_TEXT = ("Na", "Mg", "Al", "Si", "P", "S", "Cl", "K", "Ca", "Ti", "Cr", "Mn", "Fe")
SENS_DATA = (18.4, 112.6, 313.6, 687.0, 910.2, 1450.0, 2800.0, 2405.9,
             2800.0, 3292.4, 3780.0, 4000.0, 4297.3)          # cps/wt%
TEFLON_TEXT = ELEMENTS_DUST
TEFLON_DATA = ((1.89, 0.31), (5.97, 0.12), (45.96, 0.10), (2.68, 0.22),
               (16.65, 0.16), (7.39, 0.24), (1.47, 0.02), (17.99, 0.21))  # (cps, err)
BLAB_TEXT = SENS_TEXT
BLAB_DATA = (0.020, 0.473, 2.537, 19.862, 0.165, 0.022, 0.009, 0.928,
             24.061, 6.774, 0.237, 0.737, 44.175)              # 无尘 BHVO 相对计数
BFLIGHT_TEXT = SENS_TEXT
BFLIGHT_DATA = (0.029, 0.462, 2.428, 19.585, 0.211, 0.312, 0.237, 0.917,
                23.725, 6.719, 0.241, 0.737, 44.398)           # 含尘 BHVO 相对计数
START_TEXT = SENS_TEXT
START_DATA = (2.43, 7.77, 9.71, 43.69, 0.97, 7.77, 0.97, 0.49,
              6.80, 0.97, 0.49, 0.49, 17.48)                   # 初始成分猜测 wt%


def norm100(v):
    v = np.asarray(v, dtype=float)
    return 100.0 * v / v.sum()


def load_support():
    """loadData() 移植：返回 8 元素对齐的 input 字典（与 .m 同字段名）。"""
    cps_all = norm100(TEFLON_DATA and [r[0] for r in TEFLON_DATA])
    err_scale = 100.0 / sum(r[0] for r in TEFLON_DATA)
    sens_map = dict(zip(SENS_TEXT, SENS_DATA))
    blab = dict(zip(BLAB_TEXT, norm100(BLAB_DATA)))
    bflight = dict(zip(BFLIGHT_TEXT, norm100(BFLIGHT_DATA)))
    start = dict(zip(START_TEXT, norm100(START_DATA)))
    teflon_cps = dict(zip(TEFLON_TEXT, cps_all))
    teflon_err = {e: dict(zip(TEFLON_TEXT,
                  [r[1] * err_scale for r in TEFLON_DATA]))[e] for e in TEFLON_TEXT}
    inp = {}
    for e in ELEMENTS_DUST:
        inp[e] = {"Sens": sens_map.get(e, 0.0),
                  "BL_cps": blab.get(e, 0.0),
                  "BF_cps": bflight.get(e, 0.0),
                  "Start_wt": start.get(e, 0.0),
                  "PTFE_cps": teflon_cps.get(e, 0.0),
                  "PTFE_cps_err": teflon_err.get(e, 0.0)}
    return inp


def load_fractions(csv_path):
    """读采样分数表：返回 (T 数组 [g/cm²], LF 矩阵 [8 元素 × nT])。

    与 .m 一致：跳过 1 行表头；过滤 T 为 NaN 或 < 1E-6 的行。
    """
    rows = []
    with open(csv_path, encoding="utf-8", errors="replace", newline="") as f:
        for i, row in enumerate(csv.reader(f)):
            if i == 0 or not row or not row[0].strip():
                continue
            try:
                rows.append([float(x) for x in row[:9]])
            except ValueError:
                continue
    if not rows:
        raise ValueError("采样分数表为空或不可解析：%s" % csv_path)
    arr = np.array(rows, dtype=float)
    t_raw, lf_raw = arr[:, 0], arr[:, 1:9].T
    ok = np.isfinite(t_raw) & (t_raw >= 1e-6)
    return t_raw[ok], lf_raw[:, ok]


def build_input(inp_e: dict, lf: np.ndarray):
    """把 8 元素 input 字典与 LF 表组装成 perform_cal 的输入。"""
    return {
        "LF": lf,
        "Sens": np.array([inp_e[e]["Sens"] for e in ELEMENTS_DUST]),
        "PTFE_cps": np.array([inp_e[e]["PTFE_cps"] for e in ELEMENTS_DUST]),
        "PTFE_cps_err": np.array([inp_e[e]["PTFE_cps_err"] for e in ELEMENTS_DUST]),
        "BL_cps": np.array([inp_e[e]["BL_cps"] for e in ELEMENTS_DUST]),
        "BF_cps": np.array([inp_e[e]["BF_cps"] for e in ELEMENTS_DUST]),
        "Start_wt": np.array([inp_e[e]["Start_wt"] for e in ELEMENTS_DUST]),
    }


def perform_cal(comp, inp, t_grid):
    """performCal() 移植：给定 8 元素尘成分，扫 T 网格算 χ²。

    返回 (PTFE_Chi2, PTFE_T, BHVO_Chi2, BHVO_T, totChi)——T 网格上
    χ² 最小处的取值（与 .m 输出顺序一致）。
    """
    lf, sens = inp["LF"], inp["Sens"]
    rel_c = np.asarray(comp, dtype=float) * sens
    rel_c = 100.0 * rel_c / rel_c.sum()

    # PTFE：测到的 8 元素全来自尘层 → model = LF × 尘信号（按列归一）
    model_ptfe = lf * rel_c[:, None]
    model_ptfe = 100.0 * model_ptfe / model_ptfe.sum(axis=0, keepdims=True)
    penalty = ((t_grid - PRIOR_T) ** 2) / (2.0 * PRIOR_SIGMA ** 2)
    weights = 1.0 / (inp["PTFE_cps_err"] ** 2 + np.finfo(float).eps)
    sq = (inp["PTFE_cps"][:, None] - model_ptfe) ** 2
    x_ptfe = (weights[:, None] * sq).sum(axis=0) / weights.sum() + 0.1 * penalty
    i_ptfe = int(np.argmin(x_ptfe))

    # BHVO：只看 P、S（尘最灵敏的两元素）
    # model = LF × 尘信号 + (1 − LF) × 无尘 BHVO 信号
    rs = [3, 4]                                   # P、S（0 基索引；.m 的 [4,5]）
    model_b = lf * rel_c[:, None] + (1.0 - lf) * inp["BL_cps"][:, None]
    sq_b = (inp["BF_cps"][rs][:, None] - model_b[rs, :]) ** 2
    x_b = sq_b.sum(axis=0) / len(rs)
    i_b = int(np.argmin(x_b))

    tot = float(np.sqrt(x_ptfe[i_ptfe] ** 2 + x_b[i_b] ** 2))
    return (float(x_ptfe[i_ptfe]), float(t_grid[i_ptfe]),
            float(x_b[i_b]), float(t_grid[i_b]), tot)


def mcmc(inp, t_grid, numiter=100, seed=None):
    """蒙特卡洛成分搜索（.m 主循环移植）。

    每步：随机选一元素试 +tweak（归一到 100），χ² 改善则保留；否则试
    −tweak，改善则保留。tweak 阶梯 [10 … 0.01]，每级 30 步。
    返回 (Summary[niter×5], Comp_Log[8×niter])，列序同 .m：
    Summary = [totChi, PTFE_Chi2, BHVO_Chi2, PTFE_T, BHVO_T]。
    """
    rng = np.random.default_rng(seed)
    tweaks = [10.0, 5.0, 2.5, 1.0, 0.5, 0.25, 0.1, 0.05, 0.01]
    numturns = 5 * len(ELEMENTS_DUST)
    comp = np.asarray(inp["Start_wt"], dtype=float)
    if not np.isfinite(comp).all() or comp.sum() <= 0:
        raise ValueError("Start_wt 全零或含非有限值——无法归一（检查 input 组装）")
    comp = 100.0 * comp / comp.sum()
    _, _, _, _, base_chi = perform_cal(comp, inp, t_grid)
    summary, comp_log = [], []
    for _ in range(numiter):
        for tweak in tweaks:
            for _ in range(numturns):
                re = int(rng.integers(0, len(comp)))
                test = comp.copy()
                test[re] += tweak
                test = 100.0 * test / test.sum()
                _, _, _, _, tc = perform_cal(test, inp, t_grid)
                if tc < base_chi:
                    comp, base_chi = test, tc
                    continue
                test = comp.copy()
                test[re] = max(comp[re] - tweak, 0.0)
                test = 100.0 * test / test.sum()
                _, _, _, _, tc = perform_cal(test, inp, t_grid)
                if tc < base_chi:
                    comp, base_chi = test, tc
        _, ptfe_chi, bhvo_chi, ptfe_t, bhvo_t = perform_cal(comp, inp, t_grid)
        summary.append([base_chi, ptfe_chi, bhvo_chi, ptfe_t, bhvo_t])
        comp_log.append(comp.copy())
    return np.array(summary), np.vstack(comp_log)


def summarize(summary, comp_log):
    """中位结果 + 95% CI + µm 换算（.m 输出段：Summary 列序
    [totChi, PTFE_Chi2, BHVO_Chi2, PTFE_T, BHVO_T]）。"""
    med = np.median(summary, axis=0)
    comp_med = np.median(comp_log, axis=0)
    comp_lo = np.percentile(comp_log, 2.5, axis=0)
    comp_hi = np.percentile(comp_log, 97.5, axis=0)
    return {
        "tot_chi2": float(med[0]),
        "ptfe_chi2": float(med[1]), "bhvo_chi2": float(med[2]),
        "ptfe_T_gcm2": float(med[3]), "ptfe_T_um": float(T_TO_UM(med[3])),
        "bhvo_T_gcm2": float(med[4]), "bhvo_T_um": float(T_TO_UM(med[4])),
        "dust_composition": {e: {"median_wt": float(comp_med[i]),
                                 "ci95": [float(comp_lo[i]), float(comp_hi[i])]}
                             for i, e in enumerate(ELEMENTS_DUST)},
    }
