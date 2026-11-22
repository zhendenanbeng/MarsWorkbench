# -*- coding: utf-8 -*-
"""APXS 基本参数法-lite（FP-lite）：从原始 rsp 能谱自算氧化物 wt%。

对 03 号脚本标定结果的运行期引擎；物理模型（吸收修正 Sherman-lite）：
  厚靶产额  area_i ∝ C_i · u_i / (μ/ρ)_m(E_i)
  - C_i：元素质量分数（基体含化学计量 O）；
  - u_i：逐元素有效灵敏度——α 激发截面、荧光产额、探测器效率、Pu L XRF
    项折入，由 25 apb + 40 apa 自然靶标定（数据准备脚本 03）；
  - (μ/ρ)_m(E_i)：基体在元素 i Kα1 能量处的总质量衰减系数
    （xraylib CS_Total 权重和，17 元素 × 14 能量冻结表）。
求解：C_i ∝ area_i·(μ/ρ)_m^(t)(E_i)/u_i → 氧化物口径归一 Σ=100 wt% →
迭代到收敛。每靶几何/注量因子 k_t 被归一化消去——与官方 rwp 口径一致
（官方亦为归一化 wt%）。

简化与边界（诚实声明）：
  1. XRF（子体 240Pu L 线 14.28/18.30/21.42 keV）的独立入射衰减依赖折入
     u_i：对 Na-Fe 主量元素 μ_m(E_i) ≫ μ_m(E_L)，两项基体依赖形状相同；
     Zn/Ni/Br 量级能区除外（痕量元素，误差并入留出验证）。
  2. 次级荧光增强（如 Ti Kα→Ca）未建模：τ_Ca(4.5 keV)·C_Ti 量级，主量
     影响 <1%。
  3. 探测器效率随任务时间辐照退化：u_i 为标定时代（sol 629-4327 apb +
     apa 靶）时代均值；跨时代外推受限（rocknest3_rp 边界外先例见 02 号
     脚本台账，负结果留档）。
  4. α 能量依赖 σ_α(E) 在逃逸深度内近似常数：χ_i·R_α ≫ 1 对全部主量成立
     （Si Kα χR≈60，Fe Kα χR≈7）。
  5. P 无自研峰位（P Kα 2.010 未拟合），预测模式基体不含 P——P2O5 典型
     <1 wt%，对 μ_m 影响有界，留出验证覆盖。
"""
from __future__ import annotations

import json
from pathlib import Path

# 标准原子量（IUPAC 2021；氧化物换算系数由此计算，与官方 rwp 氧化物口径一致）
_AT = {"O": 15.999, "Na": 22.98977, "Mg": 24.305, "Al": 26.98154,
       "Si": 28.085, "P": 30.97376, "S": 32.06, "Cl": 35.45, "K": 39.0983,
       "Ca": 40.078, "Ti": 47.867, "Cr": 51.9961, "Mn": 54.938,
       "Fe": 55.845, "Ni": 58.6934, "Zn": 65.38, "Br": 79.904}


# 元素 → (官方 rwp 氧化物名, 元素在氧化物中的质量分数)；Cl/Ni/Zn/Br 为元素行
OXIDE_OF = {
    "Na": ("Na2O", 2 * _AT["Na"] / (2 * _AT["Na"] + _AT["O"])),
    "Mg": ("MgO", _AT["Mg"] / (_AT["Mg"] + _AT["O"])),
    "Al": ("Al2O3", 2 * _AT["Al"] / (2 * _AT["Al"] + 3 * _AT["O"])),
    "Si": ("SiO2", _AT["Si"] / (_AT["Si"] + 2 * _AT["O"])),
    "P": ("P2O5", 2 * _AT["P"] / (2 * _AT["P"] + 5 * _AT["O"])),
    "S": ("SO3", _AT["S"] / (_AT["S"] + 3 * _AT["O"])),
    "Cl": ("Cl", 1.0),
    "K": ("K2O", 2 * _AT["K"] / (2 * _AT["K"] + _AT["O"])),
    "Ca": ("CaO", _AT["Ca"] / (_AT["Ca"] + _AT["O"])),
    "Ti": ("TiO2", _AT["Ti"] / (_AT["Ti"] + 2 * _AT["O"])),
    "Cr": ("Cr2O3", 2 * _AT["Cr"] / (2 * _AT["Cr"] + 3 * _AT["O"])),
    "Mn": ("MnO", _AT["Mn"] / (_AT["Mn"] + _AT["O"])),
    "Fe": ("FeO", _AT["Fe"] / (_AT["Fe"] + _AT["O"])),
    "Ni": ("Ni", 1.0),
    "Zn": ("Zn", 1.0),
    "Br": ("Br", 1.0),
}
OXIDE_ORDER = ("SiO2", "TiO2", "Al2O3", "FeO", "MgO", "CaO", "Na2O", "K2O",
               "SO3", "Cl", "P2O5", "Cr2O3", "MnO", "Ni", "Zn", "Br")

# 预测模式可自算元素（apxs_fit.KA_KEV 全集；P 不在内）
FP_ELEMENTS = ("Na", "Mg", "Al", "Si", "S", "Cl", "K", "Ca", "Ti",
               "Cr", "Mn", "Fe", "Ni", "Zn")

_TABLES = None
_CALIB = None


def _tables():
    global _TABLES
    if _TABLES is None:
        p = Path(__file__).resolve().parent / "apxs_fp_tables.json"
        _TABLES = json.loads(p.read_text(encoding="utf-8"))
    return _TABLES


def load_calibration():
    """冻结标定（u_i + 留出界）→ dict；不在位返回 None。"""
    global _CALIB
    if _CALIB is None:
        p = Path(__file__).resolve().parent / "apxs_fp_calibration.json"
        if not p.is_file():
            return None
        _CALIB = json.loads(p.read_text(encoding="utf-8"))
    return _CALIB


def matrix_mass_atten(comp, e_kev):
    """基体总质量衰减系数 (μ/ρ)_m(E)（cm²/g）。comp: {元素: 质量分数}。

    e_kev 必须在冻结能量网格内（14 个 Kα1 能量）；网格外取最近邻——
    只允许 O/P/S/Cl 等基体行使用，待拟合元素能量恒在网格内。
    """
    t = _tables()
    energies = t["energies_kev"]
    idx = min(range(len(energies)), key=lambda i: abs(energies[i] - e_kev))
    total, wsum = 0.0, 0.0
    for el, w in comp.items():
        if el not in t["cs_total"]:
            raise KeyError("基体元素 %s 不在物理表（17 元素集）" % el)
        w = float(w)
        if w < 0:
            raise ValueError("基体负质量分数 %s=%g" % (el, w))
        total += w * t["cs_total"][el][idx]
        wsum += w
    if wsum <= 0:
        raise ValueError("基体为空")
    return total / wsum


def elements_from_oxides(oxides):
    """官方 rwp 口径氧化物 wt → ({元素: wt}, O wt)。仅用于标定/对账。"""
    els, o_total = {}, 0.0
    for el, (ox, f) in OXIDE_OF.items():
        v = oxides.get(ox)
        if v is None:
            continue
        el_wt = v * f
        els[el] = el_wt
        o_total += v - el_wt
    return els, o_total


def oxides_from_elements(el_wt):
    """元素 wt（Σ 按 100 氧化物口径归一）→ {氧化物名: wt}。

    内部迭代用元素键版本见 solve()；本函数是公开输出口径。
    """
    fsum = 0.0
    for el, w in el_wt.items():
        if el not in OXIDE_OF:
            raise KeyError("未知元素 %s" % el)
        fsum += (w / OXIDE_OF[el][1]) if w > 0 else 0.0
    if fsum <= 0:
        raise ValueError("元素质量全为空")
    return {OXIDE_OF[el][0]: 100.0 * w / OXIDE_OF[el][1] / fsum
            for el, w in el_wt.items() if w > 0}


def _normalize_oxides(el_wt):
    """元素 wt → {元素: 氧化物 wt}（元素键内部口径，Σ=100）。"""
    fsum = 0.0
    for el, w in el_wt.items():
        fsum += (w / OXIDE_OF[el][1]) if w > 0 else 0.0
    if fsum <= 0:
        raise ValueError("元素质量全为空")
    return {el: 100.0 * w / OXIDE_OF[el][1] / fsum
            for el, w in el_wt.items() if w > 0}


def correct_kbeta(areas, half_width=0.12):
    """Kβ 串扰剥离：污染源 Kβ 落入目标 Kα ±half_width 窗内则减除。

    剥离比例 = Kβ/Kα 分支比（xraylib RadRate 冻结表）；污染源按 Z 序
    处理，保证其自身 Kα 窗先被更高 Z 的串扰修正过（K→Ca、Cr→Mn→Fe
    链）。负面积保留（检出判定的真实信号），由 solve 的 active 过滤。
    """
    kb = _tables()["kbeta"]
    ka = _tables()["line_energies_kev"]
    corrected = dict(areas)
    for c in sorted(kb, key=lambda el: kb[el]["energy_kev"]):
        r_c, e_kb = kb[c]["ratio"], kb[c]["energy_kev"]
        if (corrected.get(c) or 0) <= 0:
            continue
        for t, e_t in ka.items():
            if t == c or abs(e_kb - e_t) > half_width:
                continue
            if corrected.get(t) is not None and c in corrected:
                corrected[t] = corrected[t] - r_c * corrected[c]
    return corrected


def solve(areas, sens, max_iter=30, tol=1e-6):
    """Sherman-lite 迭代求解氧化物 wt%。

    areas: {元素: cps 谱净面积}（apxs_fit.spectrum_areas 的 area 字段）；
    sens: {元素: u_i}（load_calibration()["sensitivities"]）。
    返回 {"available", "reason"/"oxides", "element_wt", "o_wt",
          "iterations", "converged", "mu_by_el", "active"}；
    oxides 以官方 rwp 氧化物名为键。
    """
    active = [el for el in FP_ELEMENTS
              if (areas.get(el) or 0) > 0 and (sens.get(el) or 0) > 0]
    if "Si" not in active:
        return {"available": False,
                "reason": "Si 峰未检出或灵敏度不在册——FP 求解拒绝（不硬算）"}
    if len(active) < 5:
        return {"available": False,
                "reason": "可拟合元素不足 5 个（%d 个）——FP 求解拒绝" % len(active)}
    # Kβ 串扰剥离后再进 Sherman-lite（标定 u 与预测同口径）
    areas = correct_kbeta(areas)
    active = [el for el in active if areas.get(el, 0) > 0]
    if "Si" not in active:
        return {"available": False,
                "reason": "Kβ 剥离后 Si 净面积≤0——FP 求解拒绝（不硬算）"}
    # 初始：无吸收修正估计（C ∝ area/u）
    ox = _normalize_oxides({el: areas[el] / sens[el] for el in active})
    mu_by_el, converged, it = {}, False, 0
    for it in range(1, max_iter + 1):
        el_wt = {el: ox[el] * OXIDE_OF[el][1] for el in ox}
        o_wt = 100.0 - sum(el_wt.values())          # 化学计量 O（100 基准）
        comp = dict(el_wt)
        comp["O"] = o_wt
        wsum = sum(comp.values())
        comp = {el: v / wsum for el, v in comp.items()}
        mu_by_el = {el: matrix_mass_atten(comp, _tables()["line_energies_kev"][el])
                    for el in active}
        ox_new = _normalize_oxides({el: areas[el] * mu_by_el[el] / sens[el]
                                    for el in active})
        delta = max(abs(ox_new.get(el, 0.0) - ox.get(el, 0.0))
                    for el in set(ox) | set(ox_new))
        ox = ox_new
        if delta < tol:
            converged = True
            break
    el_wt = {el: ox[el] * OXIDE_OF[el][1] for el in ox}
    return {"available": True,
            "oxides": {OXIDE_OF[el][0]: round(v, 3) for el, v in ox.items()},
            "element_wt": {k: round(v, 3) for k, v in el_wt.items()},
            "o_wt": round(100.0 - sum(el_wt.values()), 3),
            "iterations": it, "converged": converged,
            "mu_by_el": {k: round(v, 2) for k, v in mu_by_el.items()},
            "active": active}


def predict(areas, calib=None):
    """工作台入口：净面积 dict → FP 预测 + 对账行（官方 rwp 在位时）。"""
    calib = calib or load_calibration()
    if not calib:
        return {"available": False,
                "reason": "FP 标定不在册（apxs_fp_calibration.json）"}
    res = solve(areas, calib.get("sensitivities") or {})
    if not res.get("available"):
        return res
    res["scope"] = calib.get("scope", "")
    res["calib_era"] = calib.get("calib_era", "")
    res["holdout_bounds"] = calib.get("holdout", {})
    return res
