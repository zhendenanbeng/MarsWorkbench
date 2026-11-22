#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""APXS 数据准备脚本 00：从 xraylib 提取基本参数法物理表 → 冻结 JSON。

构建期工具（需要 pip install xraylib；运行期 apxs_fp.py 只读冻结 JSON，
完全离线）。提取内容：
  * 17 个基体元素（O Na Mg Al Si P S Cl K Ca Ti Cr Mn Fe Ni Zn Br）在
    14 个待拟合元素 Kα1 能量处的总质量衰减系数 CS_Total（cm²/g，
    Chantler/Hell 数据，xraylib 4.3.0）；
  * 各元素 Kα1 线能量与 K 边（核对用）。

物理口径说明（详见 instruments/apxs_fp.py 模块头）：
  α 激发项与 Pu L XRF 项对主量元素（Na-Fe）的基体依赖同为出射吸收
  1/μ_m(E_i)，故引擎只需要 Kα 能量处的 μ 表；XRF 的独立入射衰减依赖
  仅在 Zn/Ni/Br 量级能区显著（痕量，忽略并入 u_i）。
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "instruments" / "apxs_fp_tables.json"

# 引擎待拟合元素（Kα 能量网格来源）
FIT_ELEMENTS = ("Na", "Mg", "Al", "Si", "S", "Cl", "K", "Ca", "Ti",
                "Cr", "Mn", "Fe", "Ni", "Zn")
# 基体元素（含化学计量 O 与官方表全部氧化物阳离子 + Cl/Br 阴离子）
MATRIX_ELEMENTS = ("O", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "K",
                   "Ca", "Ti", "Cr", "Mn", "Fe", "Ni", "Zn", "Br")

XRAYLIB_HINT = "构建期依赖 xraylib（pip install xraylib），运行期离线只读本表"


def main():
    try:
        import xraylib as x
    except ImportError:
        print(XRAYLIB_HINT)
        return 1
    doc = {
        "scope": "APXS 基本参数法-lite 物理表：CS_Total 总质量衰减系数（cm²/g）"
                 "@14 个待拟合元素 Kα1 能量；17 基体元素。来源 xraylib %s "
                 "(Chantler/Hell 数据)；Kα1 线能量与 K 边一并冻结核对。"
                 % x.__version__ if hasattr(x, "__version__") else "xraylib",
        "generated_by": "数据准备脚本/00_gen_fp_tables.py（构建期，xraylib 4.3.0）",
        "xraylib_version": getattr(x, "__version__", "unknown"),
        "fit_elements": list(FIT_ELEMENTS),
        "matrix_elements": list(MATRIX_ELEMENTS),
        "energies_kev": [],
        "line_energies_kev": {},
        "k_edges_kev": {},
        "cs_total": {},                     # 元素 → [14 个能量处的 cm²/g]
        "pu_l_xrf_kev": {                   # XRF 激发线（折入 u_i；文档口径）
            "La1": x.LineEnergy(94, x.LA1_LINE),
            "Lb1": x.LineEnergy(94, x.LB1_LINE),
            "Lg1": x.LineEnergy(94, x.LG1_LINE),
            "note": "244Cm 源的 XRF 激发来自子体 240Pu 的 L 线"
                    "（Gellert et al. MSL APXS；WUSTL APXS 帮助页）",
        },
    }
    energies = []
    for el in FIT_ELEMENTS:
        Z = x.SymbolToAtomicNumber(el)
        e = x.LineEnergy(Z, x.KA1_LINE)
        energies.append(round(e, 5))
        doc["line_energies_kev"][el] = round(e, 5)
        doc["k_edges_kev"][el] = round(x.EdgeEnergy(Z, x.K_SHELL), 5)
    doc["energies_kev"] = energies
    for el in MATRIX_ELEMENTS:
        Z = x.SymbolToAtomicNumber(el)
        doc["cs_total"][el] = [round(x.CS_Total(Z, e), 4) for e in energies]
    # Kβ/Kα 分支比（RadRate：Kβ1+Kβ3 / Kα1+Kα2；低 Z 无 Kβ3 时回退 Kβ1；
    # Kβ2 ~1-2% 忽略）
    doc["kbeta"] = {}
    for el in FIT_ELEMENTS:
        Z = x.SymbolToAtomicNumber(el)
        try:
            r_ka = x.RadRate(Z, x.KA1_LINE) + x.RadRate(Z, x.KA2_LINE)
            r_kb = 0.0
            for ln in (x.KB1_LINE, x.KB3_LINE):
                try:
                    r_kb += x.RadRate(Z, ln)
                except Exception:                        # noqa: BLE001
                    pass
            if r_ka > 0 and r_kb > 0:
                doc["kbeta"][el] = {"ratio": round(r_kb / r_ka, 5),
                                    "energy_kev": round(x.LineEnergy(Z, x.KB1_LINE), 5)}
        except Exception:                                # noqa: BLE001
            pass
    OUT.write_text(json.dumps(doc, ensure_ascii=False, indent=1),
                   encoding="utf-8")
    print("物理表写入", OUT)
    print("能量网格:", energies)
    return 0


if __name__ == "__main__":
    sys.exit(main())
