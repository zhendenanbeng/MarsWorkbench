# -*- coding: utf-8 -*-
"""突出度门槛扫描：prom_sigma 取不同值时，各靶点的判定怎么变。

用法::

    python verify\\prom_sweep.py            # 默认扫 0 / 1 / 2 / 3 / 4
    python verify\\prom_sweep.py 0 3        # 只比这两个

目的：sol 370（Alfalfa）有 42 个点次掉进"未归类"，怀疑是宽包上被噪声
推出来的假极大（突出度只有零点几个 σ）。这个脚本就是拿数据回答
"把突出度门槛加上去，到底是修好了还是改坏了"。
"""

from __future__ import annotations

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = str(_layout.project_layout.APP)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from sherloc import minerals                        # noqa: E402
from sherloc import qc                              # noqa: E402
from sherloc.organize import scan_folder            # noqa: E402

DATA = os.path.join(os.path.dirname(ROOT), "03_下载数据/SHERLOC/练习数据")

# 只看有点数、有矿物判定的那几次观测（小定标扫描跳过，省时间）。
TARGETS = {
    ("141", "srlc10203"): "对照 Foux",
    ("349", "srlc11420"): "Montpezat Survey(空)",
    ("162", "srlc10500"): "Máaz Guillaumes",
    ("186", "srlc10500"): "Bellegarde(100)",
    ("186", "srlc11420"): "Bellegarde survey",
    ("257", "srlc11360"): "Séítah Dourbes",
    ("293", "srlc11360"): "Séítah Quartier",
    ("349", "srlc11360"): "Montpezat(100)",
    ("370", "srlc11360"): "Alfalfa(100)",
    ("370", "srlc11420"): "Alfalfa survey",
}


def analyse_all():
    """跑一遍数据目录，返回 {(sol, seq): (最多点数的那个观测的 res, 描述)}。"""
    best = {}
    _products, observations = scan_folder(DATA)
    for ob in observations:
        if ob.spectra is None:
            continue
        key = (ob.sol, ob.sequence_id)
        if key not in TARGETS:
            continue
        wl, r1, _r2, _r3 = ob.read_regions()
        if wl is None or r1 is None:
            continue
        res = qc.analyse_spectra(wl, r1)
        if key not in best or res["n_points_total"] > best[key]["n_points_total"]:
            best[key] = res
    return best


def line(res):
    named = {}
    unassigned_pts, unassigned_groups = set(), 0
    for g in res["mineral_groups"]:
        if g["classes"] and g["classes"][0]["key"] != "unassigned":
            name = g["classes"][0]["class"]
            n, pos = named.get(name, (0, g["band"]["center"]))
            if g["n_points"] > n:
                pos = g["band"]["center"]
            named[name] = (n + g["n_points"], pos)
        else:
            unassigned_groups += 1
            unassigned_pts.update(g["pts"])
    # 只统计"有合格波段"的点，避免把 1296 点的空扫描混进来
    with_band = sum(1 for d in (res.get("point_detail") or [])
                    if d["pos"] is not None)
    txt = " | ".join("%s %d@%.0f" % (k, v[0], v[1])
                     for k, v in sorted(named.items(), key=lambda kv: -kv[1][0]))
    return ("有点 %3d | 未归类 %2d簇/%2d点 | %s"
            % (with_band, unassigned_groups, len(unassigned_pts), txt or "（无）"))


def main():
    proms = [float(v) for v in sys.argv[1:]] or [0.0, 1.0, 2.0, 3.0, 4.0]
    for prom in proms:
        minerals.PROM_SIGMA = prom
        print("=" * 100)
        print("prom_sigma = %.1f" % prom)
        best = analyse_all()
        for (sol, seq), desc in TARGETS.items():
            res = best.get((sol, seq))
            if res is None:
                print("  sol %-4s %-9s %-18s （没跑到）" % (sol, seq, desc))
                continue
            print("  sol %-4s %-9s %-18s %s" % (sol, seq, desc, line(res)))


if __name__ == "__main__":
    main()
