# -*- coding: utf-8 -*-
"""回答一个问题：那些"未归类"的簇，到底是"库里没这种矿物"，还是"区间没覆盖"？

做法：把全库每个簇都拿 62 种矿物的实测谱库打一遍分，分成两组比：

  A 已判出类别的簇   —— 谱库打分在这组上应该明显更高（否则打分方法本身不灵）
  B 未归类的簇       —— 如果 B 里出现高分候选，说明**库里有、只是峰位区间没覆盖**
                        如果 B 全是低分，说明这些簇的形状本来就对不上任何标准谱

这个区别很重要：前者是判类规则的问题（可以修），后者是数据的问题（修不了）。

用法：python verify\\library_vs_unassigned.py
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

from sherloc import library, minerals, qc            # noqa: E402
from sherloc.organize import scan_folder             # noqa: E402

DATA = os.path.join(os.path.dirname(ROOT), "03_下载数据/SHERLOC/练习数据")


def main():
    _p, observations = scan_folder(DATA)
    assigned, unassigned = [], []
    for ob in observations:
        if ob.spectra is None:
            continue
        wl, r1, _r2, _r3 = ob.read_regions()
        if wl is None or r1 is None:
            continue
        # 必须走和质控台同一条路（qc.analyse_spectra），否则死区裁剪和
        # 窗口口径都不一样，簇会跟报告里的对不上——实测踩过。
        res = qc.analyse_spectra(wl, r1)
        for g in res["mineral_groups"]:
            key = g["classes"][0]["key"] if g["classes"] else "unassigned"
            tag = (ob.sol, ob.sequence_id, round(g["pos"]), g["n_points"])
            lib = g.get("library") or []
            best = lib[0] if lib else None
            rec = (tag, key, g["band"]["center"], g["band"]["fwhm"],
                   (best["mineral"] if best else "-"),
                   (best["score"] if best else float("nan")),
                   (best["class_key"] if best else "-"),
                   [p[0] for p in (g.get("observed_peaks") or [])])
            (unassigned if key == "unassigned" else assigned).append(rec)

    print("谱库 = Razzell Hollis et al. (2021) 的 62 种矿物实测峰位（248.6 nm 激发）")
    print("得分 = 观测峰位与该矿物标准峰位的匹配度，0–1\n")

    for name, rows in (("A 已判出类别的簇", assigned), ("B 未归类的簇", unassigned)):
        sc = np.array([r[5] for r in rows if np.isfinite(r[5])])
        print("%s：%d 个" % (name, len(rows)))
        if sc.size:
            print("   谱库最高分  中位 %.2f   均值 %.2f   最高 %.2f   ≥0.5 的有 %d 个"
                  % (np.median(sc), sc.mean(), sc.max(), int((sc >= 0.5).sum())))
        print()

    print("未归类的簇逐条（按谱库得分排序）")
    print("%-6s %-9s %6s %5s %6s %6s  %-22s %s"
          % ("sol", "seq", "簇心", "点数", "质心", "FWHM", "谱库最佳", "观测峰"))
    for tag, key, center, fw, mineral, score, ck, peaks in sorted(
            unassigned, key=lambda r: -(r[5] if np.isfinite(r[5]) else -1)):
        print("%-6s %-9s %6.0f %5d %6.0f %6.0f  %-14s %.2f  %s"
              % (tag[0], tag[1], tag[2], tag[3], center, fw, mineral, score, peaks))

    print("\n对照：已判出类别的簇里谱库得分最高的 8 个")
    for tag, key, center, fw, mineral, score, ck, peaks in sorted(
            assigned, key=lambda r: -r[5])[:8]:
        print("  sol %-5s %-9s 簇心 %4.0f → 判为 %-12s 谱库 %-14s %.2f  %s"
              % (tag[0], tag[1], tag[2], minerals.LABELS.get(key, key),
                 mineral, score, peaks))
    return 0


if __name__ == "__main__":
    sys.exit(main())
