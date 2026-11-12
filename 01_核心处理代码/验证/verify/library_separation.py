# -*- coding: utf-8 -*-
"""规则判名组与未归类组的谱库分数区分度（代理 AUC）。

标签来自本软件的峰位规则，不能作为独立矿物真值。此 AUC 只度量两套
方法的内部一致性；分数再高也不证明矿物识别准确率，更不能单独验收算法。

为什么要这个东西：2026-09-15 量过一次，谱库比对给"已判出类别的组"和
"未归类的组"打出的最高分**中位数完全一样**（都是 0.33）——等于随便怎么对
都差不多，这种分数不能当证据用。所以以后每次改比对算法，都要跑这个脚本：
看这两个规则分组是否被分开（0.5 近似随机排序，1.0 完全分开）；
真实准确率需要另建未用于调参、带外部证据的盲标测试集。

用法：python verify\\library_separation.py [tol1 tol2 ...]
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

from sherloc import library, qc                       # noqa: E402
from sherloc.organize import scan_folder              # noqa: E402

DATA = os.path.join(os.path.dirname(ROOT), "03_下载数据/SHERLOC/练习数据")


def collect():
    """把全库的组收集起来：(是否判出类别, 观测峰, 类别键)。"""
    _p, observations = scan_folder(DATA)
    rows = []
    for ob in observations:
        if ob.spectra is None:
            continue
        wl, r1, _r2, _r3 = ob.read_regions()
        if wl is None or r1 is None:
            continue
        res = qc.analyse_spectra(wl, r1)
        for g in res.get("mineral_groups") or []:
            obs = g.get("observed_peaks") or []
            if not obs:
                continue
            key = g["classes"][0]["key"] if g["classes"] else "unassigned"
            rows.append((key != "unassigned", obs, key, g.get("pos")))
    return rows


def auc(pos, neg):
    a, b = np.asarray(pos, float), np.asarray(neg, float)
    if not a.size or not b.size:
        return float("nan")
    return float((a[:, None] > b[None, :]).mean()
                 + 0.5 * (a[:, None] == b[None, :]).mean())


def main():
    tols = [float(t) for t in sys.argv[1:]] or [8.0]
    rows = collect()
    print("全库 %d 个组（判出类别的 %d 个，未归类 %d 个）"
          % (len(rows), sum(1 for r in rows if r[0]), sum(1 for r in rows if not r[0])))
    print()
    print("%6s %8s %8s %8s %8s %8s" % ("tol", "已判定中位", "未归类中位", "AUC", "已判定最高", "未归类最高"))
    for tol in tols:
        pos, neg = [], []
        for named, obs, _k, _p in rows:
            m = library.match(obs, tol=tol)
            top = max((r["score"] for r in m), default=0.0)
            (pos if named else neg).append(top)
        print("%6.0f %8.3f %8.3f %8.3f %8.3f %8.3f"
              % (tol, np.median(pos), np.median(neg), auc(pos, neg),
                 max(pos), max(neg)))
    print()
    print("这是软件规则标签的代理AUC，不能解释为矿物识别准确率。")
    print("0.5约为随机排序、1.0完全分开这两个规则组；无通用的可用性验收阈值。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
