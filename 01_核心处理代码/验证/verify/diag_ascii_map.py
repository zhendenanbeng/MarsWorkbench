# -*- coding: utf-8 -*-
"""把一次扫描的逐点判定打成 10x10 的字符图（不用浏览器就能看）。

用法::

    python verify\\diag_ascii_map.py <数据文件夹> [--scan srlc11360] [--cluster]

为什么需要它：判断一个"未归类"的簇到底是真信号还是噪声，
最直接的办法就是看它的点在图上是否**连成一片**——
成片的极可能是真实矿物斑块，散点的多半是噪声。
浏览器里的那张图是同一份数据，这里只是换成终端能看的形式。
"""

from __future__ import annotations

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout

import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = str(_layout.project_layout.APP)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from sherloc import minerals, organize, qc              # noqa: E402


def render(res, by_cluster=False):
    """10x10 蛇形铺开。第 1 点在左上，第 10 点右上，第 11 点折回下一行右端。"""
    detail = {d["point"]: d for d in (res.get("point_detail") or [])}
    # 按点数从大到小给每个簇一个字母，图上就能看出"哪一片属于哪个簇"。
    groups = sorted(res["mineral_groups"], key=lambda g: (-g["n_points"], g["pos"]))
    tag = {round(g["pos"], 1): chr(65 + k) for k, g in enumerate(groups[:26])}
    legend = {chr(65 + k): g for k, g in enumerate(groups[:26])}
    grid = [[" "] * 10 for _ in range(10)]
    for i in range(min(100, res["n_points_analysed"])):
        r, c = divmod(i, 10)
        col = c if r % 2 == 0 else 9 - c
        d = detail.get(i)
        if d is None or (d["pos"] is None and not d["groups"]):
            grid[r][col] = "."
            continue
        if by_cluster:
            keys = [round(g, 1) for g in d["groups"]]
            grid[r][col] = tag.get(keys[0], "?") if keys else "?"
        else:
            grid[r][col] = minerals.SHORT.get(d["classes"][0]["key"], "?")
    return [" ".join(row) for row in grid], legend


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("folder")
    ap.add_argument("--scan", default=None)
    ap.add_argument("--cluster", action="store_true", help="按簇心百位画，而不是按类别")
    a = ap.parse_args()

    _products, observations = organize.scan_folder(a.folder)
    for ob in observations:
        if a.scan and a.scan not in ob.sequence_id:
            continue
        if ob.spectra is None:
            continue
        wl, r1, _r2, _r3 = ob.read_regions()
        if wl is None or r1 is None:
            continue
        res = qc.analyse_spectra(wl, r1)
        print("=" * 60)
        print("%s %s %s" % (ob.sol, ob.sequence_id, ob.sclk[:14]))
        print(qc.headline(res))
        for g in res["mineral_groups"]:
            print("   簇心 %6.1f 点 %2d 类别 %s"
                  % (g["pos"], g["n_points"],
                     "、".join(c["class"] for c in g["classes"])))
        print("   " + "-" * 19)
        rows, legend = render(res, a.cluster)
        for row in rows:
            print("   " + row)
        print("   " + "-" * 19)
        if a.cluster:
            for k, g in legend.items():
                print("   %s = 簇心 %6.1f（%d 点，%s）"
                      % (k, g["pos"], g["n_points"],
                         "、".join(c["class"] for c in g["classes"])))
            print("   ? = 落在簇里但簇没给字母；. = 没有合格波段")
        else:
            print("   图例 " + "  ".join("%s=%s" % (v, k)
                                        for k, v in minerals.SHORT.items())
                  + "  ?=未归类或单点")


if __name__ == "__main__":
    main()
