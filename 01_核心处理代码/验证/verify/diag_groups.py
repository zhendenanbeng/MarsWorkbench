# -*- coding: utf-8 -*-
r"""诊断工具：把某次扫描的聚类组逐条打出来，必要时按"段"拆开看。

用法::

    python verify\diag_groups.py <数据文件夹> [--scan srlc11360] [--per-part]

为什么需要它：质控台的 HTML 只给结论，看不到"这一簇的波段位置到底散在哪、
平均谱上的次级峰长什么样"。要判断一个组是聚类合并错了、还是判类区间没覆盖，
就必须把中间量摊开看。

它走的是和质控台**完全相同**的代码路径（qc.analyse_spectra），
只是把里面 minerals.analyse_scan 的组信息原样打印。
"""

from __future__ import annotations

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout

import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = str(_layout.project_layout.APP)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from sherloc import organize, qc                        # noqa: E402


def show(tag, wl, r1):
    print("=" * 78)
    print("%s — %d 点" % (tag, r1.shape[0]))
    res = qc.analyse_spectra(wl, r1)
    print(qc.headline(res))
    print("  有合格波段的点 %d | 判定分布 %s"
          % (res.get("n_points_with_band") or 0, res["counts"]))
    for g in res["mineral_groups"]:
        b = g["band"]
        cls = "、".join(c["class"] for c in g["classes"])
        bp = np.array(g["band_pos"])
        print("  -- 簇心 %6.1f 点 %3d 质心 %6.1f FWHM %5.1f SNR %5.1f 法 %-8s"
              % (g["pos"], g["n_points"], b["center"], b["fwhm"],
                 b.get("snr", float("nan")), b["method"]))
        print("     波段位置 %s" % np.array2string(bp, precision=0, max_line_width=180))
        print("     组噪声 sigma %.1f | 谱库用观测峰 %s"
              % (g["sigma_mean"], g["observed_peaks"]))
        print("     类别 %s" % cls)
        if g.get("fallback"):
            print("     拟合退回：%s" % g["fallback"])
        for lib in (g.get("library") or [])[:3]:
            print("     谱库：%s %s" % (lib.get("name"), lib.get("score")))
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("folder")
    ap.add_argument("--scan", default=None, help="只看序列号里含这个串的扫描")
    ap.add_argument("--per-part", action="store_true", help="把多段扫描按段分开看")
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
        show("%s %s %s（整段）" % (ob.sol, ob.sequence_id, ob.sclk[:14]), wl, r1)

        if a.per_part and len(ob.parts) > 1:
            off = 0
            for p in ob.parts:
                n = p.n_spectra or 0
                if n:
                    show("    |_ 段 expID=%s sclk=%s"
                         % (p.pid.flags.get("expid", "?"), p.pid.sclk[:14]),
                         wl, r1[off:off + n])
                off += n


if __name__ == "__main__":
    main()
