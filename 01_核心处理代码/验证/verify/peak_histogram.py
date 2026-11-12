# -*- coding: utf-8 -*-
"""看 1000-1120 cm-1 这一段里，逐点峰位到底是双峰还是单峰。

要回答的问题：sol 293 的碳酸盐报了 48 个点、未定硅酸盐报 0 个点。
这是"聚类把 1050 的硅酸盐宽带和 1083 的碳酸盐峰并成一组"造成的，
还是这一段本来就是**一个** 1079 附近的宽峰？

判据很直接：
  * 逐点峰位直方图出现两个峰（~1050 和 ~1083）→ 聚类在合并，要改聚类
  * 只有一个峰（~1079）→ 不是聚类问题，是峰本来就这么宽

用法：python verify\\peak_histogram.py
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

from sherloc import qc                      # noqa: E402
from sherloc.organize import scan_folder    # noqa: E402

DATA = os.path.join(os.path.dirname(ROOT), "03_下载数据/SHERLOC/练习数据")

TARGETS = {
    "293": ("碳酸盐峰 1085 · 硅酸盐宽带 ~1050", (950.0, 1150.0)),
    "257": ("橄榄石 827 · 碳酸盐 1083 · 辉石 993（论文：46/11/9）", (780.0, 1150.0)),
    "141": ("对照组：论文说这一段没有矿物", (950.0, 1150.0)),
}


def hist(positions, lo, hi, step=10.0):
    nb = int(round((hi - lo) / step))
    bins = [0] * nb
    for p in positions:
        k = int((p - lo) // step)
        if 0 <= k < nb:
            bins[k] += 1
    return bins


def main():
    _products, observations = scan_folder(DATA)
    todo = [ob for ob in observations
            if ob.sol in TARGETS and ob.spectra is not None]
    if not todo:
        print("没找到目标观测")
        return 1

    for ob in todo:
        wl, r1, r2, r3 = ob.read_regions()
        if wl is None or r1 is None:
            continue
        res = qc.analyse_spectra(wl, r1, r2, r3)
        lo, hi = TARGETS[ob.sol][1]
        det = res.get("point_detail") or []
        if not det:
            print("sol %s %s：没有逐点结果" % (ob.sol, ob.sequence_id))
            continue

        print("=" * 78)
        print("sol %s  %s  %s" % (ob.sol, ob.sequence_id, ob.sclk))
        print("  参考：" + TARGETS[ob.sol][0])

        sel = [d for d in det if d["pos"] is not None and lo <= d["pos"] <= hi]
        print("  落在 %.0f-%.0f cm-1 的点共 %d 个" % (lo, hi, len(sel)))
        if not sel:
            continue

        for tag, sub in (("全部（含 SNR<3）", sel),
                         ("只看 SNR>=3", [d for d in sel if d["snr"] >= 3.0])):
            if not sub:
                continue
            bins = hist([d["pos"] for d in sub], lo, hi)
            print("  %s，n=%d" % (tag, len(sub)))
            top = max(bins) or 1
            for k, n in enumerate(bins):
                if n == 0:
                    continue
                left = lo + k * 10
                bar = "#" * max(1, int(round(n / top * 34)))
                print("    %4.0f-%4.0f |%-34s| %d" % (left, left + 10, bar, n))

        # 现在被并进同一组的点，各自的峰位分布
        groups = res.get("mineral_groups") or []
        big = groups
        for g in big:
            bp = g.get("band_pos") or []
            name = g["classes"][0]["class"] if g["classes"] else "未归类"
            if not bp:
                continue
            print("  组：%s，聚类峰位 %.0f，%d 个点 / %d 个波段，"
                  "进组波段跨 %.0f cm-1（%.0f ~ %.0f），平均谱拟合 "
                  "峰心 %.0f / FWHM %.0f / SNR %.1f（%s）"
                  % (name, g["pos"], g["n_points"], len(bp),
                     max(bp) - min(bp), min(bp), max(bp),
                     g["band"]["center"], g["band"]["fwhm"], g["band"]["snr"],
                     g["band"].get("method", "?")
                     + ("" if "fallback" not in g["band"]
                        else " 退回原因:" + g["band"]["fallback"])))
            oh = g.get("oh")
            print("      O-H 带：%s" % (
                "没跑" if not oh else
                ("否掉了——" + (oh.get("reject") or "原因不明")) if not oh.get("ok") else
                "峰位 %.0f / FWHM %.0f / SNR %.1f / %s"
                % (oh["pos"], oh["fwhm"], oh["snr"], oh["method"])))
            st = g.get("split_test")
            if st is None:
                print("      拆分检验：没跑")
            elif st.get("ok"):
                print("      拆分检验：成立 → %.0f 与 %.0f"
                      % (st["comps"][0]["center"], st["comps"][1]["center"]))
            else:
                print("      拆分检验：不成立——%s" % st.get("reject"))
            # 把这一组的波段位置也打出来，看是不是两撮
            bins = hist(bp, lo, hi)
            row = " ".join("%.0f:%d" % (lo + k * 10, n)
                           for k, n in enumerate(bins) if n)
            print("      进组波段分布 " + row)
    return 0


if __name__ == "__main__":
    sys.exit(main())
