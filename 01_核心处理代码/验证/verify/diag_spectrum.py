# -*- coding: utf-8 -*-
"""把谱画成终端里的字符图 —— 没有浏览器、没有 matplotlib 也能看形状。

用法::

    python verify\\diag_spectrum.py <数据文件夹> --scan srlc11360 --mean
    python verify\\diag_spectrum.py <数据文件夹> --scan srlc11360 --point 7
    python verify\\diag_spectrum.py <数据文件夹> --scan srlc11360 --group 976.6

看的是什么：扣完基线之后的残差谱。判断"这一簇到底是一个宽包、还是几个真峰"，
最直接的证据就是形状本身。
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

from sherloc import minerals, organize, qc              # noqa: E402


def ascii_plot(x, y, lo, hi, width=78, height=17, marks=()):
    sel = (x >= lo) & (x <= hi)
    xs, ys = x[sel], y[sel]
    if xs.size == 0:
        print("   （这个区间没有数据）")
        return
    # 按像素宽度分箱，取每箱中位数，避免长谱挤成一团。
    idx = np.linspace(0, xs.size, width + 1).round().astype(int)
    bx, by = [], []
    for a, b in zip(idx[:-1], idx[1:]):
        if b > a:
            bx.append(float(np.median(xs[a:b])))
            by.append(float(np.median(ys[a:b])))
    by = np.array(by)
    top, bot = float(by.max()), float(by.min())
    rng = max(top - bot, 1e-9)
    grid = [[" "] * len(bx) for _ in range(height)]
    for j, v in enumerate(by):
        row = height - 1 - int(round((v - bot) / rng * (height - 1)))
        grid[row][j] = "*"
    print("   %10.1f ┐%s" % (top, "".join(grid[0])))
    for r in range(1, height):
        print("   %10s │%s" % ("", "".join(grid[r])))
    print("   %10.1f └%s" % (bot, "─" * len(bx)))
    # 参考线画在下面一行，避免把谱形挡掉
    tick = [" "] * len(bx)
    for m in marks:
        if not (bx[0] <= m <= bx[-1]):
            continue
        k = int(np.argmin(np.abs(np.array(bx) - m)))
        lab = "|%g" % m
        for t, ch in enumerate(lab):
            if k + t < len(bx):
                tick[k + t] = ch
    print("   %10s  %s" % ("", "".join(tick)))
    print("   %10s  %g cm-1 ~ %g cm-1" % ("", bx[0], bx[-1]))


def binned(x, y, step):
    """把谱按 step cm-1 分箱取中位数，打印成表。比字符图好读，也不受噪声影响。"""
    lo = float(np.floor(x.min() / step) * step)
    hi = float(x.max())
    rows = []
    b = lo
    while b < hi:
        sel = (x >= b) & (x < b + step)
        if sel.any():
            rows.append((0.5 * (b + (b + step)), float(np.median(y[sel])),
                         float(np.max(y[sel]))))
        b += step
    return rows


def resid_of(x, y):
    return y - minerals.baseline(x, y)


def prominence(y, j, span):
    """峰 j 的突出度：峰高减去「到左右第一个更高点之间最低谷」的较大者。

    这是判断一个"局部极大"到底是不是真峰的关键量。宽包上被噪声推出来的
    假极大，突出度往往只有几十计数（和噪声同量级）；真峰则有几百上千。
    """
    lo_end = max(0, j - span)
    hi_end = min(len(y) - 1, j + span)
    left = j
    lmin = y[j]
    while left > lo_end and y[left] <= y[j]:
        lmin = min(lmin, y[left])
        left -= 1
    right = j
    rmin = y[j]
    while right < hi_end and y[right] <= y[j]:
        rmin = min(rmin, y[right])
        right += 1
    return float(y[j] - max(lmin, rmin))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("folder")
    ap.add_argument("--scan", default=None)
    ap.add_argument("--lo", type=float, default=800.0)
    ap.add_argument("--hi", type=float, default=1500.0)
    ap.add_argument("--point", type=int, default=None)
    ap.add_argument("--group", type=float, default=None, help="簇心，画这一簇所有点的平均谱")
    ap.add_argument("--mean", action="store_true", help="画全扫描平均谱")
    ap.add_argument("--all-bands", action="store_true", help="列出所有点的全部合格波段")
    ap.add_argument("--profile", type=float, default=None,
                    help="按这个 cm-1 步长打印分箱中位数表")
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
        xs = res["x"]                      # 800-1500 波段轴
        sub = res["sub"]                   # 逐点扣基线后的残差
        detail = res["point_detail"]
        marks = (947, 970, 1010, 1050, 1085, 1134)

        print("=" * 78)
        print("%s %s %s  —  %s" % (ob.sol, ob.sequence_id, ob.sclk[:14],
                                   qc.headline(res)))

        if a.mean:
            print("--- 全扫描 %d 点平均谱（扣基线后） ---" % sub.shape[0])
            ascii_plot(xs, sub.mean(axis=0), a.lo, a.hi, marks=marks)

        if a.point is not None and 0 <= a.point < sub.shape[0]:
            d = detail[a.point]
            print("--- 第 %d 点（行 %d 列 %d）sigma %.1f ---"
                  % (a.point, a.point // 10 + 1,
                     (a.point % 10 if (a.point // 10) % 2 == 0 else 9 - a.point % 10) + 1,
                     d["sigma"]))
            ascii_plot(xs, sub[a.point], a.lo, a.hi, marks=marks)

        if a.group is not None:
            for g in res["mineral_groups"]:
                if abs(g["pos"] - a.group) > 5:
                    continue
                print("--- 簇心 %.1f 的 %d 个点平均谱（扣基线后） ---"
                      % (g["pos"], g["n_points"]))
                ascii_plot(xs, sub[g["pts"]].mean(axis=0), a.lo, a.hi, marks=marks)
                print("     波段位置 %s" % np.array2string(
                    np.array(g["band_pos"]), precision=0, max_line_width=180))
                print("     类别 %s" % "、".join(c["class"] for c in g["classes"]))

        if a.all_bands:
            print("--- 每个点的全部合格波段：pos / fwhm / snr / 突出度σ ---")
            n_weak = 0
            for d in detail:
                if not d["bands"]:
                    continue
                items = []
                for b in d["bands"]:
                    j = int(np.argmin(np.abs(xs - b["pos"])))
                    prom = prominence(sub[d["point"]], j, 60)
                    psig = prom / d["sigma"] if d["sigma"] else float("nan")
                    if psig < 1.0:
                        n_weak += 1
                    items.append("%.0f/%.0f/%.1f/%.1f"
                                 % (b["pos"], b["fwhm"], b["snr"], psig))
                txt = "  ".join(items)
                print("   点 %3d：%s" % (d["point"], txt))
            print("   突出度 < 1σ 的波段数：%d" % n_weak)

        if a.profile:
            def table(title, y):
                print("--- %s（每 %.0f cm-1 一箱：中位数 / 最大值） ---"
                      % (title, a.profile))
                line = []
                for c, med, mx in binned(xs, y, a.profile):
                    line.append("%.0f:%+.0f/%.0f" % (c, med, mx))
                for k in range(0, len(line), 4):
                    print("   " + "   ".join(line[k:k + 4]))

            table("全扫描平均谱", sub.mean(axis=0))
            if a.point is not None and 0 <= a.point < sub.shape[0]:
                table("第 %d 点残差谱" % a.point, sub[a.point])
            if a.group is not None:
                for g in res["mineral_groups"]:
                    if abs(g["pos"] - a.group) <= 5:
                        table("簇心 %.1f（%d 点）平均谱" % (g["pos"], g["n_points"]),
                              sub[g["pts"]].mean(axis=0))


if __name__ == "__main__":
    main()
