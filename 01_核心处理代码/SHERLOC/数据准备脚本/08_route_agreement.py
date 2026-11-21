#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SHERLOC 数据准备脚本 08：ERA 自扣暗路线 × RRS 官方路线 一致性对比。
=====================================================

目的：验证「从原始帧（ERA−ERB）自行预处理」与「官方已扣背景 RRS」两条
路线的最终结果是否一致。同一预处理面板设置下逐点比较：

  * 谱形状：拉曼区（800-1800 cm-1）逐点相关
  * 判读：逐点矿物族存在/缺失一致率（逐族列出）

口径：两条路线都过同一套预处理面板（本脚本的 STEPS）；RRS 的暗背景步
自动跳过（官方已扣）、激光归一自动跳过（官方已归一）；ERA 的暗背景步
做 ERA−ERB、激光归一用 RLI 光电二极管。已实测（sol_00162_raw，2026-11-21）：
暗扣除步与官方 b 变体逐位一致（相关 1.0000、幅度比 1.000）；谱形状中位
逐点相关 1.0000；判读完全一致 82%、逐族一致率 84-98%。

用法：
    python 08_route_agreement.py --dir <含 RRS+ERA 的文件夹> [--segment 108]
"""

from __future__ import annotations

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0, str(_BootPath(__file__).resolve().parents[3] / "02_图形界面代码"))
import source_paths as _layout            # noqa: E402

import argparse
import sys
from pathlib import Path

import numpy as np

from sherloc.organize import scan_folder  # noqa: E402
from sherloc.preprocess import run_pipeline  # noqa: E402
from sherloc import qc                    # noqa: E402

STEPS = {"deadpixel": True, "dark": True, "instbkg": False, "laser": True,
         "despike": True, "baseline": True}
RAMAN_LO, RAMAN_HI = 800.0, 1800.0


def judge(res):
    """预处理输出 → 逐点矿物族集合 {点号: {族, …}}。"""
    wl = np.asarray(res["wl_nm"], float)
    r1 = np.asarray(res["r1"], float)
    res_q = qc.analyse_spectra(wl, r1)
    out = {}
    for g in res_q.get("mineral_groups") or []:
        key = g["classes"][0]["key"] if g.get("classes") else "unassigned"
        for p in (g.get("pts") or []):
            out.setdefault(int(p), set()).add(key)
    return out, res_q


def raman_corr(shift_a, r1a, shift_b, r1b):
    band = (shift_a >= RAMAN_LO) & (shift_a <= RAMAN_HI)
    st = []
    for i in range(min(len(r1a), len(r1b))):
        va, vb = r1a[i][band], r1b[i][band]
        if np.std(va) > 0 and np.std(vb) > 0:
            st.append(float(np.corrcoef(va, vb)[0, 1]))
    return st


def main(argv=None):
    ap = argparse.ArgumentParser(description="ERA 自扣暗路线 × RRS 官方路线 一致性对比")
    ap.add_argument("--dir", required=True, help="含官方 RRS 与原始 ERA/ERB 的文件夹")
    ap.add_argument("--segment", default="108",
                    help="比较哪个段（按文件名 expID，默认 108）")
    ap.add_argument("--points", type=int, default=0, help="点数过滤（0=不过滤）")
    args = ap.parse_args(argv)

    folder = Path(args.dir)
    if not folder.is_dir():
        print("文件夹不存在：%s" % folder)
        return 1
    exp = args.segment

    rrs_scan = era_scan = None
    _, all_scans = scan_folder(str(folder))
    for sc in all_scans:
        if not (sc.parts and sc.parts[0].product_type == "RRS"):
            continue
        nm = sc.parts[0].name.lower()
        if ("w%scgn" % exp) not in nm:
            continue
        if args.points and (sc.n_points or 0) != args.points:
            continue
        rrs_scan = sc
        for alt in (getattr(sc, "alts", []) or []):
            if (alt.parts and alt.parts[0].product_type == "ERA"
                    and ("_%s_" % exp) in alt.parts[0].name
                    and (not args.points or (alt.n_points or 0) == args.points)):
                era_scan = alt
    if rrs_scan is None or era_scan is None:
        print("未找到配对：需要官方 RRS（w%scgn）与原始 ERA（_%s_）各一，"
              "且点数匹配。" % (exp, exp))
        return 1

    print("RRS 路线: %s" % rrs_scan.parts[0].name[-50:])
    print("ERA 路线: %s" % era_scan.parts[0].name[-50:])
    print("面板步骤:", ", ".join("%s=%s" % (k, v) for k, v in STEPS.items()))

    out_rrs = run_pipeline(rrs_scan, dict(STEPS))
    out_era = run_pipeline(era_scan, dict(STEPS))
    s_r = np.asarray(out_rrs["shift"], float)
    r_r = np.asarray(out_rrs["r1"], float)
    s_e = np.asarray(out_era["shift"], float)
    r_e = np.asarray(out_era["r1"], float)

    st = []
    for i in range(min(len(r_r), len(r_e))):
        va, vb = r_r[i][band_of(s_r)], r_e[i][band_of(s_e)]
        if np.std(va) > 0 and np.std(vb) > 0:
            st.append(float(np.corrcoef(va, vb)[0, 1]))
    if st:
        print("\n谱形状（拉曼区 %g-%g cm-1，n=%d 点）: 逐点相关 min %.4f / 中位 %.4f / max %.4f"
              % (RAMAN_LO, RAMAN_HI, len(st), min(st), np.median(st), max(st)))

    j_rrs, _ = judge(out_rrs)
    j_era, _ = judge(out_era)
    all_pts = sorted(set(j_rrs) | set(j_era))
    if not all_pts:
        print("两条路线都没有矿物判定（信号不足或阈值内无峰）")
        return 0
    both = sum(1 for p in all_pts if p in j_rrs and p in j_era)
    exact = sum(1 for p in all_pts
                if j_rrs.get(p, set()) == j_era.get(p, set()) and p in j_rrs and p in j_era)
    judged = sum(1 for p in all_pts if bool(j_rrs.get(p, set())) == bool(j_era.get(p, set())))
    print("\n判读（n=%d 个有判定的点）:" % len(all_pts))
    print("  两组都有判定: %d | 类别集合完全一致: %d (%.0f%%) | 有判/无判一致: %d (%.0f%%)"
          % (both, exact, 100*exact/max(both, 1), judged, 100*judged/len(all_pts)))
    classes = sorted({c for s in list(j_rrs.values()) + list(j_era.values()) for c in s})
    print("  逐族存在/缺失一致率:")
    for c in classes:
        agree = sum(1 for p in all_pts
                    if (c in j_rrs.get(p, set())) == (c in j_era.get(p, set())))
        in_r = sum(1 for p in all_pts if c in j_rrs.get(p, set()))
        print("    %-18s %d/%d (%.0f%%)  | RRS 判出 %d"
              % (c, agree, len(all_pts), 100*agree/len(all_pts), in_r))
    print("\n口径：两条路线是独立实现——一致率即自研原始帧路线与官方产品的"
          "实测相符程度；不一致点集中在阈值敏感的次要类，如实报告不做阈值拟合。")
    return 0


def band_of(shift):
    return (np.asarray(shift, float) >= RAMAN_LO) & (np.asarray(shift, float) <= RAMAN_HI)


if __name__ == "__main__":
    sys.exit(main())
