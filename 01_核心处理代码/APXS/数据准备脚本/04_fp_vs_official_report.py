#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""APXS 数据准备脚本 04：FP 层 vs 官方 rwp 全量对账报告。

对 03 号脚本冻结标定的 69 靶逐靶跑 FP 预测，与官方 rwp（含 1σ 统计
精度列）逐氧化物对账：
  * 散点矩阵（FP vs 官方，拟合/留出分色，1:1 线）→ 缓存/fp_vs_official_comparison.png
  * 逐氧化物偏差统计（中位/均值相对误差、官方 ±2σ 吻合率）→
    缓存/fp_vs_official_summary.md
确定性：面积/灵敏度均读冻结台账，不重新拟合——任何机器复跑同数字。

用法：python 04_fp_vs_official_report.py [--no-write]
"""
from __future__ import annotations

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0, str(_BootPath(__file__).resolve().parents[3] / "02_图形界面代码"))
import source_paths as _layout            # noqa: E402

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

import project_layout                      # noqa: E402
from instruments import apxs_fp            # noqa: E402

HERE = Path(__file__).resolve().parent
LEDGER_FALLBACK = (_BootPath(__file__).resolve().parents[2]
                   / "验证/verify/fixtures/apxs/fp_calibration_ledger.json")


def load_official_with_sigma(rwp_path):
    """rwp 全 16 行 → {氧化物: (wt, 1σ 统计精度)}。"""
    rows = list(csv.reader(open(rwp_path, encoding="utf-8", errors="replace")))
    out = {}
    for r in rows[1:]:
        if len(r) >= 3:
            try:
                out[r[0].strip()] = (float(r[1]), float(r[2]))
            except ValueError:
                continue
    return out


def find_rwp(sol, prefix, dr):
    base = dr / "APXS" / "下载" / ("calib_targets" if prefix == "apb"
                                   else "calib_targets_apa")
    sol4 = sol[-4:]                       # 文件名内嵌 4 位 sol（目录名可能 5 位）
    hits = sorted(base.glob("ap?_*rsp%s*.csv" % sol4))
    return hits[0].with_name(hits[0].name.replace("rsp", "rwp", 1)) if hits else None


def main(argv=None):
    ap = argparse.ArgumentParser(description="FP vs 官方全量对账")
    ap.add_argument("--no-write", action="store_true")
    args = ap.parse_args(argv)

    dr = Path(project_layout.data_root())
    led_p = dr / "APXS" / "缓存" / "apxs_fp_calib.json"
    if not led_p.is_file():
        led_p = LEDGER_FALLBACK
    led = json.loads(led_p.read_text(encoding="utf-8"))
    cal = apxs_fp.load_calibration()
    if not cal:
        print("运行期标定不在册")
        return 1
    sens = cal["sensitivities"]

    recs = []
    miss = 0
    for rec in led["dataset"]:
        rwp = find_rwp(rec["sol"], rec["prefix"], dr)
        if rwp is None or not rwp.is_file():
            miss += 1
            continue
        official = load_official_with_sigma(rwp)
        pred = apxs_fp.solve(rec["areas"], sens)
        if not pred.get("available"):
            continue
        for ox, pv in pred["oxides"].items():
            ov = official.get(ox)
            if not ov or ov[0] <= 0.2:            # 0.2 wt% 以下不计
                continue
            recs.append({"sol": rec["sol"], "prefix": rec["prefix"],
                         "oxide": ox, "fp": pv, "official": ov[0],
                         "sigma": ov[1],
                         "is_val": (led["dataset"].index(rec) % 3 == 2)})
    print("对账行：%d（缺 rwp 的靶 %d 个）" % (len(recs), miss))

    # ---- 逐氧化物统计 ----
    order = [ox for ox in apxs_fp.OXIDE_ORDER
             if any(r["oxide"] == ox for r in recs)]
    stats = {}
    for ox in order:
        rows = [r for r in recs if r["oxide"] == ox]
        val = [r for r in rows if r["is_val"]]
        fit = [r for r in rows if not r["is_val"]]
        def agg(rs):
            if not rs:
                return {}
            rel = [abs(r["fp"] - r["official"]) / r["official"] for r in rs]
            in2s = sum(1 for r in rs
                       if abs(r["fp"] - r["official"]) <= 2 * max(r["sigma"], 1e-6))
            return {"n": len(rs),
                    "median_rel_err": float(np.median(rel)),
                    "mean_rel_err": float(np.mean(rel)),
                    "in_official_2sigma": in2s / len(rs)}
        stats[ox] = {"all": agg(rows), "fit": agg(fit), "holdout": agg(val)}

    lines = ["# APXS FP 层 vs 官方 rwp 全量对账（%d 靶，%d 行）" %
             (led["n_targets"], len(recs)), "",
             "FP = 冻结标定自算 wt%；官方 = 配对 rwp（ver 2，含 1σ 统计精度）。"
             "±2σ 吻合率 = |FP−官方| ≤ 2σ 的靶占比（σ 为官方统计精度，"
             "只含计数噪声，不含系统差——吻合率是保守下界口径）。", "",
             "| 氧化物 | n | 全量中位误差 | 全量±2σ吻合 | 留出中位误差 | 留出±2σ吻合 |",
             "|---|---|---|---|---|---|"]
    for ox in order:
        s = stats[ox]
        a, h = s["all"], s.get("holdout") or {}
        lines.append("| %s | %d | %.1f%% | %.0f%% | %s | %s |" % (
            ox, a["n"], 100 * a["median_rel_err"], 100 * a["in_official_2sigma"],
            ("%.1f%%" % (100 * h["median_rel_err"])) if h else "—",
            ("%.0f%%" % (100 * h["in_official_2sigma"])) if h else "—"))
    lines += ["", "拟合/留出切分与口径见 instruments/apxs_fp.py 与 03 号脚本台账；"
              "边界靶（rocknest3_rp sol 102）不在本数据集内，见台账 boundary_targets。"]
    summary = "\n".join(lines)

    # ---- 散点矩阵 ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
    plt.rcParams["axes.unicode_minus"] = False

    C_FIT, C_VAL, C_LN = "#1e5ad2", "#c82828", "#666666"
    n_ox = len(order)
    ncol = 4
    nrow = (n_ox + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.0 * ncol, 2.8 * nrow))
    axes = np.atleast_2d(axes)
    for k, ox in enumerate(order):
        ax = axes[k // ncol][k % ncol]
        rows = [r for r in recs if r["oxide"] == ox]
        for pop, color, label in ((False, C_FIT, "拟合靶"), (True, C_VAL, "留出靶")):
            rs = [r for r in rows if r["is_val"] == pop]
            if rs:
                ax.scatter([r["official"] for r in rs], [r["fp"] for r in rs],
                           s=14, alpha=.65, color=color, label=label if k == 0 else None)
        lims = [0, max(max(r["official"] for r in rows),
                       max(r["fp"] for r in rows)) * 1.15]
        ax.plot(lims, lims, ls="--", lw=.8, color=C_LN)
        ax.set_xlim(lims); ax.set_ylim(lims)
        s = stats[ox]["all"]
        ax.set_title("%s（中位误差 %.1f%%，n=%d）"
                     % (ox, 100 * s["median_rel_err"], s["n"]), fontsize=9)
        ax.tick_params(labelsize=7)
        if k // ncol == nrow - 1:
            ax.set_xlabel("官方 wt%", fontsize=8)
        if k % ncol == 0:
            ax.set_ylabel("FP 自算 wt%", fontsize=8)
    for k in range(n_ox, nrow * ncol):
        axes[k // ncol][k % ncol].axis("off")
    fig.legend(loc="upper right", fontsize=9, frameon=False)
    fig.suptitle("APXS 基本参数法-lite 自算 wt%% vs 官方 rwp（%d 靶逐靶对账，虚线=1:1）"
                 % led["n_targets"], fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))

    out_png = dr / "APXS" / "缓存" / "fp_vs_official_comparison.png"
    out_md = dr / "APXS" / "缓存" / "fp_vs_official_summary.md"
    if not args.no_write:
        fig.savefig(out_png, dpi=150)
        out_md.write_text(summary, encoding="utf-8")
        print("图 →", out_png)
        print("表 →", out_md)
    print()
    print(summary)
    return 0


if __name__ == "__main__":
    sys.exit(main())
