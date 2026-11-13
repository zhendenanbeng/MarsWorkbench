#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""APXS 数据准备脚本 02：Si 内标比值标定层。
==========================================

背景（对 01_multitarget_calib 结论的修正）：
  绝对峰面积的经验标定斜率≈0——火星自然靶 25 靶的跨靶几何/接触距离/
  停留时间因子主导绝对强度，成分信号被淹没。改用 **Si 内标比值**
  （元素/Si 峰面积比 vs 官方 wt 比）后几何因子相消，Ca/Fe/K/Mg 的
  跨靶 log-log 相关 0.89-0.98，幂律指数稳定（0.72-0.93，自吸收所致）。

输入：03_下载数据/APXS/缓存/apxs_calib.json（01 号脚本的 25 靶面积+官方 wt）
输出：03_下载数据/APXS/缓存/apxs_ratio_calib.json
  * 每元素幂律标定 wt比 = A × (面积比)^b，70% 拟合 / 30% 留出验证
  * 比值预测是半定量信号（元素/SiO2 比值），不是绝对浓度；
    绝对浓度需基本参数法（GUAPX/PIQUANT 类），不在本层。

用法：python 02_ratio_calibration.py [--no-write]
"""

from __future__ import annotations

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0, str(_BootPath(__file__).resolve().parents[3] / "02_图形界面代码"))
import source_paths as _layout            # noqa: E402

import argparse
import datetime
import json
import sys
from pathlib import Path

import numpy as np

DR = _layout.project_layout.data_root()
CAL_P = DR / "APXS" / "缓存" / "apxs_calib.json"
OUT_P = DR / "APXS" / "缓存" / "apxs_ratio_calib.json"
ELEMENTS = ("Na", "Mg", "Al", "K", "Ca", "Ti", "Fe")


def main(argv=None):
    ap = argparse.ArgumentParser(description="APXS Si 内标比值标定（25 靶）")
    ap.add_argument("--no-write", action="store_true")
    args = ap.parse_args(argv)
    if not CAL_P.is_file():
        print("输入不在位：%s\n先运行 01_multitarget_calib.py" % CAL_P)
        return 1
    doc = json.loads(CAL_P.read_text(encoding="utf-8"))
    ds = doc.get("dataset") or []

    # 确定性切分：按 sol 排序，idx%3==2 留出（与 01 号脚本 17/8 惯例一致）
    targets = sorted(ds, key=lambda t: str(t.get("sol")))
    fit_idx = [i for i in range(len(targets)) if i % 3 != 2]
    val_idx = [i for i in range(len(targets)) if i % 3 == 2]
    print("靶 %d 个：拟合 %d / 留出 %d" % (len(targets), len(fit_idx), len(val_idx)))

    out = {"scope": "APXS Si 内标比值标定（元素/Si 峰面积比 vs 官方 wt 比，幂律双参，"
                    "70%% 拟合/30%% 留出）。比值预测是半定量信号，不是绝对浓度；"
                    "绝对浓度需基本参数法（GUAPX/PIQUANT 类），不在本层。"
                    "背景：绝对面积标定斜率≈0（跨靶几何因子主导），比值层修复之。",
           "generated_utc": datetime.datetime.now(datetime.timezone.utc)
           .strftime("%Y-%m-%dT%H:%M:%SZ"),
           "n_targets": len(targets), "n_fit": len(fit_idx), "n_val": len(val_idx),
           "elements": {}}

    print()
    print("=== 幂律标定 log(wt比) = logA + b×log(面积比)，留出验证 ===")
    for e in ELEMENTS:
        rows = []
        for i, t in enumerate(targets):
            x = (t.get("elements") or {}).get(e) or {}
            si = (t.get("elements") or {}).get("Si") or {}
            if not (x.get("area") and x.get("official_wt")
                    and si.get("area") and si.get("official_wt")):
                continue
            if min(x["area"], si["area"], x["official_wt"], si["official_wt"]) <= 0:
                continue
            rows.append((i, np.log(x["area"] / si["area"]),
                         np.log(x["official_wt"] / si["official_wt"])))
        if len(rows) < 10:
            print("  %-3s 有效靶不足（%d）" % (e, len(rows)))
            continue
        la = np.array([r[1] for r in rows])
        lw = np.array([r[2] for r in rows])
        idxs = [r[0] for r in rows]
        fmask = np.array([i in fit_idx for i in idxs])
        vmask = ~fmask
        # 全样本相关 + 幂律参数（拟合子集）
        c_all = float(np.corrcoef(la, lw)[0, 1])
        b, log_a = np.polyfit(la[fmask], lw[fmask], 1)
        # 留出验证：预测 wt 比 → 相对误差
        errs = []
        for j in np.where(vmask)[0]:
            pred_ratio = np.exp(log_a + b * la[j])
            true_ratio = np.exp(lw[j])
            if true_ratio > 0:
                errs.append(abs(pred_ratio / true_ratio - 1))
        med_rel = float(np.median(errs)) if errs else None
        out["elements"][e] = {
            "corr_all": c_all, "power_b": float(b), "power_A": float(np.exp(log_a)),
            "n_rows": len(rows), "val_rel_err_median": med_rel}
        print("  %-3s 全样本相关 %+.4f  b=%.3f A=%.3f  留出中位相对误差 %s（n=%d）"
              % (e, c_all, b, float(np.exp(log_a)),
                 "%.1f%%" % (med_rel * 100) if med_rel is not None else "-", len(rows)))

    strong = {e: v for e, v in out["elements"].items() if v["corr_all"] >= 0.88}
    print()
    print("比值层可用元素（相关 ≥0.88）: %s" % ", ".join(sorted(strong)))

    # 外部留出预测：rocknest3_rp（sol 102）不在 25 靶内（靶范围 sol 629-4327）——
    # 冻结幂律预测其元素/SiO2 比值，与官方 rwp 对比（独立靶 held-out 验证）
    fix = _layout.project_layout.VERIFY / "fixtures/apxs/sol00102"
    rsp = fix / "apa_406545731rsp01020050178_______p1.csv"
    rwp = fix / "apa_406545731rwp01020050178_______p1.csv"
    holdout = None
    if rsp.is_file() and rwp.is_file():
        from instruments import apxs as apxs_backend
        from instruments import apxs_fit
        data = apxs_backend.pds3.read_table(str(rsp), skiprows=1)
        energy = data[:, 0] / 1000.0
        counts = data[:, 1] if data.shape[1] == 2 else np.median(data[:, 1:-1], axis=1)
        re = apxs_fit.spectrum_areas(energy[1:-1], counts[1:-1])
        rwp_ox = {}
        for line in rwp.read_text(encoding="utf-8").splitlines()[1:]:
            parts = [p.strip() for p in line.split(",")]
            if len(parts) >= 2 and parts[0]:
                try:
                    rwp_ox[parts[0]] = float(parts[1])
                except ValueError:
                    continue
        si_a = (re.get("Si") or {}).get("area")
        si_w = rwp_ox.get("SiO2")
        if si_a and si_w:
            holdout = {}
            print()
            print("=== 外部留出预测：rocknest3_rp（sol 102，不在 25 靶内）===")
            for e in ELEMENTS:
                cal = out["elements"].get(e)
                r = re.get(e)
                if not cal or not r or r.get("area") is None:
                    continue
                ox = {"Na": "Na2O", "Mg": "MgO", "Al": "Al2O3", "K": "K2O",
                      "Ca": "CaO", "Ti": "TiO2", "Fe": "FeO"}.get(e)
                true_ratio = (rwp_ox.get(ox) / si_w) if (ox and rwp_ox.get(ox) and si_w) else None
                if not true_ratio:
                    continue
                pred_ratio = cal["power_A"] * (r["area"] / si_a) ** cal["power_b"]
                rel = pred_ratio / true_ratio - 1
                holdout[e] = {"pred_ratio": pred_ratio, "true_ratio": true_ratio,
                              "rel_dev": rel}
                print("  %-3s 预测 wt比 %.3f vs 官方 %.3f（%+.1f%%）"
                      % (e, pred_ratio, true_ratio, rel * 100))
    if holdout:
        out["holdout_rocknest3_rp"] = {
            "sol": 102, "verdict": "边界外——预测偏差 -85%~+60%，不可用",
            "note": "sol 102 早于标定时代（25 靶 sol 629-4327）。APXS 探测器随任务"
                    "时间辐照退化，效率漂移使幂律不可跨时代外推；比值层适用域="
                    "标定时代内的靶。负结果如实记录，不作成功声明。",
            "ratios": holdout}

    if args.no_write:
        print("--no-write：未写缓存")
        return 0
    OUT_P.parent.mkdir(parents=True, exist_ok=True)
    OUT_P.write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    print("写入", OUT_P)
    return 0


if __name__ == "__main__":
    sys.exit(main())
