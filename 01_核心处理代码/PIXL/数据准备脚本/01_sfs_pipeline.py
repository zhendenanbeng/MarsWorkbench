#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""PIXL 数据准备脚本 01：SFS 逐 PMC 自研定量管线（正式版）。
=====================================================

输入：
  * PDS data_raw_spectroscopy 的 SFS csv 目录（逐 PMC 4096 通道 XRF 计数，
    ENA/ENB 双探测器；下载用工作台"原始谱"页或 pixl_download）
  * 可选：同扫描 RQC 成分表（官方衍生产品，只作对账基准）

输出（写 03_下载数据/PIXL/缓存/）：
  * sfs_energy_axis.json —— 能量轴自证（v4 多点线系）+ revisions
  * sfs_vs_rqc.json —— 逐 PMC Kα 峰面积（A/B 探测器低者）+ RQC wt
  * sfs_ratio_matrix.json —— 全元素两两 log-log 比值相关矩阵
  * 控制台摘要

自算/对账边界声明：
  自研变换 = cps 归一 → 能量轴多点自证 → 元素 Kα 窗净面积 → 逐 PMC
  A/B 探测器择低（Schmidt 2025 Table S4 口径）→ 元素比值信号。
  对账基准 = RQC（官方 wt%）逐 PMC 相关与比值矩阵；跨扫描绝对浓度
  标定不在本层（缺地形/衍射/死时间修正，属 PIQUANT 基本参数层，
  见 published_calibrations 的 maaz_crossscan_calibration_boundary）。

用法：
    python 01_sfs_pipeline.py                    # 默认 Bellegarde 全扫描集
    python 01_sfs_pipeline.py --dir <目录> --rqc <csv> [--no-write]
"""

from __future__ import annotations

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0, str(_BootPath(__file__).resolve().parents[3] / "02_图形界面代码"))
import source_paths as _layout            # noqa: E402

import argparse
import csv
import datetime
import json
import sys
from pathlib import Path

import numpy as np

from instruments import sfs_spectra       # noqa: E402

DR = _layout.project_layout.data_root()
CACHE = DR / "PIXL" / "缓存"
RATIO_ELEMENTS = ("Na", "Mg", "Al", "Si", "P", "K", "Ca", "Ti", "Mn", "Fe", "S", "Cl")
RQC_COL = {"Na": "Na2O", "Mg": "MgO", "Al": "Al2O3", "Si": "SiO2", "P": "P2O5",
           "K": "K2O", "Ca": "CaO", "Ti": "TiO2", "Mn": "MnO", "Fe": "FeO-T",
           "S": "SO3", "Cl": "Cl"}
UTC_NOW = lambda: datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")  # noqa: E731


def load_scan(scan_dir: Path):
    """读双探测器 SFS：{pmc: {det: cps 谱}}（行级 live_time 归一）。"""
    spec_by_pmc = {}
    lt_by_pmc = {}
    for f in sorted(scan_dir.glob("*.csv")):
        up = f.name.upper()
        det = "ENA" if "ENA" in up else ("ENB" if "ENB" in up else None)
        if det is None:
            continue
        sp, meta = sfs_spectra.read_sfs(f)
        if sp is None:
            continue
        sp = np.asarray(sp, dtype=float)
        for i, m in enumerate(meta):
            p = int(m.get("PMC") or 0)
            if p <= 0:
                continue
            lt = float(m.get("live_time") or 0)
            spec_by_pmc.setdefault(p, {})[det] = sp[i]
            lt_by_pmc.setdefault(p, {})[det] = lt if lt > 0 else np.nan
    return spec_by_pmc, lt_by_pmc


def patch_areas(pmc, spec_by_pmc, lt_by_pmc, axis, elements):
    """逐 PMC 峰面积：A/B 探测器取低者（缺失探测器件取另一件）。"""
    out = {el: None for el in elements}
    for el in elements:
        vals = []
        for det in ("ENA", "ENB"):
            if det not in spec_by_pmc[pmc]:
                continue
            pa = sfs_spectra.peak_areas(spec_by_pmc[pmc][det] / lt_by_pmc[pmc][det],
                                        axis, elements=(el,))
            if pa[el]["area"] is not None:
                vals.append(pa[el]["area"])
        out[el] = min(vals) if vals else None
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description="PIXL SFS 自研定量管线（能量自证+峰面积+比值矩阵）")
    ap.add_argument("--dir", default=None, help="SFS 下载目录（默认 sfs_sol0186）")
    ap.add_argument("--rqc", default=None, help="RQC csv 路径（默认随包示例）")
    ap.add_argument("--no-write", action="store_true", help="只打印不写缓存")
    args = ap.parse_args(argv)

    scan_dir = Path(args.dir) if args.dir else DR / "PIXL" / "下载" / "sfs_sol0186"
    if not scan_dir.is_dir():
        print("SFS 目录不存在：%s\n先用工作台下载 SFS（或 pixl_download）" % scan_dir)
        return 1
    rqc_p = Path(args.rqc) if args.rqc else next((DR / "PIXL" / "随包示例").glob("*rqc*.csv"), None)

    spec_by_pmc, lt_by_pmc = load_scan(scan_dir)
    pmcs = sorted(spec_by_pmc)
    if len(pmcs) < 3:
        print("有效 PMC 不足 3（现 %d）" % len(pmcs))
        return 1
    print("PMC 总数 %d（%d-%d）" % (len(pmcs), pmcs[0], pmcs[-1]))

    # 1) 轴自证（v4 多点线系）
    det_means = []
    for det in ("ENA", "ENB"):
        rows = [spec_by_pmc[p][det] / lt_by_pmc[p][det] for p in pmcs
                if det in spec_by_pmc[p] and np.isfinite(lt_by_pmc[p][det])]
        det_means.append(np.mean(np.vstack(rows), axis=0))
    axis = sfs_spectra.calibrate_energy(np.mean(det_means, axis=0))
    print("能量自证: g=%.6f b=%+.4f（%s）" % (axis["keV_per_ch"], axis["zero_keV"], axis["mode"]))

    # 2) 逐 PMC 峰面积（A/B 择低）
    elements = list(RATIO_ELEMENTS)
    areas = {el: [] for el in elements}
    for p in pmcs:
        pa = patch_areas(p, spec_by_pmc, lt_by_pmc, axis, elements)
        for el in elements:
            areas[el].append(pa[el])

    # 3) RQC 对账 + 比值矩阵
    ca_area = fe_area = ca_wt = fe_wt = None
    ca_corr = fe_corr = ratio_corr = None
    ratios = []
    if rqc_p and Path(rqc_p).is_file():
        rqc = {}
        for r in csv.DictReader(open(rqc_p, encoding="utf-8-sig", newline="")):
            try:
                rqc[int(float(r["PMC"]))] = r
            except (ValueError, KeyError):
                continue
        common = [p for p in pmcs if p in rqc]

        def wt_of(p, el):
            try:
                return float(rqc[p].get(RQC_COL[el] + "_wt%", "") or "nan")
            except ValueError:
                return float("nan")

        ca_a = np.array([areas["Ca"][i] for i, p in enumerate(pmcs) if p in rqc])
        fe_a = np.array([areas["Fe"][i] for i, p in enumerate(pmcs) if p in rqc])
        ca_wt = np.array([wt_of(p, "Ca") for p in common])
        fe_wt = np.array([wt_of(p, "Fe") for p in common])
        good = np.isfinite(ca_wt) & np.isfinite(fe_wt)
        ca_corr = float(np.corrcoef(ca_a[good], ca_wt[good])[0, 1])
        fe_corr = float(np.corrcoef(fe_a[good], fe_wt[good])[0, 1])
        ratio_corr = float(np.corrcoef(np.log(ca_a[good] / np.maximum(fe_a[good], 1e-9)),
                                       np.log(ca_wt[good] / np.maximum(fe_wt[good], 1e-9)))[0, 1])
        print("逐PMC相关: Ca %.4f | Fe %.4f | Ca/Fe 比值 log-log %.4f（n=%d）"
              % (ca_corr, fe_corr, ratio_corr, len(common)))
        ca_area = ca_a.tolist()
        fe_area = fe_a.tolist()
        ca_wt = ca_wt.tolist()
        fe_wt = fe_wt.tolist()

        # 全元素两两 log-log 比值矩阵
        wts = {el: np.array([wt_of(p, el) for p in common]) for el in elements}
        a_by_el = {el: np.array([areas[el][i] for i, p in enumerate(pmcs) if p in rqc])
                   for el in elements}
        for i, e1 in enumerate(elements):
            for e2 in elements[i + 1:]:
                ok = (np.isfinite(wts[e1]) & np.isfinite(wts[e2])
                      & (a_by_el[e1] > 0) & (a_by_el[e2] > 0)
                      & (wts[e1] > 0) & (wts[e2] > 0))
                if ok.sum() < 200:
                    continue
                c = float(np.corrcoef(np.log(a_by_el[e1][ok] / a_by_el[e2][ok]),
                                      np.log(wts[e1][ok] / wts[e2][ok]))[0, 1])
                ratios.append({"pair": "%s/%s" % (e1, e2), "loglog_corr": c,
                               "n": int(ok.sum())})
        ratios.sort(key=lambda r: -r["loglog_corr"])
        top = ratios[:6]
        print("比值矩阵（%d 对，前 6）: %s" % (
            len(ratios), "; ".join("%s %.4f" % (r["pair"], r["loglog_corr"]) for r in top)))
    else:
        print("RQC 不在位——跳过对账与比值矩阵")

    if args.no_write:
        print("--no-write：缓存未改动")
        return 0

    # 4) 写缓存（revisions 追加，旧值保留）
    CACHE.mkdir(parents=True, exist_ok=True)
    stamp = UTC_NOW()
    axis_p = CACHE / "sfs_energy_axis.json"
    old = json.loads(axis_p.read_text(encoding="utf-8")) if axis_p.is_file() else {}
    rev = old.get("revisions", [])
    if old.get("axis") and abs(old["axis"].get("keV_per_ch", 0) - axis["keV_per_ch"]) > 1e-12:
        rev = rev + [{"date": stamp,
                      "note": "轴自证更新（扫描集 %s）" % scan_dir.name,
                      "old_axis": old.get("axis"),
                      "old_ca_corr": old.get("ca_corr"), "old_fe_corr": old.get("fe_corr"),
                      "old_n_pmc": old.get("n_pmc")}]
    axis_p.write_text(json.dumps({
        "scope": "PIXL SFS 能量轴自证（v4 多点线系）。轴为扫描集平均 cps 谱的自证值，"
                 "可由 01_sfs_pipeline.py 复现；标定文件到位后应替换。",
        "axis": {"keV_per_ch": axis["keV_per_ch"], "zero_keV": axis["zero_keV"]},
        "mode": axis.get("mode"), "ca_corr": ca_corr, "fe_corr": fe_corr,
        "n_pmc": len(pmcs), "pmc_range": [pmcs[0], pmcs[-1]],
        "scan_dir": scan_dir.name,
        "source_sha256": {f.name: __import__("hashlib").sha256(f.read_bytes()).hexdigest()
                          for f in sorted(scan_dir.glob("*.csv"))},
        "generated_utc": stamp, "revisions": rev,
    }, ensure_ascii=False, indent=1), encoding="utf-8")

    if rqc_p and ca_area is not None:
        (CACHE / "sfs_vs_rqc.json").write_text(json.dumps({
            "pmc": common, "ca_area": ca_area, "fe_area": fe_area,
            "ca_wt": ca_wt, "fe_wt": fe_wt,
            "axis_keV_per_ch": axis["keV_per_ch"], "axis_zero_keV": axis["zero_keV"],
            "protocol": "逐 PMC Kα 峰面积取 A/B 探测器低者（Schmidt 2025 Table S4 口径）",
            "generated_utc": stamp}, ensure_ascii=False), encoding="utf-8")
        (CACHE / "sfs_ratio_matrix.json").write_text(json.dumps({
            "scope": "全元素两两 log-log 比值相关矩阵（面积比 vs RQC wt 比，逐 PMC）。"
                     "比值信号是半定量对账层，不是浓度。",
            "axis_keV_per_ch": axis["keV_per_ch"], "n_pmc": len(common),
            "ratios": ratios, "generated_utc": stamp}, ensure_ascii=False, indent=1),
            encoding="utf-8")
    print("缓存写入 %s" % CACHE)
    return 0


if __name__ == "__main__":
    sys.exit(main())
