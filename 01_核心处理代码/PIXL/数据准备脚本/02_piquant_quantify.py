#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""PIXL 数据准备脚本 02：PIQUANT（官方开源定量）× 官方 RQC 逐 PMC 对账。
=====================================================

这是三层对账的第二层（PIQUANT 层）：
  * 第一层（01_sfs_pipeline.py）：自研峰面积 vs RQC——半定量（比值/扫描内标定）
  * 第二层（本脚本）：PIQUANT 基本参数法 vs RQC——绝对 wt% 对账
  * PIQUANT 是 PIXL 官方开源定量代码（github.com/pixlise/piquant，
    Caltech/UW 许可），本层用其 Windows 构建对同一份 SFS 原始谱做
    官方算法复现，输出与本仓库自研层口径分开。

探测器口径：FM 标定文件（Calibration_PIXL_FM_SurfaceOps_5minECFs_Rev1_
Jul2021.csv，随 PIQUANT 仓库分发）官方 COMMENT 明示"Use Detector B
spectra"——本脚本量化 ENB 谱；RQC wt% 是官方 A/B 择低后逐 PMC 值
（Schmidt 2025 Table S4 口径），两边口径差异计入对账偏差，不作合并。

流程：SFS 双探测器读取（复用 01 的 load_scan）→ 能量轴自证（v4）→
按步长抽样 PMC → 逐 PMC 导出 MSA（PIQUANT 原生谱格式，带本扫描
自证轴）→ Piquant quant（FM 配置+FM 标定）→ 解析 wt% → 与 RQC
逐元素偏差统计 → 写 缓存/piquant_vs_rqc.json + 报告 md。

用法：
    python 02_piquant_quantify.py                 # 默认 sol0186，步长 10
    python 02_piquant_quantify.py --stride 5 --limit 50
    python 02_piquant_quantify.py --dir <目录> --rqc <csv> --no-write
"""

from __future__ import annotations

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0, str(_BootPath(__file__).resolve().parents[3] / "02_图形界面代码"))
import source_paths as _layout            # noqa: E402

import argparse
import csv
import datetime
import importlib.util
import json
import sys
from pathlib import Path

from instruments import piquant_bridge    # noqa: E402

DR = _layout.project_layout.data_root()
CACHE = DR / "PIXL" / "缓存"
RQC_COL = {"Na2O": "Na2O", "MgO": "MgO", "Al2O3": "Al2O3", "SiO2": "SiO2",
           "P2O5": "P2O5", "K2O": "K2O", "CaO": "CaO", "TiO2": "TiO2",
           "MnO": "MnO", "FeO-T": "FeO-T", "SO3": "SO3", "Cl": "Cl"}
# FM 标定覆盖的元素（缺 P/Mn/Cr/Sr 的 ECF 条目）→ PIQUANT 可对账的氧化物
PIQUANT_OXIDES = ("Na2O", "MgO", "Al2O3", "SiO2", "K2O", "CaO", "TiO2",
                  "FeO-T", "SO3", "Cl")
ELEMENT_LINES = "Na_K Mg_K Al_K Si_K K_K Ca_K Ti_K Fe_K S_K Cl_K Ar_I"
UTC_NOW = lambda: datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")  # noqa: E731


def _load_pipeline():
    """导入 01_sfs_pipeline（数字开头模块名，需 importlib 按路径加载）。"""
    p = Path(__file__).resolve().parent / "01_sfs_pipeline.py"
    spec = importlib.util.spec_from_file_location("sfs_pipeline_01", p)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["sfs_pipeline_01"] = mod
    spec.loader.exec_module(mod)
    return mod


def load_rqc(rqc_p: Path):
    rows = {}
    for r in csv.DictReader(open(rqc_p, encoding="utf-8-sig", newline="")):
        try:
            rows[int(float(r["PMC"]))] = r
        except (ValueError, KeyError):
            continue
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description="PIQUANT 官方开源定量 × RQC 逐 PMC 对账")
    ap.add_argument("--dir", default=None, help="SFS 下载目录（默认 sfs_sol0186）")
    ap.add_argument("--rqc", default=None, help="RQC csv 路径（默认随包示例）")
    ap.add_argument("--stride", type=int, default=10,
                    help="PMC 抽样步长（默认 10；每谱 ~3s，步长 10 ≈ 280 点）")
    ap.add_argument("--limit", type=int, default=0, help="最多量化多少个 PMC（0=不限）")
    ap.add_argument("--detector", choices=("ENB", "ENA"), default="ENB",
                    help="量化探测器（默认 ENB=FM 标定官方口径；ENA 未标定仅作参考）")
    ap.add_argument("--no-write", action="store_true", help="只打印不写缓存")
    args = ap.parse_args(argv)

    st = piquant_bridge.status()
    if not st["available"]:
        print("PIQUANT 不可用：%s" % st["reason"])
        return 1
    print("PIQUANT: %s" % st["bin"])
    print("FM 配置: %s" % st["config"])
    print("FM 标定: %s" % st["calibration"])

    pipe = _load_pipeline()
    scan_dir = Path(args.dir) if args.dir else DR / "PIXL" / "下载" / "sfs_sol0186"
    if not scan_dir.is_dir():
        print("SFS 目录不存在：%s" % scan_dir)
        return 1
    rqc_p = Path(args.rqc) if args.rqc else next((DR / "PIXL" / "随包示例").glob("*rqc*.csv"), None)

    spec_by_pmc, lt_by_pmc = pipe.load_scan(scan_dir)
    pmcs = sorted(spec_by_pmc)
    if len(pmcs) < 3:
        print("有效 PMC 不足 3（现 %d）" % len(pmcs))
        return 1
    print("PMC 总数 %d（%d-%d）" % (len(pmcs), pmcs[0], pmcs[-1]))

    # 1) 轴自证（与 01 相同的两探测器均值 v4 流程）
    det_means = []
    for det in ("ENA", "ENB"):
        rows_d = [spec_by_pmc[p][det] / lt_by_pmc[p][det] for p in pmcs
                  if det in spec_by_pmc[p] and lt_by_pmc[p][det] > 0
                  and str(lt_by_pmc[p][det]) != "nan"]
        if rows_d:
            det_means.append(__import__("numpy").mean(__import__("numpy").vstack(rows_d), axis=0))
    import numpy as np
    axis = pipe.sfs_spectra.calibrate_energy(np.mean(det_means, axis=0))
    print("能量自证: g=%.6f b=%+.4f（%s）" % (axis["keV_per_ch"], axis["zero_keV"], axis["mode"]))

    # 2) RQC 覆盖率门（PMC 是扫描内局部计数，双向覆盖 ≥95% 才算同扫描）
    if not (rqc_p and Path(rqc_p).is_file()):
        print("RQC 不在位——无对账基准，退出（量化结果不落盘）")
        return 1
    rqc = load_rqc(Path(rqc_p))
    common = [p for p in pmcs if p in rqc]
    cover_sfs = len(common) / max(len(pmcs), 1)
    cover_rqc = len(common) / max(len(rqc), 1)
    if cover_sfs < 0.95 or cover_rqc < 0.95:
        print("对账跳过：PMC 覆盖率不足（SFS 侧 %.0f%% / RQC 侧 %.0f%%）——扫描集与该 RQC 不匹配"
              % (cover_sfs * 100, cover_rqc * 100))
        return 1
    print("RQC: %s（%d 行，覆盖率 %.0f%%/%.0f%%）" % (Path(rqc_p).name, len(rqc),
                                                     cover_sfs * 100, cover_rqc * 100))

    # 3) 抽样逐 PMC 量化
    sampled = [p for i, p in enumerate(pmcs) if i % max(args.stride, 1) == 0]
    if args.limit > 0:
        sampled = sampled[:args.limit]
    print("抽样 %d/%d PMC（步长 %d，探测器 %s）——每谱约 3s，预计 ~%.0f 分钟"
          % (len(sampled), len(pmcs), args.stride, args.detector,
             len(sampled) * 3.0 / 60))
    work = CACHE / "piquant_work"
    rows, fails = [], 0
    for n, p in enumerate(sampled, 1):
        det = args.detector
        if det not in spec_by_pmc[p]:
            det = "ENA" if det == "ENB" else "ENB"
            if det not in spec_by_pmc[p]:
                continue
        raw = spec_by_pmc[p][det]
        lt = lt_by_pmc[p][det]
        if not lt or lt <= 0 or str(lt) == "nan":
            continue
        try:
            msa = piquant_bridge.sfs_to_msa(
                raw, lt, axis, work / ("pmc%05d_%s.msa" % (p, det)),
                title="%s PMC%d %s" % (scan_dir.name, p, det))
        except (OSError, ValueError) as exc:
            fails += 1
            print("  [%d/%d] PMC %d MSA 导出失败：%s" % (n, len(sampled), p, exc))
            continue
        q, why = piquant_bridge.run_quant(msa, work, ELEMENT_LINES)
        if q is None:
            fails += 1
            print("  [%d/%d] PMC %d quant 失败：%s" % (n, len(sampled), p, why))
            continue
        r = rqc[p]

        def wt_of(col):
            try:
                return float(r.get(col + "_wt%", "") or "nan")
            except ValueError:
                return float("nan")

        row = {"pmc": p, "detector": det, "live_time_s": lt}
        for ox in PIQUANT_OXIDES:
            row["piquant_" + ox] = q["oxides"].get(ox)
            row["rqc_" + ox] = wt_of(ox)
        rows.append(row)
        if n % 20 == 0 or n == len(sampled):
            s = q["oxides"]
            print("  [%d/%d] PMC %d: SiO2 %.2f FeO-T %.2f CaO %.2f（sum %.1f%%）"
                  % (n, len(sampled), p, s.get("SiO2", -1), s.get("FeO-T", -1),
                     s.get("CaO", -1), q.get("element_sum_pct", -1)))
    if fails:
        print("失败 PMC：%d/%d（量化中断或导出失败的点不进对账）" % (fails, len(sampled)))
    ok_rows = [r for r in rows
               if all(r.get("piquant_" + ox) is not None for ox in ("SiO2", "FeO-T"))]
    if len(ok_rows) < 10:
        print("有效量化 <10 点——不足以对账")
        return 1

    # 4) 逐元素偏差统计（PIQUANT vs RQC，同 PMC）
    import numpy as np
    stats = {}
    for ox in PIQUANT_OXIDES:
        pv = np.array([r["piquant_" + ox] for r in ok_rows], dtype=float)
        rv = np.array([r["rqc_" + ox] for r in ok_rows], dtype=float)
        good = np.isfinite(pv) & np.isfinite(rv) & (rv > 0)
        if good.sum() < 10:
            continue
        rel = (pv[good] - rv[good]) / rv[good]
        corr = float(np.corrcoef(pv[good], rv[good])[0, 1])
        stats[ox] = {"n": int(good.sum()),
                     "median_rel_dev_pct": float(np.median(rel) * 100),
                     "mean_abs_rel_dev_pct": float(np.mean(np.abs(rel)) * 100),
                     "p90_abs_rel_dev_pct": float(np.percentile(np.abs(rel), 90) * 100),
                     "corr": corr,
                     "piquant_median_wt": float(np.median(pv[good])),
                     "rqc_median_wt": float(np.median(rv[good]))}
    print("\nPIQUANT(ENB, FM标定) vs 官方 RQC（同 PMC，正偏差=PIQUANT 偏高）:")
    for ox, s in sorted(stats.items(), key=lambda kv: kv[1]["mean_abs_rel_dev_pct"]):
        print("  %-5s n=%-4d 中位偏差 %+6.1f%%  平均|偏差| %5.1f%%  P90 %5.1f%%  相关 %.4f"
              % (ox, s["n"], s["median_rel_dev_pct"], s["mean_abs_rel_dev_pct"],
                 s["p90_abs_rel_dev_pct"], s["corr"]))

    if args.no_write:
        print("--no-write：缓存未改动")
        return 0

    # 5) 落盘
    stamp = UTC_NOW()
    # 中间文件清理：MSA 输入与 plot 曲线可由本脚本重新生成（全部结果行
    # 已进 piquant_vs_rqc.json），成功收尾后删除，避免缓存目录膨胀。
    import shutil
    shutil.rmtree(work, ignore_errors=True)
    doc = {
        "scope": "PIQUANT（PIXL 官方开源定量，github.com/pixlise/piquant 本地构建）对"
                 "同一份 SFS 原始谱做官方算法复现：ENB 谱 × FM 配置/FM 标定"
                 "（SurfaceOps Rev1 Jul2021，官方 COMMENT 指定 B 探测器）→ wt% 与"
                 "官方 RQC（A/B 择低逐 PMC 值）对账。两边口径差异（探测器选择、"
                 "标定代次）计入偏差，不代表任一方错误。与本仓库自研峰面积层"
                 "（sfs_vs_rqc.json）口径分开。",
        "piquant_bin_sha256": __import__("hashlib").sha256(
            Path(piquant_bridge.piquant_bin()).read_bytes()).hexdigest(),
        "scan_dir": scan_dir.name,
        "rqc": Path(rqc_p).name,
        "detector": args.detector,
        "stride": args.stride,
        "axis_keV_per_ch": axis["keV_per_ch"], "axis_zero_keV": axis["zero_keV"],
        "n_pmc_total": len(pmcs), "n_quantified": len(ok_rows), "n_failed": fails,
        "element_lines": ELEMENT_LINES,
        "stats": stats,
        "rows": ok_rows,
        "generated_utc": stamp,
    }
    CACHE.mkdir(parents=True, exist_ok=True)
    (CACHE / "piquant_vs_rqc.json").write_text(
        json.dumps(doc, ensure_ascii=False), encoding="utf-8")
    rep = ["# PIQUANT × RQC 对账（%s，%s）" % (scan_dir.name, stamp),
           "",
           "PIQUANT 官方开源定量（ENB × FM 标定）vs 官方 RQC，同 PMC 抽样 %d 点"
           "（步长 %d）。正偏差 = PIQUANT 偏高。" % (len(ok_rows), args.stride),
           "",
           "| 氧化物 | n | 中位偏差 | 平均\\|偏差\\| | P90 | 相关 |",
           "|---|---|---|---|---|---|"]
    for ox, s in sorted(stats.items(), key=lambda kv: kv[1]["mean_abs_rel_dev_pct"]):
        rep.append("| %s | %d | %+0.1f%% | %.1f%% | %.1f%% | %.4f |"
                   % (ox, s["n"], s["median_rel_dev_pct"], s["mean_abs_rel_dev_pct"],
                      s["p90_abs_rel_dev_pct"], s["corr"]))
    rep.append("")
    rep.append("FM 标定不覆盖 P2O5/MnO/Cr2O3/Sr（无 ECF 条目），该四项不在本层对账内；"
               "Na2O 属低能已知难点（表面效应+标定代次差异），偏差单独解读。")
    (CACHE / "piquant_vs_rqc_report.md").write_text("\n".join(rep), encoding="utf-8")
    print("\n缓存写入 %s（piquant_vs_rqc.json + 报告 md）" % CACHE)
    return 0


if __name__ == "__main__":
    sys.exit(main())
