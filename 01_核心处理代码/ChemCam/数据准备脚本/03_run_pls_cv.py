#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ChemCam 数据准备脚本 03：自研 PLS 留一靶交叉验证 + 对账
==========================================================
输入：
  * 02 号脚本的 training_manifest.json（CCS 观测谱）
  * CCS 靶文献成分表（ccs_composition.json，从 Fabre 2011/Vaniman 2012
    人工转录；结构见 --make-template）
输出：
  * <数据根>/ChemCam/缓存/chemcam_pls_cv.json —— 逐氧化物 RMSE/bias/预测表
  * 与 gold_intervals.json 的 PLS_loo_rmse_ground 区间对照结论

口径：这是独立复现的冻结测量，不是方法学裁决；精度预期弱于官方
（官方含数百地面标样 + ICA + 距离校正）。数字只对当时代码/数据指纹负责。

用法：
    python 03_run_pls_cv.py --make-template   # 生成成分转录模板（填数值）
    python 03_run_pls_cv.py                   # 跑 CV
"""

from __future__ import annotations

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0, str(_BootPath(__file__).resolve().parents[3] / "02_图形界面代码"))
import source_paths as _layout            # noqa: E402

import argparse
import datetime
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

import project_layout                      # noqa: E402
from instruments import chemcam, libs_quant, base as cc_base   # noqa: E402

B = _layout.project_layout.VERIFY / "fixtures/历史数值基准"


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_point_spectra(csv_path):
    """读一个 CCS CSV：返回 (波长 nm 数组, 发谱矩阵 [通道 x 发])。

    与 chemcam.py load_point 的读取完全同构：第 0 列波长、末两列
    MEDIAN/AVERAGE 丢弃、中间每列一发。
    """
    data = chemcam.pds3.read_table(str(csv_path))
    if data.ndim != 2 or data.shape[1] < 3:
        raise ValueError("CSV 结构不对（%d 列）" % (data.shape[1] if data.ndim == 2 else -1))
    wl = data[:, 0]
    shots = data[:, 1:-2]
    return wl, shots


def preprocess_spectra(wl, shots):
    """与 chemcam.py _process 完全同构：除尘→去尖峰→逐发基线。

    返回 [通道 x 发] 矩阵（发数已除尘减少）。"""
    shots = cc_base.drop_dust(shots, 5)
    shots = cc_base.despike_shots(shots)
    cols = [cc_base.poly_baseline(wl, shots[:, j]) for j in range(shots.shape[1])]
    return np.column_stack(cols)


def main(argv=None):
    ap = argparse.ArgumentParser(description="ChemCam 自研 PLS 留一靶 CV")
    ap.add_argument("--make-template", action="store_true")
    ap.add_argument("--components", type=int, default=3)
    args = ap.parse_args(argv)

    data_root = Path(project_layout.data_root())
    manifest_path = data_root / "ChemCam" / "下载" / "ccs_training" / "training_manifest.json"
    comp_path = data_root / "ChemCam" / "缓存" / "ccs_composition.json"

    if args.make_template:
        template = {
            "scope": "CCS 靶文献成分转录表。数值逐靶逐氧化物填入；来源 Fabre 2011（玻璃）/"
                     "Vaniman 2012（陶瓷），页码/表号写进 locator。占位 null 表示未转录。",
            "targets": {name: {"composition": {ox: None for ox in libs_quant.OXIDES},
                               "locator": "Fabre 2011 Table ? / Vaniman 2012 Table ?"}
                        for name in sorted(json.loads(manifest_path.read_text(encoding="utf-8"))
                                           ["targets"]) if manifest_path.is_file()}
            if manifest_path.is_file() else {},
        }
        comp_path.parent.mkdir(parents=True, exist_ok=True)
        comp_path.write_text(json.dumps(template, ensure_ascii=False, indent=1), encoding="utf-8")
        print("成分模板写入", comp_path, "——逐靶逐氧化物填数值后重跑本脚本")
        return 0

    if not manifest_path.is_file():
        print("训练清单不存在：%s\n先运行 02_download_ccs_training.py" % manifest_path)
        return 1
    if not comp_path.is_file():
        print("成分表不存在：%s\n先 --make-template 生成模板并填数值" % comp_path)
        return 1
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    composition = json.loads(comp_path.read_text(encoding="utf-8"))["targets"]

    # --- 读谱、预处理、按靶聚合 ---
    spectra_by_target = {}
    missing_comp = []
    skipped_official = []
    include_excluded = "--include-excluded" in sys.argv
    for target, bucket in manifest["targets"].items():
        comp_entry = composition.get(target, {})
        comp = comp_entry.get("composition")
        if not comp or any(comp.get(ox) is None for ox in libs_quant.OXIDES):
            missing_comp.append(target)
            continue
        if comp_entry.get("excluded_by_official") and not include_excluded:
            skipped_official.append(target)
            continue
        rows = []
        for row in bucket.get("downloaded") or bucket.get("selected", []):
            try:
                wl, shots = load_point_spectra(Path(data_root / "ChemCam" / "下载" / "ccs_training") / row["file"])
                clean = preprocess_spectra(wl, shots)
                rows.append(clean.mean(axis=1))
            except Exception as exc:                     # noqa: BLE001
                print("跳过 %s：%s" % (row["file"], str(exc)[:60]))
        if rows:
            spectra_by_target[target] = np.vstack(rows)
    if missing_comp:
        print("成分未转录的靶（已跳过）：%s" % ", ".join(missing_comp))
    if skipped_official:
        print("官方校准排除的靶（--include-excluded 可纳入）：%s"
              % ", ".join(skipped_official))
    if len(spectra_by_target) < 3:
        print("可用靶不足 3 个，无法留一靶 CV；先补全成分表。")
        return 1

    print("参与 CV 的靶：%s（%d 个）" % (", ".join(sorted(spectra_by_target)), len(spectra_by_target)))
    cv_rows = {}
    for ox in libs_quant.OXIDES:
        try:
            cv_rows[ox] = libs_quant.loo_target_cv(spectra_by_target, {n: v["composition"] for n, v in composition.items()}, ox,
                                                   components=args.components)
            print("  %-6s RMSE %.3f  bias %+.3f  (n=%d)"
                  % (ox, cv_rows[ox]["rmse"], cv_rows[ox]["bias"], cv_rows[ox]["n_targets"]))
        except ValueError as exc:
            print("  %-6s 跳过：%s" % (ox, exc))

    # --- 与 gold_intervals 对照 ---
    gi = json.loads((B / "gold_intervals.json").read_text(encoding="utf-8"))
    ref = gi["anchors"]["chemcam"]["PLS_loo_rmse_ground"]["oxides"]
    verdict = {}
    for ox, r in cv_rows.items():
        band = ref.get(ox)
        verdict[ox] = (None if band is None else
                       {"in_band": bool(band["lower"] <= r["rmse"] <= band["upper"]),
                        "band": [band["lower"], band["upper"]],
                        "reference_rmse": band["reference"]})
        if verdict[ox]:
            print("  对照 %s: RMSE %.3f vs Wiens2012 %.2f（区间 %s）-> %s"
                  % (ox, r["rmse"], verdict[ox]["reference_rmse"], verdict[ox]["band"],
                     "在内" if verdict[ox]["in_band"] else "超出"))

    out = data_root / "ChemCam" / "缓存" / "chemcam_pls_cv.json"
    doc = {
        "scope": "自研 PLS 留一靶交叉验证的冻结测量。独立复现（在轨靶最小训练集），"
                 "精度预期弱于官方（Clegg 2017 含地面数据库）；不是方法学裁决，"
                 "不构成任何官方成分替代。数字只对本文件指纹负责。",
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "components": args.components,
        "n_targets": len(spectra_by_target),
        "targets": sorted(spectra_by_target),
        "cv": {ox: {"rmse": r["rmse"], "bias": r["bias"],
                    "predictions": r["predictions"], "true": r["true"]}
               for ox, r in cv_rows.items()},
        "gold_intervals_verdict": verdict,
        "training_manifest_sha256": sha256_of(manifest_path),
        "composition_sha256": sha256_of(comp_path),
        "code_sha256": project_layout.source_hashes(include_verification=False),
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(doc, ensure_ascii=False, indent=1), encoding="utf-8")
    print("冻结 CV 数字写入", out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
