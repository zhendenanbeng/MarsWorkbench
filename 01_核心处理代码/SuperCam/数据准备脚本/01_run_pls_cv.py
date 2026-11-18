# -*- coding: utf-8 -*-
"""SuperCam 数据准备脚本 01：SCCT 自研 PLS 留一靶 CV + 最终模型（Phase 2）。

输入：03_下载数据/SuperCam/缓存/scct_composition.json（Cousin 2022 转录）
      03_下载数据/SuperCam/下载/scct_training/training_manifest.json（校准谱）
输出：scam_pls_cv.json（冻结 CV 数字）+ --fit-model 生成 pls_model.json
口径：与 ChemCam 03 完全同构；对账锚点 = Ytsma 2023 SuperCam Mars-based RMSE。
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
from instruments import libs_quant         # noqa: E402
from instruments import supercam  # noqa: E402

B = _layout.project_layout.VERIFY / "fixtures/历史数值基准"


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def mean_spectrum(fits_path, backend):
    """读一份校准谱 FITS：全点平均后的预处理谱（与画图口径同构）。"""
    class _P:                                        # 最小 PointFile 桩
        files = {"cl1": str(fits_path)}
        sclk = ""
        meta = {}
    pt = _P()
    d = backend.load_point(pt, None, {"dust": True, "despike": True, "baseline": True})
    return np.asarray(d["y"], dtype=float), np.asarray(d["x"], dtype=float)


def main(argv=None):
    ap = argparse.ArgumentParser(description="SuperCam 自研 PLS 留一靶 CV")
    ap.add_argument("--components", type=int, default=3)
    ap.add_argument("--auto-components", action="store_true",
                    default=True, help="每氧化物按留一靶 CV 自动选成分数")
    ap.add_argument("--fit-model", action="store_true")
    args = ap.parse_args(argv)

    data_root = Path(project_layout.data_root())
    comp = json.loads((data_root / "SuperCam" / "缓存" / "scct_composition.json")
                      .read_text(encoding="utf-8"))["targets"]
    manifest = json.loads((data_root / "SuperCam" / "下载" / "scct_training"
                           / "training_manifest.json").read_text(encoding="utf-8"))
    dest = data_root / "SuperCam" / "下载" / "scct_training"

    backend = supercam.SuperCamBackend()
    spectra = {}
    for target, bucket in manifest["targets"].items():
        if target not in comp:
            continue
        rows = []
        for entry in bucket["downloaded"]:
            try:
                y, _x = mean_spectrum(dest / entry["file"], backend)
                rows.append(y)
            except Exception as exc:                     # noqa: BLE001
                print("跳过 %s：%s" % (entry["file"], str(exc)[:60]))
        if rows:
            spectra[target] = np.vstack(rows)
            print("%-12s %d 条平均谱（%d 通道）" % (target, len(rows), len(rows[0])))
    if len(spectra) < 3:
        print("可用靶不足 3 个。")
        return 1

    # 成分不完整的靶按氧化物排除（如 Calcite 无 Na2O/K2O）
    cv_rows = {}
    chosen = {}
    for ox in libs_quant.OXIDES:
        spectra_ox = {n: s for n, s in spectra.items()
                      if comp.get(n, {}).get("composition", {}).get(ox) is not None}
        if len(spectra_ox) < 3:
            print("  %-6s 跳过：完整成分靶不足 3 个" % ox)
            continue
        try:
            if args.auto_components:
                chosen_sc = libs_quant.auto_components(
                    spectra_ox,
                    {n: v["composition"] for n, v in comp.items()}, ox,
                    max_components=min(4, len(spectra_ox) - 2))[0]
            else:
                chosen_sc = args.components
            chosen[ox] = chosen_sc
            cv_rows[ox] = libs_quant.loo_target_cv(
                spectra_ox, {n: v["composition"] for n, v in comp.items()}, ox,
                components=chosen_sc)
            print("  %-6s RMSE %.3f  bias %+.3f  (n=%d)"
                  % (ox, cv_rows[ox]["rmse"], cv_rows[ox]["bias"], cv_rows[ox]["n_targets"]))
        except ValueError as exc:
            print("  %-6s 跳过：%s" % (ox, exc))

    gi = json.loads((B / "gold_intervals.json").read_text(encoding="utf-8"))
    ytsma = gi["anchors"].get("supercam", {}).get("moc_rmsep", {}).get("oxides", {})
    verdict = {}
    for ox, r in cv_rows.items():
        band = ytsma.get(ox)
        if band:
            verdict[ox] = {"in_band": bool(band["lower"] <= r["rmse"] <= band["upper"]),
                           "band": [band["lower"], band["upper"]],
                           "reference": band["reference"]}
            print("  对照 %s: %.3f vs 官方MOC RMSEP %.1f（区间 %s）"
                  % (ox, r["rmse"], band["reference"], verdict[ox]["band"]))
    n_in = sum(1 for v in verdict.values() if v["in_band"])
    n_all = len(verdict)

    if args.fit_model and cv_rows:
        names = sorted(spectra)
        final = {}
        for ox in cv_rows:
            ok_names = [n for n in names
                        if comp.get(n, {}).get("composition", {}).get(ox) is not None]
            Xo = np.vstack([np.mean(np.asarray(spectra[n], dtype=float), axis=0)
                            for n in ok_names])
            yo = np.array([float(comp[n]["composition"][ox]) for n in ok_names])
            # 与 CV 同口径：每氧化物用其 auto-chosen 成分数训练
            final[ox] = libs_quant.PLS1Model(
                components=chosen.get(ox, args.components)).fit(Xo, yo)
        comp_path = data_root / "SuperCam" / "缓存" / "scct_composition.json"
        meta = {"components": dict(chosen), "trained_targets": names,
                "composition_sha256": sha256_of(comp_path),
                "cv_rmse": {ox: r["rmse"] for ox, r in cv_rows.items()},
                "trained_utc": datetime.datetime.now(datetime.timezone.utc)
                .strftime("%Y-%m-%dT%H:%M:%SZ")}
        mpath = data_root / "SuperCam" / "缓存" / "pls_model.json"
        libs_quant.save_model(final, mpath, meta)
        print("最终模型保存:", mpath)

    out = data_root / "SuperCam" / "缓存" / "scam_pls_cv.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    doc = {"scope": "SuperCam 自研 PLS 留一靶 CV 冻结测量（在轨 SCCT 最小训练集，"
                    "独立复现口径；精度预期与官方 Anderson 2022 334 标样不可比）。",
           "generated_utc": datetime.datetime.now(datetime.timezone.utc)
           .strftime("%Y-%m-%dT%H:%M:%SZ"),
           "components": dict(chosen), "n_targets": len(spectra),
           "targets": sorted(spectra),
           "cv": {ox: {"rmse": r["rmse"], "bias": r["bias"],
                       "predictions": r["predictions"], "true": r["true"]}
                  for ox, r in cv_rows.items()},
           "moc_rmsep_verdict": verdict,
           "moc_rmsep_verdict_note": "区间判定 %d/%d 落带内（锚点=Ytsma 2023 Table 1 转录，"
                                     "SuperCam PDS 列=Anderson 2022 Table 5 RMSEP，±60%% 带）；"
                                     "带外氧化物如实记录——%d 靶训练集精度弱于官方 334 标样为预期，"
                                     "不作官方替代。" % (n_in, n_all, len(spectra)),
           "composition_sha256": sha256_of(data_root / "SuperCam" / "缓存"
                                           / "scct_composition.json"),
           "training_manifest_sha256": sha256_of(data_root / "SuperCam" / "下载"
                                                 / "scct_training"
                                                 / "training_manifest.json"),
           "code_sha256": project_layout.source_hashes(include_verification=False)}
    out.write_text(json.dumps(doc, ensure_ascii=False, indent=1), encoding="utf-8")
    print("冻结:", out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
