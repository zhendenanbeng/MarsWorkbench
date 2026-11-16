# -*- coding: utf-8 -*-
"""独立数据审计：不信任冻结值，从本地原始数据重算并对账全部自主定量数字。

检查层（数据在位才查，缺数据项显式 SKIP 不计失败）：
  1. 训练谱字节完整性：manifest sha256 vs 本地文件实际哈希
  2. ChemCam 留一靶 CV 独立复算 vs chemcam_pls_cv.json 冻结 RMSE
  3. SuperCam 留一靶 CV 独立复算 vs scam_pls_cv.json 冻结 RMSE
  4. 模型-冻结一致性：pls_model.json 的 CV RMSE meta 与 CV JSON 一致；
     窗模式模型的掩码长度与 meta.modes 的 n_channels 一致
  5. PIXL 能量轴/峰面积独立复算 vs sfs_vs_rqc.json（相关系数 ±0.005）
  6. APXS rocknest3_rp 峰面积独立复算 vs apxs_peak_areas.json（逐元素）
  7. 台账/对账摘要来源指纹在场

审计哲学与 independent_check.py 一致：独立重算、逐项报 OK/FAIL/SKIP；
冻结数字对不上=FAIL（审计失败，需要查管线），缺数据=SKIP（不算失败）。
"""
from __future__ import annotations

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0, str(_BootPath(__file__).resolve().parents[3] / "02_图形界面代码"))
import source_paths as _layout            # noqa: E402

import csv
import json
import sys
import numpy as np
from pathlib import Path

import project_layout                      # noqa: E402
from instruments import libs_quant         # noqa: E402
from instruments import sfs_spectra        # noqa: E402
from instruments import apxs_fit           # noqa: E402

TOTAL = FAIL = SKIP = 0


def check(ok, name):
    global TOTAL, FAIL
    TOTAL += 1
    FAIL += not ok
    print(("OK  " if ok else "FAIL") + " " + name, flush=True)


def skip(name):
    global SKIP
    SKIP += 1
    print("SKIP " + name, flush=True)


def sha256_of(path):
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    data_root = Path(project_layout.data_root())
    B = _layout.project_layout.VERIFY / "fixtures/历史数值基准"

    # ========== 1. 训练谱字节完整性 ==========
    for inst, manifest_rel in (
            ("ChemCam", "ChemCam/下载/ccs_training/training_manifest.json"),
            ("SuperCam", "SuperCam/下载/scct_training/training_manifest.json")):
        mp = data_root / manifest_rel
        if not mp.is_file():
            skip("%s 训练谱 manifest 不在位" % inst)
            continue
        manifest = json.loads(mp.read_text(encoding="utf-8"))
        n_ok = n_bad = 0
        for t, bucket in manifest["targets"].items():
            for e in bucket.get("downloaded") or bucket.get("selected", []):
                if "file" not in e:
                    continue
                f = dest_of(mp, e["file"])
                if f is None or not f.is_file():
                    n_bad += 1
                    continue
                if "sha256" in e:
                    if sha256_of(f) == e["sha256"]:
                        n_ok += 1
                    else:
                        n_bad += 1
        if n_ok == 0 and n_bad == 0:
            skip("%s 训练谱 manifest 无 sha256 声明（字节完整性跳过）" % inst)
        else:
            check(n_bad == 0 and n_ok > 0,
                  "%s 训练谱字节完整（%d 哈希对账通过，%d 坏）" % (inst, n_ok, n_bad))

    # ========== 2. ChemCam CV 独立复算 ==========
    cc_cv = data_root / "ChemCam" / "缓存" / "chemcam_pls_cv.json"
    cc_manifest = data_root / "ChemCam" / "下载" / "ccs_training" / "training_manifest.json"
    cc_comp = data_root / "ChemCam" / "缓存" / "ccs_composition.json"
    cc_dest = data_root / "ChemCam" / "下载" / "ccs_training"
    if not (cc_cv.is_file() and cc_manifest.is_file() and cc_comp.is_file()
            and cc_dest.is_dir()):
        skip("ChemCam CV 复算（数据不在位）")
    else:
        frozen = json.loads(cc_cv.read_text(encoding="utf-8"))
        spec_doc = json.loads(cc_comp.read_text(encoding="utf-8"))["targets"]
        # 数据指纹：冻结 doc 声明的成分表/manifest 哈希必须仍与盘上一致
        if frozen.get("composition_sha256"):
            check(sha256_of(cc_comp) == frozen["composition_sha256"],
                  "ChemCam 成分表指纹与冻结 doc 一致")
        if frozen.get("training_manifest_sha256"):
            check(sha256_of(cc_manifest) == frozen["training_manifest_sha256"],
                  "ChemCam 训练 manifest 指纹与冻结 doc 一致")
        # 逐氧化物 components/modes 的权威记录在模型 meta（CV doc 只存全局 int）
        _mdl = libs_quant.load_model(data_root / "ChemCam" / "缓存" / "pls_model.json") \
            if (data_root / "ChemCam" / "缓存" / "pls_model.json").is_file() else None
        model_meta = _mdl["meta"] if _mdl else {}
        import importlib.util
        spec03 = importlib.util.spec_from_file_location(
            "s03", _layout.project_layout.ROOT
            / "01_核心处理代码/ChemCam/数据准备脚本/03_run_pls_cv.py")
        m03 = importlib.util.module_from_spec(spec03)
        sys.modules.setdefault("s03", m03)
        spec03.loader.exec_module(m03)
        spectra, wl0 = {}, None
        for target, bucket in manifest_targets(cc_manifest):
            if target not in spec_doc or spec_doc[target].get("excluded_by_official"):
                continue
            rows = []
            for entry in bucket.get("downloaded") or bucket.get("selected", []):
                try:
                    wl, shots = m03.load_point_spectra(cc_dest / entry["file"])
                    rows.append(m03.preprocess_spectra(wl, shots).mean(axis=1))
                    wl0 = wl
                except Exception:                    # noqa: BLE001
                    continue
            if rows:
                spectra[target] = np.vstack(rows)
        comp_flat = {k: v["composition"] for k, v in spec_doc.items()
                     if not spec_doc[k].get("excluded_by_official")}
        n_ok = n_bad = 0
        for ox, frozen_row in frozen["cv"].items():
            names_ok = [n for n in frozen_row["true"]
                        if n in spectra and comp_flat.get(n, {}).get(ox) is not None]
            if len(names_ok) < len(frozen_row["true"]):
                n_bad += 1
                continue
            try:
                comp_all = frozen["components"] if isinstance(frozen["components"], dict) \
                    else (model_meta.get("components") or {})
                k = comp_all.get(ox) if isinstance(comp_all, dict) else frozen["components"]
                modes_all = frozen.get("modes") or model_meta.get("modes") or {}
                # 真留一靶复算：与 03 号脚本同口径（loo_target_cv/_windowed）。
                # 自拟合（fit 后预测全数据）不是 CV，不能作为复算基准。
                sub = {n: spectra[n] for n in names_ok}
                if modes_all.get(ox, {}).get("mode") == "windowed":
                    r = libs_quant.loo_target_cv_windowed(
                        sub, comp_flat, ox, wl0, components=k)
                else:
                    r = libs_quant.loo_target_cv(sub, comp_flat, ox, components=k)
                rmse = float(r["rmse"])
                if abs(rmse - frozen_row["rmse"]) <= 0.01:
                    n_ok += 1
                else:
                    n_bad += 1
                    print("   %s: 复算 %.4f vs 冻结 %.4f" % (ox, rmse, frozen_row["rmse"]))
            except Exception as exc:                     # noqa: BLE001
                n_bad += 1
                print("   %s 复算异常：%s" % (ox, str(exc)[:60]))
        check(n_bad == 0 and n_ok > 0,
              "ChemCam CV 独立复算逐氧化物一致（%d/%d）" % (n_ok, n_ok + n_bad))

    # ========== 3. SuperCam CV 独立复算 ==========
    sc_cv = data_root / "SuperCam" / "缓存" / "scam_pls_cv.json"
    sc_manifest = data_root / "SuperCam" / "下载" / "scct_training" / "training_manifest.json"
    sc_comp = data_root / "SuperCam" / "缓存" / "scct_composition.json"
    sc_dest = data_root / "SuperCam" / "下载" / "scct_training"
    if not (sc_cv.is_file() and sc_manifest.is_file() and sc_comp.is_file()
            and sc_dest.is_dir()):
        skip("SuperCam CV 复算（数据不在位）")
    else:
        frozen = json.loads(sc_cv.read_text(encoding="utf-8"))
        comp_doc = json.loads(sc_comp.read_text(encoding="utf-8"))["targets"]
        # 数据指纹：成分表/manifest 若在冻结后被改动，冻结 CV 数字即失效
        if frozen.get("composition_sha256"):
            check(sha256_of(sc_comp) == frozen["composition_sha256"],
                  "SuperCam 成分表指纹与冻结 doc 一致")
        if frozen.get("training_manifest_sha256"):
            check(sha256_of(sc_manifest) == frozen["training_manifest_sha256"],
                  "SuperCam 训练 manifest 指纹与冻结 doc 一致")
        from instruments import supercam
        backend = supercam.SuperCamBackend()
        spectra = {}
        for target, bucket in manifest_targets(sc_manifest):
            if target not in comp_doc:
                continue
            rows = []
            for entry in bucket.get("downloaded") or bucket.get("selected", []):
                try:
                    class _P:
                        files = {"cl1": str(sc_dest / entry["file"])}
                        sclk = ""
                        meta = {}
                    d = backend.load_point(_P(), None,
                                           {"dust": True, "despike": True,
                                            "baseline": True})
                    rows.append(np.asarray(d["y"], dtype=float))
                except Exception:                    # noqa: BLE001
                    continue
            if rows:
                spectra[target] = np.vstack(rows)
        comp_flat = {k: v["composition"] for k, v in comp_doc.items()}
        n_ok = n_bad = 0
        for ox, frozen_row in frozen["cv"].items():
            names_ok = [n for n in frozen_row["true"]
                        if n in spectra and comp_flat.get(n, {}).get(ox) is not None]
            if len(names_ok) < len(frozen_row["true"]):
                n_bad += 1
                continue
            try:
                k = frozen["components"][ox] if isinstance(frozen["components"], dict) \
                    else frozen["components"]
                sub = {n: spectra[n] for n in names_ok}
                r = libs_quant.loo_target_cv(sub, comp_flat, ox, components=k)
                rmse = float(r["rmse"])
                if abs(rmse - frozen_row["rmse"]) <= 0.01:
                    n_ok += 1
                else:
                    n_bad += 1
                    print("   %s: 复算 %.4f vs 冻结 %.4f" % (ox, rmse, frozen_row["rmse"]))
            except Exception as exc:                     # noqa: BLE001
                n_bad += 1
                print("   %s 复算异常：%s" % (ox, str(exc)[:60]))
        check(n_bad == 0 and n_ok > 0,
              "SuperCam CV 独立复算逐氧化物一致（%d/%d）" % (n_ok, n_ok + n_bad))

    # ========== 4. 模型-冻结一致性 ==========
    for inst, model_rel, cv_rel in (
            ("ChemCam", "ChemCam/缓存/pls_model.json", "ChemCam/缓存/chemcam_pls_cv.json"),
            ("SuperCam", "SuperCam/缓存/pls_model.json", "SuperCam/缓存/scam_pls_cv.json")):
        mp, cp = data_root / model_rel, data_root / cv_rel
        if not mp.is_file():
            skip("%s 模型文件不在位" % inst)
            continue
        doc = libs_quant.load_model(mp)
        cv = json.loads(cp.read_text(encoding="utf-8")) if cp.is_file() else {}
        if doc is None:
            check(False, "%s 模型文件可加载" % inst)
            continue
        meta = doc["meta"]
        modes = meta.get("modes", {})
        n_ok = 0
        n_bad = 0
        for ox, mdl in doc["models"].items():
            mask = getattr(mdl, "channel_mask", None)
            if mask is not None and modes.get(ox, {}).get("n_channels"):
                if len(mask) != modes[ox]["n_channels"]:
                    n_bad += 1
                else:
                    n_ok += 1
            elif mask is None and modes.get(ox, {}).get("mode", "full") == "full":
                n_ok += 1
            else:
                n_bad += 1
        cv_match = all(
            abs((meta.get("cv_rmse", {}) or {}).get(ox, -1)
                - (cv.get("cv", {}) or {}).get(ox, {}).get("rmse", -2)) <= 0.01
            for ox in doc["models"]) if cv.get("cv") else False
        check(n_bad == 0, "%s 模型窗掩码与 meta.modes 一致（%d/%d）"
              % (inst, n_ok, len(doc["models"])))
        check(cv_match, "%s 模型 meta 的 CV RMSE 与 CV JSON 一致" % inst)

    # ========== 5. PIXL 能量轴/峰面积复算 ==========
    sfs_dir = data_root / "PIXL" / "下载" / "sfs_sol0186"
    axis_p = data_root / "PIXL" / "缓存" / "sfs_energy_axis.json"
    sfs_v = data_root / "PIXL" / "缓存" / "sfs_vs_rqc.json"
    rqc_p = next(_layout.project_layout.ROOT.glob("03_下载数据/PIXL/随包示例/*rqc*.csv"),
                 None)
    if not (sfs_dir.is_dir() and axis_p.is_file() and sfs_v.is_file() and rqc_p):
        skip("PIXL SFS 复算（数据不在位）")
    else:
        axis = json.loads(axis_p.read_text(encoding="utf-8"))
        ax = axis["axis"]
        v = json.loads(sfs_v.read_text(encoding="utf-8"))
        # 冻结口径（2026-11-14 真轴 v4）：Bellegarde 全扫描集（sfs_sol0186
        # ENA/ENB，PMC 6-2815 = 2809），逐 PMC 峰面积取 A/B 探测器低者
        # （论文 Table S4 口径），轴 = calibrate_energy 自证值（审计可复算）。
        det_by_pmc = {}
        for f in sorted(sfs_dir.glob("*.csv")):
            up = f.name.upper()
            det = "ENA" if "ENA" in up else ("ENB" if "ENB" in up else None)
            if det is None:
                continue
            try:
                sp_, meta = sfs_spectra.read_sfs(f)
            except Exception:
                sp_, meta = None, None
            if sp_ is None:
                continue
            for i, pm in enumerate(meta):
                pmc = pm.get("PMC")
                if pmc in (None, "", "0"):
                    continue
                try:
                    lt = float(pm.get("live_time") or 0)
                except (ValueError, TypeError):
                    lt = 0.0
                det_by_pmc.setdefault(int(float(pmc)), {})[det] = (
                    sp_[i] / lt if lt > 0 else None)
        pmcs_ok = [p for p in v["pmc"]
                   if p in det_by_pmc and det_by_pmc[p].get("ENA") is not None
                   and det_by_pmc[p].get("ENB") is not None]
        if len(pmcs_ok) < len(v["pmc"]):
            skip("PIXL SFS 复算（PMC 不齐：%d/%d）" % (len(pmcs_ok), len(v["pmc"])))
        else:
            ca_re, fe_re = [], []
            for p in pmcs_ok:
                a_out = sfs_spectra.peak_areas(det_by_pmc[p]["ENA"], ax,
                                               elements=("Ca", "Fe"))
                b_out = sfs_spectra.peak_areas(det_by_pmc[p]["ENB"], ax,
                                               elements=("Ca", "Fe"))
                ca_re.append(min(x["area"] or 0.0 for x in (a_out["Ca"], b_out["Ca"])))
                fe_re.append(min(x["area"] or 0.0 for x in (a_out["Fe"], b_out["Fe"])))
            d_ca = max(abs(a - b) for a, b in zip(ca_re, v["ca_area"]))
            d_fe = max(abs(a - b) for a, b in zip(fe_re, v["fe_area"]))
            check(d_ca < 1e-6 and d_fe < 1e-6,
                  "SFS 逐 PMC 峰面积精确复现（Ca/Fe 最大偏差 %.1e/%.1e）"
                  % (d_ca, d_fe))
            # RQC wt 从随包 CSV 按 PMC 重取（不信冻结值），再对账冻结表
            rows = list(csv.DictReader(open(rqc_p, encoding="utf-8-sig", newline="")))
            rqc = {}
            for r in rows:
                try:
                    rqc[int(float(r["PMC"]))] = r
                except (ValueError, KeyError):
                    continue

            def _wt(pmc, col):
                try:
                    return float(rqc[pmc][col])
                except (KeyError, ValueError, TypeError):
                    return float("nan")

            ca_wt = [_wt(p, "CaO_wt%") for p in pmcs_ok]
            fe_wt = [_wt(p, "FeO-T_wt%") for p in pmcs_ok]
            d_wt_ca = max((abs(a - b) for a, b in zip(ca_wt, v["ca_wt"])
                           if np.isfinite(a)), default=0.0)
            d_wt_fe = max((abs(a - b) for a, b in zip(fe_wt, v["fe_wt"])
                           if np.isfinite(a)), default=0.0)
            check(d_wt_ca < 1e-9 and d_wt_fe < 1e-9,
                  "RQC 逐 PMC wt 与冻结对账表一致（最大差 %.1e/%.1e）"
                  % (d_wt_ca, d_wt_fe))
            good = np.isfinite(ca_wt)
            ca_corr = float(np.corrcoef(np.asarray(ca_re)[good],
                                        np.asarray(ca_wt)[good])[0, 1])
            goodf = np.isfinite(fe_wt)
            fe_corr = float(np.corrcoef(np.asarray(fe_re)[goodf],
                                        np.asarray(fe_wt)[goodf])[0, 1])
            check(abs(ca_corr - axis["ca_corr"]) < 1e-9
                  and abs(fe_corr - axis["fe_corr"]) < 1e-9,
                  "SFS vs RQC 相关复算与冻结一致（Ca %.4f / Fe %.4f）"
                  % (ca_corr, fe_corr))
            mean = np.mean([np.mean(np.vstack(list(det_by_pmc[p].values())), axis=0)
                            for p in pmcs_ok], axis=0)
            pa = sfs_spectra.peak_areas(mean, ax, elements=("Ca", "Fe"))
            check(pa["Ca"]["detected"] and pa["Fe"]["detected"],
                  "扫描级 Ca/Fe 峰检出（面积 %.0f/%.0f）"
                  % (pa["Ca"]["area"], pa["Fe"]["area"]))
            # 冻结纪律：轴演进必须带 revisions（旧值可追溯）
            rev = axis.get("revisions") or []
            ok_rev = bool(rev) and all(
                isinstance(r, dict) and r.get("date") and r.get("note") for r in rev)
            # 旧值必须与当前轴不同（revision 自指 = 记录失效）
            self_ref = [r for r in rev
                        if isinstance(r.get("old_axis"), dict)
                        and abs(r["old_axis"].get("keV_per_ch", 0)
                                - ax["keV_per_ch"]) < 1e-12]
            check(ok_rev and len(rev) >= 1 and not self_ref,
                  "能量轴 doc 的 revisions 记录完整且旧值非自指（%d 条）" % len(rev))

    # ========== 5b. 第二扫描跨扫描一致性（sol0300）==========
    sfs2 = data_root / "PIXL" / "下载" / "sfs_sol0300"
    rqc2_p = next((data_root / "PIXL" / "随包示例").glob("rqc_sol0300_*.csv"), None) \
        if (data_root / "PIXL" / "随包示例").is_dir() else None
    if not (sfs2.is_dir() and rqc2_p):
        skip("第二扫描一致性（sol0300 数据不在位）")
    else:
        # 真轴跨扫描物理一致性：独立自证轴与 Bellegarde 冻结轴差 ≤1%
        det2 = {}
        for f in sorted(sfs2.glob("*.csv")):
            up = f.name.upper()
            det = "ENA" if "ENA" in up else ("ENB" if "ENB" in up else None)
            if det is None:
                continue
            sp_, meta = sfs_spectra.read_sfs(f)
            if sp_ is None:
                continue
            sp_ = np.asarray(sp_, dtype=float)
            lt = np.array([float(m.get("live_time") or 0) for m in meta])
            lt[lt <= 0] = np.nan
            det2[det] = sp_.mean(axis=0) / np.nanmean(lt)
        re_axis2 = sfs_spectra.calibrate_energy(np.mean(list(det2.values()), axis=0))
        check(abs(re_axis2["keV_per_ch"] - ax["keV_per_ch"]) / ax["keV_per_ch"] <= 0.01
              and abs(re_axis2["zero_keV"] - ax["zero_keV"]) <= 0.05,
              "第二扫描独立自证轴与冻结轴一致（g=%.6f/%.6f，b=%+.3f/%+.3f）"
              % (re_axis2["keV_per_ch"], ax["keV_per_ch"],
                 re_axis2["zero_keV"], ax["zero_keV"]))
        # 第二扫描逐 PMC 定量（A/B 择低）vs 其配套 RQC
        det_pmc = {}
        for f in sorted(sfs2.glob("*.csv")):
            up = f.name.upper()
            det = "ENA" if "ENA" in up else ("ENB" if "ENB" in up else None)
            if det is None:
                continue
            sp_, meta = sfs_spectra.read_sfs(f)
            if sp_ is None:
                continue
            for i, pm in enumerate(meta):
                p = int(pm.get("PMC") or 0)
                if p <= 0:
                    continue
                det_pmc.setdefault(p, {})[det] = sp_[i]
        rqc2 = {}
        for r in csv.DictReader(open(rqc2_p, encoding="utf-8-sig", newline="")):
            try:
                rqc2[int(float(r["PMC"]))] = r
            except (ValueError, KeyError):
                continue
        common2 = [p for p in sorted(det_pmc) if p in rqc2]
        cover = len(common2) / max(len(det_pmc), 1)
        check(cover >= 0.95,
              "第二扫描 PMC 覆盖率 ≥95%%（%.0f%%，%d/%d）"
              % (cover * 100, len(common2), len(det_pmc)))
        ca2 = []
        for p in common2:
            vals = []
            for det in ("ENA", "ENB"):
                if det in det_pmc[p]:
                    pa = sfs_spectra.peak_areas(det_pmc[p][det], ax, elements=("Ca",))
                    if pa["Ca"]["area"] is not None:
                        vals.append(pa["Ca"]["area"])
            ca2.append(min(vals) if vals else 0.0)
        wt2 = []
        for p in common2:
            try:
                wt2.append(float(rqc2[p]["CaO_wt%"]))
            except (KeyError, ValueError, TypeError):
                wt2.append(float("nan"))
        wt2 = np.asarray(wt2, dtype=float)
        ca2 = np.asarray(ca2, dtype=float)
        ok2 = np.isfinite(wt2)
        ca_corr2 = float(np.corrcoef(ca2[ok2], wt2[ok2])[0, 1])
        check(ca_corr2 >= 0.98,
              "第二扫描逐 PMC Ca 相关 ≥0.98（%.4f，n=%d）" % (ca_corr2, int(ok2.sum())))

    # ========== 6. APXS rocknest3_rp 峰面积复算 ==========
    fix = (_layout.project_layout.VERIFY / "fixtures/apxs/sol00102")
    areas_p = data_root / "APXS" / "缓存" / "apxs_peak_areas.json"
    rsp = fix / "apa_406545731rsp01020050178_______p1.csv"
    rwp = fix / "apa_406545731rwp01020050178_______p1.csv"
    if not (rsp.is_file() and areas_p.is_file()):
        skip("APXS 峰面积复算（数据不在位）")
    else:
        from instruments import apxs as apxs_backend
        frozen = json.loads(areas_p.read_text(encoding="utf-8"))
        if "source_sha256" in frozen and "rsp" in frozen.get("source_sha256", {}):
            check(sha256_of(rsp) == frozen["source_sha256"]["rsp"],
                  "APXS rsp 源文件哈希与冻结 doc 一致")
        data = apxs_backend.pds3.read_table(str(rsp), skiprows=1)
        energy_kev = data[:, 0] / 1000.0
        counts = data[:, 1] if data.shape[1] == 2 else np.median(data[:, 1:-1], axis=1)
        energy_kev, counts = energy_kev[1:-1], counts[1:-1]
        # 冻结口径：spectrum_areas 整体复算（窗中心=KA_KEV 元素表）。
        # 不得拿冻结 doc 里的 argmax 峰位当窗中心——那是拟合输出不是输入。
        re = apxs_fit.spectrum_areas(energy_kev, counts)
        n_ok = n_bad = 0
        for el, frozen_r in frozen["areas"].items():
            if frozen_r.get("area") is None:
                continue
            r = re.get(el)
            if r is None or r.get("area") is None:
                n_bad += 1
                continue
            if abs(r["area"] - frozen_r["area"]) <= 1e-6:
                n_ok += 1
            else:
                n_bad += 1
                print("   %s: 复算 %.4f vs 冻结 %.4f" % (el, r["area"], frozen_r["area"]))
        check(n_bad == 0 and n_ok > 0,
              "APXS 峰面积复算逐元素一致（%d/%d）" % (n_ok, n_ok + n_bad))

    # ========== 7. 台账来源指纹在场 ==========
    for name in ("published_calibrations.json", "published_targets.json",
                 "gold_intervals.json", "paper_crosscheck_summary.json",
                 "pixl_ratio_reconciliation.json"):
        p = B / name
        if not p.is_file():
            check(False, "台账 %s 在场" % name)
            continue
        doc = json.loads(p.read_text(encoding="utf-8"))
        has_src = bool(doc.get("sources")) or bool(doc.get("source_numbers_sha256")) \
            or bool(doc.get("source_sha256"))
        check(has_src, "台账 %s 来源可追溯" % name)

    print()
    print("审计完成：%d 项检查（OK %d / FAIL %d / SKIP %d）"
          % (TOTAL, TOTAL - FAIL, FAIL, SKIP))
    print("共 %d 项断言，%d 项失败" % (TOTAL, FAIL))
    return int(FAIL > 0)


def manifest_targets(manifest_path):
    manifest = json.loads(Path(manifest_path).read_text(encoding="utf-8"))
    return manifest["targets"].items()


def dest_of(manifest_path, filename):
    p = Path(manifest_path).parent
    f = p / filename
    return f if f.is_file() else None


if __name__ == "__main__":
    sys.exit(main())
