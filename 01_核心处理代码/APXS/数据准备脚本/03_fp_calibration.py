#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""APXS 数据准备脚本 03：基本参数法-lite 标定（25 apb + 40 apa 自然靶）。

流程：
  1. 数据集构建：本地 calib_targets 的 25 个 apb 靶（01 号脚本遗产，滤
     live-time 异常）+ PDS 均匀采样下载的 apa 常规科学靶（默认 40 个，
     覆盖任务时代）——每靶 rsp 能谱（cps 归一）+ 配对 rwp 官方 16 行表。
  2. 逐靶自研峰面积（apxs_fit.spectrum_areas 全 14 元素）。
  3. 灵敏度 u_i：u_i = area_i·(μ/ρ)_m(E_i)/C_el_i（基体含化学计量 O，
     自算元素口径）；逐靶计算后取中位数（稳健）。70% 拟合 / 30% 留出。
  4. 留出验证：Sherman-lite 迭代预测氧化物 wt% vs 官方，逐氧化物中位
     相对误差；同靶无吸收修正基线（μ≡1）对照，量化吸收修正增量。
  5. --freeze 时把 u_i + 留出界写进 instruments/apxs_fp_calibration.json
     （运行期引擎的离线标定）+ 缓存/apxs_fp_calib.json（全量台账）。

用法：python 03_fp_calibration.py [--n-apa 40] [--no-write] [--freeze]
"""
from __future__ import annotations

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0, str(_BootPath(__file__).resolve().parents[3] / "02_图形界面代码"))
import source_paths as _layout            # noqa: E402

import argparse
import csv
import datetime
import hashlib
import json
import os
import re
import sys
import urllib.request
from pathlib import Path

import numpy as np

import project_layout                      # noqa: E402
from instruments import apxs_fit           # noqa: E402
from instruments import pds3               # noqa: E402
from instruments import apxs_fp            # noqa: E402

BASE = "https://pds-geosciences.wustl.edu/msl/msl-m-apxs-4_5-rdr-v1/mslapx_1xxx/data/"
UA = {"User-Agent": "sherloc-workbench/1.0"}

# 峰面积质量过滤：live time 与 SiO2 合理域（APXS 积分几十秒~几小时）
LIVE_MIN, LIVE_MAX = 10.0, 1e6
SIO2_MIN, SIO2_MAX = 20.0, 80.0


def http_get(url, retries=3):
    last = None
    for _ in range(retries):
        try:
            return urllib.request.urlopen(urllib.request.Request(url, headers=UA),
                                          timeout=180).read()
        except Exception as exc:                         # noqa: BLE001
            last = exc
    raise last


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_spectrum(rsp_path):
    """rsp → (energy_kev, cps 谱, live_time)；首通道=live time 剔除。"""
    data = pds3.read_table(str(rsp_path), skiprows=1)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError("谱表结构异常")
    live = float(data[0, 1])
    energy = data[:, 0] / 1000.0
    counts = data[:, 1] if data.shape[1] == 2 else np.median(data[:, 1:-1], axis=1)
    energy, counts = energy[1:-1], counts[1:-1]
    if not (LIVE_MIN <= live <= LIVE_MAX):
        raise ValueError("live time 异常 %.1f s" % live)
    return energy, counts / live, live


def load_official(rwp_path):
    """rwp 全 16 行 → {氧化物: wt}（官方 ver2 口径，含统计误差列忽略）。"""
    rows = list(csv.reader(open(rwp_path, encoding="utf-8", errors="replace")))
    oxides = {}
    for r in rows[1:]:
        if len(r) >= 2:
            try:
                oxides[r[0].strip()] = float(r[1])
            except ValueError:
                continue
    return oxides


def _pair_identity(fname):
    """与后端 apxs._pair_identity 同构：去 rsp/rwp 类型段与版本段。"""
    stem = fname.rsplit("/", 1)[-1][:-4]
    s = re.sub(r"r(?:sp|wp)", "", stem.lower(), count=1)
    return re.sub(r"_p\d+$", "", s)


def _file_version(fname):
    m = re.search(r"_p(\d+)\.csv$", fname.lower())
    return int(m.group(1)) if m else 0


def scan_all_apa_sols(sols):
    """并发扫描全部 sol 目录 → {sol: [(rsp_url, rwp_url), ...]}。

    wustl 按下传 sol 组织目录且 bundle 根无索引，逐目录并发扫是唯一
    确定性全量口径；只记录有 apa rsp+rwp 配对的目录。
    """
    from concurrent.futures import ThreadPoolExecutor

    def probe(sol):
        try:
            h2 = http_get(BASE + "sol%s/" % sol)
        except Exception:                                # noqa: BLE001
            return None
        h2 = h2.decode("utf-8", "replace")
        rsp_list = sorted(set(re.findall(r'href="(/msl[^"]*apa[^"]*rsp[^"]*\.csv)"',
                                         h2, re.I)))
        if not rsp_list:
            return None
        rwp_map = {}
        for u in set(re.findall(r'href="(/msl[^"]*apa[^"]*rwp[^"]*\.csv)"',
                                h2, re.I)):
            key = _pair_identity(u)
            if key not in rwp_map or _file_version(u) > _file_version(rwp_map[key]):
                rwp_map[key] = u
        pairs = []
        for u in rsp_list:
            key = _pair_identity(u)
            if key in rwp_map:
                pairs.append((u, rwp_map[key]))
        return pairs or None

    out = {}
    with ThreadPoolExecutor(max_workers=16) as pool:
        for sol, pairs in zip(sols, pool.map(probe, sols)):
            if pairs:
                out[sol] = pairs
    return out


def download_apa_targets(n_apa, dest):
    """全量扫描 apa 配对 → 均匀采样下载（确定性、断点续传）。

    每 sol 取第一个配对（与 01 号脚本"每 sol 一靶"惯例一致），
    均匀采样覆盖任务时代（探测器退化随时代漂移的覆盖要求）。
    """
    dest.mkdir(parents=True, exist_ok=True)
    h = http_get(BASE).decode("utf-8", "replace")
    sols = sorted(set(re.findall(r'sol(\d{5})', h)))
    if not sols:
        raise RuntimeError("PDS 目录列表为空")
    apa_map = scan_all_apa_sols(sols)
    apa_sols = sorted(apa_map)
    print("全量扫描：%d 个 sol 目录含 apa 配对" % len(apa_sols))
    if not apa_sols:
        return []
    idx = np.unique(np.linspace(0, len(apa_sols) - 1, n_apa).astype(int))
    got = []
    for i in idx:
        if len(got) >= n_apa:
            break
        sol = apa_sols[int(i)]
        pair = apa_map[sol][0]
        rsp_f = pair[0].split("/")[-1]
        rwp_f = pair[1].split("/")[-1]
        local_rsp, local_rwp = dest / rsp_f, dest / rwp_f
        for url_tail, local in (pair[0], local_rsp), (pair[1], local_rwp):
            if not (local.is_file() and local.stat().st_size > 0):
                try:
                    data = http_get("https://pds-geosciences.wustl.edu" + url_tail)
                    open(str(local) + ".part", "wb").write(data)
                    os.replace(str(local) + ".part", str(local))
                except Exception as exc:                 # noqa: BLE001
                    print("  下载失败 %s：%s" % (local.name, str(exc)[:60]))
        if local_rsp.is_file() and local_rwp.is_file():
            got.append({"sol": sol, "prefix": "apa",
                        "rsp": local_rsp, "rwp": local_rwp})
    return got


def sol_from_filename(fname):
    """文件名第二数字段前 4 位 = 观测 sol（如 apb_453330389rsp06290311330
    → 0629；apa_597995182rsp22590731206 → 2259；与目录名核验一致）。"""
    m = re.match(r"ap[ab]_\d{9}rsp(\d{4})", fname.lower())
    return m.group(1) if m else None


def collect_targets(n_apa):
    """数据集 = 本地 apb（滤异常）+ 下载 apa；返回 [{sol,prefix,rsp,rwp}]。"""
    dr = Path(project_layout.data_root())
    targets = []
    apb_dir = dr / "APXS" / "下载" / "calib_targets"
    if apb_dir.is_dir():
        rwp_map = {}
        for p in apb_dir.glob("*rwp*.csv"):
            rwp_map[p.name.replace("rwp", "rsp", 1)] = p
        for rsp in sorted(apb_dir.glob("*rsp*.csv")):
            rwp = rwp_map.get(rsp.name)
            if rwp is None:
                continue
            sol = sol_from_filename(rsp.name)
            if sol is None:
                continue
            targets.append({"sol": sol, "prefix": rsp.name[:3],
                            "rsp": rsp, "rwp": rwp})
    apa_dir = dr / "APXS" / "下载" / "calib_targets_apa"
    targets += download_apa_targets(n_apa, apa_dir)
    return targets


def build_dataset(targets):
    """逐靶峰面积 + 官方表 + 质量过滤 → dataset 行列表 + 拒收台账。"""
    dataset, rejected = [], []
    for t in targets:
        try:
            energy, cps, live = load_spectrum(t["rsp"])
            official = load_official(t["rwp"])
        except Exception as exc:                         # noqa: BLE001
            rejected.append({"file": t["rsp"].name, "reason": str(exc)[:60]})
            continue
        sio2 = official.get("SiO2")
        if sio2 is None or not (SIO2_MIN <= sio2 <= SIO2_MAX):
            rejected.append({"file": t["rsp"].name,
                             "reason": "SiO2 越域 %s" % sio2})
            continue
        areas_raw = apxs_fit.spectrum_areas(energy, cps,
                                            elements=apxs_fp.FP_ELEMENTS)
        areas = {el: r["area"] for el, r in areas_raw.items()
                 if r["area"] is not None}
        if "Si" not in areas:
            rejected.append({"file": t["rsp"].name, "reason": "Si 峰未检出"})
            continue
        dataset.append({"sol": t["sol"], "prefix": t["prefix"],
                        "live_s": round(live, 1),
                        "sha256": sha256_of(str(t["rsp"])),
                        "areas": {k: round(v, 4) for k, v in areas.items()},
                        "official": official})
    return dataset, rejected


def matrix_from_official(official):
    """官方氧化物 → 基体（元素+O 质量分数）；缺行按 0 处理（自算口径对齐）。"""
    els, o_wt = apxs_fp.elements_from_oxides(official)
    comp = {k: v / 100.0 for k, v in els.items()}
    comp["O"] = o_wt / 100.0
    return comp, els


def fit_sensitivities(dataset, fit_idx):
    """u_i = median_t[area_i·μ_m,t(E_i)/C_el_i,t]（官方基体，逐靶）。

    同时给出无吸收修正基线灵敏度 sens0（μ≡1 口径）用于增量对照。
    """
    sens, sens0, per_target = {}, {}, {}
    for el in apxs_fp.FP_ELEMENTS:
        e_kev = apxs_fp._tables()["line_energies_kev"][el]
        uvals, u0vals = [], []
        for i in fit_idx:
            rec = dataset[i]
            ox = rec["official"].get(apxs_fp.OXIDE_OF[el][0])
            areas_c = apxs_fp.correct_kbeta(rec["areas"])   # 与 solve 同口径
            a = areas_c.get(el)
            if not ox or not a or a <= 0:
                continue
            comp, els = matrix_from_official(rec["official"])
            if els.get(el, 0) <= 0:
                continue
            uvals.append(a * apxs_fp.matrix_mass_atten(comp, e_kev) / els[el])
            u0vals.append(a / els[el])
        if len(uvals) < 8:
            continue
        med = float(np.median(uvals))
        sens[el] = med
        sens0[el] = float(np.median(u0vals))
        per_target[el] = {"n": len(uvals),
                          "iqr_rel": float((np.percentile(uvals, 75)
                                            - np.percentile(uvals, 25)) / med)}
    return sens, sens0, per_target


def predict_record(rec, sens, sens0, absorb=True):
    """单靶预测氧化物。absorb=True：Sherman-lite 迭代；False：无修正基线。"""
    if absorb:
        return apxs_fp.solve(rec["areas"], sens)
    w = {el: rec["areas"][el] / sens0[el] for el in sens0
         if rec["areas"].get(el, 0) > 0}
    if "Si" not in w:
        return {"available": False, "reason": "Si 缺"}
    return {"available": True, "oxides": apxs_fp.oxides_from_elements(w),
            "iterations": 0, "converged": True}


def evaluate(dataset, val_idx, sens, sens0):
    """留出验证：FP（Sherman-lite）vs 无修正基线，逐氧化物相对误差。"""
    def err_rows(idx_list, absorb):
        errs, iter_log = {}, []
        for i in idx_list:
            rec = dataset[i]
            pred = predict_record(rec, sens, sens0, absorb)
            if not pred.get("available"):
                continue
            iter_log.append((rec["sol"], pred.get("iterations"),
                             pred.get("converged")))
            for ox, pv in pred["oxides"].items():
                tv = rec["official"].get(ox)
                if tv and tv > 0.2:                    # 0.2 wt% 以下不计误差
                    errs.setdefault(ox, []).append(abs(pv - tv) / tv)
        return ({ox: {"median_rel_err": float(np.median(v)),
                      "n": len(v)} for ox, v in errs.items()}, iter_log)

    fp, iter_log = err_rows(val_idx, True)
    base, _ = err_rows(val_idx, False)
    return fp, base, iter_log


def main(argv=None):
    ap = argparse.ArgumentParser(description="APXS FP-lite 标定与留出验证")
    ap.add_argument("--n-apa", type=int, default=40)
    ap.add_argument("--no-write", action="store_true")
    ap.add_argument("--freeze", action="store_true",
                    help="把 u_i+留出界写进 instruments/apxs_fp_calibration.json")
    args = ap.parse_args(argv)

    print("== 1. 数据集构建 ==")
    targets = collect_targets(args.n_apa)
    print("靶文件配对数:", len(targets))
    dataset, rejected = build_dataset(targets)
    print("有效靶:", len(dataset), "| 拒收:", len(rejected))
    for r in rejected[:8]:
        print("  拒收:", r["file"], "←", r["reason"])

    # 确定性切分：按 sol 排序 idx%3==2 留出（与 01/02 号脚本惯例一致）
    dataset.sort(key=lambda r: (r["sol"], r["prefix"]))
    n = len(dataset)
    fit_idx = [i for i in range(n) if i % 3 != 2]
    val_idx = [i for i in range(n) if i % 3 == 2]
    print("== 2. 灵敏度拟合 ==")
    print("靶 %d：拟合 %d / 留出 %d" % (n, len(fit_idx), len(val_idx)))
    sens, sens0, per_target = fit_sensitivities(dataset, fit_idx)
    for el, u in sorted(sens.items()):
        print("  %-3s u=%10.3f  IQR/中位=%.2f (n=%d)"
              % (el, u, per_target[el]["iqr_rel"], per_target[el]["n"]))

    print("== 3. 留出验证 ==")
    fp_err, base_err, iter_log = evaluate(dataset, val_idx, sens, sens0)
    bad_iter = [t for t in iter_log if not t[2]]
    print("  Sherman 迭代：%d 留出靶，未收敛 %d（最大迭代轮数见台账）"
          % (len(iter_log), len(bad_iter)))
    order = apxs_fp.OXIDE_ORDER
    print("  %-7s %12s %12s" % ("氧化物", "FP 中位误差", "无修正基线"))
    for ox in order:
        f_ = fp_err.get(ox, {}).get("median_rel_err")
        b_ = base_err.get(ox, {}).get("median_rel_err")
        if f_ is None and b_ is None:
            continue
        print("  %-7s %11.1f%%%12s"
              % (ox, 100 * f_ if f_ is not None else -1,
                 ("%11.1f%%" % (100 * b_)) if b_ is not None else "—"))

    # 比值层对照（同留出靶；氧化物 wt 比 / SiO2 口径，来自 02 号脚本冻结）
    rc_p = Path(project_layout.data_root()) / "APXS" / "缓存" / "apxs_ratio_calib.json"
    ratio_err = {}
    if rc_p.is_file():
        rc = json.loads(rc_p.read_text(encoding="utf-8"))
        cal = rc.get("elements") or {}
        errs = {}
        for i in val_idx:
            rec = dataset[i]
            a_si = rec["areas"].get("Si")
            si_ox = rec["official"].get("SiO2")
            if not a_si or not si_ox:
                continue
            for el, c in cal.items():
                a = rec["areas"].get(el)
                if not a:
                    continue
                pred = c["power_A"] * (a / a_si) ** c["power_b"]
                true_ox = rec["official"].get(apxs_fp.OXIDE_OF[el][0])
                if true_ox and true_ox > 0.2:
                    errs.setdefault(el, []).append(
                        abs(pred - true_ox / si_ox) / (true_ox / si_ox))
        ratio_err = {el: float(np.median(v)) for el, v in errs.items() if len(v) >= 3}

    # 边界靶诚实记录：rocknest3_rp（sol 102）在标定时代边缘稀疏区，
    # FP 偏差大——与 02 号脚本比值层同靶负结果一致，入档不隐藏
    boundary = []
    fix_dir = (_layout.project_layout.VERIFY / "fixtures/apxs/sol00102")
    fix_rsp = fix_dir / "apa_406545731rsp01020050178_______p1.csv"
    fix_rwp = fix_dir / "apa_406545731rwp01020050178_______p1.csv"
    if fix_rsp.is_file() and fix_rwp.is_file() and sens:
        try:
            energy_f, cps_f, live_f = load_spectrum(fix_rsp)
            areas_f = {el: r["area"] for el, r in
                       apxs_fit.spectrum_areas(energy_f, cps_f,
                                               elements=apxs_fp.FP_ELEMENTS).items()
                       if r["area"] is not None}
            pred_f = apxs_fp.solve(areas_f, sens)
            off_f = load_official(fix_rwp)
            if pred_f.get("available"):
                boundary.append({
                    "target": "rocknest3_rp", "sol": 102,
                    "note": "标定时代（sol 46-4216）边缘稀疏区靶；02 号脚本"
                            "比值层同靶亦为负结果（探测器效率 sol<600 漂移）。"
                            "边界声明：该靶 FP 预测不可引用。",
                    "fp_oxides": pred_f["oxides"], "official": off_f,
                    "rel_dev": {ox: round((pred_f["oxides"][ox] - off_f[ox]) / off_f[ox], 3)
                                for ox in pred_f["oxides"]
                                if off_f.get(ox) and off_f[ox] > 0.2},
                })
        except Exception as exc:                         # noqa: BLE001
            boundary.append({"target": "rocknest3_rp",
                             "note": "边界靶计算失败：%s" % str(exc)[:60]})

    out = {
        "scope": "APXS 基本参数法-lite 标定（吸收修正 Sherman-lite，u_i 逐元素灵敏度，"
                 "25 apb + %d apa 自然靶，70/30 切分）。预测为自算氧化物 wt%%（归一化 "
                 "100，含化学计量 O 基体）；对官方 rwp 留出验证。简化与边界见 "
                 "instruments/apxs_fp.py 模块头。" % sum(1 for r in dataset if r["prefix"] == "apa"),
        "generated_utc": datetime.datetime.now(datetime.timezone.utc)
                                 .strftime("%Y-%m-%dT%H:%M:%SZ"),
        "n_targets": n,
        "n_apb": sum(1 for r in dataset if r["prefix"] == "apb"),
        "n_apa": sum(1 for r in dataset if r["prefix"] == "apa"),
        "rejected": rejected,
        "boundary_targets": boundary,
        "sensitivities": sens,
        "sens_per_target": per_target,
        "holdout_fp": fp_err,
        "holdout_baseline_no_absorption": base_err,
        "holdout_ratio_layer": ratio_err,
        "n_fit": len(fit_idx), "n_val": len(val_idx),
        "dataset": dataset,
        "code_sha256": project_layout.source_hashes(include_verification=False),
    }

    dr = Path(project_layout.data_root())
    out_p = dr / "APXS" / "缓存" / "apxs_fp_calib.json"
    if not args.no_write:
        out_p.parent.mkdir(parents=True, exist_ok=True)
        tmp = str(out_p) + ".part"
        open(tmp, "w", encoding="utf-8").write(json.dumps(out, ensure_ascii=False, indent=1))
        os.replace(tmp, out_p)
        print("台账写入", out_p)
    if args.freeze:
        freeze_p = _BootPath(__file__).resolve().parents[1] / "instruments" / "apxs_fp_calibration.json"
        era_sols = sorted(int(r["sol"]) for r in dataset)
        frozen = {
            "scope": out["scope"],
            "generated_utc": out["generated_utc"],
            "calib_era": "sol %d-%d（%d apb + %d apa 靶）"
                         % (era_sols[0], era_sols[-1], out["n_apb"], out["n_apa"]),
            "sensitivities": sens,
            "holdout": {ox: {"median_rel_err": fp_err[ox]["median_rel_err"],
                             "n": fp_err[ox]["n"]} for ox in fp_err},
            "n_targets": n,
        }
        tmp = str(freeze_p) + ".part"
        open(tmp, "w", encoding="utf-8").write(json.dumps(frozen, ensure_ascii=False, indent=1))
        os.replace(tmp, freeze_p)
        print("运行期标定冻结 →", freeze_p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
