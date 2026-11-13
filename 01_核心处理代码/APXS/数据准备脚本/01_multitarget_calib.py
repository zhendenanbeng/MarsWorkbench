# -*- coding: utf-8 -*-
"""APXS 数据准备脚本 01：多靶下载 → 峰面积-浓度经验标定（Phase 4 后续阶段）。

流程：
  1. PDS 列 sol 目录 → 随机采样 --n-sols 个
  2. 每 sol 下载 rsp+rwp 对（断点续传；rwp 缺失的 sol 跳过）
  3. 逐靶峰面积（apxs_fit.spectrum_areas）+ rwp 官方元素 wt
  4. 元素级经验标定：log(峰面积) ~ log(元素 wt) 线性拟合（留 30% 靶验证）
  5. 输出标定曲线 JSON（拟合参数+验证误差+全部数据点+指纹）
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
import random
import re
import sys
import urllib.request
from pathlib import Path

import numpy as np

import project_layout                      # noqa: E402
from instruments import apxs_fit           # noqa: E402
from instruments import apxs as apxs_backend  # noqa: E402

BASE = "https://pds-geosciences.wustl.edu/msl/msl-m-apxs-4_5-rdr-v1/mslapx_1xxx/data/"
UA = {"User-Agent": "sherloc-workbench/1.0"}
ELEM_OF = {"Na": ("Na2O", 0.742), "Mg": ("MgO", 0.603), "Al": ("Al2O3", 0.529),
           "Si": ("SiO2", 0.4667), "K": ("K2O", 0.830), "Ca": ("CaO", 0.7147),
           "Ti": ("TiO2", 0.5995), "Mn": ("MnO", 0.7745), "Fe": ("FeO", 1.0)}
ELEMENTS = tuple(ELEM_OF)


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


def load_rsp_rwp(rsp_path, rwp_path):
    data = apxs_backend.pds3.read_table(str(rsp_path), skiprows=1)
    energy_kev = data[:, 0] / 1000.0
    counts = data[:, 1] if data.shape[1] == 2 else np.median(data[:, 1:-1], axis=1)
    live_time_s = float(counts[0])          # 首通道=总积分时间（秒）
    energy_kev, counts = energy_kev[1:-1], counts[1:-1]
    counts = counts / max(live_time_s, 1e-6)   # 计数率 cps——跨靶可比的前提
    rows = list(csv.reader(open(rwp_path, encoding="utf-8", errors="replace")))
    oxides = {}
    for r in rows[1:]:
        if len(r) >= 2:
            try:
                oxides[r[0].strip()] = float(r[1])
            except ValueError:
                continue
    return energy_kev, counts, oxides


def main(argv=None):
    ap = argparse.ArgumentParser(description="APXS 多靶峰面积-浓度经验标定")
    ap.add_argument("--n-sols", type=int, default=25)
    ap.add_argument("--dest", default=None)
    args = ap.parse_args(argv)

    data_root = Path(project_layout.data_root())
    dest = Path(args.dest) if args.dest else data_root / "APXS" / "下载" / "calib_targets"
    dest.mkdir(parents=True, exist_ok=True)

    h = http_get(BASE).decode("utf-8", "replace")
    sols = sorted(set(re.findall(r'sol(\d{5})', h)))
    print("sol 目录:", len(sols))
    random.seed(7)
    sample = random.sample(sols, min(args.n_sols * 3, len(sols)))

    targets = []
    for sol in sample:
        if len(targets) >= args.n_sols:
            break
        try:
            h2 = http_get(BASE + "sol%s/" % sol).decode("utf-8", "replace")
        except Exception:                                # noqa: BLE001
            continue
        rsp = sorted(set(re.findall(r'href="(/msl[^"]*rsp[^"]*\.csv)"', h2, re.I)))
        rwp = sorted(set(re.findall(r'href="(/msl[^"]*rwp[^"]*\.csv)"', h2, re.I)))
        if not rsp or not rwp:
            continue
        rsp_f = rsp[0].split("/")[-1]
        rwp_f = rwp[0].split("/")[-1]
        local_rsp, local_rwp = dest / rsp_f, dest / rwp_f
        for url, local in (("https://pds-geosciences.wustl.edu" + rsp[0], local_rsp),
                           ("https://pds-geosciences.wustl.edu" + rwp[0], local_rwp)):
            if not (local.is_file() and local.stat().st_size > 0):
                try:
                    data = http_get(url)
                    open(str(local) + ".part", "wb").write(data)
                    import os as _os
                    _os.replace(str(local) + ".part", str(local))
                except Exception as exc:                 # noqa: BLE001
                    print("下载失败 %s：%s" % (local.name, str(exc)[:50]))
        if local_rsp.is_file() and local_rwp.is_file():
            targets.append((sol, local_rsp, local_rwp))
    print("收集 rsp+rwp 靶:", len(targets))

    # 逐靶峰面积 + 官方 wt
    dataset = []
    for sol, rsp_path, rwp_path in targets:
        try:
            energy_kev, counts, oxides = load_rsp_rwp(rsp_path, rwp_path)
        except Exception as exc:                         # noqa: BLE001
            print("跳过 sol%s：%s" % (sol, str(exc)[:50]))
            continue
        areas = apxs_fit.spectrum_areas(energy_kev, counts, elements=ELEMENTS)
        rec = {"sol": sol, "sha256": sha256_of(str(rsp_path)), "elements": {}}
        for el in ELEMENTS:
            ox, factor = ELEM_OF[el]
            a = areas[el].get("amp") or areas[el].get("area")
            if a is None or ox not in oxides:
                continue
            rec["elements"][el] = {"amp": areas[el].get("amp"),
                                   "area": areas[el].get("area"),
                                   "official_wt": oxides[ox]}
        if len(rec["elements"]) >= 6:
            dataset.append(rec)
    print("有效靶:", len(dataset))

    # 元素级经验标定：log(amp) ~ log(element wt)，70% 拟合 / 30% 验证
    out = {"scope": "APXS 峰面积-浓度经验标定（25 个自然靶，log-log 线性，留出验证）。"
                    "实测结论（2026-11-13）：cps 归一后各元素 log-log 斜率≈0——火星自然靶"
                    "组成范围太窄（SiO2 主量集中于 40-50 wt%），峰面积-浓度经验标定不可分辨；"
                    "官方定量依赖 408 实验室标样+基本参数法。面积层的有效用途是检出判定与"
                    "元素序对账（rocknest3_rp 已验证），不是浓度标定。基本参数法不在本层。",
           "generated_utc": datetime.datetime.now(datetime.timezone.utc)
           .strftime("%Y-%m-%dT%H:%M:%SZ"),
           "n_targets": len(dataset),
           "calibrations": {}, "dataset": dataset,
           "code_sha256": project_layout.source_hashes(include_verification=False)}
    for el in ELEMENTS:
        pts = [(np.log10(r["elements"][el]["amp"]),
                np.log10(r["elements"][el]["official_wt"]))
               for r in dataset
               if el in r["elements"] and r["elements"][el]["amp"] > 0
               and r["elements"][el]["official_wt"] > 0]
        if len(pts) < 8:
            continue
        rng = random.Random(11)
        pts_shuffled = pts[:]
        rng.shuffle(pts_shuffled)
        n_fit = int(len(pts_shuffled) * 0.7)
        fit_pts, val_pts = pts_shuffled[:n_fit], pts_shuffled[n_fit:]
        xs, ys = zip(*fit_pts)
        slope, intercept = np.polyfit(xs, ys, 1)
        preds = [10 ** (slope * x + intercept) for x, _ in val_pts]
        trues = [10 ** y for _, y in val_pts]
        rmses = [float(np.sqrt(np.mean([(p - t) / t ** 2 * t for p, t in
                                        zip(preds, trues)]) ** 2))] if val_pts else []
        rel_err = [abs(p - t) / t for p, t in zip(preds, trues)]
        out["calibrations"][el] = {
            "slope": float(slope), "intercept": float(intercept),
            "n_fit": len(fit_pts), "n_val": len(val_pts),
            "val_rel_err_median": float(np.median(rel_err)) if rel_err else None,
        }
        print("  %-3s 标定 log-log 斜率 %.2f；验证中位相对误差 %.0f%%（n=%d）"
              % (el, slope, 100 * (out["calibrations"][el]["val_rel_err_median"] or 0),
                 len(val_pts)))

    out_path = data_root / "APXS" / "缓存" / "apxs_calib.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = str(out_path) + ".part"
    open(tmp, "w", encoding="utf-8").write(json.dumps(out, ensure_ascii=False, indent=1))
    import os as _os
    _os.replace(tmp, str(out_path))
    print("标定曲线写入", out_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
