#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
研究脚本：拉曼基线方法对比——现行对称多项式 vs ALS vs airPLS
============================================================
背景（Carey et al., LPSC 2015 #2464）：六种自动基线法在纯矿物谱上的
矿物匹配实验里，airPLS（Zhang 2010）种级正确率最高（48.4%），ALS 次之
（46.7%），不去基线只有 8.2%。本项目管线/分类用的是对称稳健多项式
（degree 5, sigma-clip），从未与 Whittaker 系方法做过同数据对比。

做什么：
  1. 纯 numpy 实现 ALS（Eilers 2005 非对称最小二乘）与 airPLS
     （Zhang 2010 自适应迭代重加权惩罚最小二乘），不引入新依赖
  2. 合成谱（已知真值）定量：基线保真度 RMS、噪声窗偏移/平整度、
     噪声窗假峰数、真实峰位/峰高恢复——5 个噪声种子取均值±标准差
  3. 随包 sol293 真实谱：预处理（baseline 关）后的中位谱上三种基线
     各自扣减，列指纹区检出带清单与噪声窗统计——只记录差异，
     **不判优劣、不改管线默认**
  4. 产出 numbers.json（scope 免责 + 代码/数据 SHA 指纹）+ report.md

用法：
    python baseline_compare.py                 # 默认随包示例
    python baseline_compare.py --sample-dir <目录>

注意：
  * 这是一次冻结测量，不是方法学裁决；换基线属于分类输入变更，
    须另行经金标准台账评估后单独立项。
"""

from __future__ import annotations

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0, str(_BootPath(__file__).resolve().parents[2] / "02_图形界面代码"))
import source_paths as _layout            # noqa: E402

import argparse
import datetime
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

import project_layout                      # noqa: E402
from sherloc import minerals, organize, preprocess   # noqa: E402

LASER_NM = 248.6
NOISE_CM = (2500.0, 3000.0)
FINGERPRINT = (800.0, 1800.0)
SEEDS = (11, 23, 37, 41, 53)


# ---------- 第三方方法的纯 numpy 实现 ----------

def _diff_matrix(n, order):
    """n 阶单位矩阵的 k 次差分（Whittaker 平滑的惩罚核）。"""
    d = np.eye(n)
    for _ in range(order):
        d = np.diff(d, axis=0)
    return d


def als_baseline(y, lam=1.0e6, p=0.01, n_iter=10, order=2):
    """ALS：Eilers & Boelens 2005。非对称权重 + 二阶差分惩罚。"""
    d = _diff_matrix(y.size, order)
    h = lam * (d.T @ d)
    w = np.ones(y.size)
    z = y.astype(float)
    for _ in range(n_iter):
        z = np.linalg.solve(np.diag(w) + h, w * y)
        wn = np.where(y > z, p, 1.0 - p)
        if np.array_equal(wn, w):
            break
        w = wn
    return z


def airpls_baseline(y, lam=1.0e6, n_iter=15, order=2):
    """airPLS：Zhang et al. 2010 原版口径——二阶差分惩罚、λ=1e6、
    负残差指数权重 w=exp(i·|d|/T−1)（T 取首轮负残差绝对值之和，
    起阻尼作用）、迭代上限默认 15（原版默认）。

    实测敏感性（Phase 3 冻结测量）：本轴上 15 轮欠拟合（噪声窗偏移
    +16）、26 轮过冲（z 掉到数据下方 −80）。脚本另附轮数敏感性扫描。
    """
    d = _diff_matrix(y.size, order)
    h = lam * (d.T @ d)
    w = np.ones(y.size)
    z = y.astype(float)
    t = None
    for i in range(1, n_iter + 1):
        z = np.linalg.solve(np.diag(w) + h, w * y)
        resid = y - z
        neg = resid[resid < 0]
        if neg.size == 0:
            break
        if t is None:
            t = float(np.abs(neg).sum())
        w = np.where(resid < 0, np.exp(i * np.abs(resid) / t - 1.0), 0.0)
    return z


# ---------- 评测 ----------

def peak_curves(x, peaks):
    return sum(a * np.exp(-4.0 * np.log(2.0) * ((x - c) / f) ** 2)
               for c, a, f in peaks)


def synthetic_round(seed, als_lam):
    rng = np.random.default_rng(seed)
    x = np.linspace(214.0, 4124.0, 2048)
    peaks = [(831.0, 60.0, 24.0), (1029.0, 45.0, 30.0), (1084.0, 55.0, 28.0),
             (1600.0, 40.0, 40.0), (2497.0, 50.0, 60.0)]
    truth = peak_curves(x, peaks)
    baseline = (60.0 * np.exp(-0.5 * ((x - 1400.0) / 800.0) ** 2)
                + 0.012 * x + 20.0)
    noise = rng.normal(0.0, 3.0, x.size)
    y = truth + baseline + noise
    sigma = float(np.std(y[(x >= NOISE_CM[0]) & (x <= NOISE_CM[1])]))
    out = {}
    methods = {
        "poly5_symmetric": minerals.baseline(x, y),
        "als": als_baseline(y, lam=als_lam),
        "airpls": airpls_baseline(y),
    }
    noise_w = (x >= NOISE_CM[0]) & (x <= NOISE_CM[1])
    for name, base in methods.items():
        corrected = y - base
        finger = (x >= FINGERPRINT[0]) & (x <= FINGERPRINT[1])
        fidelity = float(np.sqrt(np.mean((corrected[finger] - truth[finger]) ** 2)))
        # 真实峰位/峰高恢复
        recover = []
        for c, a, _f in peaks:
            win = (x >= c - 15.0) & (x <= c + 15.0)
            found = float(x[win][int(np.argmax(corrected[win]))])
            amp = float(corrected[win].max())
            recover.append({"center_true": c, "center_found": found,
                            "pos_err": found - c, "amp_true": a, "amp_found": amp})
        # 噪声窗假峰：>5σ 的局部极大
        nw = corrected[noise_w]
        local = (nw[1:-1] > nw[:-2]) & (nw[1:-1] > nw[2:]) & (nw[1:-1] > 5 * sigma)
        out[name] = {
            "baseline_fidelity_rms": fidelity,
            "noise_offset": float(np.median(corrected[noise_w])),
            "noise_std": float(np.std(corrected[noise_w])),
            "false_bands_noise_window": int(local.sum()),
            "peak_recovery": recover,
        }
    return out, sigma


def summarize_synthetic(runs):
    keys = list(runs[0].keys())
    summary = {}
    for method in keys:
        rows = [r[method] for r in runs]
        def ms(field):
            vals = [row[field] for row in rows]
            return {"mean": float(np.mean(vals)), "std": float(np.std(vals))}
        recs = [row["peak_recovery"] for row in rows]
        pos_err = [abs(rec["pos_err"]) for row in recs for rec in row]
        summary[method] = {
            "baseline_fidelity_rms": ms("baseline_fidelity_rms"),
            "noise_offset": ms("noise_offset"),
            "noise_std": ms("noise_std"),
            "false_bands_noise_window": ms("false_bands_noise_window"),
            "peak_pos_abs_err_max": {"mean": float(np.mean(pos_err)),
                                     "std": float(np.std(pos_err))},
        }
    return summary


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def real_sample(sample, als_lam):
    steps = {"deadpixel": True, "dark": True, "laser": True,
             "despike": False, "baseline": False}      # 基线关：本脚本自己扣
    _, scans = organize.scan_folder(str(sample))
    scans = [s for s in scans if s.parts]
    out = []
    for sc in scans:
        res = preprocess.run_pipeline(sc, dict(steps))
        shift = np.asarray(res["shift"], dtype=float)
        spec = np.asarray(res["r1"], dtype=float)
        if spec.size == 0:
            continue
        y = np.median(spec, axis=0)
        sigma = float(np.std(y[(shift >= NOISE_CM[0]) & (shift <= NOISE_CM[1])]))
        noise_w = (shift >= NOISE_CM[0]) & (shift <= NOISE_CM[1])
        finger = (shift >= FINGERPRINT[0]) & (shift <= FINGERPRINT[1])
        methods = {
            "poly5_symmetric": minerals.baseline(shift, y),
            "als": als_baseline(y, lam=als_lam),
            "airpls": airpls_baseline(y),
        }
        per = {}
        for name, base in methods.items():
            corrected = y - base
            bands = minerals.find_bands(shift[finger], corrected[finger], sigma)
            centers = [float(b["pos"]) for b in (bands or [])]
            per[name] = {
                "bands_fingerprint": [round(v, 1) for v in centers],
                "n_bands": len(centers),
                "noise_offset": float(np.median(corrected[noise_w])),
                "noise_std": float(np.std(corrected[noise_w])),
            }
        out.append({"scan": sc.key, "points": int(spec.shape[0]), "methods": per})
        print(sc.key, json.dumps(per, ensure_ascii=False)[:200])
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description="基线方法对比（冻结测量，不做裁决）")
    ap.add_argument("--sample-dir", default=None,
                    help="示例数据目录；默认 <数据根>/SHERLOC/工作示例（退回随包示例）")
    ap.add_argument("--out", default=None, help="输出目录")
    args = ap.parse_args(argv)

    data_root = project_layout.data_root()
    sample = Path(args.sample_dir) if args.sample_dir else None
    if sample is None:
        for name in ("工作示例", "随包示例"):
            cand = data_root / "SHERLOC" / name
            if cand.is_dir() and any(cand.glob("*.csv")):
                sample = cand
                break
    if sample is None:
        print("找不到示例数据；先加载示例或用 --sample-dir 指定。")
        return 1

    # 公平起见给 ALS 一档"轻调参"：在种子 11 上扫 λ，按保真 RMS 选优；
    # 现行 poly5 是在役参数，不动。airPLS 用原版固定口径（λ=1e6、15 轮），
    # 另测迭代轮数敏感性（15/20/26 轮）——原版默认在本轴上欠拟合且敏感，
    # 这本身就是要冻结的测量结论，避免稻草人也避免藏问题。
    lam_grid = [1.0e4, 1.0e5, 1.0e6]
    probe = 214.0 + np.arange(2048) * (4124.0 - 214.0) / 2047.0
    peaks_true = [(831.0, 60.0, 24.0), (1029.0, 45.0, 30.0), (1084.0, 55.0, 28.0),
                  (1600.0, 40.0, 40.0), (2497.0, 50.0, 60.0)]
    truth_probe = peak_curves(probe, peaks_true)
    base_probe = (60.0 * np.exp(-0.5 * ((probe - 1400.0) / 800.0) ** 2)
                  + 0.012 * probe + 20.0)
    y_probe = (truth_probe + base_probe
               + np.random.default_rng(11).normal(0.0, 3.0, probe.size))
    finger_probe = (probe >= FINGERPRINT[0]) & (probe <= FINGERPRINT[1])
    scored = []
    for lam in lam_grid:
        corr = y_probe - als_baseline(y_probe, lam=lam)
        scored.append((float(np.sqrt(np.mean((corr[finger_probe]
                                              - truth_probe[finger_probe]) ** 2))), lam))
    scored.sort()
    tuning = {"als": {"lam": scored[0][1], "grid": lam_grid,
                      "probe_rms_by_lam": {str(lam): r for r, lam in scored}}}
    print("轻调参 als：λ=%g（探针 RMS %.1f）" % (scored[0][1], scored[0][0]))

    airpls_lam, airpls_niter = 1.0e6, 15
    sens = {}
    nw_probe = (probe >= NOISE_CM[0]) & (probe <= NOISE_CM[1])
    for niter in (15, 20, 26):
        corr = y_probe - airpls_baseline(y_probe, lam=airpls_lam, n_iter=niter)
        sens[str(niter)] = {
            "baseline_fidelity_rms": float(np.sqrt(np.mean(
                (corr[finger_probe] - truth_probe[finger_probe]) ** 2))),
            "noise_offset": float(np.median(corr[nw_probe])),
        }
        print("airPLS 迭代敏感性：%2d 轮 -> 保真RMS %.1f，噪声窗偏移 %+.1f"
              % (niter, sens[str(niter)]["baseline_fidelity_rms"],
                 sens[str(niter)]["noise_offset"]))

    print("合成谱评测（%d 个噪声种子）..." % len(SEEDS))
    runs = [synthetic_round(seed, tuning["als"]["lam"])[0] for seed in SEEDS]
    synthetic = summarize_synthetic(runs)
    for name, row in synthetic.items():
        print("  %-16s 保真RMS %.1f±%.1f  噪声窗偏移 %.2f  假峰 %.1f"
              % (name, row["baseline_fidelity_rms"]["mean"],
                 row["baseline_fidelity_rms"]["std"], row["noise_offset"]["mean"],
                 row["false_bands_noise_window"]["mean"]))

    print("真实样本评测：%s" % sample)
    real = real_sample(sample, tuning["als"]["lam"])

    stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = Path(args.out) if args.out else Path(__file__).resolve().parent / "out" / "baseline_compare" / stamp
    out.mkdir(parents=True, exist_ok=True)
    sample_files = sorted(p for p in sample.rglob("*.csv") if "rrs" in p.name.lower())
    doc = {
        "scope": "基线方法的冻结测量，不是方法学裁决，不构成换基线的决定；"
                 "换基线属于分类输入变更，须经金标准台账评估后单独立项。"
                 "数字只对本文件记录的代码/数据指纹负责。",
        "reference": "Carey et al., LPSC 2015 #2464; Zhang et al. 2010 Analyst 135:1138; "
                     "Eilers & Boelens 2005",
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "methods": {
            "poly5_symmetric": "现行：degree5 + sigma-clip（minerals.baseline / 管线同口径）",
            "als": "ALS 二阶差分, p=0.01, ≤10 轮；λ=%.0e（种子11轻调参）" % tuning["als"]["lam"],
            "airpls": "airPLS 原版口径：二阶差分, λ=1e6, 15 轮（见迭代敏感性）",
        },
        "lambda_tuning": tuning,
        "airpls_iteration_sensitivity": {
            "note": "原版默认 15 轮在本轴欠拟合；轮数增大先改善后过冲（z 掉到数据下方）。"
                    "这是冻结的稳健性测量，不是调参建议。",
            "probe": sens,
        },
        "synthetic": {"seeds": list(SEEDS), "peaks_true": [[831, 60, 24], [1029, 45, 30],
                        [1084, 55, 28], [1600, 40, 40], [2497, 50, 60]],
                        "summary": synthetic},
        "real": {"sample_dir": str(sample),
                 "sample_sha256": {p.name: sha256_of(p) for p in sample_files},
                 "per_scan": real},
        "code_sha256": project_layout.source_hashes(include_verification=False),
    }
    (out / "numbers.json").write_text(json.dumps(doc, ensure_ascii=False, indent=1),
                                      encoding="utf-8")

    lines = ["# 基线方法对比（冻结测量）", "",
             "- 参考：" + doc["reference"],
             "- 方法：" + "；".join(doc["methods"].values()), "",
             "## 合成谱（真值已知，%d 种子均值±std）" % len(SEEDS), "",
             "| 方法 | 保真RMS | 噪声窗偏移 | 噪声窗假峰 | 峰位最大误差 |",
             "|---|---|---|---|---|"]
    for name, row in synthetic.items():
        lines.append("| %s | %.1f±%.1f | %.2f±%.2f | %.1f±%.1f | %.1f±%.1f |"
                     % (name, row["baseline_fidelity_rms"]["mean"],
                        row["baseline_fidelity_rms"]["std"],
                        row["noise_offset"]["mean"], row["noise_offset"]["std"],
                        row["false_bands_noise_window"]["mean"],
                        row["false_bands_noise_window"]["std"],
                        row["peak_pos_abs_err_max"]["mean"],
                        row["peak_pos_abs_err_max"]["std"]))
    lines += ["", "## 随包 sol293（中位谱，指纹区检出带）", ""]
    for row in real:
        lines.append("### %s（%d 点）" % (row["scan"], row["points"]))
        for name, per in row["methods"].items():
            lines.append("- %s：%d 带 %s；噪声窗偏移 %.2f、std %.1f"
                         % (name, per["n_bands"], per["bands_fingerprint"],
                            per["noise_offset"], per["noise_std"]))
        lines.append("")
    (out / "report.md").write_text("\n".join(lines), encoding="utf-8")
    print("报告与数字已写入 %s" % out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
