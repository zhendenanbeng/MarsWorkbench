#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
研究脚本：量化 SHERLOC 仪器本底（熔石英 490/800/1055/1190 cm-1）在
现行预处理输出里的残留幅度，以及"缩放相减"能拿回多少。
======================================================================
依据 Jakubek et al. 2024：光学件本底在单谱里 SNR<1 不可见，但逐点
中位/平均时会显形；橄榄石 815-865 cm-1 窗口与 800 本底带直接重叠。

做什么：
  1. 对随包 sol293 示例逐扫描跑现行默认预处理（deadpixel/dark/laser/baseline）
  2. 把收纳臂本底谱插值到各扫描位移轴，按四本底带窗中位比估计缩放 M
  3. 分别量化：四带窗口中位残留（减前/减后）、橄榄石窗峰形、
     噪声窗 2500-3000 cm-1 的标准差（减除不应扰动无本底区）
  4. 产出 report.md + numbers.json（带 scope 免责与全部 SHA 指纹）

用法：
    python bkg_residual_probe.py                       # 用默认示例与本底
    python bkg_residual_probe.py --bg <stowed_background.json>

注意：
  * 这是一次冻结测量，不是方法学裁决；数字只对当时代码/数据指纹负责。
  * 不修改任何预处理或分类默认行为。
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
from sherloc import organize, preprocess   # noqa: E402

LASER_NM = 248.6
BANDS_CM = {"b490": 490.0, "b800": 800.0, "b1055": 1055.0, "b1190": 1190.0}
BAND_HALF_WIDTH = 25.0
OLIVINE_CM = (815.0, 865.0)                # 与 minerals.py RANGES["olivine"] 同一口径
NOISE_CM = (2500.0, 3000.0)
DEFAULT_STEPS = {"deadpixel": True, "dark": True, "laser": True,
                 "despike": False, "baseline": True}          # 与工作台默认勾选一致
M_MIN, M_MAX = 0.0, 2.0


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def window_median(spec, shift, lo, hi):
    m = (shift >= lo) & (shift <= hi)
    return None if not m.any() else float(np.median(spec[m]))


def main(argv=None):
    ap = argparse.ArgumentParser(description="量化仪器本底在预处理输出中的残留")
    ap.add_argument("--sample-dir", default=None,
                    help="示例数据目录；默认 <数据根>/SHERLOC/工作示例（退回随包示例）")
    ap.add_argument("--bg", default=None,
                    help="本底 JSON；默认 <数据根>/SHERLOC/缓存/stowed_background.json")
    ap.add_argument("--out", default=None, help="输出目录；默认 研究脚本/out/bkg_residual/<时间戳>")
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
        print("找不到 sol293 示例；先在工作台点『加载示例数据』，"
              "或用 --sample-dir 指定目录。")
        return 1
    bg_path = Path(args.bg) if args.bg else data_root / "SHERLOC" / "缓存" / "stowed_background.json"
    if not bg_path.is_file():
        print("本底谱不存在：%s\n先运行 数据准备脚本/04 与 05 生成。" % bg_path)
        return 1
    bg = json.loads(bg_path.read_text(encoding="utf-8"))
    bg_shift = np.asarray(bg["shift_cm1"], dtype=float)
    bg_spec = np.asarray(bg["background"], dtype=float)

    _, scans = organize.scan_folder(str(sample))
    scans = [s for s in scans if s.parts]
    if not scans:
        print("示例目录里没有可用扫描：%s" % sample)
        return 1

    stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = Path(args.out) if args.out else Path(__file__).resolve().parent / "out" / "bkg_residual" / stamp
    out.mkdir(parents=True, exist_ok=True)

    per_scan = []
    for sc in scans:
        res = preprocess.run_pipeline(sc, dict(DEFAULT_STEPS))
        shift, spec = res["shift"], res["r1"]
        if shift is None or spec is None or spec.size == 0:
            continue
        spec = np.asarray(spec, dtype=float)          # [点 x 通道]
        shift = np.asarray(shift, dtype=float)
        if shift[0] > shift[-1]:                      # 保险：升序
            shift, spec = shift[::-1], spec[:, ::-1]
        bg_i = np.interp(shift, bg_shift, bg_spec)

        w = np.zeros_like(shift, dtype=bool)
        for c in BANDS_CM.values():
            w |= (shift >= c - BAND_HALF_WIDTH) & (shift <= c + BAND_HALF_WIDTH)
        med_scan = float(np.median(spec[:, w]))
        med_bg = float(np.median(bg_i[w]))
        M = float(np.clip(med_scan / med_bg if med_bg else 0.0, M_MIN, M_MAX))
        # 论文口径的 M 是"调到主带刚好抹平"；四窗联合中位比会被次强带拖偏
        # （Phase 1 实测：800 带过校正）。所以同时给一组只用 800 主带估的 M。
        w8 = (shift >= BANDS_CM["b800"] - BAND_HALF_WIDTH) & (shift <= BANDS_CM["b800"] + BAND_HALF_WIDTH)
        M8 = float(np.clip(float(np.median(spec[:, w8])) / float(np.median(bg_i[w8]))
                           if np.median(bg_i[w8]) else 0.0, M_MIN, M_MAX))
        subtracted = spec - M * bg_i
        subtracted8 = spec - M8 * bg_i

        bands = {}
        for name, c in BANDS_CM.items():
            lo, hi = c - BAND_HALF_WIDTH, c + BAND_HALF_WIDTH
            m = (shift >= lo) & (shift <= hi)
            bands[name] = {
                "before_median": float(np.median(spec[:, m])),
                "after_median": float(np.median(subtracted[:, m])),
                "after_median_M800": float(np.median(subtracted8[:, m])),
            }
        ol = (shift >= OLIVINE_CM[0]) & (shift <= OLIVINE_CM[1])
        ol_pre, ol_post = spec[:, ol], subtracted[:, ol]
        n_zap = int((np.abs(ol_pre) > 1e-12).any(axis=0).sum())
        bands["olivine_window"] = {
            "before_median": float(np.median(ol_pre)),
            "after_median": float(np.median(ol_post)),
            "after_median_M800": float(np.median(spec[:, ol] - M8 * bg_i[ol][None, :])),
            "before_peak_max": float(np.max(np.median(ol_pre, axis=0))) if ol_pre.size else 0.0,
            "after_peak_max": float(np.max(np.median(ol_post, axis=0))) if ol_post.size else 0.0,
            "channels": int(n_zap),
        }
        nz = (shift >= NOISE_CM[0]) & (shift <= NOISE_CM[1])
        bands["noise_window"] = {
            "before_std": float(np.std(spec[:, nz])),
            "after_std": float(np.std(subtracted[:, nz])),
            "after_std_M800": float(np.std(subtracted8[:, nz])),
        }
        per_scan.append({
            "scan": sc.key, "points": int(spec.shape[0]), "channels": int(spec.shape[1]),
            "scale_M_union": M, "scale_M_800": M8, "bands": bands,
        })
        print("%s  M并=%.3f M800=%.3f  b800 %.1f->%.1f(M并)/%.1f(M800)  橄榄石窗 %.1f->%.1f"
              % (sc.key, M, M8,
                 bands["b800"]["before_median"], bands["b800"]["after_median"],
                 bands["b800"]["after_median_M800"],
                 bands["olivine_window"]["before_median"], bands["olivine_window"]["after_median"]))

    if not per_scan:
        print("没有可分析的扫描。")
        return 1

    hashes = project_layout.source_hashes(include_verification=False)
    sample_files = sorted(p for p in sample.rglob("*.csv") if "rrs" in p.name.lower())
    doc = {
        "scope": "仪器本底残留的冻结测量。不是方法学裁决，不构成矿物鉴定准确率"
                 "或任何分类判定的证据；数字只对本文件记录的代码/数据指纹负责。",
        "reference": "Jakubek et al., LPSC 2024 #1174 / Applied Spectroscopy 2024",
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "sample_dir": str(sample),
        "background_file": str(bg_path),
        "background_sha256": sha256_of(bg_path),
        "background_scans_used": bg.get("scans_used", []),
        "sample_sha256": {p.name: sha256_of(p) for p in sample_files},
        "steps": dict(DEFAULT_STEPS),
        "bands_cm": {k: [v - BAND_HALF_WIDTH, v + BAND_HALF_WIDTH]
                     for k, v in BANDS_CM.items()},
        "olivine_window_cm": list(OLIVINE_CM),
        "noise_window_cm": list(NOISE_CM),
        "code_sha256": hashes,
        "per_scan": per_scan,
    }
    (out / "numbers.json").write_text(json.dumps(doc, ensure_ascii=False, indent=1),
                                      encoding="utf-8")

    lines = ["# 仪器本底残留量化（sol293 随包示例）", "",
             "- 参考：%s" % doc["reference"],
             "- 本底谱：%s（%d 个收纳臂扫描）" % (bg_path.name, len(bg.get("scans_used", []))),
             "- 预处理开关：%s" % doc["steps"],
             "- M并=四本底带窗联合中位比；M800=仅 800 主带中位比（Phase 1 实测联合口径会过校正）", ""]
    for r in per_scan:
        lines.append("## 扫描 %s（M并=%.3f，M800=%.3f）" % (r["scan"], r["scale_M_union"], r["scale_M_800"]))
        lines.append("")
        lines.append("| 窗口 | 减前中位 | 减后(M并) | 减后(M800) |")
        lines.append("|---|---|---|---|")
        for name in list(BANDS_CM) + ["olivine_window"]:
            b = r["bands"][name]
            lines.append("| %s | %.1f | %.1f | %.1f |"
                         % (name, b["before_median"], b["after_median"], b["after_median_M800"]))
        b = r["bands"]["noise_window"]
        lines.append("| 噪声窗std | %.1f | %.1f | %.1f |"
                     % (b["before_std"], b["after_std"], b["after_std_M800"]))
        lines.append("")
    (out / "report.md").write_text("\n".join(lines), encoding="utf-8")
    print("报告与数字已写入 %s" % out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
