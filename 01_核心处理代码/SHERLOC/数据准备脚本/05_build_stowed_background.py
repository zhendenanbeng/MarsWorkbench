#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SHERLOC 数据准备脚本 05：把收纳臂扫描聚合成仪器本底谱
======================================================
输入：04 号脚本下载的收纳臂扫描（RRS）。
输出：<数据根>/SHERLOC/缓存/stowed_background.json —— 拉曼位移轴 + 本底
强度谱（全部收纳臂光谱逐通道中位数，论文口径：用最高 SNR 的集体谱做减数）。

方法（对应 Jakubek et al. 2024）：
  1. 每条收纳臂光谱：裁死区（50:2098）→ 换拉曼位移轴 → 有配对暗帧则扣
  2. 逐扫描取 Region 1 全点的中位数谱（单谱 SNR<1，必须集体平均才显形）
  3. 跨扫描稳健剔除：按 800 cm-1 带窗中位数偏离中位数的倍数剔除离群扫描
     （防止同 sol 混入的常规观测污染本底）
  4. 保留扫描逐通道中位数 = 本底谱

绝对强度不做激光归一——减除时用缩放因子 M 对齐量纲（见工作台
「仪器本底扣除」开关），所以本底谱只保证"形状+相对幅度"。

用法：
    python 05_build_stowed_background.py
    python 05_build_stowed_background.py --scan "0413|SERIAL|SCLK"   # 只要这些扫描
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

import numpy as np

import project_layout                      # noqa: E402
from sherloc import organize, preprocess   # noqa: E402

LASER_NM = 248.6
# 跨扫描一致性检验窗：800 cm-1 是本底最强带（Jakubek Fig.1）
BAND_CHECK_CM = (775.0, 825.0)
BAND_REJECT_MAD = 5.0


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def scan_median_spectrum(scan, region=1):
    """一条扫描的 Region 谱矩阵：裁死区+位移轴+（可选）扣暗帧。

    返回 (shift 升序数组, 谱矩阵 [点 x 通道])；无光谱/无波长表返回 (None, None)。
    """
    raw = preprocess.load_scan(scan)
    arr = raw.get("r%d" % region)
    if arr is None or raw["wl_nm"] is None:
        return None, None
    sel = slice(preprocess.DEAD_LO, preprocess.DEAD_HI)
    wl = np.asarray(raw["wl_nm"])[sel]
    R = arr[:, sel].astype(float)
    dark = raw.get("d%d" % region)
    if raw.get("has_dark") and dark is not None and dark.shape[0] == R.shape[0]:
        R = R - dark[:, sel].astype(float)
    shift = 1.0e7 / LASER_NM - 1.0e7 / wl
    order = np.argsort(shift)
    return shift[order], R[:, order]


def main(argv=None):
    ap = argparse.ArgumentParser(description="聚合收纳臂扫描为仪器本底谱")
    ap.add_argument("--stowed-dir", default=None,
                    help="收纳臂数据根目录；默认 <数据根>/SHERLOC/下载/stowed")
    ap.add_argument("--out", default=None,
                    help="输出 JSON；默认 <数据根>/SHERLOC/缓存/stowed_background.json")
    ap.add_argument("--scan", action="append", default=[],
                    help="只用这个扫描（key = sol|序列|SCLK，可重复）")
    ap.add_argument("--region", type=int, default=1, choices=(1, 2, 3))
    args = ap.parse_args(argv)

    stowed = Path(args.stowed_dir) if args.stowed_dir else (
        project_layout.data_root() / "SHERLOC" / "下载" / "stowed")
    if not stowed.is_dir():
        print("收纳臂目录不存在：%s\n先运行 04_fetch_stowed_background.py 下载数据。" % stowed)
        return 1
    only = set(args.scan)

    per_scan = []          # (key, sol, median 谱, n_spectra)
    shift_ref = None
    for sol_dir in sorted(p for p in stowed.iterdir() if p.is_dir()):
        _, scans = organize.scan_folder(str(sol_dir))
        for sc in scans:
            if not sc.parts:
                continue
            key = sc.key
            if only and key not in only:
                continue
            shift, R = scan_median_spectrum(sc, args.region)
            if R is None or not np.isfinite(R).all():
                print("跳过 %s（无光谱或含非有限值）" % key)
                continue
            if shift_ref is None:
                shift_ref = shift
            elif len(shift) != len(shift_ref) or not np.allclose(shift, shift_ref, atol=1e-4):
                print("跳过 %s（位移轴与首个扫描不一致）" % key)
                continue
            per_scan.append((key, sc.sol, np.median(R, axis=0), int(R.shape[0])))
            print("收下 %s（%d 点）" % (key, R.shape[0]))

    if len(per_scan) < 2:
        print("可用扫描不足 2 个（实际 %d 个），本底谱会太脆；先补数据。" % len(per_scan))
        return 1

    # --- 跨扫描稳健剔除：800 cm-1 带窗中位数偏离倍数 ---
    lo, hi = BAND_CHECK_CM
    band = (shift_ref >= lo) & (shift_ref <= hi)
    amps = np.array([np.median(s[band]) for _, _, s, _ in per_scan])
    center = float(np.median(amps))
    mad = float(1.4826 * np.median(np.abs(amps - center))) or 1.0
    keep, rejected = [], []          # (元数据 dict, 本底谱数组) 并行两路
    keep_spec = []
    for (key, sol, spec, n), amp in zip(per_scan, amps):
        row = {"scan": key, "sol": sol, "n_spectra": n, "band_800_median": float(amp)}
        if abs(amp - center) > BAND_REJECT_MAD * mad:
            rejected.append(row)
        else:
            keep.append(row)
            keep_spec.append(spec)
    if len(keep) < 2:
        print("一致性剔除后可用扫描不足 2 个；数据可能混入了非收纳臂观测，"
              "请人工核对清单后用 --scan 指定。")
        return 1
    for r in rejected:
        print("剔除离群扫描 %s（800 带中位 %.1f vs 集体中位 %.1f）"
              % (r["scan"], r["band_800_median"], center))

    background = np.median(np.vstack(keep_spec), axis=0)

    hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
              for p in sorted((_layout.project_layout.APP).glob("*.py"))}
    out = Path(args.out) if args.out else (
        project_layout.data_root() / "SHERLOC" / "缓存" / "stowed_background.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    doc = {
        "scope": "SHERLOC 仪器本底拉曼谱（收纳臂/暗天空口径）。形状与相对幅度可信，"
                 "绝对强度未做激光归一；减除时必须按窗口中位比缩放（M）。"
                 "不是矿物真值，不参与任何分类判定。",
        "method": "逐扫描 Region 全点中位 -> 800 带一致性剔除 -> 保留扫描逐通道中位",
        "reference": "Jakubek et al., LPSC 2024 #1174; Applied Spectroscopy 2024",
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "laser_nm": LASER_NM,
        "region": args.region,
        "n_channels": int(background.size),
        "shift_cm1": [round(float(v), 4) for v in shift_ref],
        "background": [float(v) for v in background],
        "scans_used": [r for r in keep],
        "scans_rejected": rejected,
        "band_check_cm": list(BAND_CHECK_CM),
        "band_reject_mad": BAND_REJECT_MAD,
        "code_sha256": hashes,
    }
    out.write_text(json.dumps(doc, ensure_ascii=False), encoding="utf-8")
    print("本底谱写入 %s（%d 通道，来自 %d 个扫描 / %d 条光谱）"
          % (out, background.size, len(keep), sum(r["n_spectra"] for r in keep)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
