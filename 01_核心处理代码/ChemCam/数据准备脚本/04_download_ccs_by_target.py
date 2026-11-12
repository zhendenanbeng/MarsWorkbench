#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ChemCam 数据准备脚本 04：按靶名从官方 MOC 筛观测并下载 CCS 训练谱
====================================================================
取代 01/02 的 census 路线：官方 MOC 表每行给出观测文件名 + 靶名 +
官方成分预测（PDS mslccm_1xxx/data/moc/，Ytsma 2023 同源数据）。

流程：
  1. 下载全部 MOC 分段（sol 0000-3000+，42 个小 csv）
  2. 解析数据行：文件名 / 靶名 / 官方 9 氧化物预测±误差
  3. 用本地 libsindex.tab 把文件名映射到 sol 路径
  4. 每靶均匀采样 --per-target 个观测，下载 CCS csv（完整 3MB 文件）
  5. 输出 training_manifest.json：{靶名: [{file, sha256, pds_url, official: {氧化物: [值, err]}}]}

用法：
    python 04_download_ccs_by_target.py                    # 每靶 6 个
    python 04_download_ccs_by_target.py --per-target 10
"""

from __future__ import annotations

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0, str(_BootPath(__file__).resolve().parents[3] / "02_图形界面代码"))
import source_paths as _layout            # noqa: E402

import argparse
import collections
import csv
import hashlib
import json
import os
import re
import sys
import urllib.request
from pathlib import Path

BASE = "https://pds-geosciences.wustl.edu/msl/msl-m-chemcam-libs-4_5-rdr-v1/mslccm_1xxx/"
MOC_DIR = BASE + "data/moc/"
UA = {"User-Agent": "sherloc-workbench/1.0"}
# 靶名与官方 MOC 表完全一致（col1 原样小写化匹配）
TARGETS = ("coronation", "norite", "picrite", "shergottite", "macusanite",
           "nau-2_hi-s", "nau-2_med-s", "nau-2_low-s", "kga-d_med-s")
OXIDE_ORDER = ("SiO2", "TiO2", "Al2O3", "FeOT", "MgO", "CaO", "Na2O", "K2O", "MnO")


def http_get(url, head=None, retries=4):
    headers = dict(UA)
    if head:
        headers["Range"] = "bytes=0-%d" % (head - 1)
    last = None
    for _ in range(retries):
        try:
            return urllib.request.urlopen(urllib.request.Request(url, headers=headers),
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


def main(argv=None):
    ap = argparse.ArgumentParser(description="按靶名下载 CCS 训练观测")
    ap.add_argument("--per-target", type=int, default=6)
    ap.add_argument("--dest", default=None)
    args = ap.parse_args(argv)

    try:
        import project_layout
        data_root = Path(project_layout.data_root())
    except Exception:                                    # noqa: BLE001
        data_root = Path(_layout.project_layout.ROOT or ".") / "03_下载数据"
    cache = data_root / "ChemCam" / "缓存"
    cache.mkdir(parents=True, exist_ok=True)

    # --- 1/2. 下载并解析 MOC ---
    obs = collections.defaultdict(list)                  # 靶名 -> [(fname, official)]
    n_moc = 0
    h = http_get(MOC_DIR).decode("utf-8", "replace")
    moc_files = sorted(set(re.findall(r'href="(/msl[^"]+moc[^"]*\.csv)"', h, re.I)))
    for rel in moc_files:
        data = http_get("https://pds-geosciences.wustl.edu" + rel).decode("utf-8", "replace")
        for line in data.splitlines():
            parts = [p.strip() for p in line.split(",")]
            if len(parts) < 30 or not parts[1]:
                continue
            target = parts[1]
            if target.lower() not in TARGETS:
                continue
            fname = parts[0]
            try:
                # MOC 列布局（表头每氧化物 4 列）：[值, '+/-', err, precision]
                vals = {ox: [float(parts[2 + 4 * i]), float(parts[4 + 4 * i])]
                        for i, ox in enumerate(OXIDE_ORDER)}
                obs[target].append((fname, vals))
                n_moc += 1
            except (ValueError, IndexError):
                continue
        n_moc += 0
    print("MOC 在轨靶观测总数:", n_moc, "（%d 个靶）" % len(obs))
    for t, rows in sorted(obs.items()):
        print("  %-12s %5d 个观测" % (t, len(rows)))

    # --- 3. libsindex 文件名 -> sol 路径 ---
    index_path = cache / "libsindex.tab"
    if not index_path.is_file():
        print("下载官方索引 libsindex.tab ...")
        http_get(BASE + "index/libsindex.tab", head=None)
        open(str(index_path) + ".part", "wb").write(http_get(BASE + "index/libsindex.tab"))
        os.replace(str(index_path) + ".part", str(index_path))
    name_to_path = {}
    for row in csv.reader(index_path.read_text(encoding="utf-8", errors="replace").splitlines()):
        if len(row) >= 3 and row[2].strip().upper().endswith(".CSV"):
            name_to_path[row[2].strip().upper()] = row[1].strip().replace("\\", "/") + row[2].strip()

    # --- 4. 每靶采样下载 ---
    dest = Path(args.dest) if args.dest else data_root / "ChemCam" / "下载" / "ccs_training"
    dest.mkdir(parents=True, exist_ok=True)
    manifest = {"scope": "CCS 在轨标定靶训练观测（谱级原始训练集）。official=官方 MOC 成分预测"
                         "（仅作对账基准，绝不作为训练真值；训练真值见 ccs_composition.json）。"
                if False else "CCS 在轨标定靶训练观测清单。official 字段是官方 MOC 成分预测"
                              "（对账基准，不是训练真值；训练真值在 ccs_composition.json）。",
                "per_target": args.per_target, "targets": {}}
    failed = 0
    for t, rows in sorted(obs.items()):
        rows = sorted(rows)                              # 按文件名/SCLK 稳定排序
        step = max(1, len(rows) // args.per_target)
        picked = rows[::step][:args.per_target]
        entries = []
        for fname, official in picked:
            rel = name_to_path.get(fname.upper())
            if not rel:
                print("  索引缺失：%s" % fname)
                continue
            local = dest / fname
            if not (local.is_file() and local.stat().st_size > 0):
                try:
                    data = http_get(BASE + rel)
                    tmp = str(local) + ".part"
                    open(tmp, "wb").write(data)
                    os.replace(tmp, str(local))
                except Exception as exc:                 # noqa: BLE001
                    print("  下载失败 %s：%s" % (fname, str(exc)[:60]))
                    failed += 1
                    continue
            entries.append({"file": fname, "bytes": local.stat().st_size,
                            "sha256": sha256_of(str(local)),
                            "pds_url": BASE + rel, "official": official})
        manifest["targets"][t] = {"n_available": len(rows), "downloaded": entries}
        print("  %-12s 采样 %d/%d" % (t, len(entries), len(rows)))

    out = dest / "training_manifest.json"
    tmp = str(out) + ".part"
    open(tmp, "w", encoding="utf-8").write(json.dumps(manifest, ensure_ascii=False, indent=1))
    os.replace(tmp, str(out))
    total = sum(len(v["downloaded"]) for v in manifest["targets"].values())
    print("manifest 写入 %s（%d 文件，%d 失败）" % (out, total, failed))
    return 0 if total else 1


if __name__ == "__main__":
    sys.exit(main())
