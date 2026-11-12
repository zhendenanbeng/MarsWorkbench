#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ChemCam 数据准备脚本 02：下载 CCS 标定靶训练观测（PLS 训练集）
================================================================
输入：01 号脚本的 census（ccs_target_census.json）。
输出：<数据根>/ChemCam/下载/ccs_training/ 下的完整 CCS CSV +
      training_manifest.json（{靶名: [{file, sha256, bytes, pds_url}]}）。

采样策略：每靶最多 --per-target 个观测，跨 sol 均匀取样（时间分布代表性）。
断点续传：已存在且大小一致的文件跳过。

用法：
    python 02_download_ccs_training.py                    # 每靶 6 个观测
    python 02_download_ccs_training.py --per-target 12
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import urllib.request
from pathlib import Path

BASE = "https://pds-geosciences.wustl.edu/msl/msl-m-chemcam-libs-4_5-rdr-v1/"
UA = {"User-Agent": "sherloc-workbench/1.0"}


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main(argv=None):
    ap = argparse.ArgumentParser(description="下载 CCS 训练观测")
    ap.add_argument("--per-target", type=int, default=6)
    ap.add_argument("--census", default=None, help="census JSON 路径")
    ap.add_argument("--dest", default=None, help="下载目录")
    args = ap.parse_args(argv)

    try:
        import project_layout
        data_root = Path(project_layout.data_root())
    except Exception:                                    # noqa: BLE001
        from pathlib import Path as P
        data_root = P(".").resolve()
    census_path = Path(args.census) if args.census else data_root / "ChemCam" / "缓存" / "ccs_target_census.json"
    census = json.loads(census_path.read_text(encoding="utf-8"))
    targets = census["targets"]
    cc = {"macusanite", "norite", "picrite", "shergottite", "graphite",
          "kga-2med-s", "nau-2lo-s", "nau-2med-s", "nau-2hi-s", "ti"}
    buckets = {t: v for t, v in targets.items() if t.lower() in cc}
    if not buckets:
        print("census 里没有 CCS 标定靶观测；先运行 01_ccs_target_census.py")
        return 1

    dest = Path(args.dest) if args.dest else data_root / "ChemCam" / "下载" / "ccs_training"
    dest.mkdir(parents=True, exist_ok=True)
    manifest = {"scope": "CCS 在轨标定靶训练观测清单（谱级原始训练集）。成分对账值在"
                         "published_calibrations.json，本清单只含谱文件与指纹。",
                "per_target": args.per_target, "targets": {}}
    failed = 0
    for t, bucket in sorted(buckets.items()):
        paths = bucket.get("paths") or [e["path"] + e["file"] for e in bucket.get("examples", [])]
        n = len(paths)
        step = max(1, n // args.per_target)
        picked = paths[::step][:args.per_target]
        rows = []
        for rel in picked:
            fname = rel.split("/")[-1]
            local = dest / fname
            url = BASE + rel
            if not (local.is_file() and local.stat().st_size > 0):
                try:
                    req = urllib.request.Request(url, headers=UA)
                    data = urllib.request.urlopen(req, timeout=180).read()
                    tmp = str(local) + ".part"
                    open(tmp, "wb").write(data)
                    os.replace(tmp, str(local))
                    print("  %s <- %s (%.1f MB)" % (t, fname, len(data) / 1e6))
                except Exception as exc:                 # noqa: BLE001
                    print("  失败 %s：%s" % (fname, str(exc)[:60]))
                    failed += 1
                    continue
            else:
                print("  已有 %s" % fname)
            rows.append({"file": fname, "bytes": local.stat().st_size,
                         "sha256": sha256_of(str(local)), "pds_url": url})
        manifest["targets"][t] = {"n_available": n, "downloaded": rows}
        print("靶 %-16s 可用 %5d，下载 %d" % (t, n, len(rows)))

    out = dest / "training_manifest.json"
    tmp = str(out) + ".part"
    open(tmp, "w", encoding="utf-8").write(json.dumps(manifest, ensure_ascii=False, indent=1))
    os.replace(tmp, str(out))
    total = sum(len(v["downloaded"]) for v in manifest["targets"].values())
    print("manifest 写入 %s（%d 个文件，%d 个失败）" % (out, total, failed))
    return 1 if failed and total == 0 else 0


if __name__ == "__main__":
    sys.exit(main())
