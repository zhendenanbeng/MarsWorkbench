#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ChemCam 数据准备脚本 01：CCS 在轨标定靶观测普查
==================================================
目的（Phase 1：自主 PLS 定量）：找出现在轨打在 9+1 个标定靶上的 CCS 观测，
按 TARGET 分桶——它们是自研 PLS 的训练集（谱）×文献成分（对账值）的谱来源。

方法：
  1. 下载 PDS 官方索引 libsindex.tab（全量 12 万行，缓存本地）
  2. 筛 PRODUCT_TYPE=CCS 的行
  3. 对每个 CCS 用 HTTP Range 只拉文件头部（16KB），读 "# TARGET = ..." 行
     ——PDS 服务器支持 Range（206），普查 2188 个产品只需 ~35MB 流量
  4. 按 TARGET 分桶，输出 census JSON（观测数/样例文件清单 per 靶）

产物：<数据根>/ChemCam/缓存/ccs_target_census.json

用法：
    python 01_ccs_target_census.py                # 全量普查
    python 01_ccs_target_census.py --limit 50     # 只普查前 N 个（调试）

注意：联网走环境变量代理；TARGET 为 UNK/空的行保留但不计入靶桶。
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
import urllib.request

BASE = "https://pds-geosciences.wustl.edu/msl/msl-m-chemcam-libs-4_5-rdr-v1/mslccm_1xxx/"
INDEX_URL = BASE + "index/libsindex.tab"
UA = {"User-Agent": "sherloc-workbench/1.0"}
HEAD_BYTES = 16384
CCS_TARGETS = {"macusanite", "norite", "picrite", "shergottite", "graphite",
               "kga-2med-s", "nau-2lo-s", "nau-2med-s", "nau-2hi-s", "ti"}


def http_get(url, dest=None, head=None):
    req = urllib.request.Request(url, headers=dict(UA, **({"Range": "bytes=0-%d" % (head - 1)} if head else {})))
    with urllib.request.urlopen(req, timeout=120) as r:
        data = r.read()
    if dest:
        tmp = dest + ".part"
        open(tmp, "wb").write(data)
        os.replace(tmp, dest)
    return data


def main(argv=None):
    ap = argparse.ArgumentParser(description="CCS 在轨标定靶观测普查")
    ap.add_argument("--limit", type=int, default=0, help="只普查前 N 个 CCS（调试）")
    ap.add_argument("--out", default=None, help="census JSON 输出路径")
    args = ap.parse_args(argv)

    import tempfile
    from pathlib import Path
    try:
        import project_layout
        cache = project_layout.data_root() / "ChemCam" / "缓存"
    except Exception:                                    # noqa: BLE001
        cache = Path(tempfile.gettempdir()) / "marsworkbench_ccs"
    cache.mkdir(parents=True, exist_ok=True)

    index_path = cache / "libsindex.tab"
    if not index_path.is_file():
        print("下载官方索引 libsindex.tab（约 32MB，仅一次）...")
        http_get(INDEX_URL, str(index_path))
    lines = index_path.read_text(encoding="utf-8", errors="replace").splitlines()
    print("索引行数:", len(lines))

    ccs = []
    seen_sclk = set()
    for row in csv.reader(lines):
        if len(row) < 6:
            continue
        fname = row[2].strip()
        ptype = row[5].strip()
        if "CCS" in fname.upper() and "CCS" in ptype.upper():
            # 同一观测（SCLK 相同，文件名首段）多版本只普查一次
            m = re.match(r"([A-Z0-9]+?)_(\d+)(ccs)", fname.upper())
            key = m.group(2) if m else fname
            if key in seen_sclk:
                continue
            seen_sclk.add(key)
            ccs.append((row[1].strip(), fname))
    print("CCS 观测（按 SCLK 去重）:", len(ccs))
    if args.limit:
        ccs = ccs[:args.limit]

    targets = {}
    unk = 0
    path_of = {fname: path for path, fname in ccs}

    def probe(job):
        path, fname = job
        url = BASE + path + fname
        try:
            head = http_get(url, head=HEAD_BYTES).decode("utf-8", "replace")
            m = re.search(r"#\s*TARGET\s*=\s*\"?([^,\r\n\"]+)", head, re.I)
            return fname, (m.group(1).strip() if m else "")
        except Exception as exc:                         # noqa: BLE001
            return fname, "ERR:" + str(exc)[:40]

    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=12) as pool:
        for i, (fname, target) in enumerate(pool.map(probe, ccs)):
            if target.startswith("ERR:"):
                if (i + 1) % 200 == 0:
                    print("  [%d/%d] 失败样例 %s: %s" % (i + 1, len(ccs), fname, target[4:44]))
                unk += 1
                continue
            if not target or target.upper() == "UNK":
                unk += 1
                continue
            bucket = targets.setdefault(target, {"n": 0, "examples": [], "paths": []})
            bucket["n"] += 1
            bucket["paths"].append(path + fname)
            if len(bucket["examples"]) < 3:
                bucket["examples"].append({"path": path_of[fname], "file": fname})
            if (i + 1) % 200 == 0:
                print("  [%d/%d] 已识别 %d 个靶、%d 个 UNK" % (i + 1, len(ccs), len(targets), unk))

    print("普查完成：识别 %d 个 TARGET（%d 个 UNK/空）" % (len(targets), unk))
    for t in sorted(targets):
        mark = " <-- CCS 标定靶" if t.lower() in CCS_TARGETS else ""
        print("  %-20s %5d 个观测%s" % (t, targets[t]["n"], mark))

    doc = {
        "scope": "CCS 在轨标定靶观测普查结果（谱级原始训练集清单）。TARGET 取自产品 CSV 头部官方字段；"
                 "本文件只列观测清单，不含任何自研科学量。",
        "index_source": INDEX_URL,
        "n_ccs_products": len(ccs),
        "n_unknown_target": unk,
        "ccs_targets_reference": sorted(CCS_TARGETS),
        "targets": targets,
    }
    out = Path(args.out) if args.out else cache / "ccs_target_census.json"
    tmp = str(out) + ".part"
    open(tmp, "w", encoding="utf-8").write(json.dumps(doc, ensure_ascii=False, indent=1))
    os.replace(tmp, str(out))
    print("census 写入", out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
