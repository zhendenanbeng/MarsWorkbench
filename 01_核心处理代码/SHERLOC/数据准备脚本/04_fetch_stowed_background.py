#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SHERLOC 数据准备脚本 04：下载收纳臂（stowed-arm）扫描，供仪器本底谱构建
==========================================================================
依据 Jakubek et al.（LPSC 2024, #1174；Applied Spectroscopy 2024）：
SHERLOC 光学件（无镀膜 SiO2）在 490/800/1055/1190 cm-1 产生固定本底拉曼。
把机械臂收纳、视场指向火星暗天空时采的扫描就是仪器本底——论文用了
sol 413 / 642 / 678 / 842 / 867 五次收纳臂扫描。

做什么：
  1. 从 PDS 逐 sol 列目录并下载 RRS+RLS+RLI（与 fetch_sol.py 同一套骨架）
  2. 已存在且字节数一致的文件跳过（可断点续跑）
  3. 下载后用 organize.scan_folder 列出每个 sol 的扫描清单
     （sol / 序列 / 点数 / 时间 / 产品类型），供人工圈定哪些是收纳臂扫描

用法：
    python 04_fetch_stowed_background.py                 # 全部五个 sol
    python 04_fetch_stowed_background.py --sol sol_0413  # 只下一个
    python 04_fetch_stowed_background.py --list-only     # 只列目录不下载

依赖：
    需要 numpy 环境（工作台 .venv 即可）。

注意：
  * 联网走环境变量代理：bash 里 `export HTTPS_PROXY=http://127.0.0.1:7890`
    再运行；cmd 里 `set HTTPS_PROXY=...`。
  * 论文只给了 sol 号；同一 sol 若有多个活动，未必全是收纳臂扫描。
    构建本底谱的 05 号脚本会按 800 cm-1 带一致性做稳健剔除，但
    下载后的清单仍值得人工过目一遍。
"""

from __future__ import annotations

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0, str(_BootPath(__file__).resolve().parents[3] / "02_图形界面代码"))
import source_paths as _layout            # noqa: E402

import argparse
import os
import sys
import urllib.request

ROOT = _layout.project_layout.APP
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import view                                # noqa: E402
import project_layout                      # noqa: E402
from sherloc import organize               # noqa: E402

# Jakubek et al. 2024 使用的五次收纳臂扫描所在 sol（PDS 目录用 5 位编号）
STOWED_SOLS = ("sol_00413", "sol_00642", "sol_00678", "sol_00842", "sol_00867")


def download(url, path, chunk=1 << 16):
    req = urllib.request.Request(url, headers={"User-Agent": "sherloc-workbench/1.0"})
    tmp = path + ".part"
    with urllib.request.urlopen(req, timeout=180) as r, open(tmp, "wb") as f:
        total = 0
        while True:
            buf = r.read(chunk)
            if not buf:
                break
            f.write(buf)
            total += len(buf)
    os.replace(tmp, path)
    return total


def list_scans(sol_dir):
    """打印一个 sol 目录里的扫描清单，返回 (扫描数, 总点数)。"""
    _, scans = organize.scan_folder(sol_dir)
    rows = organize.summary_rows(scans)
    for r in rows:
        print("    扫描 %s | %s | 点数 %s | %s | %s"
              % (r["sequence_id"], r["sclk"], r["n_points"], r["time"], r["types"]))
    return len(rows), sum(r["n_points"] or 0 for r in rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description="下载 SHERLOC 收纳臂扫描（仪器本底）")
    ap.add_argument("--sol", action="append", default=[],
                    help="只下这个 sol（可重复）；默认论文五连 sol")
    ap.add_argument("--types", default="RRS,RLS,RLI",
                    help="逗号分隔的产品类型，默认 RRS,RLS,RLI")
    ap.add_argument("--dest", default=None,
                    help="下载根目录；默认 <数据根>/SHERLOC/下载/stowed")
    ap.add_argument("--limit-mb", type=float, default=60.0,
                    help="单个文件超过这个大小就跳过")
    ap.add_argument("--list-only", action="store_true",
                    help="只列 PDS 目录和本地扫描清单，不下载")
    args = ap.parse_args(argv)

    sols = tuple(args.sol) if args.sol else STOWED_SOLS
    want = {t.strip().upper() for t in args.types.split(",") if t.strip()}
    base = Path(args.dest) if args.dest else project_layout.data_root() / "SHERLOC" / "下载" / "stowed"
    base.mkdir(parents=True, exist_ok=True)

    failed = 0
    for sol in sols:
        sol_dir = base / sol
        sol_dir.mkdir(parents=True, exist_ok=True)
        listing = view.list_sol_files("processed", sol)
        if "error" in listing:
            print("[ %s ] 列目录失败：%s" % (sol, listing["error"]))
            failed += 1
            continue
        print("[ %s ] PDS 目录共 %d 个文件；要下 %s"
              % (sol, listing["total_files"], "、".join(sorted(want))))
        if not args.list_only:
            for grp in listing["groups"]:
                if grp["key"] not in want:
                    continue
                for e in grp["files"]:
                    path = str(sol_dir / e["name"])
                    if os.path.exists(path) and os.path.getsize(path) == e["size"]:
                        print("    已有 %s" % e["name"])
                        continue
                    if e["size"] / 1e6 > args.limit_mb:
                        print("    跳过 %s（%.1f MB 超过上限）" % (e["name"], e["size"] / 1e6))
                        continue
                    try:
                        n = download(e["url"], path)
                        print("    下载 %s  %.2f MB" % (e["name"], n / 1e6))
                    except Exception as exc:            # noqa: BLE001
                        print("    失败 %s：%s" % (e["name"], exc))
                        failed += 1
        try:
            n_scans, n_pts = list_scans(str(sol_dir))
            print("    -> 本地 %d 个扫描，共 %s 点" % (n_scans, n_pts))
        except Exception as exc:                        # noqa: BLE001
            print("    扫描清单失败：%s" % exc)

    print("完成%s" % ("（有 %d 个失败项）" % failed if failed else "。全部就绪"))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
