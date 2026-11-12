# -*- coding: utf-8 -*-
"""按 sol 从 PDS 下载质控台需要的那几类产品。

只下需要的：RRS（光谱，含波长表）、RLS（激光点坐标）、RLI（激光能量）。
原始产品 ERA/ERB/EPA/ESP 不在默认范围内——那要另外走路线 B。

用法：
    python verify\\fetch_sol.py sol_0257 --dest "C:\\...\\raw\\sol_0257"
    python verify\\fetch_sol.py sol_0257 --kind raw --types ERA,ERB,EPA,ESP
"""
from __future__ import annotations

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout

import argparse
import os
import sys
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = str(_layout.project_layout.APP)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import view                                 # noqa: E402


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


def main(argv=None):
    ap = argparse.ArgumentParser(description="按 sol 下载 SHERLOC 产品")
    ap.add_argument("sol", help="例如 sol_0257")
    ap.add_argument("--kind", default="processed",
                    choices=sorted(view.PDS_KINDS.keys()))
    ap.add_argument("--types", default="RRS,RLS,RLI",
                    help="逗号分隔的产品类型，例如 RRS,RLS,RLI")
    ap.add_argument("--dest", default=None, help="下载到哪个目录")
    ap.add_argument("--limit-mb", type=float, default=60.0,
                    help="单个文件超过这个大小就跳过（防止误下大文件）")
    args = ap.parse_args(argv)

    want = {t.strip().upper() for t in args.types.split(",") if t.strip()}
    dest = args.dest or os.path.join(os.path.dirname(ROOT),
                                     "03_下载数据/SHERLOC/练习数据", "raw", args.sol)
    os.makedirs(dest, exist_ok=True)

    listing = view.list_sol_files(args.kind, args.sol)
    if "error" in listing:
        print("列目录失败：%s" % listing["error"])
        return 1

    print("PDS 目录 %s，共 %d 个文件；要下的是 %s"
          % (listing["url"], listing["total_files"], "、".join(sorted(want))))
    got = skipped = 0
    for grp in listing["groups"]:
        if grp["key"] not in want:
            continue
        for e in grp["files"]:
            mb = e["size"] / 1e6
            path = os.path.join(dest, e["name"])
            if os.path.exists(path) and os.path.getsize(path) == e["size"]:
                print("  已有 %s" % e["name"])
                skipped += 1
                continue
            if mb > args.limit_mb:
                print("  跳过 %s（%.1f MB 超过上限）" % (e["name"], mb))
                skipped += 1
                continue
            try:
                n = download(e["url"], path)
                print("  下载 %s  %.2f MB" % (e["name"], n / 1e6))
                got += 1
            except Exception as exc:                # noqa: BLE001
                print("  失败 %s：%s" % (e["name"], exc))
                return 1
    print("完成：下载 %d 个，跳过 %d 个 -> %s" % (got, skipped, dest))
    return 0


if __name__ == "__main__":
    sys.exit(main())
