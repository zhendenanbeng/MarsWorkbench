# -*- coding: utf-8 -*-
"""命令行连续特征地图与光谱分型；每次运行使用独立产物目录。"""
from __future__ import annotations

import source_paths  # source roots; no-op for bundled imports

import argparse
import json
import os

from sherloc.organize import scan_folder
from sherloc import explore


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("folder", help="数据文件夹")
    ap.add_argument("--scan", type=int, default=0, help="扫描索引，从0开始；先用--list查看")
    ap.add_argument("--list", action="store_true", help="列出扫描并退出")
    ap.add_argument("-o", "--out", help="输出根目录，默认数据目录下_analysis")
    ap.add_argument("--band-a", type=float, nargs=2, default=explore.DEFAULT_BANDS[0])
    ap.add_argument("--band-b", type=float, nargs=2, default=explore.DEFAULT_BANDS[1])
    ap.add_argument("--min-snr", type=float, default=3.)
    ap.add_argument("-k", type=int, default=3, help="请求的光谱分型数，2–8")
    ap.add_argument("--scope", choices=("unknown", "all"), default="unknown")
    args = ap.parse_args(argv)
    _, scans = scan_folder(args.folder)
    if args.list:
        for i, sc in enumerate(scans):
            print(i, "Sol", sc.sol, sc.sequence_id, sc.n_points, "点", sc.all_sclks)
        return 0
    if not 0 <= args.scan < len(scans):
        ap.error("扫描索引超出范围；可先用--list查看")
    try:
        result = explore.run(scans[args.scan], args.out or os.path.join(args.folder,"_analysis"),
                             (args.band_a,args.band_b), args.min_snr, args.k, args.scope,
                             progress=lambda d,t,m: print("[%d/%d] %s"%(d,t,m),flush=True))
    except (ValueError, OSError) as exc:
        ap.error(str(exc))
    print(json.dumps(result,ensure_ascii=False,indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
