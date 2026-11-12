# -*- coding: utf-8 -*-
"""抽样扫描 SHERLOC 数据目录，估计每次扫描的点数，找出 100 → 25 的转折点。

原理：RRS 的 CSV 每个点一行，点数可以从文件大小反推。
标定（都用 label 里的 <records> 核过）：
    sol 0293 srlc11360w108cgn   100 点（2 段各 50）→ 5432726 B
    sol 1480 srlc11470w104cgn   100 点 → 5535923 B
    sol 1640 srlc12114w104cgn    25 点 → 1476203 B
    sol 1853 srlc12122w104cgn    25 点 → 1464664 B
每点约 54–59 KB（取 55 KB）。这个估计有 ±10% 误差，
只用来区分 25 / 50 / 100 / 1296 几档，不当精确值用。

用法：python 06_算法实验与依据\\scan_points.py [取样间隔] [起始sol] [结束sol]
"""
import re
import sys
import time
import urllib.request

BASE = ("https://pds-geosciences.wustl.edu/m2020/"
        "urn-nasa-pds-mars2020_sherloc/data_processed")

BYTES_PER_POINT = 55000.0

ROW = re.compile(r"(\d+)\s+<A HREF=.([^>]+).>([^<]+)</A>")


def get(url):
    req = urllib.request.Request(url, headers={"User-Agent": "sherloc-workbench/1.0"})
    return urllib.request.urlopen(req, timeout=45).read().decode("utf-8", "replace")


def main(argv):
    step = int(argv[0]) if argv else 8
    lo = int(argv[1]) if len(argv) > 1 else 0
    hi = int(argv[2]) if len(argv) > 2 else 99999

    sols = sorted(set(re.findall(r"sol_(\d+)", get(BASE + "/"))))
    span = [s for s in sols if lo <= int(s) <= hi]      # 先按范围筛，再按步长取
    picks = span[::step]
    print("共 %d 个 sol，取样 %d 个（每 %d 个取一个）" % (len(sols), len(picks), step))
    print("%-7s %s" % ("sol", "RRS 文件的估计点数"))
    for s in picks:
        try:
            h = get("%s/sol_%s/" % (BASE, s))
        except Exception as exc:                     # noqa: BLE001
            print("%-7s FAIL %s" % (s, str(exc)[:40]))
            time.sleep(1.2)
            continue
        est = set()
        for sz, _href, name in ROW.findall(h):
            if name.lower().endswith(".csv") and "rrs" in name.lower():
                est.add(int(round(int(sz) / BYTES_PER_POINT)))
        print("%-7s %s" % (s, "; ".join("~%d 点" % n for n in sorted(est))
                           if est else "（无 RRS）"))
        time.sleep(1.2)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
