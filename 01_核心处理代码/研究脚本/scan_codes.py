# -*- coding: utf-8 -*-
"""从 SHERLOC 的产品总清单里，统计每种扫描代号（srlcXXXXX）出现在哪些 sol。

为什么要这样做：抽样扫目录只能看到取到的那些天。而清单 CSV 里列了**全部**
11195 个产品的 LID，一次请求就能拿到所有 sol 的产品名，
于是每种代号"最早出现在哪一天"就一目了然了。

用法：python 06_算法实验与依据\\scan_codes.py
"""
import collections
import re
import urllib.request

INV = ("https://pds-geosciences.wustl.edu/m2020/urn-nasa-pds-mars2020_sherloc/"
       "data_processed/collection_data_processed_inventory.csv")

LINE = re.compile(r":ss__(\d{4})_.*?rrs__\w+?(srlc\d+)")


def main():
    req = urllib.request.Request(INV, headers={"User-Agent": "sherloc-workbench/1.0"})
    text = urllib.request.urlopen(req, timeout=180).read().decode("utf-8", "replace")
    first = {}
    last = {}
    count = collections.Counter()
    lines = text.replace("\r", "").split("\n")
    for ln in lines:
        m = LINE.search(ln.lower())
        if not m:
            continue
        sol, code = m.group(1), m.group(2)
        count[code] += 1
        if code not in first or sol < first[code]:
            first[code] = sol
        if code not in last or sol > last[code]:
            last[code] = sol
    print("清单共 %d 行，识别出 RRS 产品 %d 条，%d 种扫描代号"
          % (len(lines), sum(count.values()), len(count)))
    print()
    print("%-12s %8s %8s %7s" % ("代号", "最早 sol", "最晚 sol", "产品数"))
    for code, n in sorted(count.items(), key=lambda kv: first[kv[0]]):
        print("%-12s %8s %8s %7d" % (code, first[code], last[code], n))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
