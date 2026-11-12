# -*- coding: utf-8 -*-
"""找并下载 ACI 灰度图。

为什么要爬目录：ACI 图像属于 `mars2020_imgops` 数据集，它挂在 JPL 影像节点上，
页面上 USGS 镜像那条链接是被注释掉的，而直接的 `/data/` 路径实测超时。
与其猜路径，不如从几个可能的根目录往下爬，按文件名找。

用法：
    python verify\\fetch_aci.py --dest "03_下载数据/SHERLOC/练习数据\\_aci"
    python verify\\fetch_aci.py --dest ... --only SC3_0257   # 只找某一批

注意：请求之间会刻意停顿，避免被对方限流（之前请求太密被封过一次）。
"""
from __future__ import annotations

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout

import argparse
import os
import re
import sys
import time
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = str(_layout.project_layout.APP)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import view                                    # noqa: E402

# 几个可能的入口。脚本会逐个试着往下列目录，不假设哪个对。
ROOTS = [
    "https://pds-imaging.jpl.nasa.gov/data/mars2020/mars2020_imgops/",
    "https://pds-imaging.jpl.nasa.gov/data/mars2020_imgops/",
    "https://pds-imaging.jpl.nasa.gov/data/mars2020/",
    "https://pdsimage2.wr.usgs.gov/archive/mars2020/mars2020_imgops/",
    "https://pdsimage.wr.usgs.gov/archive/mars2020/mars2020_imgops/",
]

# 五个科学靶点需要的图（只要主干名，扩展名和版本号后面自己配）
WANTED = [
    "SC3_0293_0692979624_464ECM_N0090000SRLC11360_0000LMJ",
    "SC3_0293_0692980965_546ECM_N0090000SRLC11360_0000LMJ",
    "SC3_0257_0689783886_964ECM_N0080000SRLC11360_0000LMJ",
    "SC3_0257_0689785319_285ECM_N0080000SRLC11421_0000LMJ",
    "SC3_0207_0685346149_136ECM_N0071836SRLC10501_0000LMJ",
    "SC3_0207_0685347492_097ECM_N0071836SRLC10501_0000LMJ",
    "SC3_0207_0685348920_929ECM_N0071836SRLC11421_0000LMJ",
    "SC3_0162_0681350725_253ECM_N0060000SRLC10500_0000LMJ",
    "SC3_0162_0681351873_164ECM_N0060000SRLC11420_0000LMJ",
    "SC3_0141_0679457520_796ECM_N0051812SRLC10203_0000LUJ",
    "SC3_0141_0679457866_796ECM_N0051812SRLC10203_0000LUJ",
    "SC3_0141_0679458211_687ECM_N0051812SRLC10203_0000LUJ",
]

EXTS = (".png", ".PNG", ".tif", ".TIF", ".jpg", ".JPG", ".img", ".IMG")
PAUSE = 1.5                     # 每次请求之间停这么久，别把对方惹毛


def fetch(url, timeout=40):
    req = urllib.request.Request(url, headers={"User-Agent": "sherloc-workbench/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def crawl(root, wanted, dest, max_pages=80, max_depth=4, log=print):
    """从 root 往下爬目录，找到 wanted 里的文件就下载。返回 (找到的, 访问过的页数)。"""
    hit, seen, queue, pages = {}, set(), [(root, 0)], 0
    while queue and pages < max_pages:
        url, depth = queue.pop(0)
        if url in seen:
            continue
        seen.add(url)
        try:
            html = fetch(url).decode("utf-8", "replace")
        except Exception as exc:                       # noqa: BLE001
            log("    ! %s  %s" % (url[:70], str(exc)[:45]))
            time.sleep(PAUSE)
            continue
        pages += 1
        links = re.findall(r'HREF="([^"]+)"', html, re.I)
        for href in links:
            if href.startswith("?") or "Parent Directory" in href:
                continue
            full = urllib.parse.urljoin(url, href)
            if not full.startswith(root):
                continue
            name = urllib.parse.unquote(os.path.basename(full.rstrip("/")))
            if href.endswith("/"):
                if depth < max_depth:
                    queue.append((full, depth + 1))
                continue
            stem = name.rsplit("::", 1)[0]
            for w in wanted:
                if stem.upper().startswith(w.upper()):
                    hit.setdefault(w, []).append(full)
        time.sleep(PAUSE)
    return hit, pages


def main(argv=None):
    ap = argparse.ArgumentParser(description="找并下载 ACI 灰度图")
    ap.add_argument("--dest", default=None, help="下载到哪（默认 03_下载数据/SHERLOC/练习数据\\_aci）")
    ap.add_argument("--only", default=None, help="只在名字含这个串的目标里找，例如 SC3_0257")
    ap.add_argument("--max-pages", type=int, default=80)
    args = ap.parse_args(argv)

    dest = args.dest or os.path.join(os.path.dirname(ROOT),
                                     "03_下载数据/SHERLOC/练习数据", "_aci")
    os.makedirs(dest, exist_ok=True)
    wanted = [w for w in WANTED if not args.only or args.only.upper() in w.upper()]
    print("要找 %d 张图，下载到 %s" % (len(wanted), dest))

    all_hits = {}
    for root in ROOTS:
        print("== 试入口 %s" % root)
        try:
            hits, pages = crawl(root, wanted, dest, max_pages=args.max_pages,
                                log=lambda s: print(s))
        except Exception as exc:                       # noqa: BLE001
            print("   入口不可用：%s" % str(exc)[:60])
            continue
        print("   访问 %d 页，命中 %d 个目标" % (pages, len(hits)))
        for k, v in hits.items():
            all_hits.setdefault(k, []).extend(v)
        if all_hits:
            break

    if not all_hits:
        print("\n没找到。可能原因：入口被限流、路径与预期不同、或该数据集不提供目录列表。")
        print("退路：用浏览器从分析员笔记本手工下载，清单见 _aci\\要下载哪些图.md")
        return 1

    got = 0
    for want, urls in sorted(all_hits.items()):
        url = urls[0]
        name = os.path.basename(urllib.parse.unquote(url))
        path = os.path.join(dest, name)
        try:
            data = fetch(url, timeout=120)
            with open(path, "wb") as f:
                f.write(data)
            print("  下载 %s  %.1f MB" % (name, len(data) / 1e6))
            got += 1
        except Exception as exc:                       # noqa: BLE001
            print("  失败 %s：%s" % (name, str(exc)[:50]))
        time.sleep(PAUSE)
    print("完成：%d/%d 张 -> %s" % (got, len(wanted), dest))
    return 0 if got else 1


if __name__ == "__main__":
    sys.exit(main())
