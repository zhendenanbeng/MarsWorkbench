# -*- coding: utf-8 -*-
"""PDS 目录索引抓取：列 sol 目录、列 sol 里的文件。

和 view.py 里的同名能力同一套路（IIS 目录索引 + 正则），但这里按
「数据包根目录」参数化，供多个仪器共用。
"""

from __future__ import annotations

import re
import time
import urllib.request

LISTING_ROW = re.compile(r'(\d+)\s+<A HREF="([^"]+)">([^<]+)</A>', re.I)
SOL_DIR = re.compile(r'href="[^"]*?/(sol_?\d+)/"', re.I)

USER_AGENT = "mars-consoles/1.0"


def http_get(url, timeout=45):
    """GET，带三次重试。PDS 这几台主机偶发 SSL 断流，不重试一次抖动就失败。"""
    last = None
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception as exc:                    # noqa: BLE001
            last = exc
            time.sleep(1.5 * (attempt + 1))
    raise last


def list_sols(root):
    """列出某数据包目录下所有 sol 子目录名（如 'sol00013' / 'sol_00700'）。"""
    html = http_get(root.rstrip("/") + "/").decode("utf-8", "replace")
    sols = SOL_DIR.findall(html)
    return sorted(set(s.lower() for s in sols))


def list_files(root, sol):
    """列出 sol 目录里的文件。返回 [{name, size, url}]；目录穿越在调用方校验。"""
    sol = re.sub(r"[^a-z0-9_]", "", (sol or "").lower())
    if not sol:
        return []
    url = "%s/%s/" % (root.rstrip("/"), sol)
    html = http_get(url).decode("utf-8", "replace")
    out = []
    for size, href, name in LISTING_ROW.findall(html):
        if name.lower().endswith("/"):             # 子目录
            continue
        link = href if href.startswith("http") else \
            "https://pds-geosciences.wustl.edu" + href
        out.append({"name": name, "size": int(size), "url": link})
    return out
