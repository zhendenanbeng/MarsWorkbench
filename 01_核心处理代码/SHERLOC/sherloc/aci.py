# -*- coding: utf-8 -*-
"""ACI 灰度图的查找与下载。

思路：不猜路径。PDS 的检索 API 能按产品编号（lid）精确查到一张图，
返回的记录里有 `ops:Data_File_Info.ops:file_ref`——那就是数据文件本身的直链。

两个容易踩的坑：

1. **要用 http，不是 https。** 实测 `pds.nasa.gov` 与 `pds-imaging.jpl.nasa.gov`
   的 HTTPS 会断（SSL EOF / 连接超时），但 HTTP 正常。
2. **记录里有两个 file_ref。** `Label_File_Info` 那个是 `.xml` 标签，
   `Data_File_Info` 才是 `.IMG` 图像。

下载下来的 `.IMG` 是 VICAR 格式，用 `sherloc.vicar` 读。
"""

from __future__ import annotations

import json
import os
import re
import urllib.error
import urllib.parse
import urllib.request

API = "http://pds.nasa.gov/api/search/1/products"
LID_PREFIX = "urn:nasa:pds:mars2020_imgops:data_aci_imgops:"
UA = {"User-Agent": "sherloc-workbench/1.0"}


def normalize_stem(name):
    """把图像名归一成 lid 用的主干：去目录、去扩展名、去尾部版本号、转小写。

    RLS 里写 `..._0000LMJ03.IMG`，磁盘上的文件可能是 `..._0000LMJ04.IMG`，
    lid 里则是 `..._0000lmj`——剥掉版本号之后三者一致。
    """
    s = os.path.splitext(os.path.basename(name))[0]
    return re.sub(r"\d+$", "", s).lower()


_CACHE = {}


def query(stem, timeout=60, use_cache=True):
    """按主干查一张 ACI 图。返回 {lid, name, url, size} 或 None。"""
    stem = normalize_stem(stem)
    if use_cache and stem in _CACHE:
        return _CACHE[stem]

    q = urllib.parse.quote('lid eq "%s%s"' % (LID_PREFIX, stem))
    url = "%s?limit=1&q=%s" % (API, q)
    try:
        req = urllib.request.Request(url, headers=dict(UA, Accept="application/json"))
        with urllib.request.urlopen(req, timeout=timeout) as r:
            d = json.loads(r.read().decode("utf-8", "replace"))
    except Exception:                                    # noqa: BLE001
        return None

    if not (d.get("summary") or {}).get("hits"):
        return None
    p = (d.get("data") or [None])[0]
    if not p:
        return None
    props = p.get("properties") or {}

    url = ""
    for key in ("ops:Data_File_Info.ops:file_ref",     # 数据文件（.IMG）← 要这个
                "ops:Label_File_Info.ops:file_ref"):   # 标签（.xml）← 兜底
        v = props.get(key)
        if v:
            url = v[0] if isinstance(v, list) else v
            break
    if not url:
        return None
    names = props.get("ops:Data_File_Info.ops:file_name") or []
    sizes = props.get("ops:Data_File_Info.ops:file_size") or []
    out = {"lid": p.get("id") or "",
           "name": (names[0] if names else url.rsplit("/", 1)[1]),
           "url": url,
           "size": int(sizes[0]) if sizes else 0}
    _CACHE[stem] = out
    return out


def needed(folder):
    """扫一个文件夹，列出它需要的 ACI 图（从各次观测的 RLS 坐标表读出来）。

    返回 [{'stem','display','sol','seq'}]，按 sol 排序、同一张图只留一条。
    """
    from .organize import scan_folder

    seen = {}
    try:
        _products, observations = scan_folder(folder)
    except Exception:                                    # noqa: BLE001
        return []
    for ob in observations:
        try:
            pts = ob.read_shot_positions()
        except Exception:                                # noqa: BLE001
            continue
        for p in pts or []:
            img = p.get("image")
            if not img:
                continue
            stem = normalize_stem(img)
            if stem not in seen:
                seen[stem] = {"stem": stem, "display": img, "sol": ob.sol,
                              "seq": ob.sequence_id}
    return sorted(seen.values(), key=lambda d: (d["sol"], d["stem"]))


def have(aci_dir):
    """已经下到这个目录里的图，返回归一化主干的集合。"""
    out = set()
    if not aci_dir or not os.path.isdir(aci_dir):
        return out
    for root, _dirs, files in os.walk(aci_dir):
        for f in files:
            if f.lower().endswith((".img", ".png", ".tif", ".tiff", ".jpg")):
                out.add(normalize_stem(f))
    return out


def download(url, path, chunk=262144, on_bytes=None, cancelled=None):
    """下载一个文件到 path，先写 .part 再改名。返回字节数。

    断流会自动重试（PDS 主机偶发 SSL EOF）。
    """
    tmp = path + ".part"
    last = None
    for attempt in range(3):
        if cancelled and cancelled():
            break
        try:
            req = urllib.request.Request(url, headers=UA)
            n = 0
            with urllib.request.urlopen(req, timeout=180) as r, open(tmp, "wb") as f:
                while True:
                    if cancelled and cancelled():
                        break
                    buf = r.read(chunk)
                    if not buf:
                        break
                    f.write(buf)
                    n += len(buf)
                    if on_bytes:
                        on_bytes(len(buf))
            if cancelled and cancelled():
                break
            os.replace(tmp, path)
            return n
        except Exception as exc:                        # noqa: BLE001
            last = exc
            import time
            time.sleep(2.0 * (attempt + 1))
    if os.path.exists(tmp):
        os.remove(tmp)
    if cancelled and cancelled():
        return 0
    raise last if last else RuntimeError("下载失败")
