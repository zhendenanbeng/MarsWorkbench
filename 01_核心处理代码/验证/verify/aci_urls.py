# -*- coding: utf-8 -*-
"""用 PDS 检索 API 查出 ACI 图像的确切下载地址。

绕过分析员笔记本的关键：PDS 的检索 API 能按 lid 精确查到一个产品，
返回的 metadata.label_url 就是它在服务器上的真实路径。有了这个，
就不用去笔记本里点来点去，直接把地址粘进浏览器就能下。

注意：API 走 http（https 从本机被 TLS 层挡掉了），主机 pds.nasa.gov 可达。

用法：
    python verify\\aci_urls.py              # 打印 12 张图的地址
    python verify\\aci_urls.py --out 清单.md  # 同时写一份 markdown
"""
from __future__ import annotations

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout

import argparse
import json
import os
import sys
import urllib.error
import urllib.parse
import urllib.request

API = "http://pds.nasa.gov/api/search/1/products"
PREFIX = "urn:nasa:pds:mars2020_imgops:data_aci_imgops:"

# (sol, 靶点, 图像名主干)  —— 主干来自各次观测的 RLS 坐标表
WANTED = [
    ("293", "Quartier · Séítah", "sc3_0293_0692979624_464ecm_n0090000srlc11360_0000lmj"),
    ("293", "Quartier · Séítah", "sc3_0293_0692980965_546ecm_n0090000srlc11360_0000lmj"),
    ("257", "Dourbes · Séítah", "sc3_0257_0689783886_964ecm_n0080000srlc11360_0000lmj"),
    ("257", "Dourbes · Séítah", "sc3_0257_0689785319_285ecm_n0080000srlc11421_0000lmj"),
    ("207", "Garde · Séítah", "sc3_0207_0685346149_136ecm_n0071836srlc10501_0000lmj"),
    ("207", "Garde · Séítah", "sc3_0207_0685347492_097ecm_n0071836srlc10501_0000lmj"),
    ("207", "Garde · Séítah", "sc3_0207_0685348920_929ecm_n0071836srlc11421_0000lmj"),
    ("162", "Guillaumes · Máaz", "sc3_0162_0681350725_253ecm_n0060000srlc10500_0000lmj"),
    ("162", "Guillaumes · Máaz", "sc3_0162_0681351873_164ecm_n0060000srlc11420_0000lmj"),
    ("141", "Foux · Máaz（天然面）", "sc3_0141_0679457520_796ecm_n0051812srlc10203_0000luj"),
    ("141", "Foux · Máaz（天然面）", "sc3_0141_0679457866_796ecm_n0051812srlc10203_0000luj"),
    ("141", "Foux · Máaz（天然面）", "sc3_0141_0679458211_687ecm_n0051812srlc10203_0000luj"),
]


def query(lid, timeout=60):
    """按 lid 查一个产品，返回它的 metadata。"""
    q = urllib.parse.quote('lid eq "%s"' % lid)
    req = urllib.request.Request(API + "?limit=1&q=" + q,
                                 headers={"User-Agent": "sherloc-workbench/1.0",
                                          "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        d = json.loads(r.read().decode("utf-8", "replace"))
    hits = (d.get("summary") or {}).get("hits") or 0
    if not hits:
        return None
    return (d.get("data") or [None])[0]


def label_to_data(label_url, ext=".png"):
    """label 地址换成同名数据文件地址。"""
    base, name = label_url.rsplit("/", 1)
    return base + "/" + name.rsplit(".", 1)[0] + ext


def main(argv=None):
    ap = argparse.ArgumentParser(description="查 ACI 图像的确切下载地址")
    ap.add_argument("--out", default=None, help="同时写一份 markdown 清单")
    ap.add_argument("--download", default=None,
                    help="顺便下载到这个目录（PDS 的原始 .IMG，VICAR 格式）")
    args = ap.parse_args(argv)

    rows, missing = [], []
    for sol, target, stem in WANTED:
        lid = PREFIX + stem
        try:
            p = query(lid)
        except Exception as exc:                       # noqa: BLE001
            missing.append((stem, "查询失败 %s" % str(exc)[:50]))
            continue
        if not p:
            missing.append((stem, "没查到"))
            continue
        props = p.get("properties") or {}
        # 陷阱：记录里有**两个** file_ref——
        #   ops:Label_File_Info.ops:file_ref -> 标签 .xml
        #   ops:Data_File_Info.ops:file_ref  -> 数据 .IMG  ← 我们要的是这个
        # 按尾串匹配会撞上先出现的 Label 那个（踩过一次，下了 12 个 xml）。
        url = ""
        for key in ("ops:Data_File_Info.ops:file_ref",
                    "ops:Label_File_Info.ops:file_ref"):
            v = props.get(key)
            if v:
                url = v[0] if isinstance(v, list) else v
                break
        if not url:
            url = (p.get("metadata") or {}).get("label_url") or ""
        if not url:
            missing.append((stem, "记录里没有 label_url"))
            continue
        rows.append({"sol": sol, "target": target, "stem": stem,
                     "name": url.rsplit("/", 1)[1],
                     "data": url})
        print("sol %-4s %s" % (sol, url.rsplit("/", 1)[1]))

    print("\n共 %d/%d 张查到地址" % (len(rows), len(WANTED)))
    if missing:
        print("没查到的：")
        for stem, why in missing:
            print("   %s  %s" % (stem, why))

    if rows:
        print("\n========== 直接用浏览器打开这些地址就能下载 ==========")
        for r in rows:
            print(r["data"])

    if args.out and rows:
        with open(args.out, "w", encoding="utf-8") as f:
            f.write("# ACI 图像直链（用 PDS 检索 API 查出来的）\n\n")
            f.write("这些是 PDS 自己的元数据里给的 `file_ref`，**直接就是数据文件的地址**，\n"
                    "粘进浏览器地址栏回车即可下载——不需要打开分析员笔记本，也不需要点任何菜单。\n\n")
            f.write("文件是 `.IMG`，**VICAR 格式**（ASCII 标签 + 二进制栅格），"
                    "PIL 直接打不开，但本项目的 `sherloc/vicar.py` 能读，"
                    "`make_mineral_map.py` 已经接上了。\n\n")
            f.write("| sol | 靶点 · 单元 | 直链 |\n| --- | --- | --- |\n")
            for r in rows:
                f.write("| %s | %s | [%s](%s) |\n"
                        % (r["sol"], r["target"], r["name"], r["data"]))
            f.write("\n## 全部地址（方便复制）\n\n```\n")
            for r in rows:
                f.write(r["data"] + "\n")
            f.write("```\n")
        print("\n已写出清单：%s" % args.out)

    if args.download and rows:
        os.makedirs(args.download, exist_ok=True)
        ok = 0
        for r in rows:
            name = r["data"].rsplit("/", 1)[1]
            path = os.path.join(args.download, name)
            if os.path.exists(path):
                print("  已有 %s" % name)
                ok += 1
                continue
            try:
                req = urllib.request.Request(r["data"], headers={"User-Agent": "sherloc-workbench/1.0"})
                with urllib.request.urlopen(req, timeout=300) as resp:
                    data = resp.read()
                with open(path, "wb") as f:
                    f.write(data)
                print("  下载 %s  %.2f MB" % (name, len(data) / 1e6))
                ok += 1
            except Exception as exc:                   # noqa: BLE001
                print("  失败 %s：%s" % (name, str(exc)[:60]))
        print("下载完成：%d/%d -> %s" % (ok, len(rows), args.download))
    return 0 if rows else 1


if __name__ == "__main__":
    sys.exit(main())
