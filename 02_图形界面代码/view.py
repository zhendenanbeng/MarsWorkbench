#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SHERLOC 光谱工作台。

整个应用就是一个网页，用 Edge 的应用模式开成独立窗口（无地址栏、无标签页）：

    python view.py [起始文件夹] [--port 8765] [--browser] [--no-window]

页面里能做：选数据文件夹、开关预处理步骤、看光谱图、导出 PNG、生成质控台。
本版只做读取 / 预处理 / 作图，不做任何矿物识别。
"""

from __future__ import annotations

import source_paths  # source roots; no-op for bundled imports

import argparse
import html
import io
import json
import os
import re
import sys
import threading
import time
import urllib.parse
import urllib.request
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import appenv
import runtime_http

import numpy as np                                            # noqa: E402
from instruments import notes as paper_notes
from instruments import spectrum_export
from instruments import history
from instruments import download_queue


def alert(title, message):
    """弹一个原生提示框。不需要任何 GUI 库，用 Win32 的 MessageBox。"""
    try:
        import ctypes
        ctypes.windll.user32.MessageBoxW(0, message, title, 0x10)
    except Exception:                                   # noqa: BLE001
        try:
            import tkinter.messagebox as mb
            mb.showerror(title, message)
        except Exception:                               # noqa: BLE001
            print("%s\n%s" % (title, message))


try:
    from sherloc import aci, minerals, preprocess, qc, render, reference  # noqa: E402
    from sherloc.desktop import find_edge, open_app_window, remember_instrument  # noqa: E402
    from sherloc.organize import scan_folder                  # noqa: E402
except ImportError as _exc:                                   # noqa: BLE001
    alert("SHERLOC 光谱工作台 · 启动失败",
          "缺少依赖：%s\n\n请在命令行里执行一次：\n\n"
          "    pip install numpy pandas pillow\n\n然后再双击启动。" % _exc)
    raise SystemExit(1)

COLORS = {"r1": (200, 40, 40), "r2": (30, 90, 210), "r3": (20, 150, 90),
          "sum": (110, 60, 160)}
SETTINGS = appenv.settings('sherloc')
PING_TIMEOUT = 180.0         # 窗口关掉后多久自动退出（秒；浏览器后台会限制定时器频率，留宽一点）

# PDS 上的 SHERLOC 数据根目录。两个来源：原始(product)与处理后(processed)。
PDS_BASE = "https://pds-geosciences.wustl.edu/m2020/urn-nasa-pds-mars2020_sherloc"
PDS_KINDS = {"processed": "data_processed", "raw": "data_raw"}

PRODUCT_LABELS = {
    "RRS": "RRS · 已扣背景的科学光谱", "RRA": "RRA · 处理后主动帧",
    "RRB": "RRB · 处理后暗背景帧", "RCS": "RCS · 定标靶扣背景光谱",
    "RCA": "RCA · 定标靶主动帧", "RCB": "RCB · 定标靶暗背景帧",
    "RCC": "RCC · AlGaN 定标拟合",
    "RLS": "RLS · 激光点坐标", "RLI": "RLI · 激光能量",
    "RMO": "RMO · 位置与强度汇总", "RLP": "RLP · 扫描镜位置表",
    "RM": "RM* · 波段积分强度图",
    "ERA": "ERA · 科学光谱 主动帧（原始）", "ERB": "ERB · 科学光谱 暗背景（原始）",
    "ECA": "ECA · 定标靶 主动帧（原始）", "ECB": "ECB · 定标靶 暗背景（原始）",
    "ESH": "ESH · 仪器健康（科学模式）", "ECH": "ECH · 仪器健康（定标模式）",
    "EPA": "EPA · 激光光电二极管", "ESP": "ESP · 扫描镜位置 (mm)",
    "EXH": "EXH · 扫描镜位置 (像素)", "ERP": "ERP · 车载处理光谱",
}

# 快捷选择：看图够用的最小集合
PRESET_VIEW = {"processed": ["RRS", "RLS", "RLI"], "raw": ["ERA", "ERB", "EPA", "ESH"]}


# --------------------------------------------------------------------------
# 设置
# --------------------------------------------------------------------------

def load_settings():
    try:
        return history.read(SETTINGS)
    except Exception:                                   # noqa: BLE001
        return {}


def save_settings(data):
    try:
        history.update(SETTINGS,lambda current:dict(current,**data))
    except (OSError,ValueError):
        appenv.logger().exception('could not save settings')


# --------------------------------------------------------------------------
# 数据
# --------------------------------------------------------------------------

def scan_listing(folder):
    products, scans = scan_folder(folder)
    # 同 sol 同序列号的多次独立扫描（sol 141 有 3 次）标签完全相同没法挑，
    # 重复的带上 SCLK 尾巴区分
    counts = {}
    for sc in scans:
        k = (sc.sol, sc.sequence_id)
        counts[k] = counts.get(k, 0) + 1
    out = []
    for i, sc in enumerate(scans):
        s = sc.spectra
        dup = counts[(sc.sol, sc.sequence_id)] > 1
        out.append({
            "index": i,
            "sol": sc.sol,
            "sequence_id": sc.sequence_id,
            "sclk": sc.sclk,
            "points": sc.n_points or 0,
            "product": s.product_type if s else "?",
            "split": ("%d 段" % len(sc.parts)) if sc.is_split else "单段",
            "label": "Sol %s · %s · %s 点 · %s%s%s"
                     % (sc.sol, sc.sequence_id, sc.n_points or 0,
                        s.product_type if s else "?",
                        "（%d 段）" % len(sc.parts) if sc.is_split else "",
                        " · %s" % sc.sclk[-6:] if dup else ""),
        })
    # 同一次观测的原始版本（ERA 主动帧）作为独立条目列出——用户可从
    # 暗背景扣除开始自行逐步预处理（organize 把 RAW 版挂在扫描的 alts 上）
    for sc in list(scans):
        for alt in (getattr(sc, "alts", []) or []):
            if not (alt.parts and alt.parts[0].product_type in ("ERA",)):
                continue
            s = alt.spectra
            scans.append(alt)
            out.append({
                "index": len(out),
                "sol": alt.sol,
                "sequence_id": alt.sequence_id,
                "sclk": alt.sclk,
                "points": alt.n_points or 0,
                "product": s.product_type if s else "?",
                "split": ("%d 段" % len(alt.parts)) if alt.is_split else "单段",
                "label": "Sol %s · %s · %s 点 · %s · 原始帧（可从暗背景扣除开始逐步预处理）"
                         % (alt.sol, alt.sequence_id, alt.n_points or 0,
                            s.product_type if s else "?"),
            })
    # 自扣暗输出子文件夹自动并入（07_raw_dark_subtract.py 的产物）：
    # 打开父文件夹即可看到「自扣暗」标记的扫描，不必切换目录。
    extra_dir = os.path.join(folder, "自扣暗_输出")
    if os.path.isdir(extra_dir):
        _, extra_scans = scan_folder(extra_dir)
        for sc in extra_scans:
            s = sc.spectra
            scans.append(sc)
            out.append({
                "index": len(out),
                "sol": sc.sol,
                "sequence_id": sc.sequence_id,
                "sclk": sc.sclk,
                "points": sc.n_points or 0,
                "product": (s.product_type if s else "?") + "·自扣暗",
                "split": ("%d 段" % len(sc.parts)) if sc.is_split else "单段",
                "label": "Sol %s · %s · %s 点 · %s · 自扣暗"
                         % (sc.sol, sc.sequence_id, sc.n_points or 0,
                            s.product_type if s else "?"),
            })
    return products, scans, out


def _count_data(full, files, depth=3):
    """数一个目录里的 PDS4 数据文件，最多往下找 depth 层。"""
    x = sum(1 for f in files if f.lower().endswith(".xml"))
    c = sum(1 for f in files if f.lower().endswith(".csv"))
    if (x and c) or depth <= 0:
        return x, c
    for sub in files[:120]:
        sp = os.path.join(full, sub)
        if not os.path.isdir(sp):
            continue
        try:
            sf = os.listdir(sp)
        except Exception:                               # noqa: BLE001
            continue
        sx, sc = _count_data(sp, sf, depth - 1)
        x += sx
        c += sc
        if x and c:
            break
    return x, c


def browse(path):
    """列出子文件夹，并标出哪些装着 SHERLOC 数据。"""
    path = os.path.abspath(path or os.path.expanduser("~"))
    try:
        names = sorted(os.listdir(path), key=lambda s: s.lower())
    except Exception as exc:                            # noqa: BLE001
        return {"path": path, "parent": None, "dirs": [], "error": str(exc)}
    items = []
    for n in names:
        full = os.path.join(path, n)
        if not os.path.isdir(full) or n.startswith((".", "_", "$")):
            continue
        try:
            files = os.listdir(full)
        except Exception:                               # noqa: BLE001
            continue
        x, c = _count_data(full, files, depth=3)
        items.append({"name": n, "path": full, "xml": x, "csv": c,
                      "data": x > 0 and c > 0})
    items.sort(key=lambda v: (not v["data"], v["name"].lower()))
    return {"path": path, "parent": os.path.dirname(path.rstrip("\\/")) or None,
            "dirs": items, "error": None}


def _slice(arrays, shift):
    if shift is not None:
        mask = (shift >= 200) & (shift <= 4200)
        return shift, mask, "拉曼位移 (cm-1)"
    n = arrays["r1"].shape[1] if arrays["r1"] is not None else 0
    x = np.arange(n)
    return x, np.ones(n, dtype=bool), "SCCD 通道号"


def _series(res, mask, mode):
    out = []
    if mode == "sum":
        # 官方 RM* 产品的口径：完整光谱 = 三个 Region 表相加（它们在波长轴上互补）
        total = None
        for k in ("r1", "r2", "r3"):
            if res[k] is not None:
                total = res[k] if total is None else total + res[k]
        if total is not None:
            out.append({"key": "sum", "label": "完整光谱 R1+R2+R3", "y": total[0][mask]})
        return out
    if mode in ("all", "r1") and res["r1"] is not None:
        out.append({"key": "r1", "label": "Region 1 拉曼", "y": res["r1"][0][mask]})
    if mode in ("all", "fluor"):
        for k, tag in (("r2", "Region 2 荧光"), ("r3", "Region 3 荧光")):
            if res[k] is not None:
                out.append({"key": k, "label": tag, "y": res[k][0][mask]})
    return out


def series_for(scan, steps, point, mode):
    res = preprocess.run_pipeline(scan, steps, points=[point])
    # 点号超出实际可读的谱行数时，管线会静默跳过选行、_series 拿第 0 行
    # ——用户以为在看第 N 点，其实在看第 0 点。这里显式报错
    for k in ("r1", "r2", "r3"):
        if res[k] is not None and res[k].shape[0] != 1:
            raise ValueError("点号 %s 超出该扫描实际可读的 %d 条光谱"
                             % (point, res[k].shape[0]))
    x, mask, xlabel = _slice(res, res["shift"])
    return {"x": x[mask].tolist(), "xlabel": xlabel,
            "series": [{"key": s["key"], "label": s["label"],
                        "y": [None if not np.isfinite(v) else round(float(v), 2)
                              for v in s["y"]]} for s in _series(res, mask, mode)],
            "log": res["log"], "product": res["product"],
            "wl_source": res["wl_source"], "n_points": scan.n_points or 0}


def scan_minerals(state, scan):
    """整个扫描的矿物类别判定（带缓存，拖点号时不重算）。"""
    key = scan.key
    if key in state["mineral_cache"]:
        return state["mineral_cache"][key]
    out = None
    try:
        wl, r1, r2, r3 = scan.read_regions()
        if wl is not None and r1 is not None:
            # 与质控报告、分布图共用入口，避免整谱基线与拉曼窗口口径分叉。
            result = qc.analyse_spectra(wl, r1, r2, r3)
            if result["n_points_with_band"] is not None:
                out = {
                    "groups": result["mineral_groups"],
                    "point_classes": [sorted({c["class"] for c in p["classes"]})
                                      for p in result["point_detail"]],
                    "n_points_with_band": result["n_points_with_band"],
                    "n_points_in_groups": result['n_points_in_groups'],
                    "n_points_with_candidate_bands": result['n_points_with_candidate_bands'],
                    "n_points_with_named_classes": result['n_points_with_named_classes'],
                    "sigma_scan": result["noise_sigma_median"],
                }
    except Exception as exc:
        # 分析失败不等于未检出。失败不缓存，以便文件修复/临时问题消失后可重试。
        return {'error':type(exc).__name__+': '+str(exc)}
    state["mineral_cache"][key] = out
    return out


def mineral_summary(mr, point=None):
    """整理成界面友好的结构。"""
    if not mr:
        return {"ok": False, "classes": [], "point": []}
    if mr.get('error'):
        return {'ok':False,'classes':[],'point':[],'error':mr['error']}
    classes = []
    for g in mr["groups"]:
        if not g["classes"]:
            continue
        c = g["classes"][0]
        classes.append({"class": c["class"], "key": c["key"],
                        "center": round(float(g["band"]["center"]), 1),
                        "snr": round(float(g["band"]["snr"]), 1),
                        "n_points": int(g["n_points"]),
                        "requires_review":bool(g.get('requires_review')),
                        "selection_basis":g.get('selection_basis',''),
                        "evidence": c["evidence"]})
    mine = []
    if point is not None and 0 <= point < len(mr["point_classes"]):
        mine = mr["point_classes"][point]
    return {"ok": True, "classes": classes, "point": mine,
            "n_with_band": mr["n_points_with_band"],
            "n_in_groups": mr.get('n_points_in_groups',mr['n_points_with_band']),
            "n_candidate_points": mr.get('n_points_with_candidate_bands'),
            "n_named_points": mr.get('n_points_with_named_classes'),
            "sigma": round(float(mr["sigma_scan"]), 1)}


def render_scan_png(scan, steps, point, mode, out_path):
    res = preprocess.run_pipeline(scan, steps, points=[point])
    x, mask, xlabel = _slice(res, res["shift"])
    series = [(s["label"], s["y"], COLORS[s["key"]]) for s in _series(res, mask, mode)]
    render.render_png(out_path, series, x[mask], xlabel,
                      "SHERLOC 光谱 · Sol %s · %s · 第 %d 点"
                      % (scan.sol, scan.sequence_id, point),
                      "产品 %s" % res["product"],
                      notes=["%s | %s" % (l["step"], l["text"]) for l in res["log"]])
    return out_path


# --------------------------------------------------------------------------
# 从 PDS 下载
# --------------------------------------------------------------------------

LISTING_ROW = re.compile(r'(\d+)\s+<A HREF="([^"]+)">([^<]+)</A>')


def export_current_csv(scan,steps,point,mode,out):
    if not 0<=point<(scan.n_points or 0):raise ValueError('点号无效')
    res=preprocess.run_pipeline(scan,steps,points=[point])
    if any(res[k] is not None and res[k].shape[0]!=1 for k in ('r1','r2','r3')):raise ValueError('点号超出实际光谱行数')
    x,mask,xlabel=_slice(res,res['shift']);curves=_series(res,mask,mode)
    metadata={'instrument':'sherloc','sol':scan.sol,'sequence':scan.sequence_id,'sclk':scan.sclk,
              'point_index':point,'product':res['product'],'preprocessing':steps,'mode':mode,'xlabel':xlabel,
              'ylabel':'产品处理后的强度（见预处理日志）','log':res['log'],
              'sources':spectrum_export.sources(p.data_path for p in scan.parts)}
    return spectrum_export.write(out,'sherloc_sol%s_%s_%s_pt%d'%(scan.sol,scan.sequence_id,scan.sclk,point),x[mask],[(s['key'],s['y']) for s in curves],metadata)

def _http_get(url, timeout=45):
    # PDS 这几个主机经常偶发 SSL 断流（实测几次里就撞一次），
    # 不重试的话一次抖动就会让整个操作失败。
    last = None
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "sherloc-workbench/1.0"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception as exc:                        # noqa: BLE001
            last = exc
            time.sleep(1.5 * (attempt + 1))
    raise last


def list_sol_dirs(kind):
    """列出 PDS 上某个集合里所有 sol 目录。"""
    html = _http_get("%s/%s/" % (PDS_BASE, PDS_KINDS[kind])).decode("utf-8", "replace")
    sols = re.findall(r'href="[^"]*?/(sol_\d+)/"', html, re.I)
    return sorted(set(s.lower() for s in sols))


def list_sol_files(kind, sol):
    """列出某个 sol 目录里的文件，按产品类型归类。"""
    sol = os.path.basename(sol)                       # 防目录穿越
    if not re.fullmatch(r"sol_\d+", sol):
        return {"error": "sol 名字不合法：%s" % sol}
    url = "%s/%s/%s/" % (PDS_BASE, PDS_KINDS[kind], sol)
    html = _http_get(url).decode("utf-8", "replace")
    entries = []
    for size, href, name in LISTING_ROW.findall(html):
        ext = name.rsplit(".", 1)[-1].lower()
        if ext not in ("csv", "xml"):
            continue
        entries.append({"size": int(size), "url": "https://pds-geosciences.wustl.edu" + href,
                        "name": name, "ext": ext})
    if not entries:
        return {"error": "这个目录里没有 csv/xml 文件（可能还没发布）"}

    groups = {}
    for e in entries:
        n = e["name"].upper()
        key = n[23:26] if len(n) >= 26 else "???"
        if key.startswith("RM") and key[2].isdigit():
            key = "RM"
        g = groups.setdefault(key, {"key": key, "label": PRODUCT_LABELS.get(key, key + " · 其它"),
                                    "files": [], "count": 0, "bytes": 0})
        g["files"].append(e)
        g["count"] += 1
        g["bytes"] += e["size"]
    out = sorted(groups.values(), key=lambda g: (-g["bytes"], g["key"]))
    has_sci = any(k in groups for k in ("RRS", "ERA"))
    preset = [k for k in PRESET_VIEW.get(kind, []) if k in groups] if has_sci else []
    return {"sol": sol, "kind": kind, "url": url, "groups": out,
            "total_files": len(entries),
            "total_bytes": sum(e["size"] for e in entries),
            "preset": preset, "has_science": has_sci}


def _aci_worker(state, folder, dest, stems):
    """后台线程：把缺的 ACI 灰度图从 PDS 拉下来。

    先查一遍每张图的直链（走 PDS 检索 API），再逐个下载。
    """
    st = state["aci"]
    try:
        os.makedirs(dest, exist_ok=True)
        st.update(total=len(stems), done=0, failed=0, out=dest, log=[])
        for stem, display in stems:
            if st.get("cancel"):
                st["log"].append("已取消")
                break
            rec = aci.query(stem)
            if not rec:
                st["failed"] += 1
                st["log"].append("查不到 %s" % display)
                continue
            path = os.path.join(dest, rec["name"])
            if os.path.exists(path) and (not rec["size"] or os.path.getsize(path) == rec["size"]):
                st["done"] += 1
                continue
            try:
                n = aci.download(rec["url"], path,
                                 on_bytes=lambda k: st.__setitem__("bytes_done", st.get("bytes_done", 0) + k),
                                 cancelled=lambda: st.get("cancel"))
                if n:
                    st["done"] += 1
                    st["log"].append("%s  %.1f MB" % (rec["name"], n / 1e6))
                else:
                    break
            except Exception as exc:                    # noqa: BLE001
                st["failed"] += 1
                st["log"].append("失败 %s：%s" % (rec["name"], str(exc)[:70]))
    except Exception as exc:                            # noqa: BLE001
        st["error"] = "%s: %s" % (type(exc).__name__, exc)
    finally:
        st["running"] = False


def _dl_worker(state, kind, sol, types, dest):
    types=frozenset(types)
    download_queue.run(state['dl'],sol,dest,lambda day,out:_dl_one(state,kind,day,types,out))

def _dl_one(state, kind, sol, types, dest):
    dl = state["dl"]
    try:
        info = list_sol_files(kind, sol)
        if "error" in info:
            dl["error"] = info["error"]
            return
        picked = [g for g in info["groups"] if g["key"] in types]
        files = []
        for g in picked:
            files.extend(g["files"])
        # csv 必须带上同名 xml；xml 则单独勾选就够
        picked_names = {f["name"].rsplit(".", 1)[0].lower() for f in files}
        all_entries = [e for g in info["groups"] for e in g["files"]]
        for e in all_entries:
            stem, ext = e["name"].rsplit(".", 1)
            if ext.lower() == "xml" and stem.lower() in picked_names and e not in files:
                files.append(e)

        if not files:dl['error']='没有匹配所选类型的文件';return
        os.makedirs(dest, exist_ok=True)
        files.sort(key=lambda e: (e["ext"].lower() != "xml", e["name"]))
        dl.update(total=len(files), done=0, bytes_total=sum(f["size"] for f in files),
                  bytes_done=0, out=dest, sol=sol, kind=kind, log=[])
        for e in files:
            if dl.get("cancel"):
                dl["log"].append("已取消")
                break
            path = download_queue.target_path(dest,e['name'])
            if os.path.exists(path) and os.path.getsize(path) == e["size"]:
                dl["done"] += 1
                dl["bytes_done"] += e["size"]
                continue
            tmp = path + ".part"
            # 断流就重试同一个文件。PDS 那几台主机偶发 SSL EOF，
            # 不重试的话一次抖动就让整个下载停在这个文件上。
            err = None
            for attempt in range(3):
                if dl.get("cancel"):
                    break
                try:
                    req = urllib.request.Request(
                        e["url"], headers={"User-Agent": "sherloc-workbench/1.0"})
                    got = 0
                    with urllib.request.urlopen(req, timeout=180) as r, open(tmp, "wb") as f:
                        while True:
                            if dl.get("cancel"):
                                break
                            chunk = r.read(262144)
                            if not chunk:
                                break
                            f.write(chunk)
                            got += len(chunk)
                    if dl.get("cancel"):
                        break
                    download_queue.check_length(r,got)
                    dl["bytes_done"] += got
                    os.replace(tmp, path)
                    dl["done"] += 1
                    err = None
                    break
                except Exception as exc:              # noqa: BLE001
                    err = exc
                    dl["log"].append("第 %d 次失败，重试 %s：%s"
                                     % (attempt + 1, e["name"], str(exc)[:50]))
                    time.sleep(2.0 * (attempt + 1))
            if dl.get("cancel"):
                if os.path.exists(tmp):
                    os.remove(tmp)
                break
            if err is not None:
                if os.path.exists(tmp):
                    try:
                        os.remove(tmp)
                    except OSError:
                        pass
                dl["log"].append("放弃 %s：%s" % (e["name"], str(err)[:70]))
        if not dl.get('cancel') and dl.get('done',0)<dl.get('total',0):dl['error']='部分文件未完成'
        if not dl.get("cancel") and not dl.get("error"):
            dl["log"].append("完成：%d 个文件，%.1f MB" % (dl["done"], dl["bytes_done"] / 1048576))
    except Exception as exc:                          # noqa: BLE001
        dl["error"] = "%s: %s" % (type(exc).__name__, exc)


# --------------------------------------------------------------------------
# 页面
# --------------------------------------------------------------------------

PAGE = """<!doctype html><html lang="zh"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>SHERLOC 光谱工作台</title><style>
:root{--bg:#f1f3f6;--panel:#fff;--panel2:#f7f8fa;--ink:#14171b;--muted:#5f6773;
--line:#e0e4ea;--line2:#eef1f5;--accent:#1a63c8;--accent-ink:#fff;
--ok:#12794a;--warn:#a86400;--err:#b3352f;--radius:10px}
@media (prefers-color-scheme:dark){:root{--bg:#0f1113;--panel:#1b1e22;--panel2:#22262b;
--ink:#e9ecf0;--muted:#98a1ac;--line:#31363d;--line2:#262b31;--accent:#78a9ee;
--accent-ink:#0f1113;--ok:#5cc08b;--warn:#e0a63c;--err:#e0796f}}
/* 手动主题切换（页头 ◐ 按钮）：显式指定时覆盖系统偏好 */
html[data-theme="dark"]{--bg:#0f1113;--panel:#1b1e22;--panel2:#22262b;
--ink:#e9ecf0;--muted:#98a1ac;--line:#31363d;--line2:#262b31;--accent:#78a9ee;
--accent-ink:#0f1113;--ok:#5cc08b;--warn:#e0a63c;--err:#e0796f}
html[data-theme="light"]{--bg:#f1f3f6;--panel:#fff;--panel2:#f7f8fa;
--ink:#14171b;--muted:#5f6773;--line:#e0e4ea;--line2:#eef1f5;--accent:#1a63c8;
--accent-ink:#fff;--ok:#12794a;--warn:#a86400;--err:#b3352f}
svg text.tick{fill:var(--muted);font-size:11px}
svg .gridline{stroke:var(--line)}
svg .axl{fill:var(--ink);font-size:12px}
svg .peakline{opacity:.5;stroke-dasharray:4 4}
svg .peaklabel{font-size:10px;font-weight:600}
svg .xhair{stroke:var(--muted);stroke-dasharray:2 3}
.charttip{position:absolute;pointer-events:none;background:var(--panel);
border:1px solid var(--line);border-radius:6px;padding:4px 9px;font-size:11px;
box-shadow:0 2px 10px rgba(0,0,0,.14);white-space:nowrap;z-index:6;color:var(--ink)}
.charttip b{font-weight:600}
.iconbtn{border:1px solid var(--line);background:var(--panel);border-radius:7px;
padding:5px 11px;cursor:pointer;font:inherit;color:var(--ink)}
.iconbtn:hover{border-color:var(--muted)}
button, .iconbtn{transition:border-color .15s,color .15s,background .15s}
button:focus-visible,.iconbtn:focus-visible{outline:2px solid var(--accent);outline-offset:1px}
.lbl{font-size:11px;font-weight:600;margin-bottom:3px}
button.primary{box-shadow:0 1px 2px rgba(26,99,200,.30)}
button.primary:hover{filter:brightness(1.06)}
*{box-sizing:border-box}
html,body{margin:0;height:100%}
body{background:var(--bg);color:var(--ink);
font:14px/1.6 ui-sans-serif,-apple-system,"Segoe UI","Microsoft YaHei",sans-serif}
.top{padding:18px 22px 10px;display:flex;align-items:baseline;gap:14px;flex-wrap:wrap}
.top h1{font-size:19px;font-weight:600;margin:0}
.top .sub{color:var(--muted);font-size:12px}
main{max-width:1280px;margin:0 auto;padding:0 18px 40px}
.card{background:var(--panel);border:1px solid var(--line);border-radius:var(--radius);
padding:12px 14px;margin-bottom:12px}
.bar{display:flex;gap:12px;align-items:center;flex-wrap:wrap}
.bar .grow{flex:1 1 220px;min-width:0}
.lbl{font-size:12px;color:var(--muted);display:block;margin-bottom:2px}
code.path{display:block;background:var(--panel2);padding:5px 9px;border-radius:6px;
font-size:12px;overflow-wrap:anywhere;word-break:break-all}
.chk{display:inline-flex;align-items:center;gap:6px;margin-right:16px;font-size:13px;
cursor:pointer;user-select:none}
button{font:inherit;border:1px solid var(--line);background:var(--panel2);color:var(--ink);
padding:7px 14px;border-radius:7px;cursor:pointer}
button:hover{border-color:var(--muted)}
button.primary{background:var(--accent);color:var(--accent-ink);border-color:transparent;font-weight:500}
button.ghost{background:transparent}
button:disabled{opacity:.45;cursor:default}
select,input[type=text]{font:inherit;padding:6px 9px;border:1px solid var(--line);
border-radius:7px;background:var(--panel);color:var(--ink);max-width:100%;min-width:0}
input[type=range]{width:100%}
#scan{width:100%}
.row{display:grid;grid-template-columns:280px minmax(0,1fr);gap:12px;align-items:start}
.row > *{min-width:0}
.log{margin:0;padding:0;list-style:none;font-size:12px;max-height:220px;overflow:auto}
.log li{display:flex;gap:8px;padding:4px 0;border-bottom:1px solid var(--line2);color:var(--muted)}
.log li b{font-weight:500;color:var(--ink);flex:0 0 76px}
.log li span{min-width:0;overflow-wrap:anywhere;word-break:break-word}
.log li.warn{color:var(--warn)}
.log li.skip{opacity:.7}
svg{display:block;width:100%;height:auto}
.legend{display:flex;flex-wrap:wrap;gap:6px 18px;font-size:12px;
color:var(--muted);margin:0 0 6px}
.legend span{display:inline-flex;align-items:center;gap:6px}
.legend i{width:16px;height:3px;border-radius:2px;display:inline-block;flex:0 0 16px}
.hint{font-size:12px;color:var(--muted);margin:8px 0 0;
overflow-wrap:anywhere;word-break:break-word}
.blist{max-height:340px;overflow:auto;margin-top:8px}
.brow{display:flex;gap:8px;align-items:center;border-bottom:1px solid var(--line2)}
.bname{flex:1;text-align:left;background:none;border:0;padding:7px 6px;border-radius:6px;
overflow-wrap:anywhere}
.bname:hover{background:var(--panel2)}
.buse{background:var(--accent);color:var(--accent-ink);border:0;padding:5px 12px;font-size:13px}
.dot{color:var(--accent);font-size:11px;margin-left:6px}
.prog{height:8px;background:var(--line2);border-radius:4px;overflow:hidden;margin-top:8px}
.prog > i{display:block;height:100%;width:0;background:var(--accent);transition:width .2s}
.msg{font-size:12px;margin-top:6px}
.msg.ok{color:var(--ok)}.msg.err{color:var(--err)}
.typelist{display:grid;grid-template-columns:repeat(auto-fill,minmax(270px,1fr));
gap:2px 18px;margin-top:6px}
.chk.ty{display:flex;align-items:center;gap:6px;margin:0;padding:3px 0}
.chk.ty .sz{margin-left:auto;color:var(--muted);font-size:11px;white-space:nowrap}
.dllog{margin:8px 0 0;padding:8px 10px;background:var(--panel2);border-radius:6px;
font-size:11px;line-height:1.5;max-height:150px;overflow:auto;white-space:pre-wrap;
overflow-wrap:anywhere;color:var(--muted)}
.minbox{margin-top:10px;border-top:1px solid var(--line2);padding-top:8px;font-size:12px}
.minbox .row2{display:flex;gap:10px;align-items:baseline;padding:3px 0;flex-wrap:wrap}
.tagm{display:inline-block;padding:1px 8px;border-radius:3px;font-size:11px;
background:color-mix(in srgb, var(--accent) 16%, transparent);color:var(--accent)}
.minbox .ev{color:var(--muted);font-size:11px}
input[type=checkbox],input[type=range]{accent-color:var(--accent)}
.top{position:sticky;top:0;z-index:20;
background:color-mix(in srgb,var(--bg) 80%,transparent);backdrop-filter:blur(8px);
border-bottom:1px solid var(--line2);padding:10px 22px 8px}
.top h1{font-size:17px;letter-spacing:.2px}
.top h1::before{content:"";display:inline-block;width:9px;height:9px;border-radius:2px;
background:var(--accent);margin-right:9px}
.card{box-shadow:0 1px 2px rgba(16,24,40,.05)}
details{margin:2px 0}
summary{cursor:pointer;font-size:12px;color:var(--muted);user-select:none;
white-space:nowrap}
summary:hover{color:var(--ink)}
.stepper{display:flex;gap:4px}
.stepper button{flex:1;padding:6px 0;background:var(--panel);font-weight:600}
.stepper button:hover{border-color:var(--accent);color:var(--accent)}
.log{max-height:170px}
@media (max-width:900px){.row{grid-template-columns:1fr}}
</style></head><body>
<div class="top">
  <h1>SHERLOC 光谱工作台</h1>
  <span class="sub">读文件 · 预处理 · 出光谱图 · 化学佐证</span>
  <button class="iconbtn" id="themeBtn" title="切换深色/浅色" style="margin-left:auto">◐ 深浅色</button>
</div>
<main>

  <div class="card">
    <div class="bar">
      <div class="grow">
        <span class="lbl">数据文件夹</span>
        <code class="path" id="folderPath">-</code>
      </div>
      <button id="browseBtn">换文件夹…</button>
      <button id="dlOpenBtn">下载数据…</button>
    </div>
    <div class="bar" style="margin-top:10px">
      <span class="lbl" style="margin:0">预处理</span>
      <span id="steps"></span>
      <details style="margin-left:auto">
        <summary>为什么有些步骤会自动跳过？</summary>
        <p class="hint" style="margin:6px 0 0;max-width:560px">
          如果这是 <b>processed</b>（LASER-NORMALIZED）数据，PDS 早就做过<b>激光归一</b>和
          <b>去宇宙射线</b>了，工具会按 label 自动跳过这两步（左边那条 log 会写明依据）。
          <b>死区裁剪</b>和<b>基线扣除</b>不算重复——PDS 从不做这两件事，不扣基线的话
          荧光背景会把几波数宽的拉曼带整个淹掉。
        </p>
      </details>
    </div>
    <div class="bar" style="margin-top:8px">
      <span class="lbl" style="margin:0">显示区间</span>
      <select id="mode">
        <option value="all" selected>全部三个 Region</option>
        <option value="sum">完整光谱（R1+R2+R3 相加）</option>
        <option value="r1">只看 Region 1（拉曼）</option>
        <option value="fluor">只看 Region 2/3（荧光）</option>
      </select>
      <details style="margin-left:auto">
        <summary>导出设置</summary>
        <div class="bar" style="margin-top:8px">
          <div class="grow">
            <span class="lbl">PNG 输出目录</span>
            <input type="text" id="outDir" style="width:100%">
          </div>
        </div>
      </details>
      <button class="primary" id="exportBtn">导出 PNG</button>
      <button id="cancelExportBtn" hidden>取消导出</button>
      <button id="qcBtn">生成质控台</button>
      <button id="mapBtn">生成矿物分布图</button>
      <button id="analysisOpenBtn">特征与分型…</button>
      <button id="openOutBtn" hidden>打开输出目录</button>
      <button class="ghost" id="quitBtn">退出</button>
    </div>
    <div class="prog" id="prog" hidden><i id="progBar"></i></div>
    <div class="msg" id="progMsg"></div>
    <div class="msg" id="qcMsg" hidden></div>
    <div class="msg" id="mapMsg" hidden></div>
  </div>

  <div class="card" id="analysisPanel" hidden>
    <h2 style="font-size:16px;margin:0 0 6px">连续特征地图与光谱分型</h2>
    <p class="hint">分析当前扫描的 Region 1。固定执行死区裁剪、必要的暗背景/激光归一和800–4000 cm⁻¹基线扣除；不跟随预览开关。积分面积不是矿物含量，类型编号不是矿物鉴定。</p>
    <div class="bar">
      <label>A窗口下限 <input id="bandALow" type="number" min="800" max="1800" value="1060" style="width:80px"></label>
      <label>上限 <input id="bandAHigh" type="number" min="800" max="1800" value="1110" style="width:80px"></label>
      <label>B窗口下限 <input id="bandBLow" type="number" min="800" max="1800" value="990" style="width:80px"></label>
      <label>上限 <input id="bandBHigh" type="number" min="800" max="1800" value="1030" style="width:80px"></label>
      <span class="hint">单位 cm⁻¹；同时输出 A/B 面积比</span>
    </div>
    <div class="bar" style="margin-top:8px">
      <label>信噪比门槛 <input id="analysisSnr" type="number" min="0.1" max="100" step="0.1" value="3" style="width:70px"></label>
      <label>请求分型数 <input id="analysisK" type="number" min="2" max="8" step="1" value="3" style="width:60px"></label>
      <label>分型范围 <select id="analysisScope"><option value="unknown">未归类且有合格波段的点</option><option value="all">全部有合格波段的点</option></select></label>
    </div>
    <div class="bar" style="margin-top:8px"><label for="analysisOut">输出根目录</label><input id="analysisOut" class="grow" type="text">
      <button id="analysisStartBtn">分析当前扫描</button><button id="analysisCancelBtn" hidden>取消分析</button>
      <button id="analysisCloseBtn" class="ghost">关闭面板</button></div>
    <p class="hint">请选择多点RRS观测进行分型。低信噪比点显示灰色；比值仅对两个积分均可靠且为正的点计算。无RLS时仍生成数据表和分型报告。每次运行写入独立目录，保留参数与输入校验值。</p>
    <div id="analysisMsg" class="msg" role="status"></div><a id="analysisResultLink" hidden target="_blank">打开分析报告</a>
  </div>

  <div class="card" id="browser" hidden></div>

  <div class="card" id="dlCard" hidden>
    <details open>
    <summary style="font-size:13px;font-weight:600;color:var(--ink)">光谱产品（按火星日下载）</summary>
    <div class="bar" style="margin-top:10px">
      <span class="lbl" style="margin:0">数据来源</span>
      <select id="dlKind">
        <option value="processed">处理后 data_processed（看图、质控用）</option>
        <option value="raw">原始 data_raw（从零做处理用）</option>
      </select>
      <span class="lbl" style="margin:0">火星日</span>
      <select id="dlSol" multiple size="6" title="Ctrl/Shift多选火星日" style="min-width:150px"></select>
      <button id="dlRefresh">刷新列表</button>
      <button class="ghost" id="dlClose" style="margin-left:auto">关闭</button>
    </div>
    <div class="msg" id="dlInfo"></div>
    <div id="dlTypes"></div>
    <div class="bar" style="margin-top:8px">
      <span class="lbl" style="margin:0">快捷选择</span>
      <button id="dlPresetView">看图够用的</button>
      <button id="dlPresetAll">全选</button>
      <button id="dlPresetNone">全不选</button>
      <span class="sub" id="dlSel" style="color:var(--muted);font-size:12px"></span>
    </div>
    <div class="bar" style="margin-top:10px">
      <div class="grow">
        <span class="lbl">下载到</span>
        <input type="text" id="dlDest" style="width:100%">
      </div>
      <button class="primary" id="dlStart">开始下载</button>
      <button id="dlCancel" hidden>取消下载</button>
    </div>
    <div class="prog" id="dlProg" hidden><i id="dlBar"></i></div>
    <div class="msg" id="dlMsg"></div>
    <button id="dlOpenFolder" hidden>打开下载目录</button>
    <pre class="dllog" id="dlLog" hidden></pre>
    </details>

    <details open style="margin-top:12px;border-top:1px solid var(--line2);padding-top:10px">
    <summary style="font-size:13px;font-weight:600;color:var(--ink)">ACI 灰度图（矿物分布图的底图）</summary>
    <div class="bar" style="margin-top:8px">
      <button id="aciScan">扫描缺哪些</button>
      <button class="primary" id="aciStart" disabled>下载缺失的</button>
      <button id="aciCancel" hidden>取消</button>
    </div>
    <div class="msg" id="aciInfo">点「扫描缺哪些」看看这个文件夹需要哪几张 ACI 图。</div>
    <div class="bar" style="margin-top:8px">
      <div class="grow">
        <span class="lbl">下载到</span>
        <input type="text" id="aciDest" style="width:100%">
      </div>
    </div>
    <div class="prog" id="aciProg" hidden><i id="aciBar"></i></div>
    <div class="msg" id="aciMsg"></div>
    <pre class="dllog" id="aciLog" hidden></pre>
    </details>

    <details style="margin-top:12px;border-top:1px solid var(--line2);padding-top:10px">
    <summary style="font-size:13px;font-weight:600;color:var(--ink)">WATSON 上色（给 ACI 底图上真彩，可选）</summary>
    <div class="bar" style="margin-top:8px">
      <button id="watStart">开始上色</button>
      <button id="watCancel" hidden>取消</button>
    </div>
    <div class="msg" id="watInfo">给 ACI 图配同定位的 WATSON 彩色图并配准上色；成败自动判定，不成的图保持灰度。需要 opencv：pip install opencv-python-headless</div>
    <div class="prog" id="watProg" hidden><i id="watBar"></i></div>
    <div class="msg" id="watMsg"></div>
    <pre class="dllog" id="watLog" hidden></pre>
    </details>
  </div>

  <div class="row">
    <div class="card">
      <span class="lbl">扫描（共 <span id="nscan">__NSCAN__</span> 次）</span>
      <select id="scan" style="width:100%"></select>
      <span class="lbl" style="margin-top:10px">点号（←/→ 翻点，Shift ±10）</span>
      <input type="range" id="pt" min="0" max="0" value="0">
      <div class="stepper" style="margin-top:6px">
        <button id="ptBack10" title="后退 10 点">« 10</button>
        <button id="ptPrev" title="上一点（←）">◀</button>
        <button id="ptNext" title="下一点（→）">▶</button>
        <button id="ptFwd10" title="前进 10 点">10 »</button>
      </div>
      <div class="sub" id="ptlabel"
        style="color:var(--muted);font-size:12px;text-align:center;margin-top:4px;
        font-variant-numeric:tabular-nums">-</div>
      <details open style="margin-top:12px">
        <summary>预处理日志</summary>
        <ul class="log" id="log"></ul>
      </details>
    </div>
    <div class="card">
      <div class="bar"><button id="lockBtn">锁定对比</button><span id="lockInfo" class="sub"></span></div>
      <div class="legend" id="legend"></div>
      <div id="plot" style="position:relative"></div>
      <p class="hint" id="hint"></p>
      <div class="minbox" id="minerals" hidden></div>
    </div>
  </div>
</main>
<script>
var STEPS=[];
var $=function(id){return document.getElementById(id);};
var scanSel=$('scan'), pt=$('pt'), plot=$('plot'), logEl=$('log'), hint=$('hint');

function stepsOn(){
  var o={}; STEPS.forEach(function(k){
    var el=document.querySelector('#steps input[data-k="'+k+'"]');
    o[k]=el?el.checked:false; });
  return o;
}
function stepParam(){
  var on=stepsOn();
  return Object.keys(on).filter(function(k){return on[k];}).join(',');
}
function qs(){
  return 'scan='+scanSel.value+'&point='+pt.value+'&mode='+$('mode').value+
         '&steps='+encodeURIComponent(stepParam());
}
function showFolder(p){ $('folderPath').textContent=p||'-';
  if(p){ $('outDir').value=p.replace(/[\\\\/]+$/,'')+'\\\\_png';
    $('analysisOut').value=p.replace(/[\\\\/]+$/,'')+'\\\\_analysis'; } }

fetch('/api/steps').then(function(r){return r.json();}).then(function(list){
  var box=$('steps');
  list.forEach(function(s){
    STEPS.push(s.key);
    var l=document.createElement('label');
    l.className='chk';
    l.innerHTML='<input type="checkbox" data-k="'+s.key+'"'+(s.on?' checked':'')+'>'+s.label;
    l.querySelector('input').addEventListener('change',load);
    box.appendChild(l);
  });
  return fetch('/api/scans').then(function(r){return r.json();});
}).then(function(d){
  showFolder(d.folder); useScans(d.scans);
  if(d.settings&&d.settings.mode){ $('mode').value=d.settings.mode; }
});

function setScans(list){
  PAPER_ID=null;PAPER=null;paperToken++;
  SCANS=list; $('nscan').textContent=list.length; scanSel.innerHTML='';
  list.forEach(function(s){
    var o=document.createElement('option');
    o.value=s.index; o.textContent=s.label; scanSel.appendChild(o);
  });
}
function useScans(list){
  setScans(list);
  if(list.length){ scanSel.value=list[0].index; onScanChange(); }
  else { plot.innerHTML='<p class="hint">这个文件夹里没有找到可用的扫描（需要成对的 .csv 与 .xml）。点上面的「换文件夹」换一个。</p>';
         logEl.innerHTML=''; $('hint').textContent=''; }
}
/* 请求序号：快速连点/连切扫描时，过期响应直接丢弃——服务端每次都跑完整
   预处理管线，先发后至会把界面打回旧状态 */
var loadToken=0;
function onScanChange(){
  VIEW=null;
  var tok=++loadToken;
  $('ptlabel').textContent='加载中…';
  $('hint').textContent='';           /* 顺带清掉上一个扫描的波长表来源 */
  fetch('/api/points?scan='+scanSel.value).then(function(r){return r.json();}).then(function(d){
    if(tok!==loadToken) return;
    pt.max=Math.max(0,d.points-1); pt.value=0; load();
  });
}
scanSel.addEventListener('change',function(){ VIEW=null; onScanChange(); });
pt.addEventListener('input',load);
function stepPt(d){
  var max=parseInt(pt.max,10)||0, v=parseInt(pt.value,10)+d;
  if(v<0)v=0; if(v>max)v=max;
  if(v!==parseInt(pt.value,10)){ pt.value=v; load(); }
}
$('ptPrev').addEventListener('click',function(){stepPt(-1);});
$('ptNext').addEventListener('click',function(){stepPt(1);});
$('ptBack10').addEventListener('click',function(){stepPt(-10);});
$('ptFwd10').addEventListener('click',function(){stepPt(10);});
document.addEventListener('keydown',function(e){
  if(e.ctrlKey||e.metaKey||e.altKey||e.isComposing)return;
  if(workflowKey(e))return;
  var t=e.target&&e.target.tagName;
  if(t==='SELECT'||t==='TEXTAREA') return;
  if(t==='INPUT'){
    /* 点号滑块聚焦时也要接管：原生 range 的 Shift+方向仍是 ±1，
       和界面上写的「Shift ±10」对不上 */
    if(e.target===pt&&(e.key==='ArrowLeft'||e.key==='ArrowRight')){
      stepPt(e.shiftKey?(e.key==='ArrowLeft'?-10:10)
                    :(e.key==='ArrowLeft'?-1:1));
      e.preventDefault();
    }
    return;
  }
  if(e.key==='ArrowLeft'){stepPt(e.shiftKey?-10:-1);e.preventDefault();}
  if(e.key==='ArrowRight'){stepPt(e.shiftKey?10:1);e.preventDefault();}
});
$('mode').addEventListener('change',load);

function load(){
  var tok=++loadToken;
  loadPaperNotes();
  fetch('/api/plot?'+qs()).then(function(r){return r.json();})
   .then(function(d){ if(tok===loadToken)draw(d); });
}
/* 类别配色/缩写从 minerals 注入（单一来源，别在前端再抄一份——磷酸盐那次
   就是两份字典漂移，点阵出透明格） */
var CLASS_COLORS=__CLASS_COLORS__;
var CLASS_SHORT=__CLASS_SHORT__;
var CHART=null, tipEl=null;
var VIEW=null;        /* {x0,x1} 当前缩放视窗，null=全谱 */
var REF=null;
var curData=null;
function resetView(){ VIEW=null; if(curData) draw(curData); }
window.addEventListener('resize',function(){if(curData)draw(curData);});
document.addEventListener('workbench-theme',function(){if(curData)draw(curData);});
$('lockBtn').addEventListener('click',function(){
  if(REF){REF=null;this.textContent='锁定对比';}
  else if(curData&&!curData.error){
    var option=scanSel.options[scanSel.selectedIndex];
    REF=freezeReference(curData,(option?option.textContent:'')+' · 点'+pt.value+' · '+stepParam());this.textContent='取消对比';
  }
  if(curData)draw(curData);
});

function draw(d){
  curData=d;
  if(d.error){
    CHART=null;$('minerals').hidden=true;if(tipEl)tipEl.style.display='none';
    plot.innerHTML='<p class="hint">'+d.error+'</p>';
    /* 出错也要把点号标签从「加载中…」恢复回来，不然永远停在加载中 */
    $('ptlabel').textContent='第 '+pt.value+' 点 / 共 '+(parseInt(pt.max,10)+1)+' 点';
    return;
  }
  var W=Math.max(360,plot.clientWidth||680), H=380, m={l:74,r:22,t:30,b:56};
  var iw=W-m.l-m.r, ih=H-m.t-m.b, n=d.x.length;
  var fullX0=d.x[0], fullX1=d.x[n-1];
  var x0=fullX0, x1=fullX1;
  if(VIEW){
    x0=Math.max(VIEW.x0,fullX0); x1=Math.min(VIEW.x1,fullX1);
    if(!(x1>x0)){ VIEW=null; x0=fullX0; x1=fullX1; }
  }
  var i0v=0, i1v=n-1;
  while(i0v<i1v&&d.x[i0v]<x0)i0v++;
  while(i1v>i0v&&d.x[i1v]>x1)i1v--;
  var all=[];
  var referenceView=referenceSelection(REF,d);
  $('lockInfo').textContent=REF?('参考：'+REF.label+'；'+referenceView.message):'';
  d.series.forEach(function(s){
    for(var vi=i0v;vi<=i1v;vi++){ var vv=s.y[vi];
      if(vv!==null&&isFinite(vv)) all.push(vv); } });
  referenceView.series.forEach(function(se){se.y.forEach(function(v,i){if(v!==null&&isFinite(v)&&REF.x[i]>=x0&&REF.x[i]<=x1)all.push(v);});});
  if(!all.length){ plot.innerHTML='<p class="hint">这个区间里没有可绘制的数据</p>'; return; }
  var lo=Math.min.apply(null,all), hi=Math.max.apply(null,all);
  if(lo===hi){lo-=1;hi+=1;}
  var pad=(hi-lo)*0.08; lo-=pad; hi+=pad;
  var fy=function(v){return m.t+(hi-v)/(hi-lo)*ih;};
  /* x 轴全程按数值线性画（曲线、刻度、峰位虚线同一套 fxw）。
     通道间距经 1/wl 变换后并不均匀（实测 10.4→7.5 cm⁻¹），
     峰位虚线如果按索引等距画会偏离谱峰 20-50 cm⁻¹ */
  var fxw=function(v){return m.l+(v-x0)/Math.max(1e-9,x1-x0)*iw;};
  var fx=function(i){return fxw(d.x[i]);};
  var COL={r1:'#c72828',r2:'#1e5ad2',r3:'#14965a',sum:'#6e3ca0'};

  var s=['<svg viewBox="0 0 '+W+' '+H+'" role="img" aria-label="光谱曲线">'];
  var main=d.series[0], mc=COL[main.key]||'#1a63c8';
  s.push('<defs><linearGradient id="ga" x1="0" y1="0" x2="0" y2="1">'
    +'<stop offset="0" stop-color="'+mc+'" stop-opacity=".20"/>'
    +'<stop offset="1" stop-color="'+mc+'" stop-opacity="0"/></linearGradient></defs>');
  var k, i, v;
  for(k=0;k<=4;k++){
    v=lo+(hi-lo)*k/4; var gy=fy(v);
    s.push('<line class="gridline" x1="'+m.l+'" y1="'+gy.toFixed(1)+'" x2="'+(m.l+iw)+'" y2="'+gy.toFixed(1)+'"/>');
    s.push('<text class="tick" x="'+(m.l-8)+'" y="'+(gy+4).toFixed(1)+'" text-anchor="end">'+v.toFixed(0)+'</text>');
  }
  /* x 刻度取“好看”的整数值（1/2/5×10^n），位置经 fxw——和曲线同一套映射 */
  var tStep=Math.pow(10,Math.floor(Math.log10((x1-x0)/8||1)));
  var tMag=(x1-x0)/8/tStep;
  if(tMag>=5) tStep*=5; else if(tMag>=2) tStep*=2;
  for(var tv=Math.ceil(x0/tStep)*tStep; tv<=x1+1e-9; tv+=tStep){
    s.push('<text class="tick" x="'+fxw(tv).toFixed(1)+'" y="'+(m.t+ih+18)+'" text-anchor="middle">'+tv.toFixed(0)+'</text>');
  }
  s.push('<text class="axl" x="'+(m.l+iw/2)+'" y="'+(m.t+ih+42)+'" text-anchor="middle">'+d.xlabel+'</text>');
  s.push('<text class="axl" x="15" y="'+(m.t+ih/2)+'" text-anchor="middle" transform="rotate(-90 15 '+(m.t+ih/2)+')">强度 (计数)</text>');

  // 扫描级判定带的峰位标记：把"判了什么"画回谱上
  var mn=d.minerals||{};
  if(mn.ok&&mn.classes){
    var prevX=-999, row=0;
    mn.classes.forEach(function(c){
      if(typeof c.center!=='number') return;
      var cx=fxw(c.center);
      if(cx<m.l+4||cx>m.l+iw-4) return;
      var col=CLASS_COLORS[c.key]||'#8a94a0';
      // 标签只留符号；相邻 40px 内自动换行，避免互相压字
      if(cx-prevX<40){ row=(row+1)%2; } else { row=0; }
      var ly=m.t+11+row*15;
      prevX=cx;
      s.push('<line class="peakline" x1="'+cx.toFixed(1)+'" y1="'+m.t+'" x2="'+cx.toFixed(1)+'" y2="'+(m.t+ih)+'" stroke="'+col+'"/>');
      s.push('<text class="peaklabel" x="'+(cx+4).toFixed(1)+'" y="'+ly+'" fill="'+col+'">'+(CLASS_SHORT[c.key]||'?')+'</text>');
    });
  }

  // 序列：主序列带渐变面积，所有线 1.8px 圆角
  var hoverDots=[];
  d.series.forEach(function(se,si){
    var p='', first=-1, last=-1;
    for(i=i0v;i<=i1v;i++){
      v=se.y[i];
      if(v===null||!isFinite(v)){ continue; }
      var X=fx(i).toFixed(1), Y=fy(v).toFixed(1);
      p+=(p?' L':'M')+X+' '+Y;
      if(first<0) first=i;
      last=i;
    }
    if(si===0&&first>=0&&last>first){
      var area='M'+fx(first).toFixed(1)+' '+(m.t+ih)+' '+p.replace(/^M/,'L')
          +' L'+fx(last).toFixed(1)+' '+(m.t+ih)+' Z';
      s.push('<path d="'+area+'" fill="url(#ga)" stroke="none"/>');
    }
    s.push('<path d="'+p+'" fill="none" stroke="'+COL[se.key]
      +'" stroke-width="1.8" stroke-linejoin="round" stroke-linecap="round"/>');
    hoverDots.push({key:se.key,color:COL[se.key],label:se.label});
  });

  referenceView.series.forEach(function(se){
    var path='',active=false;
    REF.x.forEach(function(v,i){var y=se.y[i];if(v<x0||v>x1||y===null||!isFinite(y)){active=false;return;}
      path+=(active?' L':' M')+fxw(v).toFixed(1)+' '+fy(y).toFixed(1);active=true;});
    s.push('<path class="reference-spectrum" d="'+path+'" fill="none" stroke="#87909c" stroke-width="1.6" stroke-dasharray="6 4"/>');
  });
  // 悬停十字线（每次重画重建；事件绑在 plot 容器上，只需绑一次）
  s.push('<line id="xhl" class="xhair" x1="0" y1="'+m.t+'" x2="0" y2="'+(m.t+ih)+'" visibility="hidden"/>');
  hoverDots.forEach(function(hd,si){
    s.push('<circle id="xhd'+si+'" r="3.5" fill="'+hd.color+'" stroke="var(--panel)" stroke-width="1.5" visibility="hidden"/>');
  });
  s.push('<rect x="'+m.l+'" y="'+m.t+'" width="'+iw+'" height="'+ih+'" fill="transparent"/>');
  s.push('</svg>');
  plot.innerHTML=s.join('');
  if(!tipEl){ tipEl=document.createElement('div'); tipEl.className='charttip';
    tipEl.style.display='none'; }
  plot.appendChild(tipEl);
  tipEl.style.display='none';

  CHART={d:d,fx:fx,fxw:fxw,fy:fy,m:m,iw:iw,ih:ih,n:n,x0:x0,x1:x1,i0v:i0v,i1v:i1v,hoverDots:hoverDots,referenceView:referenceView};
  hint.textContent=(VIEW?('已放大 '+x0.toFixed(0)+'–'+x1.toFixed(0)
    +' cm⁻¹（双击谱图还原）　'):'滚轮/框选放大 · 双击还原　')
    +('产品 '+d.product+(d.wl_source?('　波长表来自 '+d.wl_source):'　（没有波长表，按通道号作图）'));

  $('legend').innerHTML=d.series.map(function(se){
    return '<span><i style="background:'+COL[se.key]+'"></i>'+se.label+'</span>';
  }).join('');
  if(referenceView.series.length)$('legend').innerHTML+='<span style="color:#87909c">┄ 参考 '+referenceEscape(REF.label)+'：'+referenceView.series.map(function(se){return referenceEscape(se.label);}).join('、')+'</span>';
  logEl.innerHTML=d.log.map(function(l){
    return '<li class="'+l.kind+'"><b>'+l.step+'</b><span>'+l.text+'</span></li>';}).join('');
  $('ptlabel').textContent='第 '+pt.value+' 点 / 共 '+d.n_points+' 点';
  var box=$('minerals');
  if(mn.error){
    box.textContent='矿物分析未完成：'+mn.error; box.hidden=false;
  } else if(mn.ok && (mn.classes.length || (mn.point&&mn.point.length))){
    var h='<div style="margin-bottom:2px"><b>扫描级矿物类别候选</b>'
      + '<span class="ev">　候选波段点数：'+(mn.n_candidate_points==null?'未统计':mn.n_candidate_points)
      + '；进入峰组：'+mn.n_in_groups+'；获得命名候选：'+(mn.n_named_points==null?'未统计':mn.n_named_points)
      + '；噪声 sigma='+mn.sigma+'。组平均SNR来自筛选后的点，不能作为独立物相确认。</span></div>';
    mn.classes.forEach(function(c){
      var col=CLASS_COLORS[c.key]||'#1a63c8';
      h+='<div class="row2"><span class="tagm" style="background:color-mix(in srgb,'
        + col+' 14%,transparent);color:'+col+'">'+c.class+'</span>'
       + '<span>'+c.center+' cm⁻¹</span>'
       + '<span class="ev">'+c.n_points+' 个点 · SNR '+c.snr
       + (c.requires_review?' · 峰形/完整峰宽待复核':'')+'</span></div>';
    });
    if(mn.point&&mn.point.length)
      h+='<div class="ev" style="margin-top:6px">本点判定：'+mn.point.join('、')+'</div>';
    else
      h+='<div class="ev" style="margin-top:6px">本点没有合格波段</div>';
    box.innerHTML=h; box.hidden=false;
  } else { box.hidden=true; }
}

/* 悬停十字线：竖线 + 各序列圆点 + 浮动读数 */
plot.addEventListener('mousemove',function(e){
  if(!CHART) return;
  var r=plot.getBoundingClientRect();
  var mx=e.clientX-r.left;
  var c=CHART;
  /* 曲线按数值线性画，反查也走数值：像素 → 波数 → 二分找最近通道 */
  var val=c.x0+(mx-c.m.l)/Math.max(1,c.iw)*(c.x1-c.x0);
  var loidx=0, hiidx=c.n-1;
  while(hiidx-loidx>1){
    var mid=(loidx+hiidx)>>1;
    if(c.d.x[mid]<=val) loidx=mid; else hiidx=mid;
  }
  var i=(Math.abs(c.d.x[loidx]-val)<Math.abs(c.d.x[hiidx]-val))?loidx:hiidx;
  if(i<c.i0v)i=c.i0v; if(i>c.i1v)i=c.i1v;
  var xi=c.fxw(c.d.x[i]);
  var xhl=document.getElementById('xhl');
  if(!xhl) return;
  xhl.setAttribute('x1',xi); xhl.setAttribute('x2',xi);
  xhl.setAttribute('visibility','visible');
  var rows=[];
  c.hoverDots.forEach(function(hd,si){
    var v=c.d.series[si].y[i];
    var dot=document.getElementById('xhd'+si);
    if(v===null||!isFinite(v)){ if(dot) dot.setAttribute('visibility','hidden'); return; }
    if(dot){ dot.setAttribute('cx',xi); dot.setAttribute('cy',c.fy(v));
             dot.setAttribute('visibility','visible'); }
    rows.push('<b>'+hd.label+'</b> '+Number(v).toFixed(0));
  });
  var wv=c.d.x[i];
  if(REF&&c.referenceView)c.referenceView.series.forEach(function(se){var value=referenceRead(REF,se,wv);rows.push('参考 '+referenceEscape(se.label)+' '+(value?Number(value.y).toFixed(1)+'（x='+Number(value.x).toFixed(1)+'）':'—'));});
  tipEl.innerHTML='<b>'+Number(wv).toFixed(0)+' cm⁻¹</b><br>'+rows.join('<br>');
  tipEl.style.display='block';
  var tx=xi+14;
  if(tx+170>r.width) tx=xi-176;
  tipEl.style.left=tx+'px';
  tipEl.style.top=Math.max(4,e.clientY-r.top-16)+'px';
});
plot.addEventListener('mouseleave',function(){
  dragSel=null;
  if(tipEl) tipEl.style.display='none';
  var xhl=document.getElementById('xhl');
  if(xhl) xhl.setAttribute('visibility','hidden');
  var dots=(CHART||{}).hoverDots||[];
  dots.forEach(function(_hd,si){ var dot=document.getElementById('xhd'+si);
    if(dot) dot.setAttribute('visibility','hidden'); });
  hideZoomRect();
});

/* ---- 缩放：滚轮以光标为中心放大/缩小；框选放大；双击还原 ---- */
function hideZoomRect(){ var zr=document.getElementById('zoomrect');
  if(zr) zr.remove(); }
plot.addEventListener('wheel',function(e){
  if(!CHART) return;
  e.preventDefault();
  var c=CHART, r=plot.getBoundingClientRect();
  var mx=e.clientX-r.left;
  var frac=Math.max(0,Math.min(1,(mx-c.m.l)/Math.max(1,c.iw)));
  var val=c.x0+frac*(c.x1-c.x0);
  var lim0=Math.min.apply(null,curData.x), lim1=Math.max.apply(null,curData.x);
  var span=(c.x1-c.x0)*(e.deltaY>0?1.3:1/1.3);
  if(span>=(lim1-lim0)){ VIEW=null; draw(curData); return; }
  if(span<(lim1-lim0)*0.004) return;
  var nx0=val-frac*span, nx1=nx0+span;
  VIEW={x0:Math.max(lim0,nx0), x1:Math.min(lim1,nx1)};
  draw(curData);
},{passive:false});
var dragSel=null;
plot.addEventListener('mousedown',function(e){
  if(!CHART||e.button!==0) return;
  var r=plot.getBoundingClientRect();
  if(e.clientX-r.left<CHART.m.l||e.clientX-r.left>CHART.m.l+CHART.iw) return;
  dragSel={x0:e.clientX-r.left};
});
plot.addEventListener('mousemove',function(e){
  if(!dragSel||!CHART) return;
  var r=plot.getBoundingClientRect();
  var mx=e.clientX-r.left;
  hideZoomRect();
  var div=document.createElement('div');
  div.id='zoomrect';
  div.style.cssText='position:absolute;top:'+CHART.m.t+'px;height:'+CHART.ih
    +'px;left:'+Math.min(dragSel.x0,mx)+'px;width:'+Math.abs(mx-dragSel.x0)
    +'px;background:color-mix(in srgb,var(--accent) 12%,transparent);'
    +'border:1px solid var(--accent);pointer-events:none;';
  plot.appendChild(div);
  dragSel.x1=mx;
});
plot.addEventListener('mouseup',function(e){
  if(!dragSel||!CHART){ dragSel=null; return; }
  hideZoomRect();
  var c=CHART;
  if(dragSel.x1!=null&&Math.abs(dragSel.x1-dragSel.x0)>12){
    var xa=c.x0+(Math.min(dragSel.x0,dragSel.x1)-c.m.l)/Math.max(1,c.iw)*(c.x1-c.x0);
    var xb=c.x0+(Math.max(dragSel.x0,dragSel.x1)-c.m.l)/Math.max(1,c.iw)*(c.x1-c.x0);
    if(xb>xa) VIEW={x0:xa, x1:xb};
    draw(curData);
  }
  dragSel=null;
});
plot.addEventListener('dblclick',function(){ resetView(); });

/* 主题切换 */
var rootEl=document.documentElement;
var savedTheme=null;
try{ savedTheme=localStorage.getItem('sherloc_theme'); }catch(e){}
if(savedTheme==='dark'||savedTheme==='light')
  rootEl.setAttribute('data-theme',savedTheme);
$('themeBtn').addEventListener('click',function(){
  var cur=rootEl.getAttribute('data-theme')||
    (window.matchMedia&&matchMedia('(prefers-color-scheme:dark)').matches?'dark':'light');
  var next=(cur==='dark')?'light':'dark';
  rootEl.setAttribute('data-theme',next);
  try{ localStorage.setItem('sherloc_theme',next); }catch(e){}
});


/* ---- 文件夹浏览 ---- */
$('browseBtn').addEventListener('click',function(){ openBrowser(''); });
function openBrowser(p){
  fetch('/api/browse?path='+encodeURIComponent(p||'')).then(function(r){return r.json();}).then(function(d){
    var box=$('browser'); box.hidden=false;
    var cur=d.path.replace(/"/g,'&quot;');
    var h='<div class="bar"><span class="sub" style="color:var(--muted);font-size:12px">浏览：<code>'+d.path+'</code></span>'+
          '<button id="useCurrent" data-p="'+cur+'">用当前目录</button>'+
          '<button class="ghost" id="closeBrowse" style="margin-left:auto">关闭</button></div>'+
          '<div class="bar" style="margin-top:8px"><label for="browserPath">目录路径</label>'+
          '<input id="browserPath" type="text" class="grow" placeholder="粘贴文件夹完整路径">'+
          '<button id="browsePathBtn">前往</button><button id="usePathBtn">打开此路径</button></div>'+
          '<p id="browserError" class="msg err" role="alert" hidden></p><div class="blist">';
    if(d.parent) h+='<div class="brow"><button class="bname" data-p="'+d.parent+'">↑ 上一级</button></div>';
    if(d.error) h+='<p class="hint">打不开这个目录：'+d.error+'</p>';
    d.dirs.forEach(function(x){
      var p2=x.path.replace(/"/g,'&quot;');
      h+='<div class="brow"><button class="bname" data-p="'+p2+'">'+x.name+
         (x.data?('<span class="dot">有数据 '+x.xml+' 个标签</span>'):'')+'</button>';
      if(x.data) h+='<button class="buse" data-p="'+p2+'">用这个</button>';
      h+='</div>';
    });
    if(!d.dirs.length&&!d.error) h+='<p class="hint">这里没有子文件夹——上面的「用当前目录」可以直接选它。</p>';
    box.innerHTML=h+'</div>';
    loadRecentFolders(box);
    $('browserPath').value=d.path;
    $('browsePathBtn').addEventListener('click',function(){ openBrowser($('browserPath').value.trim()); });
    $('usePathBtn').addEventListener('click',function(){ loadFolder($('browserPath').value.trim()); });
    $('browserPath').addEventListener('keydown',function(e){
      if(e.key==='Enter'){ e.preventDefault(); openBrowser(this.value.trim()); }
    });
    $('useCurrent').addEventListener('click',function(){ loadFolder(this.dataset.p); });
    $('closeBrowse').addEventListener('click',function(){ box.hidden=true; });
    Array.prototype.forEach.call(box.querySelectorAll('.bname'),function(b){
      b.addEventListener('click',function(){ openBrowser(b.dataset.p); });
    });
    Array.prototype.forEach.call(box.querySelectorAll('.buse'),function(b){
      b.addEventListener('click',function(){ loadFolder(b.dataset.p); });
    });
  });
}
function loadFolder(p){
  fetch('/api/load?path='+encodeURIComponent(p)).then(function(r){return r.json();}).then(function(r){
    if(r.error){
      var err=$('browserError');
      if(err){ err.textContent=r.error; err.hidden=false; } else alert(r.error);
      return;
    }
    $('browser').hidden=true; showFolder(r.folder); useScans(r.scans);
    fetch('/api/settings',{method:'POST',body:JSON.stringify({folder:r.folder,mode:$('mode').value})});
  });
}

/* ---- 导出 ---- */
var timer=null, SCANS=[];
$('exportBtn').addEventListener('click',function(){
  var out=$('outDir').value.trim();
  if(!out){ alert('请先填 PNG 输出目录'); return; }
  var total=SCANS.reduce(function(a,s){return a+(s.points||0);},0);
  if(!confirm('将按当前预处理开关和显示区间，导出全部 '+SCANS.length+
              ' 次扫描，共 '+total+' 张 PNG。\\n导出过程中可点「取消导出」中止。')) return;
  $('exportBtn').disabled=true; $('cancelExportBtn').hidden=false;
  $('prog').hidden=false; $('progBar').style.width='0%';
  $('progMsg').className='msg'; $('progMsg').textContent='正在导出…';
  /* steps/mode 必须带上：导出要和预览同一套处理，否则所见非所得 */
  fetch('/api/export?out='+encodeURIComponent(out)+
        '&mode='+$('mode').value+'&steps='+encodeURIComponent(stepParam()))
   .then(function(r){return r.json();}).then(function(r){
    if(r.error){ finishExport({running:false,error:r.error}); return; }
    timer=setInterval(poll,500);
  });
});
$('cancelExportBtn').addEventListener('click',function(){
  $('cancelExportBtn').disabled=true;
  fetch('/api/export/cancel').then(function(r){return r.json();}).then(function(r){
    if(!r.cancelling){ $('cancelExportBtn').hidden=true; $('cancelExportBtn').disabled=false; }
  });
});
function poll(){ fetch('/api/status').then(function(r){return r.json();}).then(function(s){
  var e=s.export||{};
  if(e.running){
    var pct=e.total?Math.round(e.done/e.total*100):0;
    $('progBar').style.width=pct+'%';
    $('progMsg').textContent='正在导出 '+e.done+' / '+e.total+' 张…';
  } else { clearInterval(timer); timer=null; finishExport(e); }
});}
function finishExport(e){
  $('exportBtn').disabled=false; $('cancelExportBtn').hidden=true;
  $('cancelExportBtn').disabled=false;
  if(e.error){ $('progMsg').className='msg err'; $('progMsg').textContent='导出中止/失败：'+e.error; return; }
  $('progBar').style.width='100%';
  $('progMsg').className='msg ok';
  $('progMsg').textContent='导出完成，共 '+e.done+' 张：'+e.out;
  $('openOutBtn').hidden=false;
  $('openOutBtn').onclick=function(){ fetch('/api/openfolder?path='+encodeURIComponent(e.out)); };
}

/* ---- 质控台 ---- */
$('qcBtn').addEventListener('click',function(){
  $('qcBtn').disabled=true; $('qcMsg').hidden=false; $('qcMsg').className='msg';
  $('qcMsg').textContent='正在生成质控台…（要跑完全部扫描，可能几分钟）';
  fetch('/api/qc').then(function(r){return r.json();}).then(function(r){
    $('qcBtn').disabled=false;
    if(r.error){ $('qcMsg').className='msg err'; $('qcMsg').textContent='质控台失败：'+r.error; return; }
    $('qcMsg').className='msg ok';
     $('qcMsg').textContent='质控台已生成，正在新窗口打开：'+r.index;
   });
});

/* ---- 生成矿物分布图（把逐点矿物类别画到 RLS 真实坐标上）---- */
$('mapBtn').addEventListener('click',function(){
  $('mapBtn').disabled=true;
  $('mapMsg').hidden=false; $('mapMsg').className='msg';
  $('mapMsg').textContent='正在生成矿物分布图…（有 ACI 底图时会慢一点）';
  var aci=$('aciDest').value.trim();   // 下载面板里填过的 ACI 目录，没有就让后端自动找
  fetch('/api/maps?aci='+encodeURIComponent(aci))
   .then(function(r){return r.json();}).then(function(r){
     $('mapBtn').disabled=false;
     if(r.error){
       $('mapMsg').className='msg err';
       $('mapMsg').textContent='生成失败：'+r.error;
       if(r.log){ $('dlLog').hidden=false; $('dlLog').textContent=r.log; }
       return;
     }
     $('mapMsg').className='msg ok';
     $('mapMsg').textContent='矿物分布图已生成：'+r.files+' 张（'+r.count+' 次观测） → '+r.out
     + (r.aci? '（底图：'+r.aci+'）':'（无底图模式）');
     if(r.errors&&r.errors.length){
       $('mapMsg').textContent+='；'+r.errors.length+' 次观测出图失败（详见命令行日志）';
     }
     $('openOutBtn').hidden=false;
     $('openOutBtn').onclick=function(){ fetch('/api/openfolder?path='+encodeURIComponent(r.out)); };
   })
   .catch(function(e){ $('mapBtn').disabled=false;
     $('mapMsg').className='msg err'; $('mapMsg').textContent='生成失败：'+e; });
});

/* 连续特征与探索性分型：独立后台任务，不阻塞看谱。 */
var analysisTimer=null;
function renderAnalysis(s){
  if(!s) return;
  $('analysisStartBtn').disabled=!!s.running;
  $('analysisCancelBtn').hidden=!s.running;
  $('analysisCancelBtn').disabled=!!s.cancel;
  $('analysisResultLink').hidden=true;
  var msg=$('analysisMsg');
  if(s.running){ msg.className='msg'; msg.textContent=(s.cancel?'正在取消，等待当前计算段结束…':(s.message||'正在分析…'))+' ['+(s.done||0)+'/'+(s.total||6)+']'; }
  else if(s.error){ msg.className='msg err'; msg.textContent=(s.cancelled?'已取消：':'分析失败：')+s.error; }
  else if(s.result){
    msg.className='msg ok';
    var label={ok:'完成分型',insufficient:'候选不足，未分型',single_pattern:'仅一种谱形'}[s.result.type_status]||s.result.type_status;
    msg.textContent='完成：'+s.result.points+'点，'+s.result.maps+'张空间图；分型候选'+s.result.candidate_points+'点，类型数'+s.result.types+'（'+label+'）。输出：'+s.result.out;
    $('analysisResultLink').href='/analysis-result/'+s.id+'/index.html';
    $('analysisResultLink').hidden=false;
  }
  if(!s.running&&analysisTimer){clearInterval(analysisTimer);analysisTimer=null;}
}
function pollAnalysis(){fetch('/api/status').then(function(r){return r.json();}).then(function(s){renderAnalysis(s.analysis);}).catch(function(e){$('analysisMsg').textContent='无法读取任务状态：'+e.message;});}
$('analysisOpenBtn').addEventListener('click',function(){
  $('analysisPanel').hidden=false;pollAnalysis();
  if(!analysisTimer) analysisTimer=setInterval(pollAnalysis,700);
});
$('analysisCloseBtn').addEventListener('click',function(){$('analysisPanel').hidden=true;});
$('analysisStartBtn').addEventListener('click',function(){
  var payload={scan:Number(scanSel.value),bands:[[Number($('bandALow').value),Number($('bandAHigh').value)],[Number($('bandBLow').value),Number($('bandBHigh').value)]],min_snr:Number($('analysisSnr').value),k:Number($('analysisK').value),scope:$('analysisScope').value,out:$('analysisOut').value.trim()};
  $('analysisStartBtn').disabled=true;$('analysisResultLink').hidden=true;
  fetch('/api/analysis/start',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(payload)})
    .then(function(r){return r.json();}).then(function(r){
      if(r.error){$('analysisStartBtn').disabled=false;$('analysisMsg').className='msg err';$('analysisMsg').textContent=r.error;return;}
      pollAnalysis();if(!analysisTimer)analysisTimer=setInterval(pollAnalysis,700);
    }).catch(function(e){$('analysisStartBtn').disabled=false;$('analysisMsg').textContent='分析启动失败：'+e.message;});
});
$('analysisCancelBtn').addEventListener('click',function(){fetch('/api/analysis/cancel',{method:'POST'}).then(pollAnalysis);});

/* ---- 退出 ---- */
$('quitBtn').addEventListener('click',function(){
  if(!confirm('退出后需要重新双击「启动.bat」才能再用。确定退出吗？')) return;
  fetch('/api/quit').finally(function(){ document.body.innerHTML=
    '<div style="padding:60px;text-align:center;color:#5f6773;font-size:15px">'+
    '已退出，这个窗口可以关闭了。</div>'; });
});

/* ---- 从 PDS 下载数据 ---- */
var DL = {last:null, timer:null};
function fmtMB(b){ return (b/1048576).toFixed(b>10485760?0:1)+' MB'; }
var DEFAULT_DOWNLOAD_ROOT='';
fetch('/api/paths?instrument='+'sherloc').then(r=>r.json()).then(d=>{DEFAULT_DOWNLOAD_ROOT=d.download;}).catch(()=>{});
function dlDefaultDest(sol,kind){
  if(selectedSols().length>1)sol='批量';
  return DEFAULT_DOWNLOAD_ROOT ? DEFAULT_DOWNLOAD_ROOT+'/'+sol+(kind==='raw'?'_raw':'') : '';
}
$('dlOpenBtn').addEventListener('click',function(){
  $('dlCard').hidden=false;
  if($('dlSol').options.length) return;
  loadSols(false);
});
$('dlClose').addEventListener('click',function(){ $('dlCard').hidden=true; });
$('dlKind').addEventListener('change',function(){ $('dlSol').innerHTML=''; loadSols(false); });
$('dlRefresh').addEventListener('click',function(){ $('dlSol').innerHTML=''; loadSols(true); });
$('dlSol').addEventListener('change',loadTypes);

function loadSols(force){
  var kind=$('dlKind').value;
  $('dlInfo').className='msg'; $('dlInfo').textContent='正在读取 PDS 上的目录列表…';
  fetch('/api/dl/sols?kind='+kind+(force?'&refresh=1':''))
   .then(function(r){return r.json();}).then(function(d){
     if(d.error){ $('dlInfo').className='msg err'; $('dlInfo').textContent=d.error; return; }
     var sel=$('dlSol'); sel.innerHTML='';
     d.sols.forEach(function(s){
       var o=document.createElement('option'); o.value=s; o.textContent=s; sel.appendChild(o);
     });
     // 默认挑最后几个（较新的）火星日
     if(d.sols.length) sel.value=d.sols[d.sols.length-1];
     $('dlInfo').className='msg';
     $('dlInfo').textContent='PDS 上共 '+d.sols.length+' 个火星日可选';
     loadTypes();
   }).catch(function(e){ $('dlInfo').className='msg err';
     $('dlInfo').textContent='读取失败：'+e; });
}
function loadTypes(){
  var kind=$('dlKind').value, sol=$('dlSol').value;
  if(!sol) return;
  $('dlTypes').innerHTML='<span class="lbl">正在读取该火星日的文件清单…</span>';
  fetch('/api/dl/list?kind='+kind+'&sol='+sol).then(function(r){return r.json();}).then(function(d){
    if(d.error){ $('dlTypes').innerHTML='<span class="msg err">'+d.error+'</span>'; return; }
    DL.last=d;
    var preset=d.preset||[];
    var h='<div class="typelist">';
    d.groups.forEach(function(g){
      var on=preset.indexOf(g.key)>=0;
      h+='<label class="chk ty"><input type="checkbox" data-k="'+g.key+'" data-b="'+g.bytes+'" data-n="'+g.count+'"'
       + (on?' checked':'')+'>'+g.label+'<span class="sz">'+g.count+' 个 · '+fmtMB(g.bytes)+'</span></label>';
    });
    $('dlTypes').innerHTML=h+'</div>';
    Array.prototype.forEach.call($('dlTypes').querySelectorAll('input'),
      function(i){ i.addEventListener('change',updateDlSel); });
    $('dlInfo').className='msg';
    $('dlInfo').textContent = sol+'：共 '+d.total_files+' 个文件，'+fmtMB(d.total_bytes)
      + (d.has_science? '' : '　⚠ 这个火星日没有科学光谱（只有定标数据）');
    if(selectedSols().length>1)$('dlInfo').textContent+='；已选'+selectedSols().length+'天，类型清单取首个Sol，队列逐天筛选';
    $('dlDest').value=dlDefaultDest(sol, kind);
    updateDlSel();
  });
}
function dlBoxes(){ return $('dlTypes').querySelectorAll('input[type=checkbox]'); }
function updateDlSel(){
  var n=0,b=0,f=0;
  Array.prototype.forEach.call(dlBoxes(),function(i){
    if(i.checked){ n++; b+=Number(i.dataset.b); f+=Number(i.dataset.n); }
  });
  $('dlSel').textContent = (selectedSols().length>1?'首个Sol：':'')+(n? ('选中 '+n+' 类，'+f+' 个文件，'+fmtMB(b)) : '尚未选择');
}
$('dlPresetView').addEventListener('click',function(){
  var preset=(DL.last&&DL.last.preset)||[];
  Array.prototype.forEach.call(dlBoxes(),function(i){ i.checked=preset.indexOf(i.dataset.k)>=0; });
  updateDlSel();
});
$('dlPresetAll').addEventListener('click',function(){
  Array.prototype.forEach.call(dlBoxes(),function(i){ i.checked=true; }); updateDlSel();
});
$('dlPresetNone').addEventListener('click',function(){
  Array.prototype.forEach.call(dlBoxes(),function(i){ i.checked=false; }); updateDlSel();
});
$('dlStart').addEventListener('click',function(){
  var sols=selectedSols();if(!sols.length){alert('请至少选择一个火星日');return;}
  var types=[];
  Array.prototype.forEach.call(dlBoxes(),function(i){ if(i.checked) types.push(i.dataset.k); });
  if(!types.length){ alert('请至少勾选一类产品'); return; }
  var dest=$('dlDest').value.trim();
  if(!dest){ alert('请填下载目录'); return; }
  $('dlStart').disabled=true; $('dlCancel').hidden=false;
  $('dlOpenFolder').hidden=true;
  $('dlProg').hidden=false; $('dlBar').style.width='0%';
  $('dlMsg').className='msg'; $('dlMsg').textContent='正在连接 PDS…';
  fetch('/api/dl/start?kind='+$('dlKind').value+'&sols='+encodeURIComponent(sols.join(','))
        +'&types='+encodeURIComponent(types.join(','))+'&dest='+encodeURIComponent(dest))
   .then(function(r){return r.json();}).then(function(r){
     if(r.error){ $('dlStart').disabled=false; $('dlCancel').hidden=true;
       $('dlMsg').className='msg err'; $('dlMsg').textContent=r.error; return; }
     if(!DL.timer) DL.timer=setInterval(pollDl,600);
   });
});
$('dlCancel').addEventListener('click',function(){ fetch('/api/dl/cancel'); });
function pollDl(){
  fetch('/api/status').then(function(r){return r.json();}).then(function(s){
    var d=s.dl||{};
    if(d.sol_total){queueStatus(d);return;}
    if(d.running){
      var pct=d.bytes_total? Math.round(d.bytes_done/d.bytes_total*100) : 0;
      $('dlBar').style.width=pct+'%';
      $('dlMsg').className='msg';
      $('dlMsg').textContent='正在下载 '+d.done+' / '+d.total+' 个文件，'
        + fmtMB(d.bytes_done) + ' / ' + fmtMB(d.bytes_total) + '（'+pct+'%）';
      if(d.log&&d.log.length){ $('dlLog').hidden=false; $('dlLog').textContent=d.log.join('\\n'); }
    } else {
      clearInterval(DL.timer); DL.timer=null;
      $('dlStart').disabled=false; $('dlCancel').hidden=true;
      if(d.error){ $('dlMsg').className='msg err'; $('dlMsg').textContent='下载失败：'+d.error; return; }
      if(d.total){
        $('dlBar').style.width='100%';
        $('dlMsg').className='msg ok';
        $('dlMsg').textContent='下载完成：'+d.done+' 个文件 → '+d.out;
        $('openOutBtn').hidden=false;
        $('openOutBtn').onclick=function(){ fetch('/api/openfolder?path='+encodeURIComponent(d.out)); };
        $('dlMsg').textContent += '　（点上面的「打开输出目录」查看）';
      }
      if(d.log&&d.log.length){ $('dlLog').hidden=false; $('dlLog').textContent=d.log.join('\\n'); }
    }
  });
}

/* ---- 下载 ACI 灰度图（矿物分布图的底图）---- */
var ACI={timer:null};
function aciFolder(){ var f=$('folderPath').textContent.trim(); return (f&&f!=='-')?f:''; }
function aciDefaultDest(){
  var f=aciFolder().replace(/[\\\\/]+$/,'');
  return (f? f : '.') + '\\\\_aci';
}
$('aciScan').addEventListener('click',function(){
  var dest=$('aciDest').value.trim()||aciDefaultDest();
  $('aciDest').value=dest;
  $('aciInfo').className='msg'; $('aciInfo').textContent='正在扫这个文件夹的坐标表…';
  $('aciStart').disabled=true;
  fetch('/api/aci/scan?folder='+encodeURIComponent(aciFolder())+'&dest='+encodeURIComponent(dest))
   .then(function(r){return r.json();}).then(function(d){
     if(d.error){ $('aciInfo').className='msg err'; $('aciInfo').textContent=d.error; return; }
     $('aciInfo').className='msg';
     $('aciInfo').textContent='这个文件夹共需要 '+d.total+' 张 ACI 图：本地已有 '+d.have
       +' 张，缺 '+d.missing+' 张。'+(d.missing? '点「下载缺失的」补齐。' : '已经齐了。');
     if(d.missing){
       $('aciStart').disabled=false; $('aciStart').textContent='下载缺失的 '+d.missing+' 张';
       $('aciLog').hidden=false;
       $('aciLog').textContent=d.items.map(function(x){
         return 'sol '+x.sol+'  '+x.name; }).join('\\n');
     } else { $('aciStart').textContent='下载缺失的'; $('aciLog').hidden=true; }
   }).catch(function(e){ $('aciInfo').className='msg err'; $('aciInfo').textContent='扫描失败：'+e; });
});
$('aciStart').addEventListener('click',function(){
  var dest=$('aciDest').value.trim();
  if(!dest){ alert('请填下载目录'); return; }
  $('aciStart').disabled=true; $('aciCancel').hidden=false;
  $('aciProg').hidden=false; $('aciBar').style.width='0%';
  $('aciMsg').className='msg'; $('aciMsg').textContent='正在连 PDS 查地址…';
  fetch('/api/aci/start?folder='+encodeURIComponent(aciFolder())+'&dest='+encodeURIComponent(dest))
   .then(function(r){return r.json();}).then(function(r){
     if(r.error){ $('aciStart').disabled=false; $('aciCancel').hidden=true;
       $('aciMsg').className='msg err'; $('aciMsg').textContent=r.error; return; }
     if(!ACI.timer) ACI.timer=setInterval(pollAci,700);
   });
});
$('aciCancel').addEventListener('click',function(){ fetch('/api/aci/cancel'); });
function pollAci(){
  fetch('/api/status').then(function(r){return r.json();}).then(function(s){
    var d=s.aci||{};
    if(d.running){
      var pct=d.total? Math.round(d.done/d.total*100):0;
      $('aciBar').style.width=pct+'%';
      $('aciMsg').className='msg';
      $('aciMsg').textContent='正在下载 '+d.done+' / '+d.total+' 张'
        + (d.failed? '（失败 '+d.failed+'）':'') + '，已收 '+fmtMB(d.bytes_done||0);
      if(d.log&&d.log.length){ $('aciLog').hidden=false; $('aciLog').textContent=d.log.join('\\n'); }
    } else {
      clearInterval(ACI.timer); ACI.timer=null;
      $('aciStart').disabled=false; $('aciCancel').hidden=true;
      if(d.error){ $('aciMsg').className='msg err'; $('aciMsg').textContent='下载失败：'+d.error; return; }
      if(d.total){
        $('aciBar').style.width='100%';
        $('aciMsg').className='msg ok';
        $('aciMsg').textContent='ACI 下载完成：'+d.done+' 张 → '+d.out
          + (d.failed? '（'+d.failed+' 张失败）':'')
          + '　现在可以出带底图的矿物分布图了。';
        if(d.log&&d.log.length){ $('aciLog').hidden=false; $('aciLog').textContent=d.log.join('\\n'); }
      }
    }
  });
}

/* ---- WATSON 上色（给 ACI 底图上真彩，可选）---- */
var WAT={timer:null};
$('watStart').addEventListener('click',function(){
  var dir=$('aciDest').value.trim()||aciFolder()+'\\\\_aci';
  $('watStart').disabled=true; $('watCancel').hidden=false;
  $('watProg').hidden=false; $('watBar').style.width='0%';
  $('watMsg').className='msg'; $('watMsg').textContent='正在逐张配对+配准（每张要查清单、下相机模型，慢是正常的）…';
  fetch('/api/watson/start?dir='+encodeURIComponent(dir))
   .then(function(r){return r.json();}).then(function(r){
     if(r.error){ $('watStart').disabled=false; $('watCancel').hidden=true;
       $('watMsg').className='msg err'; $('watMsg').textContent=r.error; return; }
     if(!WAT.timer) WAT.timer=setInterval(pollWat,1000);
   });
});
$('watCancel').addEventListener('click',function(){ fetch('/api/watson/cancel'); });
function pollWat(){
  fetch('/api/status').then(function(r){return r.json();}).then(function(s){
    var d=s.watson||{};
    if(d.running){
      var pct=d.total? Math.round(d.done/d.total*100):0;
      $('watBar').style.width=pct+'%';
      $('watMsg').className='msg';
      $('watMsg').textContent='上色中 '+d.done+' / '+d.total
        + (d.failed? '（跳过 '+d.failed+'）':'');
      if(d.log&&d.log.length){ $('watLog').hidden=false; $('watLog').textContent=d.log.join('\\n'); }
    } else {
      clearInterval(WAT.timer); WAT.timer=null;
      $('watStart').disabled=false; $('watCancel').hidden=true;
      if(d.error){ $('watMsg').className='msg err'; $('watMsg').textContent='上色失败：'+d.error; return; }
      if(d.total){
        $('watBar').style.width='100%';
        $('watMsg').className='msg ok';
        $('watMsg').textContent='上色完成：'+d.done+' 张 → '+d.out+'\\\\colorized'
          + (d.failed? '（'+d.failed+' 张无上色，保持灰度）':'')
          + '　重新出矿物分布图即自动使用上色底图。';
        if(d.log&&d.log.length){ $('watLog').hidden=false; $('watLog').textContent=d.log.join('\\n'); }
      }
    }
  });
}

/* ---- 心跳：窗口关掉后服务端自动退出 ---- */
function ping(){ fetch('/api/ping').catch(function(){}); }
ping(); setInterval(ping, 10000);
</script></body></html>"""


def build_page(listing):
    import json as _json
    from sherloc import minerals as _min
    return (PAGE.replace('<details open style="margin-top:12px">',paper_notes.HTML+'<details open style="margin-top:12px">').replace('</script>',reference.JS+paper_notes.JS+spectrum_export.JS+history.JS.replace('__POINT_BASE__','0')+download_queue.JS+'</script>').replace('id="exportBtn">导出 PNG</button>','id="exportBtn">导出 PNG</button>'+spectrum_export.HTML).replace("__NSCAN__", str(len(listing)))
                .replace("__CLASS_COLORS__",
                         _json.dumps(_min.CLASS_COLORS, ensure_ascii=False))
                .replace("__CLASS_SHORT__",
                         _json.dumps(_min.SHORT, ensure_ascii=False)))


# --------------------------------------------------------------------------
# 服务
# --------------------------------------------------------------------------

def make_handler(state):
    page = build_page(state["listing"]).encode("utf-8")
    analysis_lock = threading.Lock()
    def paper_listing():
        listing=[dict(r) for r in state['listing']]
        try:data=paper_notes.load(state['folder'])
        except ValueError:return listing
        for scan,row in zip(state['scans'],listing):
            if paper_notes.marked(state['folder'],'sherloc',scan,data):row['label']='● '+row['label'].removeprefix('● ')
        return listing

    def analysis_worker(job, scan, out, bands, min_snr, k, scope):
        from sherloc import explore
        def update(done, total, message):
            job.update(done=done, total=total, message=message)
        try:
            job["result"] = explore.run(scan, out, bands, min_snr, k, scope,
                                        update, lambda: job.get("cancel", False))
            from pathlib import Path
            from product_theme import report_style
            report = Path(job["result"]["index"])
            report.write_text(report_style(report.read_text(encoding="utf-8")), encoding="utf-8")
            state.setdefault("analysis_results", {})[job["id"]] = job["result"]["out"]
        except explore.CancelledError as exc:
            job.update(error=str(exc) + "；输出目录：" + getattr(exc, "output_path", out), cancelled=True)
        except Exception as exc:
            job["error"] = "%s: %s" % (type(exc).__name__, exc)
        finally:
            job["running"] = False

    class Handler(BaseHTTPRequestHandler):
        def log_message(self, *a):
            pass

        def _send(self, code, ctype, body):
            self.send_response(code)
            self.send_header("Content-Type", ctype)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            try:
                self.wfile.write(body)
            except Exception:                           # noqa: BLE001
                pass

        def _json(self, obj, code=200):
            self._send(code, "application/json",
                       json.dumps(obj, ensure_ascii=False).encode("utf-8"))

        def _params(self):
            q = urllib.parse.urlparse(self.path).query
            return {k: v[0] for k, v in urllib.parse.parse_qs(q).items()}

        def _steps(self, p):
            on = set(x for x in p.get("steps", "").split(",") if x)
            return {k: (k in on) for k in preprocess.STEP_KEYS + preprocess.EXTRA_KEYS}

        def do_POST(self):
            path = urllib.parse.urlparse(self.path).path
            if path == "/api/analysis/cancel":
                job = state.get("analysis", {})
                if job.get("running"):
                    job["cancel"] = True
                return self._json({"cancelling": bool(job.get("running"))})
            try:
                n = int(self.headers.get("Content-Length") or 0)
                body = self.rfile.read(n).decode("utf-8") if n else "{}"
                data = json.loads(body or "{}")
            except Exception:                           # noqa: BLE001
                if path == "/api/analysis/start":
                    return self._json({"error": "分析请求必须是有效JSON"}, 400)
                data = {}
            if path=='/api/paper':
                try:
                    p=self._params();si=int(p.get('scan',0));pi=int(p.get('point',0))
                    if not 0<=si<len(state['scans']):raise ValueError('扫描索引无效')
                    scan=state['scans'][si]
                    if not 0<=pi<(scan.n_points or 0):raise ValueError('点索引无效')
                    return self._json(paper_notes.save(state['folder'],'sherloc',scan,pi,data.get('scope'),data.get('status'),data.get('note','')))
                except (ValueError,TypeError,AttributeError,OSError) as exc:return self._json({'error':str(exc)},400)
            if path == "/api/analysis/start":
                from sherloc import explore
                import uuid
                try:
                    index = data.get("scan", 0)
                    if isinstance(index, bool) or not isinstance(index, int) or not 0 <= index < len(state["scans"]):
                        raise ValueError("请先选择有效扫描")
                    bands = tuple(tuple(float(v) for v in band) for band in data.get("bands", explore.DEFAULT_BANDS))
                    min_snr, k, scope = float(data.get("min_snr", 3)), data.get("k", 3), data.get("scope", "unknown")
                    explore.validate_parameters(bands, min_snr, k, scope)
                    out = os.path.abspath(data.get("out") or os.path.join(state["folder"], "_analysis"))
                except (ValueError, TypeError, AttributeError) as exc:
                    return self._json({"error": str(exc)}, 400)
                with analysis_lock:
                    if state.get("analysis", {}).get("running"):
                        return self._json({"error": "已有分析任务运行，请等待或取消"}, 409)
                    job = {"id": uuid.uuid4().hex, "running": True, "done": 0, "total": 6,
                           "cancel": False, "message": "准备分析当前扫描"}
                    state["analysis"] = job
                    threading.Thread(target=analysis_worker,
                        args=(job, state["scans"][index], out, bands, min_snr, k, scope), daemon=True).start()
                return self._json({"started": True, "id": job["id"]})
            if path == "/api/settings":
                history.update(SETTINGS,lambda cfg:dict(cfg,**{k:v for k,v in data.items() if k in ('folder','mode')}))
                return self._json({"ok": True})
            return self._json({"error": "not found"}, 404)

        def do_GET(self):
            path = urllib.parse.urlparse(self.path).path
            p = self._params()
            try:
                if path.startswith("/analysis-result/"):
                    from pathlib import Path
                    import mimetypes
                    parts = path.split("/", 3)
                    job = state.get("analysis", {})
                    roots = state.get("analysis_results", {})
                    if len(parts) != 4 or parts[2] not in roots:
                        return self._json({"error": "分析报告不存在或任务尚未完成"}, 404)
                    root = Path(roots[parts[2]]).resolve()
                    file = (root / urllib.parse.unquote(parts[3])).resolve()
                    if not file.is_relative_to(root) or file.suffix.lower() not in (".html", ".png", ".csv", ".json"):
                        return self._json({"error": "无效报告路径"}, 403)
                    if not file.is_file():
                        return self._json({"error": "报告文件不存在"}, 404)
                    mime = mimetypes.guess_type(str(file))[0] or "application/octet-stream"
                    if file.suffix.lower() in (".html", ".csv", ".json"):
                        mime += "; charset=utf-8"
                    return self._send(200, mime, file.read_bytes())
                if path in ("/", "/index.html"):
                    return self._send(200, "text/html; charset=utf-8", page)
                if path == "/favicon.ico":
                    return self._send(204, "image/x-icon", b"")
                if path == "/api/ping":
                    state["last_ping"] = time.time()
                    return self._json({"ok": True})
                if path == "/api/scans":
                    return self._json({"folder": state["folder"], "scans": paper_listing(),
                                       "settings": load_settings()})
                if path=='/api/recent':return self._json({'folders':history.rows(load_settings())})
                if path == "/api/steps":
                    # 原始帧（ERA/ECA）在列表里时，despike 复选框默认勾上
                    # （对齐官方 c 步——processed 路线维持默认关）
                    raw_present = any(s.parts and s.parts[0].product_type in ("ERA", "ECA")
                                      for s in state["scans"])
                    steps_json = [{"key": k, "label": lab,
                                   "on": k in ("deadpixel", "dark", "laser", "baseline")
                                   or (k == "despike" and raw_present)}
                                  for k, lab in preprocess.STEPS]
                    steps_json += [{"key": k, "label": lab, "on": False, "extra": True}
                                   for k, lab in preprocess.EXTRA_FLAGS]
                    return self._json(steps_json)
                if path == "/api/browse":
                    return self._json(browse(p.get("path") or state["folder"]
                                             or os.path.expanduser("~")))
                if path == "/api/load":
                    target = p.get("path", "")
                    if not os.path.isdir(target):
                        return self._json({"error": "目录不存在：%s" % target})
                    products, scans, listing = scan_listing(target)
                    target=os.path.abspath(target);history.store_folder(SETTINGS,target)
                    state["folder"], state["scans"], state["listing"] = target, scans, listing
                    # 缓存键只有 sol|序列号|sclk，不含路径——换了文件夹不清缓存，
                    # 同编号扫描会显示上一个文件夹算出的矿物判定
                    state["mineral_cache"].clear()
                    return self._json({"folder": target, "scans": paper_listing()})
                if path == "/api/points":
                    sc = state["scans"][int(p.get("scan", 0))]
                    return self._json({"points": sc.n_points or 0})
                if path=='/api/paper':
                    si=int(p.get('scan',0));pi=int(p.get('point',0))
                    if not 0<=si<len(state['scans']):return self._json({'error':'扫描索引无效'},400)
                    scan=state['scans'][si]
                    if not 0<=pi<(scan.n_points or 0):return self._json({'error':'点索引无效'},400)
                    return self._json(paper_notes.read(state['folder'],'sherloc',scan,pi))
                if path == "/api/plot":
                    sc = state["scans"][int(p.get("scan", 0))]
                    point = int(p.get("point", 0))
                    d = series_for(sc, self._steps(p), point, p.get("mode", "all"))
                    d["minerals"] = mineral_summary(scan_minerals(state, sc), point)
                    remember_instrument("sherloc")
                    return self._json(d)
                if path=='/api/csv':
                    si=int(p.get('scan',0))
                    if not 0<=si<len(state['scans']):return self._json({'error':'扫描索引无效'},400)
                    scan=state['scans'][si]
                    return self._json(export_current_csv(scan,self._steps(p),int(p.get('point',0)),p.get('mode','all'),p.get('out') or os.path.join(state['folder'],'_png')))
                if path == "/api/png":
                    sc = state["scans"][int(p.get("scan", 0))]
                    point = int(p.get("point", 0))
                    buf = io.BytesIO()
                    tmp = os.path.join(state["folder"], "_png")
                    os.makedirs(tmp, exist_ok=True)
                    f = os.path.join(tmp, "_single.png")
                    render_scan_png(sc, self._steps(p), point, p.get("mode", "all"), f)
                    with open(f, "rb") as fh:
                        buf.write(fh.read())
                    return self._send(200, "image/png", buf.getvalue())
                if path == "/api/export":
                    if state["export"].get("running"):
                        return self._json({"error": "已经在导出了"})
                    out = p.get("out", "").strip()
                    if not out:
                        return self._json({"error": "没有指定输出目录"})
                    steps = self._steps(p) if p.get("steps") else \
                        {k: k in ("deadpixel", "dark", "laser", "baseline")
                         for k in preprocess.STEP_KEYS + preprocess.EXTRA_KEYS}
                    mode = p.get("mode", "all") or "all"
                    threading.Thread(target=_export_worker,
                                     args=(state, out, steps, mode), daemon=True).start()
                    return self._json({"started": True})
                if path == "/api/export/cancel":
                    if state["export"].get("running"):
                        state["export"]["cancel"] = True
                        return self._json({"ok": True, "cancelling": True})
                    return self._json({"ok": True, "cancelling": False})
                if path == "/api/status":
                    return self._json({"export": state["export"], "dl": download_queue.snapshot(state["dl"]),
                                       "aci": state["aci"], "watson": state["watson"],
                                       "analysis": state.get("analysis", {}),
                                       "build": BUILD})
                if path == "/api/aci/scan":
                    folder = p.get("folder") or state["folder"]
                    dest = p.get("dest") or os.path.join(folder, "_aci")
                    try:
                        need = aci.needed(folder)
                    except Exception as exc:            # noqa: BLE001
                        return self._json({"error": "扫描失败：%s" % exc})
                    have = aci.have(dest)
                    miss = [d for d in need if d["stem"] not in have]
                    return self._json({
                        "folder": folder, "dest": dest,
                        "total": len(need), "have": len(need) - len(miss),
                        "missing": len(miss),
                        "items": [{"stem": d["stem"], "name": d["display"],
                                   "sol": d["sol"]} for d in miss],
                    })
                if path == "/api/aci/start":
                    if state["aci"].get("running"):
                        return self._json({"error": "已经有一个 ACI 下载在进行"})
                    folder = p.get("folder") or state["folder"]
                    dest = p.get("dest") or os.path.join(folder, "_aci")
                    need = aci.needed(folder)
                    have = aci.have(dest)
                    stems = [(d["stem"], d["display"]) for d in need
                             if d["stem"] not in have]
                    if not stems:
                        return self._json({"error": "这个文件夹需要的 ACI 图都已经在本地了"})
                    state["aci"] = {"running": True, "cancel": False, "done": 0,
                                    "total": len(stems), "failed": 0,
                                    "bytes_done": 0, "error": None,
                                    "out": dest, "log": []}
                    threading.Thread(target=_aci_worker,
                                     args=(state, folder, dest, stems),
                                     daemon=True).start()
                    return self._json({"started": True, "total": len(stems)})
                if path == "/api/aci/cancel":
                    state["aci"]["cancel"] = True
                    return self._json({"ok": True})
                if path == "/api/watson/start":
                    if state["watson"].get("running"):
                        return self._json({"error": "已经有一个上色任务在进行"})
                    try:
                        import cv2                                          # noqa: F401
                        from sherloc import colorize                        # noqa: F401
                    except ImportError:
                        return self._json({"error":
                            "上色需要 opencv。在命令行执行一次："
                            "pip install opencv-python-headless，然后重启工作台"})
                    aci_dir = p.get("dir") or os.path.join(state["folder"], "_aci")
                    if not os.path.isdir(aci_dir):
                        return self._json({"error": "ACI 目录不存在：%s" % aci_dir})
                    if not os.listdir(aci_dir):
                        return self._json({"error": "ACI 目录是空的，先下载 ACI 图"})
                    state["watson"] = {"running": True, "cancel": False, "done": 0,
                                       "total": 0, "failed": 0, "error": None,
                                       "out": aci_dir, "log": []}
                    threading.Thread(target=_watson_worker,
                                     args=(state, aci_dir), daemon=True).start()
                    return self._json({"started": True})
                if path == "/api/watson/cancel":
                    state["watson"]["cancel"] = True
                    return self._json({"ok": True})
                if path == "/api/dl/sols":
                    kind = p.get("kind", "processed")
                    if kind not in PDS_KINDS:
                        return self._json({"error": "未知的数据来源：%s" % kind})
                    if p.get("refresh") or kind not in state["sol_cache"]:
                        try:
                            state["sol_cache"][kind] = list_sol_dirs(kind)
                        except Exception as exc:          # noqa: BLE001
                            return self._json({"error": "连不上 PDS：%s" % exc})
                    sols = state["sol_cache"][kind]
                    return self._json({"kind": kind, "sols": sols, "cached": not p.get("refresh")})
                if path == "/api/dl/list":
                    try:
                        return self._json(list_sol_files(p.get("kind", "processed"),
                                                         p.get("sol", "")))
                    except Exception as exc:              # noqa: BLE001
                        return self._json({"error": "读取目录失败：%s" % exc})
                if path == "/api/dl/start":
                    if state["dl"].get("running"):
                        return self._json({"error": "已经有一个下载在进行"})
                    kind = p.get("kind", "processed")
                    try:sol=download_queue.normalize(p.get('sols') or p.get('sol',''))
                    except ValueError as exc:return self._json({'error':str(exc)},400)
                    types = set(x for x in p.get("types", "").split(",") if x)
                    dest = p.get("dest", "").strip()
                    if kind not in PDS_KINDS or not types or not dest:
                        return self._json({"error": "参数不完整"})
                    state["dl"] = {"running": True, "cancel": False, "done": 0, "total": 0,
                                   "bytes_done": 0, "bytes_total": 0, "error": None,
                                   "out": dest, "sol": sol, "kind": kind, "log": []}
                    threading.Thread(target=_dl_worker,
                                     args=(state, kind, sol, types, dest), daemon=True).start()
                    return self._json({"started": True})
                if path == "/api/dl/cancel":
                    state["dl"]["cancel"] = True
                    return self._json({"ok": True})
                if path == "/api/openfolder":
                    target = p.get("path", "")
                    if os.path.isdir(target):
                        try:
                            os.startfile(target)            # noqa: S606
                        except Exception:                   # noqa: BLE001
                            pass
                    return self._json({"ok": True})
                if path == "/api/qc":
                    folder = state["folder"]
                    out = os.path.join(folder, "_qc", "index.html")
                    rc, log = _run_qc(folder)
                    if rc != 0:
                        return self._json({"error": "生成失败，请看控制台输出"})
                    try:
                        open_app_window("file:///" + out.replace("\\", "/"), True)
                    except Exception:                       # noqa: BLE001
                        pass
                    return self._json({"index": out, "log": log[-500:]})
                if path == "/api/maps":
                    folder = state["folder"]
                    out = os.path.join(folder, "_maps")
                    # 变量名不能叫 aci：Python 的函数作用域是平的，
                    # 这里一赋值就会把整个 do_GET 里的 aci 模块遮蔽成局部变量，
                    # 于是 /api/aci/scan 那几行引用 aci.needed() 直接报
                    # "cannot access local variable 'aci'"。踩过。
                    aci_dir = (p.get("aci") or "").strip()
                    try:
                        res = _run_maps(folder, aci_dir or None)
                    except Exception as exc:                # noqa: BLE001
                        return self._json({"error": "%s: %s" % (type(exc).__name__, exc)})
                    if res.get("error"):
                        return self._json(res)
                    try:
                        os.startfile(out)                   # noqa: S606
                    except Exception:                       # noqa: BLE001
                        pass
                    res["out"] = out
                    return self._json(res)
                if path == "/api/quit":
                    self._json({"ok": True})
                    threading.Thread(target=state["shutdown"], daemon=True).start()
                    return
                return self._send(404, "text/plain; charset=utf-8", b"not found")
            except Exception as exc:                        # noqa: BLE001
                if path.startswith("/api/"):
                    return self._json({"error": "%s: %s" % (type(exc).__name__, exc)}, 500)
                return self._send(500, "text/plain; charset=utf-8", str(exc).encode("utf-8"))

    return runtime_http.decorate(Handler,state)


def _export_name(sc, i):
    """导出 PNG 的文件名。必须带 SCLK：同 sol 同序列号的多次独立扫描
    （sol 141 就有 3 次）只用 sol+序列号会互相覆盖。"""
    sclk_tail = "_".join(s.replace(".", "_") for s in sc.all_sclks) or "na"
    return "sol%s_%s_%s_pt%03d.png" % (sc.sol, sc.sequence_id, sclk_tail, i)


def _export_worker(state, out, steps, mode):
    ex = state["export"]
    ex.update(running=True, done=0, total=0, out=out, error=None, cancel=False)
    try:
        os.makedirs(out, exist_ok=True)
        total = sum((sc.n_points or 0) for sc in state["scans"])
        ex["total"] = total
        for sc in state["scans"]:
            for i in range(sc.n_points or 0):
                if ex.get("cancel"):
                    ex["error"] = "已取消（已完成 %d 张）" % ex["done"]
                    return
                render_scan_png(sc, steps, i, mode,
                                os.path.join(out, _export_name(sc, i)))
                ex["done"] += 1
    except Exception as exc:                            # noqa: BLE001
        ex["error"] = "%s: %s" % (type(exc).__name__, exc)
    finally:
        ex["running"] = False


def _watson_worker(state, aci_dir):
    """后台线程：给 ACI 目录里的每张灰度图找共视的 WATSON 彩色图并上色。

    每张图独立成败：找不到配对 / 判 NO-GO / 下载失败都只记一行日志，
    不影响其他图。上色结果缓存到 <aci_dir>/colorized/，矿物分布图出图时
    会自动优先使用。
    """
    st = state["watson"]
    try:
        from sherloc import colorize as _colorize
        from sherloc import watson as _watson
        if not _colorize.HAS_CV2:
            st["error"] = "未安装 opencv：pip install opencv-python-headless"
            return
        from make_mineral_map import colorized_base, load_aci
        from PIL import Image

        imgs = sorted(f for f in os.listdir(aci_dir)
                      if f.lower().endswith((".img", ".png", ".tif", ".jpg")))
        st.update(total=len(imgs), done=0, failed=0, out=aci_dir, log=[])
        for f in imgs:
            if st.get("cancel"):
                st["log"].append("已取消")
                break
            path = os.path.join(aci_dir, f)
            col_img, why = colorized_base(path, verbose=False)
            if col_img is not None:
                col_img.save(os.path.join(aci_dir, "colorized",
                                          os.path.splitext(f)[0] + ".png"))
                st["done"] += 1
                st["log"].append("%s  ✓" % f)
            else:
                st["failed"] += 1
                st["log"].append("%s  未上色：%s" % (f, why))
    except Exception as exc:                            # noqa: BLE001
        st["error"] = "%s: %s" % (type(exc).__name__, exc)
    finally:
        st["running"] = False


def _run_qc(folder):
    if getattr(sys,'frozen',False):
        import run_qc
        import contextlib
        buffer=io.StringIO()
        try:
            with contextlib.redirect_stdout(buffer):run_qc.analyse_folder(folder,os.path.join(folder,'_qc'))
            return 0,buffer.getvalue()
        except Exception:
            appenv.logger().exception('frozen SHERLOC report failed')
            return 1,buffer.getvalue()
    import subprocess
    script = os.path.join(HERE, "run_qc.py")
    try:
        r = subprocess.run([sys.executable, script, folder], capture_output=True,
                           text=True, encoding="utf-8", errors="replace", timeout=1800)
        return r.returncode, (r.stdout or "") + (r.stderr or "")
    except Exception as exc:                            # noqa: BLE001
        return 1, str(exc)


def _run_maps(folder, aci_dir=None):
    """生成矿物分布图。直接调函数，不走子进程——这样进度和结果都能拿到结构化数据。"""
    import io
    from contextlib import redirect_stdout

    import make_mineral_map

    buf = io.StringIO()
    try:
        with redirect_stdout(buf):
            res = make_mineral_map.run(folder, aci_dir=aci_dir)
    except Exception as exc:                            # noqa: BLE001
        return {"error": "%s: %s" % (type(exc).__name__, exc),
                "log": buf.getvalue()[-2000:]}
    if not res.get("count"):
        return {"error": "这个文件夹里没有能出图的观测（需要有 RRS 和 RLS）。"
                         "只有原始件（ERA/ERB）时，请选择带「ERA · 原始帧」标记的扫描——"
                         "预处理面板从暗背景扣除开始逐步处理",
                "log": buf.getvalue()[-2000:]}
    res["log"] = buf.getvalue()[-2000:]
    return res


def build_tag():
    """整个程序的版本标记：所有 .py 的（路径, 改动时间, 大小）求和。

    不能只看 view.py——改了 sherloc/ 或 make_mineral_map.py 却检测不到，
    双击启动就会复用旧实例，新代码永远不生效（这个坑踩过一次）。
    """
    import hashlib

    h = hashlib.sha1()
    files = []
    for root in (HERE, os.path.join(HERE, "sherloc")):
        try:
            for f in sorted(os.listdir(root)):
                if f.endswith(".py"):
                    files.append(os.path.join(root, f))
        except OSError:
            pass
    for f in files:
        try:
            st = os.stat(f)
        except OSError:
            continue
        h.update(("%s|%d|%d\n" % (os.path.basename(f), int(st.st_mtime),
                                  st.st_size)).encode("utf-8"))
    return h.hexdigest()[:12] if files else "?"


BUILD = appenv.build_id()


def running_instance(port, timeout=2.0):
    """探测端口上是不是已经跑着本程序。

    为什么不用"试着绑一下端口"：Windows 上 SO_REUSEADDR 允许同一个端口被绑多次，
    所以那种探测**永远成功**，单实例判断等于没有（这就是三个实例同时占着 8765 的原因）。
    改成发一个 HTTP 请求问它，最实在。

    返回 (build 标记 或 "unknown")，没有实例则返回 None。
    故意用 /api/status 而不是 /api/ping —— 后者会刷新心跳，
    等于让一个没人看的旧实例永远不死。
    """
    url = "http://127.0.0.1:%d/api/status" % port
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "sherloc-workbench/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            d = json.loads(r.read().decode("utf-8", "replace"))
        return d.get("build") or "unknown"
    except Exception:                                   # noqa: BLE001
        return None


def _watchdog(state):
    """页面不再发心跳（窗口被关掉）就自动退出，避免留下占着端口的孤儿进程。"""
    while True:
        time.sleep(5)
        last = state.get("last_ping")
        if last and time.time() - last > PING_TIMEOUT:
            print("\n界面已经关闭，自动退出。")
            state["shutdown"]()
            return


# --------------------------------------------------------------------------
# 入口
# --------------------------------------------------------------------------

def pick_start_folder(arg):
    folder=arg or load_settings().get('folder')
    return os.path.abspath(folder) if folder and os.path.isdir(folder) else appenv.empty_data('sherloc')


def start(port,folder=None):
    runtime_http.initialize()
    folder=pick_start_folder(folder)
    scans,listing=[],[]
    if os.path.isdir(folder):_,scans,listing=scan_listing(folder)
    state={'folder':folder,'scans':scans,'listing':listing,'export':{},'dl':{},'aci':{},'watson':{},'sol_cache':{},
           'mineral_cache':{},'last_ping':None,'shutdown':lambda:None}
    server=runtime_http.bind(make_handler(state),port)
    state['shutdown']=server.shutdown
    return server,state


def main(argv=None):
    ap=argparse.ArgumentParser(description='SHERLOC光谱工作台')
    ap.add_argument('folder',nargs='?')
    ap.add_argument('--version',action='version',version=appenv.build_id())
    ap.add_argument('--port',type=int,default=8765)
    ap.add_argument('--browser',action='store_true')
    ap.add_argument('--no-window',action='store_true')
    ap.add_argument('--no-watchdog',action='store_true')
    args=ap.parse_args(argv)
    if not 0<=args.port<=65535:ap.error('端口须在0到65535之间')
    if args.port and running_instance(args.port)==BUILD:
        if not args.no_window:open_app_window('http://127.0.0.1:%d/'%args.port,not args.browser)
        return 0
    server,state=start(args.port,args.folder)
    url='http://127.0.0.1:%d/'%server.server_port
    print('MarsWorkbench '+appenv.build_id()+' · SHERLOC：'+url,flush=True)
    if not args.no_watchdog:threading.Thread(target=_watchdog,args=(state,),daemon=True).start()
    if not args.no_window:open_app_window(url,not args.browser)
    try:server.serve_forever()
    except KeyboardInterrupt:pass
    finally:server.server_close()
    return 0

if __name__=='__main__':
    try:raise SystemExit(main())
    except Exception as exc:
        print(appenv.crash(exc));raise SystemExit(1)
