#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""火星仪器操作台 · 通用光谱工作台。

和 SHERLOC 光谱工作台（view.py）同一套技术：本地 HTTP 服务 + 单页界面 +
Edge 应用窗口，零第三方服务。区别在于这边的页面与 API 是「仪器无关」的，
具体读什么数据、怎么预处理、怎么判定元素线，都由 instruments/ 里的
后端（ChemCam / SuperCam / APXS）决定：

    python console.py [--app chemcam] [--port 8766] [--browser] [--no-window]

一个服务同时带三个仪器，页面用 ?app=<id> 区分；各仪器各自记住数据文件夹。
"""

from __future__ import annotations

import source_paths  # source roots; no-op for bundled imports

import argparse
import hashlib
import json
import os
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

from instruments import REGISTRY                              # noqa: E402
from instruments import evidence as ievidence
from instruments import cross as icross
from instruments import notes as paper_notes
from instruments import spectrum_export
from instruments import absorption
from instruments import rmi, rmi_overlay, supercam_images
from instruments import filtering
from instruments import history
from instruments import download_queue
from instruments import base as ibase                         # noqa: E402
from instruments import pds_listing                           # noqa: E402
from sherloc.desktop import open_app_window, remember_instrument  # noqa: E402
from sherloc.render import render_png                        # noqa: E402
from instruments import absorption_ui as absorption_ui
from instruments import rmi_ui as rmi_ui
from instruments import filtering_ui as filtering_ui
from instruments import supercam_images_ui as supercam_images_ui
from instruments import quant_ui as quant_ui

SETTINGS = appenv.settings('console')
PING_TIMEOUT = 180.0
BACK_URL = "http://127.0.0.1:8790/"          # 启动器页面（没开也不影响使用）


def alert(title, message):
    try:
        import ctypes
        ctypes.windll.user32.MessageBoxW(0, message, title, 0x10)
    except Exception:                                       # noqa: BLE001
        print("%s\n%s" % (title, message))


def load_settings():
    try:
        return history.read(SETTINGS)
    except Exception:                                       # noqa: BLE001
        return {}


def save_settings(data):
    try:
        history.update(SETTINGS,lambda current:dict(current,**data))
    except (OSError,ValueError):
        appenv.logger().exception('could not save settings')


# --------------------------------------------------------------------------
# 扫描整理 / 取谱（每个仪器一份状态）
# --------------------------------------------------------------------------

def obs_listing(observations,folder=None,instrument=None):
    out = []
    try:note_data=paper_notes.load(folder) if folder else None
    except ValueError:note_data=None
    for i, obs in enumerate(observations):
        out.append({
            "index": i, "sol": obs.sol, "seq": obs.seq,
            "points": obs.n_points, "product": obs.product,
            "target": (getattr(obs, "target", "") or ""),
            "site": getattr(obs, "site", None),
            "drive": getattr(obs, "drive", None),
            "label": "Sol %s · %s · %d 点 · %s" % (obs.sol, obs.seq,
                                                   obs.n_points, obs.product),
        })
        if note_data and paper_notes.marked(folder,instrument,obs,note_data):out[-1]['label']='● '+out[-1]['label']
    return out


def _count_data(files):
    exts = {os.path.splitext(f)[1].lower() for f in files}
    return exts


def browse(path):
    """列子文件夹；标出哪些直接装着本程序认识的仪器数据。

    每个子目录都要真跑一遍扫描，大目录会卡——最多探 80 个子目录。
    """
    path = os.path.abspath(path or os.path.expanduser("~"))
    try:
        names = sorted(os.listdir(path), key=lambda s: s.lower())
    except Exception as exc:                                # noqa: BLE001
        return {"path": path, "parent": None, "dirs": [], "error": str(exc)}
    items = []
    for n in names[:400]:
        full = os.path.join(path, n)
        if not os.path.isdir(full) or n.startswith((".", "_", "$")):
            continue
        tag = None
        if len(items) < 80:
            for bid, b in REGISTRY.items():
                try:
                    found = b.scan_folder(full)
                except Exception:                           # noqa: BLE001
                    found = []
                if found:
                    tag = "%s %d 次观测" % (b.name, len(found))
                    break
        items.append({"name": n, "path": full, "tag": tag})
    items.sort(key=lambda v: (not v["tag"], v["name"].lower()))
    return {"path": path, "parent": os.path.dirname(path.rstrip("\\/")) or None,
            "dirs": items, "error": None}


def _steps_of(backend, p):
    on = set(x for x in p.get("steps", "").split(",") if x)
    return {k: (k in on) for k, _l, _d in backend.steps}


def _shot_list(n_shots, shot_numbers=None):
    return [{"i": -1, "label": "总谱（平均）"}] + \
           [{"i": i, "label": _shot_label(i, n_shots, shot_numbers)} for i in range(n_shots)]


def plot_for(state, backend, scan_i, point_i, shot, steps_on):
    obs = state["observations"][scan_i]
    point = obs.points[point_i]
    d = backend.load_point(point, shot, steps_on)
    passive = bool(d.get("passive"))
    if passive:
        marks, elements = [], []
    else:
        marks, elements = backend.annotate_data(d)
    out = {
        "x": [round(float(v), 4) for v in d["x"]],
        "xlabel": d.get("xlabel") or backend.xlabel,
        "ylabel": d.get("ylabel") or backend.ylabel,
        "series": [{"key": "main", "label": "全谱" if passive
                    else _spectrum_label(backend, shot, d),
                    "y": [None if not np.isfinite(v) else round(float(v), 4)
                          for v in d["y"]]}],
        "log": [{"step": a, "text": b} for a, b in d["log"]],
        "n_shots": d["n_shots"],
        "shots": ([{"i": -1, "label": "全谱"}] if (passive or backend.id == "apxs")
                  else _shot_list(d["n_shots"], d.get("shot_numbers"))),
        "meta": d["meta"],
        "passive": passive,
        "absorption":d.get('absorption'),
        "annotate": {"marks": marks, "elements": elements,"element_evidence":ievidence.element_statistics(marks,elements)},
        "n_points": obs.n_points,
        "quant": (backend.quant_payload(obs) if backend.id in ("chemcam", "supercam", "apxs") and getattr(backend, "quant_payload", None) else None),
    }
    if backend.id == "apxs":
        ox = backend.oxides(obs, point)
        if ox:
            out["oxides"] = ox
    if backend.id == "supercam" and not passive:
        out["moc"] = backend.moc_payload(point.stem)
    return out


def _shot_label(shot, n_shots, shot_numbers=None):
    if shot is None or shot < 0 or shot >= n_shots:
        return "总谱（%d 发平均）" % n_shots
    return ("原始第 %d 发" % shot_numbers[shot] if shot_numbers is not None
            else "第 %d 发" % (shot + 1))


def _spectrum_label(backend, shot, data):
    if backend.id == "apxs":
        return "全谱"
    return _shot_label(shot, data["n_shots"], data.get("shot_numbers"))


def render_point_png(backend, obs, point_i, shot, steps_on, out_path):
    d = backend.load_point(obs.points[point_i], shot, steps_on)
    x = np.asarray(d["x"], dtype=float)
    y = np.asarray(d["y"], dtype=float)
    passive = bool(d.get("passive"))
    elements = []
    if not passive:
        _marks, elements = backend.annotate_data(d)
    title = "%s · Sol %s · %s · 第 %d 点" % (backend.name, obs.sol, obs.seq,
                                             point_i + 1)
    label = "全谱" if passive else _spectrum_label(backend, shot, d)
    sub = "%s · %s" % (obs.product, label)
    notes = ["%s | %s" % (a, b) for a, b in d["log"]]
    if d.get('absorption') and d['absorption']['enabled']:
        refs=d['absorption']['references']
        notes.append('吸收参考（非矿物鉴定）：'+('；'.join('%s 深度%.1f%%'%(r['label'],100*r['depth']) for r in refs if r['detected']) or '未达浏览门槛'))
    if elements:
        notes.append("检出元素：" + " ".join(elements))
    if backend.id == "supercam" and not passive:
        moc = backend.moc_payload(obs.points[point_i].stem)
        if moc.get("matched"):
            pairs = ["%s %.1f±%s" % (r["oxide"], r["wt"], "%.1f" % r["err"] if r["err"] is not None else "未提供")
                     for r in moc["rows"] if r["wt"] is not None]
            if pairs:
                notes.append("官方成分预测 wt%：" + " ".join(pairs[:8]))
    render_png(out_path, [(label, y, _hex_color(backend.accent))],
               x, d.get("xlabel") or backend.xlabel, title, sub,
               notes=notes, ylabel=d.get("ylabel") or backend.ylabel)
    if not passive and os.path.isfile(out_path):ievidence.append_png(out_path,_marks,elements)
    return out_path


def _hex_color(s):
    s = s.lstrip("#")
    return tuple(int(s[i:i + 2], 16) for i in (0, 2, 4))


# --------------------------------------------------------------------------
# 页面（仪器无关模板 + 启动时按后端注入常量）
# --------------------------------------------------------------------------

PAGE = """<!doctype html><html lang="zh"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>__TITLE__</title><style>
:root{--bg:#f1f3f6;--panel:#fff;--panel2:#f7f8fa;--ink:#14171b;--muted:#5f6773;
--line:#e0e4ea;--line2:#eef1f5;--accent:__ACCENT__;--accent-ink:#fff;
--ok:#12794a;--warn:#a86400;--err:#b3352f;--radius:10px}
@media (prefers-color-scheme:dark){:root{--bg:#0f1113;--panel:#1b1e22;--panel2:#22262b;
--ink:#e9ecf0;--muted:#98a1ac;--line:#31363d;--line2:#262b31;--accent:__ACCENT_D__;
--accent-ink:#0f1113;--ok:#5cc08b;--warn:#e0a63c;--err:#e0796f}}
html[data-theme="dark"]{--bg:#0f1113;--panel:#1b1e22;--panel2:#22262b;
--ink:#e9ecf0;--muted:#98a1ac;--line:#31363d;--line2:#262b31;--accent:__ACCENT_D__;
--accent-ink:#0f1113;--ok:#5cc08b;--warn:#e0a63c;--err:#e0796f}
html[data-theme="light"]{--bg:#f1f3f6;--panel:#fff;--panel2:#f7f8fa;
--ink:#14171b;--muted:#5f6773;--line:#e0e4ea;--line2:#eef1f5;--accent:__ACCENT__;
--accent-ink:#fff;--ok:#12794a;--warn:#a86400;--err:#b3352f}
svg text.tick{fill:var(--muted);font-size:11px}
svg .gridline{stroke:var(--line)}
svg .axl{fill:var(--ink);font-size:12px}
svg .peakline{opacity:.55;stroke-dasharray:4 4}
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
*{box-sizing:border-box}
[hidden]{display:none!important}
html,body{margin:0;height:100%}
body{background:var(--bg);color:var(--ink);
font:14px/1.6 ui-sans-serif,-apple-system,"Segoe UI","Microsoft YaHei",sans-serif}
.top{position:sticky;top:0;z-index:20;display:flex;align-items:baseline;gap:14px;
flex-wrap:wrap;background:color-mix(in srgb,var(--bg) 80%,transparent);
backdrop-filter:blur(8px);border-bottom:1px solid var(--line2);padding:10px 22px 8px}
.top h1{font-size:17px;font-weight:600;margin:0;letter-spacing:.2px}
.top h1::before{content:"";display:inline-block;width:9px;height:9px;border-radius:2px;
background:var(--accent);margin-right:9px}
.top .sub{color:var(--muted);font-size:12px}
.top .right{margin-left:auto;display:flex;gap:8px}
main{max-width:1280px;margin:0 auto;padding:0 18px 40px}
.card{background:var(--panel);border:1px solid var(--line);border-radius:var(--radius);
padding:12px 14px;margin-bottom:12px;box-shadow:0 1px 2px rgba(16,24,40,.05)}
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
button.primary{background:var(--accent);color:var(--accent-ink);border-color:transparent;
font-weight:500;box-shadow:0 1px 2px rgba(16,24,40,.18)}
button.ghost{background:transparent}
button:disabled{opacity:.45;cursor:default}
select,input[type=text]{font:inherit;padding:6px 9px;border:1px solid var(--line);
border-radius:7px;background:var(--panel);color:var(--ink);max-width:100%;min-width:0}
input[type=range]{width:100%}
#scan{width:100%}
.row{display:grid;grid-template-columns:280px minmax(0,1fr);gap:12px;align-items:start}
.row > *{min-width:0}
.log{margin:0;padding:0;list-style:none;font-size:12px;max-height:170px;overflow:auto}
.log li{display:flex;gap:8px;padding:4px 0;border-bottom:1px solid var(--line2);color:var(--muted)}
.log li b{font-weight:500;color:var(--ink);flex:0 0 76px}
.log li span{min-width:0;overflow-wrap:anywhere;word-break:break-word}
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
overflow-wrap:anywhere;color:var(--ink);font:inherit;cursor:pointer}
.bname:hover{background:var(--panel2)}
.buse{background:var(--accent);color:var(--accent-ink);border:0;padding:5px 12px;font-size:13px}
.dot{color:var(--accent);font-size:11px;margin-left:6px}
.prog{height:8px;background:var(--line2);border-radius:4px;overflow:hidden;margin-top:8px}
.prog > i{display:block;height:100%;width:0;background:var(--accent);transition:width .2s}
.msg{font-size:12px;margin-top:6px}
.msg.ok{color:var(--ok)}.msg.err{color:var(--err)}
.typelist{display:grid;grid-template-columns:repeat(auto-fill,minmax(280px,1fr));
gap:2px 18px;margin-top:6px}
.chk.ty{display:flex;align-items:center;gap:6px;margin:0;padding:3px 0}
.chk.ty .sz{margin-left:auto;color:var(--muted);font-size:11px;white-space:nowrap}
.dllog{margin:8px 0 0;padding:8px 10px;background:var(--panel2);border-radius:6px;
font-size:11px;line-height:1.5;max-height:150px;overflow:auto;white-space:pre-wrap;
overflow-wrap:anywhere;color:var(--muted)}
details{margin:2px 0}
summary{cursor:pointer;font-size:12px;color:var(--muted);user-select:none;white-space:nowrap}
summary:hover{color:var(--ink)}
.stepper{display:flex;gap:4px}
.stepper button{flex:1;padding:6px 0;background:var(--panel);font-weight:600}
.stepper button:hover{border-color:var(--accent);color:var(--accent)}
.elebox{margin-top:10px;border-top:1px solid var(--line2);padding-top:8px;font-size:12px}
.elebox .row2{display:flex;gap:8px;align-items:baseline;padding:2px 0;flex-wrap:wrap}
.tage{display:inline-block;padding:1px 9px;border-radius:3px;font-size:12px;
background:color-mix(in srgb, var(--accent) 15%, transparent);color:var(--accent);
font-weight:600;margin:0 3px 3px 0}
.elebox .ev{color:var(--muted);font-size:11px}
.oxtable{margin-top:8px;border-top:1px solid var(--line2);padding-top:8px}
.oxrow{display:grid;grid-template-columns:64px 1fr 110px;gap:8px;align-items:center;
padding:2px 0;font-size:12px;font-variant-numeric:tabular-nums}
.oxbar{height:10px;background:var(--line2);border-radius:5px;overflow:hidden}
.oxbar > i{display:block;height:100%;background:var(--accent);opacity:.85}
input[type=checkbox],input[type=range]{accent-color:var(--accent)}
@media (max-width:900px){.row{grid-template-columns:1fr}}
</style></head><body>
<div class="top">
  <button class="ghost" id="backBtn" title="回到仪器选择页">‹ 选择仪器</button>
  <h1>__TITLE__</h1>
  <span class="sub">__ROVER__ · __SUBTITLE__</span>
  <span class="right">
    <button class="iconbtn" id="themeBtn" title="切换深色/浅色">◐ 深浅色</button>
  </span>
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
    </div>
    <div class="bar" style="margin-top:8px">
      <span id="shotLabel" class="lbl" style="margin:0" __SHOT_HIDDEN__>激光发数</span>
      <select id="shot" style="min-width:140px" __SHOT_HIDDEN__></select>
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
      <button id="qcBtn">生成质控报告</button>
      <button id="cancelExportBtn" hidden>取消导出</button>
      <button id="openOutBtn" hidden>打开输出目录</button>
      <button class="ghost" id="quitBtn">退出</button>
    </div>
    <div class="msg" id="qcMsg" hidden></div>
    <div class="prog" id="prog" hidden><i id="progBar"></i></div>
    <div class="msg" id="progMsg"></div>
  </div>

  <div class="card" id="browser" hidden></div>

  <div class="card" id="dlCard" hidden>
    <div class="bar" style="margin-top:2px">
      <span class="lbl" style="margin:0">数据来源</span>
      <span id="dlKinds" style="font-size:13px"></span>
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
  </div>

  <div class="row">
    <div class="card">
      <span class="lbl">观测（共 <span id="nscan">0</span> 次）</span>
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
      <div class="bar" style="margin-bottom:4px;gap:8px">
        <button id="lockBtn" class="ghost" style="padding:3px 12px;font-size:12px"
          title="把当前谱固定为参考，之后的每次翻点都与它叠加对比">锁定对比</button>
        <span class="ev" id="lockInfo"></span>
        <span class="ev" id="zoomInfo" style="margin-left:auto"></span>
      </div>
      <div class="legend" id="legend"></div>
      <div id="plot" style="position:relative"></div>
      <p class="hint" id="hint"></p>
      <div class="elebox" id="elements" hidden></div>
      <div class="oxtable" id="oxides" hidden></div>
      <div class="oxtable" id="mocbox" hidden></div>
      <div class="oxtable" id="crossbox" hidden></div>
    </div>
  </div>
</main>
<script>
var APP='__APP__';
var STEPS=[];
var $=function(id){return document.getElementById(id);};
var scanSel=$('scan'), pt=$('pt'), plot=$('plot'), logEl=$('log'), hint=$('hint');

function api(path){ return '/api/'+path+(path.indexOf('?')>=0?'&':'?')+'app='+APP; }
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
  return 'app='+APP+'&scan='+scanSel.value+'&point='+pt.value+
         '&shot='+$('shot').value+'&steps='+encodeURIComponent(stepParam());
}
function showFolder(p){ $('folderPath').textContent=p||'-';
  if(p){ $('outDir').value=p.replace(/[\\\\/]+$/,'')+'\\\\_png'; } }

fetch(api('steps')).then(function(r){return r.json();}).then(function(list){
  var box=$('steps');
  list.forEach(function(s){
    STEPS.push(s.key);
    var l=document.createElement('label');
    l.className='chk';
    l.innerHTML='<input type="checkbox" data-k="'+s.key+'"'+(s.on?' checked':'')+'>'+s.label;
    l.querySelector('input').addEventListener('change',function(){load();if($('scanFilter').value.trim())searchScans();});
    box.appendChild(l);
  });
  if(!list.length){ box.innerHTML='<span class="hint" style="margin:0">该仪器无需预处理</span>'; }
  return fetch(api('scans')).then(function(r){return r.json();});
}).then(function(d){
  showFolder(d.folder); useScans(d.scans);
});

var SCANS=[];
function useScans(list){
  resetScanFilter();
  PAPER_ID=null;PAPER=null;paperToken++;RMI_ID=null;RMI_TOKEN++;RMI_IMAGE_GENERATION++;
  SCANS=list; $('nscan').textContent=list.length; scanSel.innerHTML='';
  list.forEach(function(s){
    var o=document.createElement('option');
    o.value=s.index; o.textContent=s.label; scanSel.appendChild(o);
  });
  if(list.length){ scanSel.value=list[0].index; onScanChange(); }
  else { noFilteredObservation();plot.innerHTML='<p class="hint">这个文件夹里没有找到可用的观测（'
    +'ChemCam 需要 CL5 ccs/rdr 的 csv；SuperCam 需要 cl1 的 fits；'
    +'APXS 需要 rsp 的 csv）。点上面的「换文件夹」或「下载数据」。</p>';
        logEl.innerHTML=''; hint.textContent='';
        $('elements').hidden=true; $('oxides').hidden=true; $('mocbox').hidden=true; }
}
var loadToken=0;
function onScanChange(){
  VIEW=null;
  var tok=++loadToken;
  $('ptlabel').textContent='加载中…';
  fetch(api('points?scan='+scanSel.value)).then(function(r){return r.json();}).then(function(d){
    if(tok!==loadToken) return;
    pt.max=Math.max(0,d.points-1); pt.value=0; load();
  });
}
scanSel.addEventListener('change',function(){ VIEW=null; onScanChange(); });
pt.addEventListener('input',load);
$('shot').addEventListener('change',load);
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
    if(e.target===pt&&(e.key==='ArrowLeft'||e.key==='ArrowRight')){
      stepPt(e.shiftKey?(e.key==='ArrowLeft'?-10:10)
                    :(e.key==='ArrowLeft'?-1:1));
      e.preventDefault();
    }
    return;
  }
  if(e.key==='ArrowLeft'){stepPt(e.shiftKey?-10:-1);e.preventDefault();}
  if(e.key==='ArrowRight'){stepPt(e.shiftKey?10:1);e.preventDefault();}
  if(e.key==='ArrowUp'||e.key==='ArrowDown'){
    var sel=$('shot');
    if(sel&&!sel.hidden&&sel.options.length>1){
      var idx=sel.selectedIndex+(e.key==='ArrowUp'?-1:1);
      if(idx>=0&&idx<sel.options.length){ sel.selectedIndex=idx; load(); }
    }
    e.preventDefault();
  }
});
/* ---- 跨仪器对照（ChemCam × APXS / SuperCam × PIXL；换观测时刷新） ---- */
function loadCross(){
  var box=$('crossbox');
  if(APP!=='chemcam'&&APP!=='supercam'){ box.hidden=true; return; }
  box.hidden=false;
  box.innerHTML= APP==='supercam'
    ? '<div class="ev">正在查找同址 PIXL 扫描…（首次要连 PDS 取清单与成分表，可能十几秒；之后有缓存）</div>'
    : '<div class="ev">正在查找同址 APXS 观测…（首次要对观测逐点算证据，'
      +'9 点观测约需半分钟，之后有缓存）</div>';
  fetch(api('cross?scan='+scanSel.value)).then(function(r){return r.json();})
   .then(function(d){ renderCross(d); })
   .catch(function(){ box.innerHTML='<div class="ev">对照查询失败</div>'; });
}
function renderCross(d){
  var box=$('crossbox');
  if(!d||d.supported===false){ box.hidden=true; return; }
  var isPixl=APP==='supercam';
  var title=isPixl?'跨仪器对照（PIXL）':'跨仪器对照（APXS）';
  if(d.error){ box.innerHTML='<div class="ev">对照失败：'+d.error+'</div>'; return; }
  if(d.note){ box.innerHTML='<div style="display:flex;gap:10px;align-items:baseline">'
    +'<b>'+title+'</b><span class="ev">'+d.note+'</span></div>'; return; }
  if(!d.payload||!d.matches.length){
    box.innerHTML='<div style="display:flex;gap:10px;align-items:baseline">'
    +'<b>'+title+'</b><span class="ev">'+(isPixl?'PIXL 清单':'APXS 文件夹（'+(d.apxs_count||0)
    +' 次观测）')+'里没有匹配的观测：site/drive 不同、目标名不同、sol 也不相邻</span></div>';
    return;
  }
  var m=d.matches[0], pl=d.payload;
  var sideDesc= isPixl
    ? (' → PIXL Sol '+pl.pixl.sol+' · 定位 '+pl.pixl.pos+' · '+pl.pixl.n_points+' 点扫描中位')
    : (' → APXS Sol '+pl.apxs.sol+' · '+pl.apxs.seq+' · '+(pl.apxs.target||'?'));
  var h='<div style="display:flex;gap:10px;align-items:baseline;flex-wrap:wrap">'
    +'<b>'+title+'</b>'
    +'<span class="ev">'+m.detail+sideDesc+'</span>'
    +(d.matches.length>1?'<span class="ev">（另有 '+(d.matches.length-1)+' 个候选）</span>':'')
    +'</div>';
  var wmax=0;
  pl.rows.forEach(function(r){ if(r.apxs&&r.apxs.wt>wmax) wmax=r.apxs.wt; });
  pl.rows.forEach(function(r){
    var lib=r.libs?('<span style="color:var(--ok)">检出 '+r.libs.points
      +' 点 · 平均 '+Math.round(r.libs.mean_fraction*100)+'%发</span>')
      :'<span class="ev">—</span>';
    var ap;
    if(r.apxs){
      ap='<div class="oxbar" style="display:inline-block;width:120px;vertical-align:middle">'
        +'<i style="width:'+(wmax?(r.apxs.wt/wmax*100).toFixed(1):0)+'%"></i></div>'
        +' <span>'+r.apxs.oxide+' '+r.apxs.wt+'±'+r.apxs.err+' wt%</span>';
    } else { ap='<span class="ev">未测/未映射</span>'; }
    h+='<div class="oxrow"><span style="min-width:32px;font-weight:600">'+r.element
      +'</span><span style="min-width:200px">'+lib+'</span><span>'+ap+'</span></div>';
  });
  h+='<div class="ev">'+pl.note+'</div>';
  box.innerHTML=h; box.hidden=false;
}

function load(){
  if(!scanSel.options.length)return;
  var imageObs=SCANS.find(function(s){return String(s.index)===scanSel.value;});
  if(typeof resetScamImages==='function')resetScamImages(imageObs?imageObs.sol:null);
  var tok=++loadToken;
  loadPaperNotes();
  loadRMIs();
  loadCross();
  fetch('/api/plot?'+qs()).then(function(r){return r.json();})
   .then(function(d){ if(tok===loadToken)draw(d); });
}

var CHART=null, tipEl=null;
var VIEW=null;        /* {x0,x1} 当前缩放视窗，null=全谱 */
var REF=null;         /* 锁定对比的参考谱 {x,y,label} */
var curData=null;     /* 最近一次画的数据，供缩放重画 */
function resetView(){ VIEW=null; if(curData) draw(curData); }
window.addEventListener('resize',function(){if(curData)draw(curData);});
document.addEventListener('workbench-theme',function(){if(curData)draw(curData);});
function draw(d){
  if(d.error){
    CHART=null;curData=null;
    ['elements','oxides','mocbox','absorptionBox'].forEach(function(id){var box=$(id);if(box)box.hidden=true;});
    if(tipEl)tipEl.style.display='none';
    plot.innerHTML='<p class="hint">'+d.error+'</p>';
    $('ptlabel').textContent='第 '+(parseInt(pt.value,10)+1)+' 点 / 共 '+(parseInt(pt.max,10)+1)+' 点';
    return;
  }
  /* 发数下拉：后端每次带回当前点实际有的发数 */
  var sel=$('shot'), cur=sel.value;
  sel.hidden=!!d.passive||APP==='apxs';
  $('shotLabel').hidden=sel.hidden;
  document.querySelectorAll('#steps input').forEach(function(el){
    if(el.dataset.k==='dust'||el.dataset.k==='despike') el.disabled=!!d.passive;
    if(el.dataset.k==='reflectance')el.disabled=!(d.meta&&d.meta.kind==='cp2');
  });
  if(REF&&(REF.xlabel!==d.xlabel||REF.ylabel!==d.ylabel||REF.product!==d.product)){
    REF=null;$('lockBtn').textContent='锁定对比';
    $('lockInfo').textContent='产品类型或坐标/强度单位不同，已清除参考谱。';
  }
  sel.innerHTML=d.shots.map(function(s){
    return '<option value="'+s.i+'">'+s.label+'</option>'; }).join('');
  sel.value=(function(){
    for(var i=0;i<d.shots.length;i++){ if(String(d.shots[i].i)===String(cur)) return cur; }
    return '-1'; })();
  if(String(sel.value)!==String(cur)){ cur=sel.value; }

  curData=d;
  var W=Math.max(360,plot.clientWidth||680), H=380, m={l:74,r:22,t:30,b:56};
  var iw=W-m.l-m.r, ih=H-m.t-m.b, n=d.x.length;
  var fullX0=Math.min.apply(null,d.x), fullX1=Math.max.apply(null,d.x);
  var x0=fullX0, x1=fullX1;
  if(VIEW){
    x0=Math.max(VIEW.x0,fullX0); x1=Math.min(VIEW.x1,fullX1);
    if(!(x1>x0)){ VIEW=null; x0=fullX0; x1=fullX1; }
  }
  var i0v=0, i1v=n-1;
  while(i0v<i1v&&d.x[i0v]<x0)i0v++;
  while(i1v>i0v&&d.x[i1v]>x1)i1v--;
  var all=[];
  for(var vi=i0v;vi<=i1v;vi++){ var vv=d.series[0].y[vi];
    if(vv!==null&&isFinite(vv)) all.push(vv); }
  if(REF){ REF.y.forEach(function(rv,ri){
    if(rv!==null&&isFinite(rv)&&REF.x[ri]>=x0&&REF.x[ri]<=x1) all.push(rv); }); }
  if(!all.length){ plot.innerHTML='<p class="hint">这个区间里没有可绘制的数据</p>'; return; }
  var lo=Math.min.apply(null,all), hi=Math.max.apply(null,all);
  if(lo===hi){lo-=1;hi+=1;}
  var pad=(hi-lo)*0.08; lo-=pad; hi+=pad;
  var fy=function(v){return m.t+(hi-v)/(hi-lo)*ih;};
  var fxw=function(v){return m.l+(v-x0)/Math.max(1e-9,x1-x0)*iw;};
  var fx=function(i){return fxw(d.x[i]);};
  var MAIN='var(--accent)';

  var s=['<svg viewBox="0 0 '+W+' '+H+'" role="img" aria-label="光谱曲线">'];
  var k, i, v;
  for(k=0;k<=4;k++){
    v=lo+(hi-lo)*k/4; var gy=fy(v);
    s.push('<line class="gridline" x1="'+m.l+'" y1="'+gy.toFixed(1)+'" x2="'+(m.l+iw)+'" y2="'+gy.toFixed(1)+'"/>');
    var tickTxt=(Math.abs(v)>=1000)?(v/1000).toFixed(1)+'k':String(Number(v).toPrecision(3));
    s.push('<text class="tick" x="'+(m.l-8)+'" y="'+(gy+4).toFixed(1)+'" text-anchor="end">'
      +tickTxt+'</text>');
  }
  var tStep=Math.pow(10,Math.floor(Math.log10((x1-x0)/8||1)));
  var tMag=(x1-x0)/8/tStep;
  if(tMag>=5)tStep*=5; else if(tMag>=2)tStep*=2;
  for(var tv=Math.ceil(x0/tStep)*tStep; tv<=x1+1e-9; tv+=tStep){
    s.push('<text class="tick" x="'+fxw(tv).toFixed(1)+'" y="'+(m.t+ih+18)
      +'" text-anchor="middle">'+String(Number(tv.toFixed(Math.max(0,Math.min(8,Math.ceil(-Math.log10(tStep))+1)))))+'</text>');
  }
  s.push('<text class="axl" x="'+(m.l+iw/2)+'" y="'+(m.t+ih+42)+'" text-anchor="middle">'+d.xlabel+'</text>');
  s.push('<text class="axl" x="15" y="'+(m.t+ih/2)+'" text-anchor="middle" transform="rotate(-90 15 '
    +(m.t+ih/2)+')">'+d.ylabel+'</text>');

  if(d.absorption){d.absorption.references.forEach(function(b,i){
    var scale=d.absorption.axis_factor,a=b.range_um[0]*scale,z=b.range_um[1]*scale;
    if(z<x0||a>x1)return;
    var left=fxw(Math.max(a,x0)),right=fxw(Math.min(z,x1));
    s.push('<rect x="'+left+'" y="'+m.t+'" width="'+Math.max(0,right-left)+'" height="'+ih+'" fill="#728399" opacity=".07"/>');
    s.push('<text font-size="9" fill="#728399" x="'+Math.min(m.l+iw-3,(left+right)/2)+'" y="'+(m.t+9+(i%3)*12)+'" text-anchor="end">'+b.label+'</text>');
  });}
  /* 元素线标注：只画检出的，标签两行错开防压字 */
  var an=d.annotate||{};
  if(an.marks){
    var prevX=-999, row=0;
    an.marks.forEach(function(mk){
      if(!mk.detected) return;
      var cx=fxw(mk.x);
      if(cx<m.l+4||cx>m.l+iw-4) return;
      if(cx-prevX<38){ row=(row+1)%2; } else { row=0; }
      prevX=cx;
      s.push('<line class="peakline" x1="'+cx.toFixed(1)+'" y1="'+m.t+'" x2="'+cx.toFixed(1)
        +'" y2="'+(m.t+ih)+'" stroke="var(--accent)"/>');
      s.push('<text class="peaklabel" x="'+(cx+4).toFixed(1)+'" y="'+(m.t+11+row*15)
        +'" fill="var(--accent)">'+mk.elem+'</text>');
    });
  }

  if(REF&&REF.marks){
    REF.marks.forEach(function(mk){
      if(!mk.detected||mk.x<x0||mk.x>x1) return;
      var cx=fxw(mk.x);
      s.push('<line class="reference-line" x1="'+cx.toFixed(1)+'" y1="'+m.t+'" x2="'+cx.toFixed(1)+'" y2="'+(m.t+ih)+'" stroke="#8a94a0" stroke-dasharray="2 5" opacity=".5"/>');
    });
  }

  var main=d.series[0];
  var p='', first=-1, last=-1;
  for(i=i0v;i<=i1v;i++){
    v=main.y[i];
    if(v===null||!isFinite(v)){ continue; }
    p+=(p?' L':'M')+fx(i).toFixed(1)+' '+fy(v).toFixed(1);
    if(first<0) first=i;
    last=i;
  }
  if(first>=0&&last>first){
    var area='M'+fx(first).toFixed(1)+' '+(m.t+ih)+' '+p.replace(/^M/,'L')
        +' L'+fx(last).toFixed(1)+' '+(m.t+ih)+' Z';
    s.push('<path d="'+area+'" fill="var(--accent)" fill-opacity=".10" stroke="none"/>');
  }
  s.push('<path d="'+p+'" fill="none" stroke="'+MAIN
    +'" stroke-width="1.8" stroke-linejoin="round" stroke-linecap="round"/>');
  if(REF){
    var rp='';
    for(var rj=0;rj<REF.x.length;rj++){
      var rv=REF.y[rj];
      if(rv===null||!isFinite(rv)||REF.x[rj]<x0||REF.x[rj]>x1) continue;
      rp+=(rp?' L':'M')+fxw(REF.x[rj]).toFixed(1)+' '+fy(rv).toFixed(1);
    }
    if(rp) s.push('<path d="'+rp+'" fill="none" stroke="#8a94a0" stroke-width="1.6"'
      +' stroke-dasharray="6 4" stroke-linejoin="round"/>');
  }
  s.push('<line id="xhl" class="xhair" x1="0" y1="'+m.t+'" x2="0" y2="'+(m.t+ih)+'" visibility="hidden"/>');
  s.push('<circle id="xhd0" r="3.5" fill="'+MAIN+'" stroke="var(--panel)" stroke-width="1.5" visibility="hidden"/>');
  if(REF) s.push('<circle id="xhdr" r="3" fill="#8a94a0" stroke="var(--panel)" stroke-width="1.2" visibility="hidden"/>');
  s.push('<rect x="'+m.l+'" y="'+m.t+'" width="'+iw+'" height="'+ih+'" fill="transparent"/>');
  s.push('</svg>');
  plot.innerHTML=s.join('');
  if(!tipEl){ tipEl=document.createElement('div'); tipEl.className='charttip';
    tipEl.style.display='none'; }
  plot.appendChild(tipEl);
  tipEl.style.display='none';
  CHART={d:d,fx:fx,fxw:fxw,fy:fy,m:m,iw:iw,ih:ih,n:n,x0:x0,x1:x1,i0v:i0v,i1v:i1v};

  var lg=d.series.map(function(se){
    return '<span><i style="background:'+MAIN+'"></i>'+se.label+'</span>'; }).join('');
  if(REF) lg+='<span><i style="background:#8a94a0"></i>'+REF.label+'（参考）</span>';
  $('legend').innerHTML=lg;
  $('zoomInfo').innerHTML=VIEW?('已放大 '+Number(x0.toPrecision(6))+'–'+Number(x1.toPrecision(6))
    +'（双击谱图还原）'):'滚轮/框选放大 · 双击还原';
  logEl.innerHTML=d.log.map(function(l){
    return '<li><b>'+l.step+'</b><span>'+l.text+'</span></li>';}).join('');
  $('ptlabel').textContent='第 '+(parseInt(pt.value,10)+1)+' 点 / 共 '+d.n_points+' 点';
  var mt=d.meta||{};
  hint.textContent=(mt.sclk?('SCLK '+mt.sclk+' · '):'')+(mt.product||d.product||'');
  if(mt.note) hint.textContent+=' · '+mt.note;

  var eb=$('elements');
  if(d.passive){ eb.hidden=true; }
  else if(an.marks&&an.marks.length){
    var det=an.marks.filter(function(mk){return mk.detected;});
    var h='<div><b>元素线候选</b><span class="ev">　局部峰SNR≥5（噪声σ由差分MAD估计）；'
      +'LIBS需两条独立无重叠谱线支持。分数不是成分定量或检测置信概率。</span></div>';
    if(an.element_evidence&&an.element_evidence.length) h+='<div class="ev">逐发经验检出率的分母为预处理后全部保留发；SNR中位数按有局部峰的发计算（包括未达门槛发）。†为截断宽度下界，不能当完整FWHM。发次并非已验证的独立重复实验。</div>';
    if(an.elements&&an.elements.length){
      h+='<div style="margin:6px 0 4px">'+an.elements.map(function(e){
        var stats=(an.element_evidence||[]).find(function(r){return r.element===e;});
        var detail=(an.marks||[]).filter(function(m){return m.elem===e&&m.detected;}).map(function(m){return m.label+' SNR中位 '+m.shot_snr_median+' FWHM '+m.fwhm_observed+(m.fwhm_censored?'†':'')+' nm';}).join('；');
        return '<span class="tage" title="'+detail.replace(/"/g,'&quot;')+'">'+e+(stats?' '+stats.detected+'/'+stats.total+'发':'')+'</span>'; }).join('')+'</div>';
    } else {
      h+='<div class="ev" style="margin:6px 0 4px">没有达到检出阈值的线（或信噪比不足）</div>';
    }
    if(det.length){
      h+='<details style="margin-top:4px"><summary>检出的谱线（'+det.length+' 条）</summary><div>';
      det.forEach(function(mk){
        h+='<div class="row2"><span style="min-width:56px;font-weight:600">'+mk.elem
          +'</span><span>'+mk.label+' @ '+mk.x+'</span><span class="ev">SNR '
          +(mk.snr==null?'-':mk.snr)+' · '+(mk.reason||'候选')
          +(mk.shot_total==null?'':(' · '+mk.shot_detected+'/'+mk.shot_total+'发 · SNR中位 '
            +(mk.shot_snr_median==null?'—':Number(mk.shot_snr_median).toFixed(2))+'（'+mk.shot_snr_count+'发） · FWHM '
            +(mk.fwhm_observed==null?'—':Number(mk.fwhm_observed).toPrecision(3))+(mk.fwhm_censored?'†':'')+' nm'
            +(mk.shot_censored?' · '+mk.shot_censored+'发宽度截断†':'')))+'</span></div>'; });
      h+='</div></details>';
    }
    eb.innerHTML=h; eb.hidden=false;
  } else { eb.hidden=true; }
  renderAbsorption(d);
  renderQuant(d);

  var ox=$('oxides');
  if(d.oxides&&d.oxides.rows){
    var rows=d.oxides.rows, wmax=0;
    rows.forEach(function(r){ if(r.wt>wmax) wmax=r.wt; });
    var h='<div style="display:flex;align-items:baseline;gap:10px"><b>氧化物成分（APXS RWP）</b>'
      +'<span class="ev">'+d.oxides.file+' · 合计 '+d.oxides.total+' wt%</span></div>';
    rows.forEach(function(r){
      h+='<div class="oxrow"><span>'+r.oxide+'</span>'
        +'<div class="oxbar"><i style="width:'+(wmax?r.wt/wmax*100:0).toFixed(1)+'%"></i></div>'
        +'<span>'+r.wt+' ± '+r.err+' wt%</span></div>'; });
    ox.innerHTML=h; ox.hidden=false;
  } else { ox.hidden=true; }

  /* SuperCam 官方成分预测（PDS MOC 表查表展示，非本软件计算） */
  var mb=$('mocbox');
  if(d.moc){
    var m=d.moc;
    if(!m.have_table){
      mb.innerHTML='<div style="display:flex;align-items:baseline;gap:10px;flex-wrap:wrap">'
        +'<b>SuperCam 成分预测</b><span class="ev">需要官方成分表（约 2 MB，一次下载全程共用）</span>'
        +'<button id="mocDl" style="padding:3px 12px;font-size:12px">下载官方成分表</button>'
        +'<span class="ev" id="mocMsg"></span></div>';
      mb.hidden=false;
      $('mocDl').addEventListener('click',function(){
        var btn=this; btn.disabled=true;
        $('mocMsg').textContent='正在从 PDS 下载…';
        fetch('/api/moc/ensure?app=supercam').then(function(r){return r.json();}).then(function(r){
          if(r.error){ btn.disabled=false; $('mocMsg').textContent='失败：'+r.error; return; }
          $('mocMsg').textContent='完成（'+r.rows+' 个点的预测）';
          load();
        });
      });
    } else if(!m.matched){
      mb.innerHTML='<div style="display:flex;align-items:baseline;gap:10px">'
        +'<b>SuperCam 成分预测</b><span class="ev">该点在官方表里没有预测'
        +'（表共 '+m.table_rows+' 点；可能是新数据尚未收录或产品版本不同）</span></div>';
      mb.hidden=false;
    } else {
      var wmax2=0;
      m.rows.forEach(function(r){ if(r.wt!=null&&r.wt>wmax2) wmax2=r.wt; });
      var h='<div style="display:flex;align-items:baseline;gap:10px;flex-wrap:wrap">'
        +'<b>SuperCam 成分预测（PDS 官方模型）</b>'
        +'<span class="ev">合计 '+m.total+' wt%'
        +(m.distance?(' · 距离 '+Math.round(m.distance)+' mm'):'')
        +' · 表共 '+m.table_rows+' 点</span></div>';
      m.rows.forEach(function(r){
        var pct=(r.wt!=null&&wmax2)?(r.wt/wmax2*100).toFixed(1):0;
        var val=(r.wt==null)?'—':(r.wt+' ± '+(r.err==null?'?':r.err));
        h+='<div class="oxrow"><span>'+r.oxide+'</span>'
          +'<div class="oxbar"><i style="width:'+pct+'%"></i></div>'
          +'<span>'+val+' wt%</span></div>'; });
      var rm=m.rmsep||{}, rk=Object.keys(rm);
      if(rk.length){
        h+='<div class="ev" style="margin-top:4px">官方 RMSEP：'
          +rk.map(function(k){return k+' ±'+rm[k];}).join(' · ')+'</div>';
      }
      if(m.note) h+='<div class="ev">'+m.note+'</div>';
      mb.innerHTML=h; mb.hidden=false;
    }
  } else { mb.hidden=true; }
}

plot.addEventListener('mousemove',function(e){
  if(!CHART) return;
  var r=plot.getBoundingClientRect();
  var mx=e.clientX-r.left;
  var c=CHART;
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
  var v=c.d.series[0].y[i];
  var dot=document.getElementById('xhd0');
  var rows=[];
  if(v!==null&&isFinite(v)&&dot){
    dot.setAttribute('cx',xi); dot.setAttribute('cy',c.fy(v));
    dot.setAttribute('visibility','visible');
    rows.push('<b>'+c.d.series[0].label+'</b> '+Number(v).toPrecision(4));
  } else if(dot){ dot.setAttribute('visibility','hidden'); }
  var rdot=document.getElementById('xhdr');
  if(rdot) rdot.setAttribute('visibility','hidden');
  if(REF&&rdot&&c.d.x[i]>=REF.x[0]&&c.d.x[i]<=REF.x[REF.x.length-1]){
    var rlo=0, rhi=REF.x.length-1;
    while(rhi-rlo>1){ var rm2=(rlo+rhi)>>1;
      if(REF.x[rm2]<=c.d.x[i]) rlo=rm2; else rhi=rm2; }
    var ri2=(Math.abs(REF.x[rlo]-c.d.x[i])<Math.abs(REF.x[rhi]-c.d.x[i]))?rlo:rhi;
    var rv2=REF.y[ri2];
    if(rv2!==null&&isFinite(rv2)&&REF.x[ri2]>=c.x0&&REF.x[ri2]<=c.x1){
      rdot.setAttribute('cx',c.fxw(REF.x[ri2])); rdot.setAttribute('cy',c.fy(rv2));
      rdot.setAttribute('visibility','visible');
      rows.push('<b>'+REF.label+'（参考）</b> '+Number(rv2).toPrecision(4));
    } else { rdot.setAttribute('visibility','hidden'); }
  }
  tipEl.innerHTML='<b>'+Number(c.d.x[i]).toFixed(2)+' '+c.d.xlabel.replace(/.*\\((.+)\\).*/,'$1')+'</b><br>'+rows.join('<br>');
  tipEl.style.display='block';
  var tx=xi+14;
  if(tx+180>r.width) tx=xi-186;
  tipEl.style.left=tx+'px';
  tipEl.style.top=Math.max(4,e.clientY-r.top-16)+'px';
});
plot.addEventListener('mouseleave',function(){
  dragSel=null;
  var selection=document.getElementById('zoomrect');if(selection)selection.remove();
  if(tipEl) tipEl.style.display='none';
  var xhl=document.getElementById('xhl');
  if(xhl) xhl.setAttribute('visibility','hidden');
  var dot=document.getElementById('xhd0');
  if(dot) dot.setAttribute('visibility','hidden');
  var rdot=document.getElementById('xhdr');
  if(rdot) rdot.setAttribute('visibility','hidden');
});

/* ---- 缩放：滚轮以光标为中心放大/缩小；框选放大；双击还原 ---- */
plot.addEventListener('wheel',function(e){
  if(!CHART) return;
  e.preventDefault();
  var c=CHART, r=plot.getBoundingClientRect();
  var mx=e.clientX-r.left;
  var frac=(mx-c.m.l)/Math.max(1,c.iw);
  var val=c.x0+Math.max(0,Math.min(1,frac))*(c.x1-c.x0);
  var cur=curData;
  var lim0=Math.min.apply(null,cur.x), lim1=Math.max.apply(null,cur.x);
  var span=(c.x1-c.x0)*(e.deltaY>0?1.3:1/1.3);
  if(span>=(lim1-lim0)){ VIEW=null; draw(cur); return; }
  if(span<(lim1-lim0)*0.004) return;
  var nx0=val-frac*span, nx1=nx0+span;
  VIEW={x0:Math.max(lim0,nx0), x1:Math.min(lim1,nx1)};
  draw(cur);
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
  var old=document.getElementById('zoomrect');
  if(old) old.remove();
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
  var rect=document.getElementById('zoomrect');
  if(rect) rect.remove();
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

/* ---- 锁定对比：把当前谱固定为参考，之后每次翻点/换观测都与它叠加 ---- */
$('lockBtn').addEventListener('click',function(){
  if(!REF){
    if(!curData){ return; }
    var obsTxt=(scanSel.selectedIndex>=0?scanSel.options[scanSel.selectedIndex].text:'');
    REF={x:curData.x.slice(), y:curData.series[0].y.slice(),
         xlabel:curData.xlabel,ylabel:curData.ylabel,product:curData.product,
         marks:(curData.annotate||{}).marks||[],elements:(curData.annotate||{}).elements||[],
         label:obsTxt+' · 第 '+(Number(pt.value)+1)+' 点 · '+curData.series[0].label};
    this.textContent='取消对比';
    $('lockInfo').textContent='参考已锁定：'+REF.label+'（灰色虚线）；参考检出元素：'+(REF.elements.join('、')||'无/不适用');
  } else {
    REF=null;
    this.textContent='锁定对比';
    $('lockInfo').textContent='';
  }
  if(curData) draw(curData);
});

/* 主题 */
var rootEl=document.documentElement, savedTheme=null;
try{ savedTheme=localStorage.getItem('mars_consoles_theme'); }catch(e){}
if(savedTheme==='dark'||savedTheme==='light')
  rootEl.setAttribute('data-theme',savedTheme);
$('themeBtn').addEventListener('click',function(){
  var cur=rootEl.getAttribute('data-theme')||
    (window.matchMedia&&matchMedia('(prefers-color-scheme:dark)').matches?'dark':'light');
  var next=(cur==='dark')?'light':'dark';
  rootEl.setAttribute('data-theme',next);
  try{ localStorage.setItem('mars_consoles_theme',next); }catch(e){}
});
$('backBtn').addEventListener('click',function(){ location.href='__BACK_URL__'; });

/* 文件夹浏览 */
$('browseBtn').addEventListener('click',function(){ openBrowser(''); });
function openBrowser(p){
  fetch(api('browse?path='+encodeURIComponent(p||''))).then(function(r){return r.json();}).then(function(d){
    var box=$('browser'); box.hidden=false;
    var cur=d.path.replace(/"/g,'&quot;');
    var h='<div class="bar"><span class="sub" style="color:var(--muted);font-size:12px">浏览：<code>'+d.path+'</code></span>'
          +'<button id="useCurrent" data-p="'+cur+'">用当前目录</button>'
          +'<button class="ghost" id="closeBrowse" style="margin-left:auto">关闭</button></div>'
          +'<div class="bar" style="margin-top:8px"><label for="browserPath">目录路径</label>'
          +'<input id="browserPath" type="text" class="grow" placeholder="粘贴文件夹完整路径">'
          +'<button id="browsePathBtn">前往</button><button id="usePathBtn">打开此路径</button></div>'
          +'<p id="browserError" class="msg err" role="alert" hidden></p><div class="blist">';
    if(d.parent) h+='<div class="brow"><button class="bname" data-p="'+d.parent+'">↑ 上一级</button></div>';
    if(d.error) h+='<p class="hint">打不开这个目录：'+d.error+'</p>';
    d.dirs.forEach(function(x){
      var p2=x.path.replace(/"/g,'&quot;');
      h+='<div class="brow"><button class="bname" data-p="'+p2+'">'+x.name
         +(x.tag?('<span class="dot">'+x.tag+'</span>'):'')+'</button>';
      if(x.tag) h+='<button class="buse" data-p="'+p2+'">用这个</button>';
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
  fetch(api('load?path='+encodeURIComponent(p))).then(function(r){return r.json();}).then(function(r){
    if(r.error){
      var err=$('browserError');
      if(err){ err.textContent=r.error; err.hidden=false; } else alert(r.error);
      return;
    }
    $('browser').hidden=true; showFolder(r.folder); useScans(r.scans);
  });
}

/* 导出 */
var timer=null;
$('exportBtn').addEventListener('click',function(){
  var out=$('outDir').value.trim();
  if(!out){ alert('请先填 PNG 输出目录'); return; }
  var total=SCANS.reduce(function(a,s){return a+(s.points||0);},0);
  if(!confirm('将按当前预处理开关，导出全部 '+SCANS.length+
              ' 次观测、共 '+total+' 张总谱 PNG（平均谱）到\\n'+out+'\\n导出过程中可点「取消导出」中止。')) return;
  $('exportBtn').disabled=true; $('cancelExportBtn').hidden=false;
  $('prog').hidden=false; $('progBar').style.width='0%';
  $('progMsg').className='msg'; $('progMsg').textContent='正在导出…';
  fetch('/api/export?out='+encodeURIComponent(out)+
        '&steps='+encodeURIComponent(stepParam())+'&app='+APP)
   .then(function(r){return r.json();}).then(function(r){
    if(r.error){ finishExport({running:false,error:r.error}); return; }
    timer=setInterval(poll,500);
  });
});
$('cancelExportBtn').addEventListener('click',function(){
  $('cancelExportBtn').disabled=true;
  fetch(api('export/cancel')).then(function(r){return r.json();}).then(function(r){
    if(!r.cancelling){ $('cancelExportBtn').hidden=true; $('cancelExportBtn').disabled=false; }
  });
});
function poll(){ fetch(api('status')).then(function(r){return r.json();}).then(function(s){
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

/* ---- 质控报告（HTML 汇总 + 逐点表 CSV）---- */
var qcTimer=null;
$('qcBtn').addEventListener('click',function(){
  $('qcBtn').disabled=true;
  $('qcMsg').hidden=false; $('qcMsg').className='msg';
  $('qcMsg').textContent='正在生成质控报告…（逐点跑完当前文件夹的全部观测）';
  fetch(api('qc?configured=1&steps='+encodeURIComponent(stepParam()))).then(function(r){return r.json();}).then(function(r){
    if(r.error){ $('qcBtn').disabled=false; $('qcMsg').className='msg err';
      $('qcMsg').textContent='失败：'+r.error; return; }
    if(!qcTimer) qcTimer=setInterval(pollQc,600);
  });
});
function pollQc(){
  fetch(api('status')).then(function(r){return r.json();}).then(function(s){
    var q=s.qc||{};
    if(q.running){
      $('qcMsg').textContent='正在生成质控报告… '+q.done+' / '+q.total+' 点';
    } else {
      clearInterval(qcTimer); qcTimer=null; $('qcBtn').disabled=false;
      if(q.error){ $('qcMsg').className='msg err'; $('qcMsg').textContent='质控报告失败：'+q.error; return; }
      if(q.index){
        $('qcMsg').className='msg ok';
        $('qcMsg').textContent='质控报告已生成：'+q.count+' 次观测、'+q.points
          +' 个点 → '+q.index+'（CSV 逐点表在同目录）';
        fetch('/api/openreport?app='+APP+'&path='+encodeURIComponent(q.index));
      }
    }
  });
}

/* 退出 */
$('quitBtn').addEventListener('click',function(){
  if(!confirm('退出后需要重新从启动器（或启动.bat）才能再用。确定退出吗？')) return;
  fetch(api('quit')).finally(function(){ document.body.innerHTML=
    '<div style="padding:60px;text-align:center;color:#5f6773;font-size:15px">'+
    '已退出，这个窗口可以关闭了。</div>'; });
});

/* 从 PDS 下载 */
var DL={timer:null};
function fmtMB(b){ return (b/1048576).toFixed(b>10485760?0:1)+' MB'; }
var DEFAULT_DOWNLOAD_ROOT='';
fetch('/api/paths?instrument='+APP).then(r=>r.json()).then(d=>{DEFAULT_DOWNLOAD_ROOT=d.download;}).catch(()=>{});
function dlDefaultDest(sol){
  if(selectedSols().length>1)sol='批量';
  return DEFAULT_DOWNLOAD_ROOT ? DEFAULT_DOWNLOAD_ROOT+'/'+sol : '';
}
$('dlOpenBtn').addEventListener('click',function(){
  $('dlCard').hidden=false;
  if($('dlSol').options.length) return;
  loadSols(false);
});
$('dlClose').addEventListener('click',function(){ $('dlCard').hidden=true; });
$('dlRefresh').addEventListener('click',function(){ $('dlSol').innerHTML=''; loadSols(true); });
$('dlSol').addEventListener('change',loadTypes);

function loadSols(force){
  $('dlInfo').className='msg'; $('dlInfo').textContent='正在读取 PDS 上的目录列表…';
  fetch(api('dl/sols'+(force?'?refresh=1':''))).then(function(r){return r.json();}).then(function(d){
    if(d.error){ $('dlInfo').className='msg err'; $('dlInfo').textContent=d.error; return; }
    var sel=$('dlSol'); sel.innerHTML='';
    d.sols.forEach(function(s){
      var o=document.createElement('option'); o.value=s; o.textContent=s; sel.appendChild(o);
    });
    if(d.sols.length) sel.value=d.sols[d.sols.length-1];
    $('dlInfo').className='msg';
    $('dlInfo').textContent='PDS 上共 '+d.sols.length+' 个火星日可选';
    loadTypes();
  }).catch(function(e){ $('dlInfo').className='msg err';
    $('dlInfo').textContent='读取失败：'+e; });
}
function loadTypes(){
  var sol=$('dlSol').value;
  if(!sol) return;
  $('dlTypes').innerHTML='<span class="lbl">正在读取该火星日的文件清单…</span>';
  fetch(api('dl/list?sol='+sol)).then(function(r){return r.json();}).then(function(d){
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
    $('dlInfo').textContent=sol+'：共 '+d.total_files+' 个文件，'+fmtMB(d.total_bytes);
    if(selectedSols().length>1)$('dlInfo').textContent+='；已选'+selectedSols().length+'天，类型清单取首个Sol，队列逐天筛选';
    $('dlDest').value=dlDefaultDest(sol);
    updateDlSel();
  });
}
function dlBoxes(){ return $('dlTypes').querySelectorAll('input[type=checkbox]'); }
function updateDlSel(){
  var n=0,b=0,f=0;
  Array.prototype.forEach.call(dlBoxes(),function(i){
    if(i.checked){ n++; b+=Number(i.dataset.b); f+=Number(i.dataset.n); }
  });
  $('dlSel').textContent=(selectedSols().length>1?'首个Sol：':'')+(n? ('选中 '+n+' 类，'+f+' 个文件，'+fmtMB(b)) : '尚未选择');
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
  fetch('/api/dl/start?app='+APP+'&sols='+encodeURIComponent(sols.join(','))
        +'&types='+encodeURIComponent(types.join(','))+'&dest='+encodeURIComponent(dest))
   .then(function(r){return r.json();}).then(function(r){
     if(r.error){ $('dlStart').disabled=false; $('dlCancel').hidden=true;
       $('dlMsg').className='msg err'; $('dlMsg').textContent=r.error; return; }
     if(!DL.timer) DL.timer=setInterval(pollDl,600);
   });
});
$('dlCancel').addEventListener('click',function(){ fetch(api('dl/cancel')); });
function pollDl(){
  fetch(api('status')).then(function(r){return r.json();}).then(function(s){
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
        $('dlMsg').textContent='下载完成：'+d.done+' 个文件 → '+d.out+'。回到顶部点「换文件夹」选这个目录即可看图。';
      }
      if(d.log&&d.log.length){ $('dlLog').hidden=false; $('dlLog').textContent=d.log.join('\\n'); }
    }
  });
}

function ping(){ fetch('/api/ping').catch(function(){}); }
ping(); setInterval(ping, 10000);
</script></body></html>"""


def _accent_dark(hex_color):
    """深色主题下的浅一档主色（简单提亮）。"""
    s = hex_color.lstrip("#")
    r, g, b = (int(s[i:i + 2], 16) for i in (0, 2, 4))
    r, g, b = (min(255, int(c + (255 - c) * 0.45)) for c in (r, g, b))
    return "#%02x%02x%02x" % (r, g, b)


def build_page(backend):
    return (PAGE.replace('<span class="lbl">观测（共',filtering_ui.HTML+'<span class="lbl">观测（共').replace('<details open style="margin-top:12px">',paper_notes.HTML+'<details open style="margin-top:12px">').replace('<div class="elebox" id="elements" hidden></div>','<div class="elebox" id="elements" hidden></div>'+absorption_ui.HTML+quant_ui.HTML+(rmi_ui.HTML if backend.id=='chemcam' else supercam_images_ui.HTML if backend.id=='supercam' else '')).replace('</script>',paper_notes.JS+spectrum_export.JS+absorption_ui.JS+quant_ui.JS+rmi_ui.JS+supercam_images_ui.JS+filtering_ui.JS+history.JS.replace('__POINT_BASE__','1')+download_queue.JS+'</script>').replace('id="exportBtn">导出 PNG</button>','id="exportBtn">导出 PNG</button>'+spectrum_export.HTML)
            .replace("__TITLE__", backend.rover_label + " · " + backend.name)
            .replace("__ROVER__", backend.rover_label + "（" + backend.rover.upper() + "）")
            .replace("__SUBTITLE__", backend.subtitle)
            .replace("__APP__", backend.id)
            .replace("__SHOT_HIDDEN__", "hidden" if backend.id == "apxs" else "")
            .replace("__ACCENT__", backend.accent)
            .replace("__ACCENT_D__", _accent_dark(backend.accent))
            .replace("__BACK_URL__", BACK_URL))


# --------------------------------------------------------------------------
# 服务
# --------------------------------------------------------------------------

def _qc_worker(qc, backend, obs_list, folder, steps_on=None):
    """后台线程：给当前文件夹生成质控报告（HTML 汇总 + 逐点表 CSV）。"""
    from instruments import report as _report
    qc.update(running=True, done=0, total=0, error=None, out=None, index=None)
    try:
        if not obs_list:
            qc["error"] = "这个文件夹里没有观测，先生成不了报告"
            return
        out_dir = os.path.join(folder, "_qc", "consoles")

        def prog(done, total):
            qc["done"], qc["total"] = done, total

        kwargs={'progress':prog,'data_dir':folder}
        if steps_on is not None:kwargs['steps_on']=dict(steps_on)
        res = _report.build_report(backend, obs_list, out_dir, **kwargs)
        qc.update(out=res["out"], index=res["index"],
                  count=res["count"], points=res["points"])
        # 跨仪器对照页（ChemCam×APXS：记忆文件夹；SuperCam×PIXL：清单+按需下载）
        try:
            if backend.id == "chemcam":
                apxs_dir = icross.apxs_source_folder(load_settings())
                if apxs_dir:
                    cx = _report.build_cross_page(backend, obs_list, out_dir, apxs_dir)
                    if cx and _report.link_cross_in_index(res["index"], cx["pairs"]):
                        qc["cross"] = cx["path"]
            elif backend.id == "supercam":
                cx = _report.build_pixl_cross_page(backend, obs_list, out_dir)
                if cx and _report.link_cross_in_index(res["index"], cx["pairs"]):
                    qc["cross"] = cx["path"]
        except Exception as exc:                             # noqa: BLE001
            qc["cross_error"] = "%s: %s" % (type(exc).__name__, exc)
    except Exception as exc:                                # noqa: BLE001
        qc["error"] = "%s: %s" % (type(exc).__name__, exc)
    finally:
        qc["running"] = False


def export_name(obs, index):
    point = obs.points[index]
    return "sol%s_%s_%s_%s_pt%03d.png" % (obs.sol, obs.seq, obs.product.lower(), point.sclk, index + 1)


def _export_worker(state, backend, out, steps_on):
    ex = state["export"][backend.id]
    ex.update(running=True, done=0, total=0, out=out, error=None, cancel=False)
    try:
        os.makedirs(out, exist_ok=True)
        obs_list = state["app"][backend.id]["observations"]
        ex["total"] = sum(o.n_points for o in obs_list)
        for obs in obs_list:
            for i in range(obs.n_points):
                if ex.get("cancel"):
                    ex["error"] = "已取消（已完成 %d 张）" % ex["done"]
                    return
                render_point_png(backend, obs, i, None, steps_on,
                                 os.path.join(out, export_name(obs, i)))
                ex["done"] += 1
    except Exception as exc:                                # noqa: BLE001
        ex["error"] = "%s: %s" % (type(exc).__name__, exc)
    finally:
        ex["running"] = False


def _dl_worker(state, backend, sol, types, dest):
    types=frozenset(types)
    download_queue.run(state['dl'][backend.id],sol,dest,lambda day,out:_dl_one(state,backend,day,types,out))

def _dl_one(state, backend, sol, types, dest):
    dl = state["dl"][backend.id]
    try:
        root = list(backend.pds_roots.values())[0]["url"]      # v1 每仪器一个源
        entries = pds_listing.list_files(root, sol)
        if not entries:
            dl["error"] = "这个火星日在 PDS 上没有文件（可能还没发布）"
            return
        grouped = backend.parse_listing(entries)
        files = []
        for g in grouped["groups"]:
            if g["key"] in types:
                files.extend(g["files"])
        files.sort(key=lambda e: (os.path.splitext(e["name"])[1].lower()
                                  not in (".lbl", ".xml"), e["name"]))
        if not files:
            dl["error"] = "没有匹配所选类型 %s 的文件" % ",".join(sorted(types))
            return
        ibase.download_files(files, dest, dl)
    except Exception as exc:                                # noqa: BLE001
        dl["error"] = "%s: %s" % (type(exc).__name__, exc)


def make_handler(state):
    pages = {bid: build_page(b).encode("utf-8") for bid, b in REGISTRY.items()}

    class Handler(BaseHTTPRequestHandler):
        def do_POST(self):
            p=self._params();b=self._backend(p)
            if urllib.parse.urlparse(self.path).path!='/api/paper' or b is None:return self._json({'error':'未知接口'},404)
            try:
                n=int(self.headers.get('Content-Length',0))
                if not 0<n<=16000:raise ValueError('标注请求长度无效')
                data=json.loads(self.rfile.read(n).decode('utf-8'));st=state['app'][b.id]
                si=int(p.get('scan',0));pi=int(p.get('point',0))
                if not 0<=si<len(st['observations']):raise ValueError('观测索引无效')
                obs=st['observations'][si]
                if not 0<=pi<len(obs.points):raise ValueError('点索引无效')
                return self._json(paper_notes.save(st['folder'],b.id,obs,obs.points[pi],data.get('scope'),data.get('status'),data.get('note','')))
            except (ValueError,TypeError,AttributeError,OSError) as exc:return self._json({'error':str(exc)},400)

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
            except Exception:                               # noqa: BLE001
                pass

        def _json(self, obj, code=200):
            self._send(code, "application/json",
                       json.dumps(obj, ensure_ascii=False).encode("utf-8"))

        def _params(self):
            q = urllib.parse.urlparse(self.path).query
            return {k: v[0] for k, v in urllib.parse.parse_qs(q).items()}

        def _backend(self, p):
            return REGISTRY.get(p.get("app", ""))

        def do_GET(self):
            path = urllib.parse.urlparse(self.path).path
            p = self._params()
            try:
                if path in ("/", "/index.html"):
                    bid = p.get("app", "")
                    if bid not in pages:
                        return self._send(404, "text/plain; charset=utf-8",
                                          b"unknown app: %r" % bid.encode())
                    return self._send(200, "text/html; charset=utf-8", pages[bid])
                if path == "/favicon.ico":
                    return self._send(204, "image/x-icon", b"")
                if path == "/api/ping":
                    state["last_ping"] = time.time()
                    return self._json({"ok": True})
                b = self._backend(p)
                if b is None:
                    return self._json({"error": "未知仪器：%s" % p.get("app")}, 404)
                st = state["app"][b.id]
                if path=='/api/recent':return self._json({'folders':history.rows(load_settings(),b.id)})
                if path in ('/api/scamimg','/api/scamimg/image'):
                    if b.id!='supercam':return self._json({'error':'图像浏览仅支持SuperCam工作台'},400)
                    try:
                        if path=='/api/scamimg':
                            return self._json(supercam_images.listing(p.get('sol',''),p.get('refresh')=='1'))
                        if p.get('kind')=='ci':png=supercam_images.ci_png(p.get('name',''),p.get('sol',''))
                        elif p.get('kind')=='watson':png=supercam_images.watson_png(p.get('stem',''))
                        else:raise ValueError('图像类型须为ci或watson')
                        return self._send(200,'image/png',png)
                    except Exception as exc:return self._json({'error':str(exc)},400)
                if path in ('/api/rmi','/api/rmi/image','/api/rmi/map','/api/rmi/overlay'):
                    if b.id!='chemcam':return self._json({'error':'该仪器未验证RMI浏览'},400)
                    si=int(p.get('scan',0))
                    if not 0<=si<len(st['observations']):return self._json({'error':'观测索引无效'},400)
                    obs=st['observations'][si]
                    if path=='/api/rmi':return self._json(rmi.payload(obs,int(p.get('point',0))))
                    if path=='/api/rmi/map':
                        out_dir=os.path.join(st['folder'],'_maps')
                        return self._json(rmi_overlay.build_element_map(b,obs,out_dir))
                    index=int(p.get('image',0));images=getattr(obs,'rmi_images',[])
                    if not 0<=index<len(images):return self._json({'error':'图像索引无效'},400)
                    if path=='/api/rmi/overlay':
                        row=images[index]
                        stem='map_sol%05d_%s_%s'%(int(row['sol']),row['seq'],row['sclk'])
                        elem=p.get('elem')
                        if elem:
                            pj=os.path.join(st['folder'],'_maps',stem+'_points.json')
                            lbl=rmi_overlay.ensure_edr_label(row['name'],row['sol'],download=False) or ''
                            if not(os.path.exists(pj) and os.path.exists(lbl)):
                                return self._json({'error':'叠加图尚未生成，请先点「生成元素分布图」'},404)
                            try:
                                png=rmi_overlay.highlight_png(lbl,row['path'],pj,elem)
                            except Exception as exc:
                                return self._json({'error':'高亮重绘失败：%s'%exc},400)
                            return self._send(200,'image/png',png)
                        full=os.path.join(st['folder'],'_maps',stem+'.png')
                        if not os.path.exists(full):return self._json({'error':'叠加图尚未生成，请先点「生成元素分布图」'},404)
                        with open(full,'rb') as fimg:return self._send(200,'image/png',fimg.read())
                    return self._send(200,'image/png',rmi.preview(images[index]['path']))
                if path=='/api/paper':
                    si=int(p.get('scan',0));pi=int(p.get('point',0))
                    if not 0<=si<len(st['observations']):return self._json({'error':'观测索引无效'},400)
                    obs=st['observations'][si]
                    if not 0<=pi<len(obs.points):return self._json({'error':'点索引无效'},400)
                    return self._json(paper_notes.read(st['folder'],b.id,obs,obs.points[pi]))
                if path == "/api/steps":
                    return self._json([{"key": k, "label": lab, "on": on}
                                       for k, lab, on in b.steps])
                if path == "/api/scans":
                    import uuid
                    token=uuid.uuid4().hex;st['filter_token']=token
                    listing=obs_listing(st['observations'],st['folder'],b.id)
                    if p.get('filter'):
                        result=filtering.select(b,list(st['observations']),p.get('q',''),_steps_of(b,p),lambda:st.get('filter_token')!=token)
                        return self._json(dict(result,folder=st['folder'],scans=[listing[i] for i in result['indices']]))
                    return self._json({"folder": st["folder"],
                                       "scans": obs_listing(st["observations"],st['folder'],b.id)})
                if path == "/api/browse":
                    return self._json(browse(p.get("path")
                                             or st["folder"]
                                             or os.path.expanduser("~")))
                if path == "/api/load":
                    st['filter_token']=None
                    target = p.get("path", "")
                    if not os.path.isdir(target):
                        return self._json({"error": "目录不存在：%s" % target})
                    target=os.path.abspath(target);observations=b.scan_folder(target)
                    history.store_folder(SETTINGS,target,b.id)
                    st["folder"],st['observations']=target,observations
                    return self._json({"folder": st["folder"],
                                       "scans": obs_listing(st["observations"],st['folder'],b.id)})
                if path == "/api/points":
                    obs = st["observations"][int(p.get("scan", 0))]
                    return self._json({"points": obs.n_points})
                if path == "/api/plot":
                    remember_instrument(b.id, SETTINGS)
                    scan_i = int(p.get("scan", 0))
                    point_i = int(p.get("point", 0))
                    shot = int(p.get("shot", -1))
                    d = plot_for(st, b, scan_i, point_i, shot,
                                 _steps_of(b, p))
                    d["product"] = obs_listing([st["observations"][scan_i]])[0]["product"]
                    return self._json(d)
                if path=='/api/csv':
                    si=int(p.get('scan',0))
                    if not 0<=si<len(st['observations']):return self._json({'error':'观测索引无效'},400)
                    obs=st['observations'][si]
                    return self._json(spectrum_export.from_backend(b,obs,int(p.get('point',0)),int(p.get('shot',-1)),_steps_of(b,p),p.get('out') or os.path.join(st['folder'],'_png')))
                if path == "/api/export":
                    if state["export"][b.id].get("running"):
                        return self._json({"error": "已经在导出了"})
                    out = p.get("out", "").strip()
                    if not out:
                        return self._json({"error": "没有指定输出目录"})
                    threading.Thread(target=_export_worker,
                                     args=(state, b, out,
                                           _steps_of(b, p) if p.get("steps")
                                           else {k: dflt for k, _l, dflt in b.steps}),
                                     daemon=True).start()
                    return self._json({"started": True})
                if path == "/api/openreport":
                    target = os.path.abspath(p.get("path", ""))
                    folder = st["folder"]
                    from pathlib import Path
                    target = str(Path(target).resolve())
                    if not (folder and Path(target).is_relative_to((Path(folder) / "_qc" / "consoles").resolve())
                            and Path(target).suffix.lower() == ".html"):
                        return self._json({"error": "报告路径不在当前数据目录里"}, 400)
                    if not os.path.isfile(target):
                        return self._json({"error": "文件不存在：%s" % target}, 400)
                    try:
                        open_app_window("file:///" + target.replace("\\", "/"), True)
                    except Exception:                           # noqa: BLE001
                        try:
                            os.startfile(target)                # noqa: S606
                        except Exception:                       # noqa: BLE001
                            pass
                    return self._json({"ok": True})
                if path == "/api/qc":
                    qc = state["qc"][b.id]
                    if qc.get("running"):
                        return self._json({"error": "已经在生成质控报告了"})
                    state["qc"][b.id] = {"running": True}
                    threading.Thread(target=_qc_worker,
                                     args=(state["qc"][b.id], b, list(st["observations"]), st["folder"],_steps_of(b,p) if ('steps' in p or p.get('configured')) else None), daemon=True).start()
                    return self._json({"started": True})
                if path == "/api/export/cancel":
                    ex = state["export"][b.id]
                    if ex.get("running"):
                        ex["cancel"] = True
                        return self._json({"ok": True, "cancelling": True})
                    return self._json({"ok": True, "cancelling": False})
                if path == "/api/status":
                    return self._json({"export": state["export"][b.id],
                                       "dl": download_queue.snapshot(state["dl"][b.id]),
                                       "qc": state["qc"][b.id],
                                       "build": BUILD})
                if path == "/api/cross":
                    if b.id == "supercam":
                        scan_i = int(p.get("scan", 0))
                        if not 0 <= scan_i < len(st["observations"]):
                            return self._json({"error": "观测索引无效"}, 400)
                        sobs = st["observations"][scan_i]
                        if getattr(sobs, "passive", False):
                            return self._json({"supported": False,
                                               "note": "被动谱观测没有 LIBS 逐发证据，不做对照。"})
                        ck = ("cross_pixl", sobs.sol, sobs.seq)
                        if ck not in state.get("cross_cache", {}):
                            try:
                                result = icross.pixl_cross(b, sobs)
                            except Exception as exc:              # noqa: BLE001
                                result = {"supported": True,
                                          "error": "PIXL 对照失败：%s" % exc}
                            state.setdefault("cross_cache", {})[ck] = result
                        return self._json(state["cross_cache"][ck])
                    if b.id != "chemcam":
                        return self._json({"supported": False,
                                           "note": "跨仪器对照目前只在 ChemCam/SuperCam 台提供。"})
                    scan_i = int(p.get("scan", 0))
                    chem_obs = st["observations"][scan_i]
                    cfg_ = load_settings()
                    apxs_dir = icross.apxs_source_folder(cfg_)
                    ck = ("cross", apxs_dir, chem_obs.sol, chem_obs.seq)
                    if ck not in state.get("cross_cache", {}):
                        result = None
                        if apxs_dir:
                            try:
                                a_back = REGISTRY["apxs"]
                                a_obs = a_back.scan_folder(apxs_dir)
                                matches = icross.find_matches(chem_obs, a_obs)
                                result = {"supported": True, "apxs_folder": apxs_dir,
                                          "apxs_count": len(a_obs),
                                          "matches": [{k: m[k] for k in ("index", "basis", "detail")}
                                                      for m in matches],
                                          "payload": (icross.payload(b, chem_obs, a_back,
                                                                     matches[0]["obs"])
                                                      if matches else None)}
                            except Exception as exc:                  # noqa: BLE001
                                result = {"supported": True, "error": "APXS 扫描失败：%s" % exc}
                        else:
                            result = {"supported": True,
                                      "note": "还没有 APXS 数据文件夹：在 APXS 工作台里"
                                              "「换文件夹…」选一次 APXS 数据后，这里就能对照。"}
                        state.setdefault("cross_cache", {})[ck] = result
                    return self._json(state["cross_cache"][ck])
                if path == "/api/moc/ensure":
                    if b.id != "supercam":
                        return self._json({"error": "该仪器没有官方成分表"})
                    try:
                        p_ = b.ensure_moc(force=p.get("refresh"))
                        rows, _meta = b.moc_table(p_)
                        return self._json({"ok": True, "path": p_,
                                           "rows": len(rows) if rows else 0})
                    except Exception as exc:                  # noqa: BLE001
                        return self._json({"error": "下载失败：%s" % exc})
                if path == "/api/dl/sols":
                    if p.get("refresh") or "sols" not in st:
                        try:
                            root = list(b.pds_roots.values())[0]["url"]
                            st["sols"] = pds_listing.list_sols(root)
                        except Exception as exc:              # noqa: BLE001
                            return self._json({"error": "连不上 PDS：%s" % exc})
                    return self._json({"sols": st["sols"]})
                if path == "/api/dl/list":
                    try:
                        root = list(b.pds_roots.values())[0]["url"]
                        grouped = b.parse_listing(pds_listing.list_files(root,
                                                                         p.get("sol", "")))
                        if not grouped["groups"]:
                            return self._json({"error":
                                               "这个火星日没有可下载的文件（可能还没发布）"})
                        grouped["sol"] = p.get("sol", "")
                        return self._json(grouped)
                    except Exception as exc:                  # noqa: BLE001
                        return self._json({"error": "读取目录失败：%s" % exc})
                if path == "/api/dl/start":
                    dl = state["dl"][b.id]
                    if dl.get("running"):
                        return self._json({"error": "已经有一个下载在进行"})
                    try:sol=download_queue.normalize(p.get('sols') or p.get('sol',''))
                    except ValueError as exc:return self._json({'error':str(exc)},400)
                    types = set(x for x in p.get("types", "").split(",") if x)
                    dest = p.get("dest", "").strip()
                    if not sol or not types or not dest:
                        return self._json({"error": "参数不完整"})
                    state["dl"][b.id] = {"running": True, "cancel": False}
                    threading.Thread(target=_dl_worker,
                                     args=(state, b, sol, types, dest),
                                     daemon=True).start()
                    return self._json({"started": True})
                if path == "/api/dl/cancel":
                    state["dl"][b.id]["cancel"] = True
                    return self._json({"ok": True})
                if path == "/api/openfolder":
                    target = p.get("path", "")
                    if os.path.isdir(target):
                        try:
                            os.startfile(target)              # noqa: S606
                        except Exception:                     # noqa: BLE001
                            pass
                    return self._json({"ok": True})
                if path == "/api/quit":
                    self._json({"ok": True})
                    threading.Thread(target=state["shutdown"], daemon=True).start()
                    return
                return self._send(404, "text/plain; charset=utf-8", b"not found")
            except Exception as exc:                          # noqa: BLE001
                if path.startswith("/api/"):
                    return self._json({"error": "%s: %s" % (type(exc).__name__, exc)}, 500)
                return self._send(500, "text/plain; charset=utf-8",
                                  str(exc).encode("utf-8"))

    return runtime_http.decorate(Handler,state)


def build_tag():
    """console.py + instruments/*.py 的改动指纹，用于单实例新旧判断。"""
    h = hashlib.sha1()
    files = [os.path.join(HERE, "console.py"),
             os.path.join(HERE, "launch.py")]
    idir = os.path.join(HERE, "instruments")
    try:
        files += [os.path.join(idir, f) for f in sorted(os.listdir(idir))
                  if f.endswith(".py")]
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
    url = "http://127.0.0.1:%d/api/status?app=chemcam" % port
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "mars-consoles/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            d = json.loads(r.read().decode("utf-8", "replace"))
        return d.get("build") or "unknown"
    except Exception:                                       # noqa: BLE001
        return None


def _watchdog(state):
    while True:
        time.sleep(5)
        last = state.get("last_ping")
        if last and time.time() - last > PING_TIMEOUT:
            print("\n界面已经关闭，自动退出。")
            state["shutdown"]()
            return


def pick_start_folder(cfg,bid):
    folder=cfg.get('folders',{}).get(bid)
    return os.path.abspath(folder) if folder and os.path.isdir(folder) else appenv.empty_data(bid)


def start(port):
    """起服务但不进入主循环、不开窗口——给启动器（launch.py）内嵌用。

    返回 (httpd, state)；state["last_ping"] 供启动器的联合看门狗判断
    「工作台窗口还开着没有」。
    """
    runtime_http.initialize()
    cfg = load_settings()
    state = {"app": {}, "dl": {}, "export": {}, "qc": {}, "last_ping": None,
             "shutdown": lambda: None}
    for bid, b in REGISTRY.items():
        folder = pick_start_folder(cfg, bid)
        observations = []
        try:
            observations = b.scan_folder(folder)
        except Exception as exc:                             # noqa: BLE001
            print("扫描 %s 起始文件夹失败：%s" % (bid, exc))
        state["app"][bid] = {"folder": folder, "observations": observations}
        state["dl"][bid] = {}
        state["export"][bid] = {}
        state["qc"][bid] = {}
    httpd = runtime_http.bind(make_handler(state),port)
    state["shutdown"] = httpd.shutdown
    return httpd, state


def main(argv=None):
    ap=argparse.ArgumentParser(description='火星仪器操作台')
    ap.add_argument('--version',action='version',version=appenv.build_id())
    ap.add_argument('--app',choices=list(REGISTRY),default='chemcam')
    ap.add_argument('--port',type=int,default=8766)
    ap.add_argument('--browser',action='store_true')
    ap.add_argument('--no-window',action='store_true')
    ap.add_argument('--no-watchdog',action='store_true')
    args=ap.parse_args(argv)
    if not 0<=args.port<=65535:ap.error('端口须在0到65535之间')
    alive=running_instance(args.port) if args.port else None
    if alive==BUILD:
        if not args.no_window:open_app_window('http://127.0.0.1:%d/?app=%s'%(args.port,args.app),not args.browser)
        return 0
    server,state=start(args.port)
    url='http://127.0.0.1:%d/?app=%s'%(server.server_port,args.app)
    print('MarsWorkbench '+appenv.build_id()+' · 界面地址：'+url,flush=True)
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
