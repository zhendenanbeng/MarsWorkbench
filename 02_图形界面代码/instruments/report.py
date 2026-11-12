# -*- coding: utf-8 -*-
"""新仪器（SuperCam / ChemCam / APXS）的质控报告生成器。

产出与 SHERLOC 质控台同习惯：
  <数据目录>/_qc/index.html   汇总页（各观测概况 + 逐点表 + 代表谱图）
  <数据目录>/_qc/solXXXX_seq_points.csv   逐点表（检出元素、最强线、成分，
                                          留空的「论文对照」列供人工填写）

口径说明：元素线识别是快速候选筛选（局部峰SNR≥5，σ由差分MAD估计；LIBS需2条
独立无重叠线，APXS需至少1条无重叠线）；SuperCam 成分列是
PDS 官方模型预测值（查表，非本软件计算），APXS 成分列是官方 RWP 测量值。
报告里的谱图用当前仪器的默认预处理（除尘+去尖峰，不扣基线）。
"""

from __future__ import annotations

import csv
import datetime
import html
import os
import json
from pathlib import Path
import re
import urllib.parse
import appenv

import numpy as np
from . import evidence
from . import notes
from . import observation_summary
from . import observation_summary_ui


def _map_rows(folder, out_dir):
    """Only complete named map/coordinate pairs from the dataset's own _maps."""
    root = Path(folder).resolve()
    maps = root / '_maps'
    if not maps.is_dir() or not maps.resolve().is_relative_to(root):
        return []
    rows = []
    for path in sorted(maps.iterdir()):
        match = re.fullmatch(r'map_sol(\d{5})_(ccam\d{5})_(\d{9})\.png', path.name, re.I)
        if not match:
            continue
        csv_path = path.with_name(path.stem+'_points.csv')
        if not all(p.is_file() and p.resolve().is_relative_to(root) for p in (path, csv_path)):
            continue
        meta = path.with_name(path.stem+'_points.json')
        basis = '底图口径未记录'
        if meta.is_file() and meta.resolve().is_relative_to(root):
            try:
                basis = str(json.loads(meta.read_text(encoding='utf-8'))['basis']) or basis
            except (OSError, ValueError, KeyError, TypeError):
                basis = '底图口径元数据不可读，请核对坐标表'
        link = lambda p: urllib.parse.quote(os.path.relpath(p, out_dir).replace('\\', '/'), safe='/')
        rows.append({'sol': int(match[1]), 'seq': match[2].lower(), 'name': path.name,
                     'url': link(path), 'csv_url': link(csv_path), 'basis': basis})
    return rows


def _map_links(rows):
    if not rows:
        return ''
    return '<ul class="elementMapLinks">'+''.join(
        '<li><a href="%s">%s</a> · <a href="%s">逐点坐标CSV</a> · %s；元素快速候选，投影未经地面光斑核验，无比例尺</li>'
        % (html.escape(r['url'], quote=True), html.escape(r['name']), html.escape(r['csv_url'], quote=True), html.escape(r['basis']))
        for r in rows)+'</ul>'

CSS = """
:root{--bg:#f1f3f6;--panel:#fff;--panel2:#f7f8fa;--ink:#14171b;--muted:#5f6773;
--line:#e0e4ea;--line2:#eef1f5;--accent:#1a63c8;--ok:#12794a;--warn:#a86400;--radius:10px}
@media (prefers-color-scheme:dark){:root{--bg:#0f1113;--panel:#1b1e22;--panel2:#22262b;
--ink:#e9ecf0;--muted:#98a1ac;--line:#31363d;--line2:#262b31;--accent:#78a9ee;
--ok:#5cc08b;--warn:#e0a63c}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);
font:14px/1.6 ui-sans-serif,-apple-system,"Segoe UI","Microsoft YaHei",sans-serif}
main{max-width:1100px;margin:0 auto;padding:20px 18px 60px}
h1{font-size:19px;margin:0 0 2px}
.sub{color:var(--muted);font-size:12px;margin-bottom:16px}
.card{background:var(--panel);border:1px solid var(--line);border-radius:var(--radius);
padding:14px 16px;margin-bottom:14px}
h2{font-size:15px;margin:0 0 8px}
table{border-collapse:collapse;width:100%;font-size:12px}
th,td{border-bottom:1px solid var(--line2);padding:4px 8px;text-align:left;
white-space:normal;overflow-wrap:anywhere}
.card{overflow-x:auto}
th{color:var(--muted);font-weight:600;position:sticky;top:0;background:var(--panel)}
tr:hover td{background:var(--panel2)}
.tag{display:inline-block;padding:0 8px;border-radius:3px;font-size:11px;margin:0 2px 2px 0;
background:color-mix(in srgb,var(--accent) 15%,transparent);color:var(--accent);font-weight:600}
.ev{color:var(--muted);font-size:11px}
.num{font-variant-numeric:tabular-nums}
a{color:var(--accent)}
svg{display:block;width:100%;height:auto;background:var(--panel)}
"""

# 报告谱图的默认预处理 = 各仪器默认开关（与 /api/export 不带 steps 时一致）
def _default_steps(backend):
    return {k: on for k, _l, on in backend.steps}


def _fmt(v, nd=2):
    if v is None:
        return "—"
    try:
        return ("%.*f" % (nd, float(v)))
    except (TypeError, ValueError):
        return str(v)


def _spectrum_svg(backend, x, y, w=680, h=200):
    """报告内嵌代表谱图（折线 + 首末 x 刻度），下采样到 ~600 点控制体积。"""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(y)
    if not ok.any():
        return "<p class='ev'>（该点没有可绘制的数据）</p>"
    m = {  # margins
        "l": 52, "r": 10, "t": 8, "b": 20}
    iw, ih = w - m["l"] - m["r"], h - m["t"] - m["b"]
    step = max(1, len(x) // 600)
    xs, ys = x[ok][::step], y[ok][::step]
    lo, hi = float(np.min(ys)), float(np.max(ys))
    if hi <= lo:
        hi = lo + 1.0
    pad = (hi - lo) * 0.06
    lo, hi = lo - pad, hi + pad
    x0, x1 = float(x[ok][0]), float(x[ok][-1])
    fx = lambda v: m["l"] + (v - x0) / max(x1 - x0, 1e-9) * iw
    fy = lambda v: m["t"] + (hi - v) / (hi - lo) * ih
    pts = " ".join("%.1f,%.1f" % (fx(a), fy(b)) for a, b in zip(xs, ys))
    s = ['<svg viewBox="0 0 %d %d" role="img">' % (w, h)]
    for k in range(3):
        v = lo + (hi - lo) * (k + 0.5) / 3
        s.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="var(--line2)"/>'
                 % (m["l"], fy(v), w - m["r"], fy(v)))
        s.append('<text x="%d" y="%.1f" font-size="10" fill="var(--muted)" text-anchor="end">%s</text>'
                 % (m["l"] - 6, fy(v) + 3, _fmt(v, 1 if abs(v) >= 10 else 2)))
    s.append('<text x="%d" y="%d" font-size="10" fill="var(--muted)" text-anchor="middle">%s</text>'
             % (m["l"], h - 6, _fmt(x0, 1)))
    s.append('<text x="%d" y="%d" font-size="10" fill="var(--muted)" text-anchor="middle">%s</text>'
             % (w - m["r"], h - 6, _fmt(x1, 1)))
    s.append('<polyline points="%s" fill="none" stroke="var(--accent)" stroke-width="1.4"/>'
             % pts)
    s.append('</svg>')
    return "".join(s)


def moc_trend_svg(pairs, width=640, height=170):
    """SuperCam 官方成分预测趋势图（SiO2 与 FeOT 两条线，按点号）。

    pairs: [(点号, {氧化物: (wt, σ)}), …]（moc_payload 的 rows 展开形态）。
    点数 <2 或两氧化物全缺时返回 ""。只查表展示，不做任何拟合或插值。
    """
    if len(pairs) < 2:
        return ""
    series = {}
    for oxide in ("SiO2", "FeOT"):
        vals = [(no, ox.get(oxide, (None, None))[0]) for no, ox in pairs]
        if all(v is None for _n, v in vals):
            continue
        series[oxide] = vals
    if not series:
        return ""
    m = {"l": 46, "r": 12, "t": 18, "b": 26}
    iw, ih = width - m["l"] - m["r"], height - m["t"] - m["b"]
    xs = [p[0] for p in pairs]
    x0, x1 = min(xs), max(xs)
    all_vals = [v for vals in series.values() for _n, v in vals if v is not None]
    lo, hi = min(all_vals), max(all_vals)
    if hi <= lo:
        hi = lo + 1.0
    pad = (hi - lo) * 0.12
    lo, hi = lo - pad, hi + pad
    fx = lambda v: m["l"] + (v - x0) / max(x1 - x0, 1e-9) * iw
    fy = lambda v: m["t"] + (hi - v) / (hi - lo) * ih
    colors = {"SiO2": "#1a63c8", "FeOT": "#b3541e"}
    s = ['<svg viewBox="0 0 %d %d" role="img" aria-label="官方成分预测趋势">'
         % (width, height)]
    for k in range(3):
        v = lo + (hi - lo) * (k + 0.5) / 3
        s.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="var(--line2)"/>'
                 % (m["l"], fy(v), width - m["r"], fy(v)))
        s.append('<text x="%d" y="%.1f" font-size="10" fill="var(--muted)" text-anchor="end">%.0f</text>'
                 % (m["l"] - 5, fy(v) + 3, v))
    for oxide, vals in series.items():
        pts = " ".join("%.1f,%.1f" % (fx(n), fy(v)) for n, v in vals if v is not None)
        if pts:
            s.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="1.8"/>'
                     % (pts, colors.get(oxide, "#8a94a0")))
    s.append('<text x="%d" y="%d" font-size="10" fill="var(--muted)" text-anchor="middle">点号 %s</text>'
             % (m["l"] + iw // 2, height - 6, "%d–%d" % (x0, x1) if x0 != x1 else str(x0)))
    leg = " ".join('<span style="color:%s">%s</span>' % (colors.get(o, "#8a94a0"), o)
                   for o in series)
    return "".join(s) + '</svg><p class="ev">官方成分预测趋势（wt%%，查表展示非本软件计算）：%s</p>' % leg


def _top_lines(marks, n=3):
    det = [m for m in marks if m.get("detected") and not m.get("instrument")]
    det.sort(key=lambda m: -(m.get("snr") or 0))
    return [(m["label"], m["x"], m["snr"]) for m in det[:n]]


def build_report(backend, observations, out_dir, steps_on=None, progress=None, data_dir=None):
    """生成质控报告。返回 {count, points, out, index, csv}。

    progress(done, total) 每处理一个点回调一次，供界面进度条。
    """
    steps_on = _default_steps(backend) if steps_on is None else steps_on
    if data_dir is None:
        directories = [os.path.dirname(os.path.abspath(p)) for o in observations for pt in o.points for p in pt.files.values()]
        try:
            data_dir = os.path.commonpath(directories) if directories else out_dir
        except ValueError:
            data_dir = "多个数据目录"
    os.makedirs(out_dir, exist_ok=True)
    total_points = sum(o.n_points for o in observations)
    done = [0]

    summary, sections, csv_paths, line_paths = [], [], [], []
    paper_snapshot=notes.load(data_dir)['annotations']
    map_rows = _map_rows(data_dir, out_dir) if backend.id == 'chemcam' else []
    for obs in observations:
        observation_note=paper_snapshot.get(notes.obs_key(backend.id,obs))
        passive = bool(getattr(obs,'passive',False)) or (backend.id == "supercam" and obs.product in ("CP2", "CP3"))
        no_shots = passive or backend.id == "apxs"
        point_rows = []
        line_rows = []
        svg = ""
        ox_html = ""
        average=observation_summary.MeanSpectrum()
        moc_pairs=[]
        for i, pt in enumerate(obs.points):
            point_note=paper_snapshot.get(notes.point_key(backend.id,obs,pt))
            try:
                d = backend.load_point(pt, None, steps_on)
                if d.get("passive"):
                    marks, elements = [], []
                else:
                    marks, elements = backend.annotate_data(d)
            except Exception as exc:                  # noqa: BLE001
                point_rows.append({"point": i + 1, "sclk": pt.sclk, "error": str(exc)[:80],
                                   "paper_note":point_note,
                                   "elements": [], "top": [], "moc": None,
                                   "apxs": backend.oxides(obs, pt) if backend.id == "apxs" else None})
                if progress:
                    done[0] += 1
                    progress(done[0], total_points)
                continue
            average.add(d,backend)
            if not svg:
                svg = _spectrum_svg(backend, d["x"], d["y"])
                svg += "<p class='ev'>横轴：%s；纵轴：%s；谱图来自第%d点。</p>" % (
                    html.escape(d.get("xlabel") or backend.xlabel),
                    html.escape(d.get("ylabel") or backend.ylabel), i + 1)
            moc = None
            if backend.id == "supercam" and not passive:
                moc = backend.moc_payload(pt.stem)
                if moc.get("matched"):
                    moc_pairs.append((pt.point_no or (i + 1),
                                      {r["oxide"]: (r["wt"], r["err"]) for r in moc["rows"]}))
            for mark in marks:
                line_rows.append([i+1, pt.sclk, mark.get("label", ""), mark.get("x"),
                                  d.get("xlabel") or backend.xlabel, mark.get("snr"),
                                  bool(mark.get("detected")), bool(mark.get("instrument")),
                                  mark.get('peak_x'),mark.get('fwhm'),mark.get('reference_frame'),
                                  mark.get('independent_lines'),bool(mark.get('element_supported')),
                                  '；'.join(mark.get('ambiguous_with') or []),mark.get('reason',''),
                                  mark.get('shot_detected',''),mark.get('shot_total',''),mark.get('shot_available',''),
                                  mark.get('shot_snr_median',''),mark.get('shot_snr_count',''),
                                  mark.get('fwhm_observed',''),mark.get('fwhm_censored',''),mark.get('shot_censored','')])
            point_rows.append({"point": i + 1, "sclk": pt.sclk, "error": None,
                               "elements": elements, "top": _top_lines(marks),
                               "shot_evidence":evidence.summary(marks,elements),
                               "line_evidence":'；'.join(evidence.line_text(m) for m in marks if m['detected'] and 'shot_total' in m),
                               "element_evidence":evidence.element_statistics(marks,elements),
                               'all_shot_support':{m['elem']:m['element_shot_detected'] for m in marks if 'element_shot_detected' in m},
                               "paper_note":point_note,
                               'absorption':d.get('absorption'),
                               "n_shots": "" if no_shots else d.get("n_shots"), "moc": moc,
                               "apxs": backend.oxides(obs, pt) if backend.id == "apxs" else None})
            if progress:
                done[0] += 1
                progress(done[0], total_points)

        # 逐点 CSV（论文对照列留空，与 SHERLOC 质控台同习惯）
        prefix = "sol%s_%s_%s" % (obs.sol, obs.seq, obs.product.lower())
        csv_name = prefix + "_points.csv"
        csv_path = os.path.join(out_dir, csv_name)
        moc_cols = []
        if backend.id == "supercam":
            from .supercam import MOC_OXIDES
            moc_cols = list(MOC_OXIDES)
        apxs_names = list(dict.fromkeys(v["oxide"] for r in point_rows
                                       for v in (r.get("apxs") or {}).get("rows", [])))
        with open(csv_path, "w", encoding="utf-8-sig", newline="") as f:
            w = csv.writer(f)
            w.writerow(["点号", "SCLK", "发数", "检出元素", "最强谱线(SNR)","逐发元素支持","逐线经验统计","吸收参考（非鉴定）"]
                       + [name for c in moc_cols for name in ("官方预测 " + c + " wt%", "官方预测 " + c + " σ wt%")]
                       + [name for c in apxs_names for name in ("APXS " + c + " wt%", "APXS " + c + " σ wt%")]
                       + ["判定口径", "失败原因", "论文对照", "论文备注", "观测论文对照", "观测论文备注"])
            for r in point_rows:
                top = "; ".join("%s@%s(SNR%s)" % (a, b, c) for a, b, c in r["top"])
                moc_vals = []
                if moc_cols:
                    m = r["moc"] or {}
                    rec = m.get("rows") if m.get("matched") else None
                    for c in moc_cols:
                        v = next((x for x in (rec or []) if x["oxide"] == c), None)
                        moc_vals.extend(["" if not v or v["wt"] is None else v["wt"],
                                         "" if not v or v["err"] is None else v["err"]])
                components = {v["oxide"]:v for v in (r.get("apxs") or {}).get("rows", [])}
                apxs_vals = [components.get(c, {}).get(field, "") for c in apxs_names for field in ("wt", "err")]
                w.writerow([r["point"], r["sclk"], r.get("n_shots", ""),
                            " ".join(r["elements"]), top,r.get('shot_evidence',''),r.get('line_evidence',''),
                            '；'.join('%s 深度%.4g'%(a['label'],a['depth']) for a in (r.get('absorption') or {}).get('references',[]) if a['detected'])] + moc_vals + apxs_vals
                           + ["被动谱：不做发射线判定" if passive else
                              ("LIBS真空线表：需两条独立无重叠谱线支持元素候选" if backend.id!='apxs'
                               else "APXS谱线候选：排除背景与重叠；成分以配对RWP为准"), r["error"] or "",
                              notes.label(r.get('paper_note')),(r.get('paper_note') or {}).get('note',''),notes.label(observation_note),(observation_note or {}).get('note','')])
        csv_paths.append(csv_name)
        line_link = ""
        if not passive:
            line_name = prefix + "_lines.csv"
            with open(os.path.join(out_dir, line_name), "w", encoding="utf-8-sig", newline="") as f:
                w = csv.writer(f)
                w.writerow(["点号", "SCLK", "谱线", "位置", "横轴单位", "SNR", "达到门槛", "仪器线",
                            '实测局部峰位置','插值半高宽','参考坐标制','独立支持峰数',
                            '支持元素候选','重叠参考线','判定说明','检出发数','保留发数','覆盖发数',
                            '逐发SNR中位数','SNR可测发数','观测宽度或下界','宽度截断','截断发数'])
                w.writerows(line_rows)
            line_paths.append(line_name)
            line_link = ' · <a href="%s">全部谱线SNR CSV</a>' % html.escape(line_name)

        # 本节 HTML
        n_elems = sum(1 for r in point_rows if r["elements"])
        rows_html = []
        for r in point_rows:
            tags = "".join('<span class="tag">%s</span>' % html.escape(e)
                           for e in r["elements"][:10]) or '<span class="ev">—</span>'
            top = "; ".join("%s @%s (SNR %s)"
                            % (html.escape(a), b, c) for a, b, c in r["top"]) or "—"
            if r.get('shot_evidence'):top+='<br>'+html.escape(r['shot_evidence'])+'<details><summary>逐线经验统计</summary>'+html.escape(r.get('line_evidence',''))+'</details>'
            note=r.get('paper_note')
            if note:top+='<br>论文对照：'+html.escape(notes.label(note)+'；'+note['note'])
            if r.get('absorption'):
                absorption=r['absorption']
                top+='<br>吸收参考（非鉴定）：'+html.escape(absorption['column'])+'；'+('仅参考位置，当前DN未启用反射评分' if not absorption['enabled'] else html.escape(absorption['scope']))
                top+='<br>'+html.escape('；'.join('%s 深度%.3g'%(a['label'],a['depth']) for a in absorption['references'] if a['detected']))
            moc_cell = ""
            if backend.id == "supercam":
                m = r["moc"] or {}
                if passive:
                    moc_cell = '<span class="ev">不适用（被动光谱）</span>'
                elif m.get("matched"):
                    moc_cell = " · ".join(
                        "%s %s±%s" % (x["oxide"], _fmt(x["wt"], 1), _fmt(x["err"], 1))
                        for x in m["rows"] if x["wt"] is not None)
                elif m.get("have_table"):
                    moc_cell = '<span class="ev">官方表未收录该点</span>'
                else:
                    moc_cell = '<span class="ev">官方成分表未下载</span>'
            err_cell = ('<span class="ev">%s</span>' % html.escape(r["error"])
                        if r["error"] else "")
            rows_html.append(
                "<tr><td class='num'>%d</td><td class='num'>%s</td>"
                "<td class='num'>%s</td><td>%s</td><td>%s</td>%s<td>%s</td></tr>"
                % (r["point"], r["sclk"], r.get("n_shots", ""), tags, top,
                   ("<td>%s</td>" % moc_cell) if backend.id == "supercam" else "",
                   err_cell))
        moc_th = "<th>官方成分预测 wt%</th>" if backend.id == "supercam" else ""
        target = html.escape(getattr(obs, "target", "") or "")
        mean,mean_note=average.result()
        if mean is not None:svg+='<h3>观测平均谱</h3><div style="max-width:260px">'+_spectrum_svg(backend,*mean,w=260,h=85)+'</div><p class="ev">'+html.escape(mean_note)+'</p>'
        else:svg+='<p class="ev">'+html.escape(mean_note)+'</p>'
        sections.append("""
<div class="card" id="obs%(i)d">
<h2>Sol %(sol)s · %(seq)s · %(n)d 点 · %(prod)s%(tgt)s</h2>
<p class="ev">%(csv_link)s%(line_link)s · %(note)s</p>
%(svg)s
<table style="margin-top:8px"><tr><th>点</th><th>SCLK</th><th>发数</th>
<th>支持的元素候选</th><th>最强谱线候选</th>%(moc_th)s<th>备注</th></tr>
%(rows)s</table>%(map_links)s
</div>""" % {"i": len(sections), "sol": html.escape(str(obs.sol)),
             "seq": html.escape(str(obs.seq)), "n": obs.n_points,
             "prod": html.escape(str(obs.product)),
             "tgt": (" · " + target) if target else "",
             "csv_link": '<a href="%s">逐点表 CSV</a>' % html.escape(csv_name),
             "line_link": line_link,
             "note": ("被动谱不执行元素发射线识别，也不关联LIBS成分预测" if passive else "检出元素点数 %d / %d" % (n_elems, obs.n_points))+('；观测论文对照：'+html.escape(notes.label(observation_note)+'；'+observation_note['note']) if observation_note else ''),
             "ne": n_elems,
             "map_links": _map_links([r for r in map_rows if str(r['sol'])==(str(obs.sol).lstrip('0') or '0') and r['seq']==str(obs.seq).lower()]),
             "svg": svg + (moc_trend_svg(moc_pairs)
                           if backend.id == "supercam" else ""),
             "moc_th": moc_th,
             "rows": "".join(rows_html)})

        overview=observation_summary.summarize(point_rows,backend.id,passive)
        summary.append({"sol": obs.sol, "seq": obs.seq, "n": obs.n_points,
                        "product": obs.product, "ne": overview['element_count'], "i": len(sections) - 1,'overview':overview})

    # APXS：每个谱点只展示自身配对的官方氧化物表。
    if backend.id == "apxs":
        ox_sections = []
        for obs in observations:
            for pt in obs.points:
                ox = backend.oxides(obs, pt)
                if not ox:
                    continue
                rows = "".join(
                    "<tr><td>%s</td><td class='num'>%s</td><td class='num'>%s</td></tr>"
                    % (html.escape(r["oxide"]), _fmt(r["wt"], 3), _fmt(r["err"], 3))
                    for r in ox["rows"])
                ox_sections.append(
                    "<div class='card'><h2>Sol %s · %s · SCLK %s 氧化物成分（APXS RWP 官方测量值）</h2>"
                    "<p class='ev'>%s · 合计 %s wt%%</p>"
                    "<table><tr><th>氧化物</th><th>wt%%</th><th>±</th></tr>%s</table></div>"
                    % (html.escape(str(obs.sol)), html.escape(str(obs.seq)), html.escape(pt.sclk),
                       html.escape(ox["file"]), _fmt(ox.get("total"), 2), rows))
        sections.extend(ox_sections)

    official_names=('SiO2','FeOT') if backend.id=='supercam' else ('FeO','SiO2') if backend.id=='apxs' else ()
    sum_rows=[]
    for row in summary:
        overview=row['overview'];top='；'.join('%s %d/%d%s'%(r['element'],r['numerator'],r['denominator'],r['unit']) for r in overview['top']) or '—'
        cells=''
        for oxide in official_names:
            value=overview['official'][oxide]
            cells+='<td>'+('—' if value['mean'] is None else '%.3g（%d/%d点）'%(value['mean'],value['available_points'],value['total_points']))+'</td>'
        sum_rows.append('<tr data-sol="%s" data-points="%d" data-elements="%d"><td>%s</td><td>%s</td><td>%s</td><td>%d</td><td>%s</td><td>%s</td>%s<td><a href="#obs%d">打开</a></td></tr>'%(row['sol'],row['n'],-1 if row['ne'] is None else row['ne'],html.escape(str(row['sol'])),html.escape(row['seq']),html.escape(row['product']),row['n'],('未判定' if overview['applicable'] else '不适用') if row['ne'] is None else row['ne'],html.escape(top),cells,row['i']))
    sum_rows=''.join(sum_rows)
    official_headers=''.join('<th>%s 均值wt%%（覆盖点数）</th>'%name for name in official_names)

    moc_note = ""
    if backend.id == "supercam" and any(o.product == "CL1" for o in observations):
        rows_t, meta = backend.moc_table()
        if rows_t is None:
            moc_note = ("<p class='ev'>官方成分预测表未下载——工作台谱图下方点"
                        "「下载官方成分表」后重新生成报告即可带上成分列。</p>")
        else:
            rm = meta.get("RMSEP", {}) if meta else {}
            moc_note = ("<p class='ev'>成分列为 PDS 官方模型预测值（supercam_libs_moc，"
                        "共 %d 点，非本软件计算）；FeOT 为全铁。官方 RMSEP：%s</p>"
                        % (len(rows_t), " · ".join("%s ±%s" % (k, v)
                                                    for k, v in rm.items())))

    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    page = """<!doctype html><html lang="zh"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>%(title)s 质控报告</title><style>%(css)s</style></head><body><main>
<h1>%(title)s · 质控报告</h1>
<div class="sub">%(rover)s · 生成于 %(now)s · %(nobs)d 次观测、%(npts)d 个点 ·
数据目录 %(folder)s</div>
<div class="card"><h2>观测汇总</h2>
<p class="ev">点击Sol/点数/元素数排序。主要候选按全部有效点保留发数的经验支持率排序；APXS按有效点出现率。官方成分仅作有值点的描述性均值，不合并不确定度、不补零。平均谱为有效点等权，详见各节坐标对齐说明。</p>
<table id="observationSummary"><thead><tr><th data-sort="sol" tabindex="0" role="button">Sol ↕</th><th>序列</th><th>产品</th><th data-sort="points" tabindex="0" role="button">点数 ↕</th><th data-sort="elements" tabindex="0" role="button">候选元素数 ↕</th><th>主要候选前5</th>%(official_headers)s<th></th></tr></thead>
<tbody>%(sum_rows)s</tbody></table></div>
%(moc_note)s
<div class="card"><p class="ev">口径：谱线候选需局部极大、双侧背景及完整峰翼，局部峰SNR≥5（σ由差分MAD估计）。
LIBS使用真空波长，需至少两条独立无重叠线支持元素候选；APXS重叠峰须进一步解混，
仪器/大气线不计入目标元素。单线仍保留在谱线CSV。分数不是检测置信概率，也不作成分定量。
逐点表的论文对照来自数据目录_notes中的用户标注，非自动验证；未标注时留空。</p></div>
%(map_card)s%(sections)s
</main>%(sort_js)s</body></html>""" % {
        "title": html.escape(backend.rover_label + " · " + backend.name),
        "css": CSS, "rover": html.escape(backend.rover_label), "now": now,
        "nobs": len(observations), "npts": total_points,
        "folder": html.escape(str(data_dir)),
        "sum_rows": sum_rows, "official_headers":official_headers,"sort_js":observation_summary_ui.SORT_JS,"moc_note": moc_note,
        "map_card": ('<div class="card" id="elementMaps"><h2>元素分布图 %d 帧</h2>'%len(map_rows)+_map_links(map_rows)+'</div>') if map_rows else '',
        "sections": "".join(sections)}

    index_path = os.path.join(out_dir, "index.html")
    with open(index_path, "w", encoding="utf-8") as f:
        f.write(appenv.clean_report(page))
    return {"count": len(observations), "points": total_points,
            "out": out_dir, "index": index_path, "csv": csv_paths, "lines_csv": line_paths,'observations':summary}


def build_cross_page(backend, observations, out_dir, apxs_dir, max_points=9):
    """跨仪器对照页（ChemCam × APXS）：cross.html，从主 index 链接。

    只写有匹配的观测；无匹配时返回 None（主报告不加链接）。
    """
    from . import cross as icross
    from . import pds3 as _p3
    a_back = None
    for b in (backend,):
        pass
    from instruments import REGISTRY as _REG
    a_back = _REG["apxs"]
    a_obs = a_back.scan_folder(apxs_dir)
    sections = []
    for obs in observations:
        if getattr(obs, "passive", False):
            continue
        matches = icross.find_matches(obs, a_obs)
        if not matches:
            continue
        pl = icross.payload(backend, obs, a_back, matches[0]["obs"],
                            max_points=max_points)
        wmax = 0.0
        for r in pl["rows"]:
            if r["apxs"] and r["apxs"]["wt"] and r["apxs"]["wt"] > wmax:
                wmax = r["apxs"]["wt"]
        rows_html = []
        for r in pl["rows"]:
            lib = ("检出 %s 点 · 平均 %d%%发"
                   % (r["libs"]["points"], round(r["libs"]["mean_fraction"] * 100))
                   ) if r["libs"] else "—"
            if r["apxs"]:
                ap = ('<div style="display:inline-block;width:120px;height:10px;'
                      'background:var(--line2);border-radius:5px;vertical-align:middle">'
                      '<i style="display:block;height:100%%;width:%.1f%%;'
                      'background:var(--accent);opacity:.85"></i></div> %s %.2f±%.2f wt%%'
                      % (r["apxs"]["wt"] / wmax * 100 if wmax else 0,
                         html.escape(r["apxs"]["oxide"]), r["apxs"]["wt"] or 0,
                         r["apxs"]["err"] or 0))
            else:
                ap = '<span class="ev">未测/未映射</span>'
            rows_html.append(
                "<tr><td><b>%s</b></td><td>%s</td><td>%s</td></tr>"
                % (html.escape(r["element"]),
                   lib if r["libs"] else '<span class="ev">—</span>', ap))
        sections.append("""
<div class="card">
<h2>ChemCam Sol %s · %s · %s　×　APXS Sol %s · %s · %s</h2>
<p class="ev">%s（%s）</p>
<table><tr><th>元素</th><th>LIBS 检出证据（逐发口径）</th>
<th>APXS 官方 wt%%</th></tr>%s</table>
<p class="ev">%s</p>
</div>""" % (html.escape(str(obs.sol)), html.escape(str(obs.seq)),
             html.escape(getattr(obs, "target", "") or "?"),
             html.escape(str(pl["apxs"]["sol"])), html.escape(str(pl["apxs"]["seq"])),
             html.escape(pl["apxs"].get("target") or "?"),
             html.escape(matches[0]["detail"]),
             html.escape(matches[0]["basis"]),
             "".join(rows_html), html.escape(pl["note"])))
    if not sections:
        return None
    page = """<!doctype html><html lang="zh"><head><meta charset="utf-8">
<title>跨仪器对照（ChemCam × APXS）</title><style>%s</style></head><body><main>
<h1>跨仪器对照（ChemCam LIBS × APXS）</h1>
<div class="sub">匹配判据：site+drive 同停车位置 &gt; 目标名一致 &gt; sol 相邻。
<a href="index.html">返回质控汇总</a></div>
%s
</main></body></html>""" % (CSS, "".join(sections))
    path = os.path.join(out_dir, "cross.html")
    with open(path, "w", encoding="utf-8") as f:
        f.write(page)
    return {"path": path, "pairs": len(sections)}


def build_pixl_cross_page(backend, observations, out_dir, max_points=9):
    """跨仪器对照页（SuperCam × PIXL）：cross.html，从主 index 链接。

    PIXL 侧无本地文件夹——走 pixl_cross（清单缓存 + 候选按需下载）。
    只写有 site/drive 匹配的观测；无匹配返回 None。
    """
    from . import cross as icross

    sections = []
    errors = []
    for obs in observations:
        if getattr(obs, "passive", False):
            continue
        try:
            res = icross.pixl_cross(backend, obs, max_points=max_points)
        except Exception as exc:                          # noqa: BLE001
            errors.append("Sol %s: %s" % (obs.sol, str(exc)[:80]))
            continue
        pl = res.get("payload")
        if not pl:
            continue
        wmax = 0.0
        for r in pl["rows"]:
            if r["apxs"] and r["apxs"]["wt"] and r["apxs"]["wt"] > wmax:
                wmax = r["apxs"]["wt"]
        rows_html = []
        for r in pl["rows"]:
            lib = ("检出 %s 点 · 平均 %d%%发"
                   % (r["libs"]["points"], round(r["libs"]["mean_fraction"] * 100))
                   ) if r["libs"] else "—"
            if r["apxs"]:
                ap = ('<div style="display:inline-block;width:120px;height:10px;'
                      'background:var(--line2);border-radius:5px;vertical-align:middle">'
                      '<i style="display:block;height:100%%;width:%.1f%%;'
                      'background:var(--accent);opacity:.85"></i></div> %s %.2f±%.2f wt%%'
                      % (r["apxs"]["wt"] / wmax * 100 if wmax else 0,
                         html.escape(r["apxs"]["oxide"]), r["apxs"]["wt"] or 0,
                         r["apxs"]["err"] or 0))
            else:
                ap = '<span class="ev">未测/未映射</span>'
            rows_html.append(
                "<tr><td><b>%s</b></td><td>%s</td><td>%s</td></tr>"
                % (html.escape(r["element"]),
                   lib if r["libs"] else '<span class="ev">—</span>', ap))
        sections.append("""
<div class="card">
<h2>SuperCam Sol %s · %s · %s　×　PIXL Sol %s · 定位 %s · %s 点扫描中位</h2>
<p class="ev">%s（site/drive 同停车位置；rqc 为处理 sol 产品）</p>
<table><tr><th>元素</th><th>SuperCam LIBS 检出证据（逐发口径）</th>
<th>PIXL 逐点 XRF 扫描中位 wt%%</th></tr>%s</table>
<p class="ev">%s</p>
</div>""" % (html.escape(str(obs.sol)), html.escape(str(obs.seq)),
             html.escape(getattr(obs, "target", "") or "?"),
             html.escape(str(pl["pixl"]["sol"])), html.escape(str(pl["pixl"]["pos"])),
             pl["pixl"]["n_points"],
             html.escape(res["matches"][0]["detail"]),
             "".join(rows_html), html.escape(pl["note"])))
    if not sections:
        return None
    page = """<!doctype html><html lang="zh"><head><meta charset="utf-8">
<title>跨仪器对照（SuperCam × PIXL）</title><style>%s</style></head><body><main>
<h1>跨仪器对照（SuperCam LIBS × PIXL XRF）</h1>
<div class="sub">匹配判据：site+drive 同停车位置（PIXL label 无目标名）。
PIXL 侧为逐点 XRF（~120µm 斑）的扫描中位 wt%%。<a href="index.html">返回质控汇总</a></div>
%s
</main></body></html>""" % (CSS, "".join(sections))
    path = os.path.join(out_dir, "cross.html")
    with open(path, "w", encoding="utf-8") as f:
        f.write(page)
    return {"path": path, "pairs": len(sections)}


def link_cross_in_index(index_path, pairs, side="APXS"):
    """在已生成的主 index 口径卡里补一条对照页链接。"""
    with open(index_path, encoding="utf-8") as f:
        html_text = f.read()
    link = ('<div class="card"><p class="ev">跨仪器对照：%d 个观测匹配到同址/同名 '
            '%s，<a href="cross.html">打开对照页</a>。</p></div>' % (pairs, side))
    anchor = '<div class="card"><p class="ev">口径：'
    if anchor not in html_text:
        return False
    html_text = html_text.replace(anchor, link + anchor, 1)
    with open(index_path, "w", encoding="utf-8") as f:
        f.write(html_text)
    return True
