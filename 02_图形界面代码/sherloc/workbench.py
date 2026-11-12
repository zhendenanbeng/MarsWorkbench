# -*- coding: utf-8 -*-
"""生成质控台页面：把一次观测的解码结果、质控指标和逐点体检放进一个 HTML。"""

from __future__ import annotations

import html
import os
import platform

import numpy as np

from . import library, minerals, qc

# 跑报告时的解释器指纹。写进页脚，是为了让两份报告能不能直接对比一眼可辨：
# 同一份数据在不同 numpy 上跑出不同结果时，至少知道是哪一次跑的。
RUN_TAG = "Python %s · numpy %s" % (platform.python_version(), np.__version__)

CSS = """
:root{color-scheme:light dark;--bg:#eef0f3;--panel:#fff;--panel2:#f7f8fa;--ink:#15181c;
--muted:#5f6773;--line:#dfe3e9;--line2:#edf0f4;--accent:#1b62c4;--edge:#b57400;
--spike:#b3352f;--none:#c3c9d1;--ok:#1a7f4b}
@media (prefers-color-scheme:dark){:root{--bg:#101214;--panel:#1c1f23;--panel2:#23272c;
--ink:#e9ecf0;--muted:#98a1ac;--line:#343a41;--line2:#282d33;--accent:#7aa9ee;
--edge:#e0a63c;--spike:#e0796f;--none:#454b53;--ok:#5cc08b}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);
font:14px/1.6 ui-sans-serif,-apple-system,"Segoe UI","Microsoft YaHei",sans-serif}
main{max-width:1180px;margin:0 auto;padding:24px 20px 60px}
h1{font-size:20px;font-weight:500;margin:0 0 4px}
h2{font-size:15px;font-weight:500;margin:28px 0 10px;padding-bottom:6px;border-bottom:1px solid var(--line)}
.sub{color:var(--muted);font-size:12px}
.panel{background:var(--panel);border:1px solid var(--line);border-radius:8px;padding:16px}
.grid2{display:grid;grid-template-columns:minmax(0,1fr) 340px;gap:16px;align-items:start}
.verdict{display:flex;gap:12px;align-items:flex-start;border-left:3px solid var(--accent);
background:var(--panel);border-radius:6px;padding:12px 16px;margin:14px 0}
.verdict.bad{border-left-color:var(--spike)}
.verdict.good{border-left-color:var(--ok)}
.verdict b{font-weight:500}
table{width:100%;border-collapse:collapse;font-size:12px}
th{text-align:left;font-weight:400;color:var(--muted);padding:6px 10px 6px 0;border-bottom:1px solid var(--line)}
td{padding:6px 10px 6px 0;border-bottom:1px solid var(--line2);vertical-align:top;
overflow-wrap:anywhere;word-break:break-word}
td.num{font-variant-numeric:tabular-nums;white-space:nowrap}
code{font-size:12px;background:var(--panel2);padding:1px 4px;border-radius:3px;
overflow-wrap:anywhere;word-break:break-all}
main{max-width:1180px;margin:0 auto;padding:24px 20px 60px;overflow-x:hidden}
.kpi{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(150px,44%),1fr));gap:10px}
.kpi div{background:var(--panel2);border-radius:6px;padding:10px 12px}
.kpi .v{font-size:18px;font-variant-numeric:tabular-nums}
.kpi .k{font-size:11px;color:var(--muted)}
.imgwrap{position:relative;width:100%;aspect-ratio:4/3;background:var(--panel2);
border:1px solid var(--line);border-radius:6px;overflow:hidden}
.pt{position:absolute;width:18px;height:18px;margin:0;padding:0;transform:translate(-50%,-50%);
background:transparent;border:0;border-radius:50%;cursor:pointer;display:grid;place-items:center}
.pt i{display:block;width:6px;height:6px;border-radius:50%}
.pt[data-v=none] i{box-shadow:0 0 0 1px var(--none)}
.pt[data-v=spike] i{background:var(--spike)}
.pt[data-v=edge] i{background:var(--edge)}
.pt[data-v=band] i{background:var(--accent);box-shadow:0 0 0 2px color-mix(in srgb,var(--accent) 40%,transparent)}
.pt.on{box-shadow:0 0 0 2px var(--ink)}
.legend{display:flex;flex-wrap:wrap;gap:12px;font-size:11px;color:var(--muted);margin-top:8px}
.legend span{display:inline-flex;align-items:center;gap:5px}
.legend i{width:7px;height:7px;border-radius:50%}
.tag{display:inline-block;padding:1px 7px;border-radius:3px;font-size:11px}
.tag.band{background:color-mix(in srgb,var(--accent) 18%,transparent);color:var(--accent)}
.tag.edge{background:color-mix(in srgb,var(--edge) 22%,transparent);color:var(--edge)}
.tag.spike{background:color-mix(in srgb,var(--spike) 18%,transparent);color:var(--spike)}
svg{display:block;width:100%;height:auto}
.note{font-size:12px;color:var(--muted);margin:8px 0 0}
/* 矿物分布图：点按 RLS 真实坐标摆放，颜色按矿物类别 */
.mmap{position:relative;width:100%;aspect-ratio:4/3;margin-top:10px;
background:var(--panel2);border:1px solid var(--line);border-radius:6px;overflow:hidden}
/* 点上标的矿物字母用黑色（2026-09-17 改）。
   原来的白色在这个底色（--panel2 接近白）上根本读不出来：字母串比这个 20px 的圆
   宽得多（例如 "NaP+S+C+?"），大部分笔画其实落在圆**外面**的白底上，白字就没了。
   加一圈浅色描边是为了压在深色圆上时也读得清——不想要的话删掉 text-shadow 那行即可。 */
/* 圆点大小跟着点的疏密走（2026-09-17 修）。
   --dot 是「容器宽度的百分比」，由 _point_map() 按最近邻间距算好、写在 .mmap 上。
   原来写死 20px，对火山口底那种「100 点铺满 7×7 mm」的稀疏扫描没问题，
   但 1296 点的 Survey（点间距只有零点几个百分点）会重叠成一坨。
   高度乘 4/3 是因为容器是 4:3——这样圆才是正圆。 */
.mdot{position:absolute;border-radius:50%;
width:calc(var(--dot,1.6) * 1%);height:calc(var(--dot,1.6) * 1.3333%);
transform:translate(-50%,-50%);display:grid;place-items:center;
font-size:9px;font-weight:600;color:#111;line-height:1;
text-shadow:0 0 2px rgba(255,255,255,.95),0 0 1px rgba(255,255,255,.95);
border:1px solid rgba(255,255,255,.6);font-variant-numeric:tabular-nums}
.mdot.pm-empty{background:transparent;border:1px dashed var(--muted);color:var(--muted);
width:calc(var(--dot,1.6) * 0.55%);height:calc(var(--dot,1.6) * 0.7333%)}
.mdot.pm-faint{opacity:.45}
.plegend{display:flex;flex-wrap:wrap;gap:10px;font-size:11px;color:var(--muted);
margin-top:10px}
.plegend span{display:inline-flex;align-items:center;gap:5px}
.plegend i{width:11px;height:11px;border-radius:3px;display:inline-block}
.arrow{color:var(--muted);font-size:11px}
@media (max-width:900px){.grid2{grid-template-columns:1fr}}
"""

# 类别配色从 minerals.CLASS_COLORS 生成，网页和出图脚本共用同一份
CSS += "\n" + "\n".join(".pm-%s{background:%s}" % (k, v)
                        for k, v in minerals.CLASS_COLORS.items())

JS = """
(function(){
  var AXIS=/*__AXIS__*/;
  var SPEC=/*__SPEC__*/;
  var POS=/*__POS__*/, WID=/*__WID__*/, SNR=/*__SNR__*/, HGT=/*__HGT__*/;
  var GRID=/*__GRID__*/;
  var MINCH=/*__MINCH__*/;
  var state={pt:/*__FIRST__*/};
  var plot=document.getElementById('plot');
  var rows=document.getElementById('rows');
  function verdict(w){ if(w>=MINCH) return 'band'; if(w>=3) return 'edge'; if(w>=1) return 'spike'; return 'none';}
  var TAG={band:'宽峰候选',edge:'临界',spike:'窄尖刺',none:'无候选峰'};
  function build(){
    var wrap=document.getElementById('imgwrap');
    GRID.forEach(function(g,i){
      var b=document.createElement('button');
      b.type='button'; b.className='pt'; b.setAttribute('data-i',i);
      b.setAttribute('data-v',verdict(WID[i]));
      b.setAttribute('aria-label','点 '+i+'，'+TAG[verdict(WID[i])]);
      b.style.left=g[0]+'%'; b.style.top=g[1]+'%';
      b.innerHTML='<i></i>';
      b.addEventListener('click',function(){select(i);});
      wrap.appendChild(b);
    });
  }
  function draw(){
    var w=Math.max(320,plot.clientWidth||520),h=250,m={l:64,r:14,t:14,b:54};
    var iw=w-m.l-m.r, ih=h-m.t-m.b, n=AXIS.length, y=SPEC[state.pt];
    var lo=-320, hi=420;
    var fx=function(i){return m.l+i*iw/(n-1);};
    var fy=function(v){return m.t+(hi-v)*ih/(hi-lo);};
    var pxc=iw/(n-1);
    var s=['<svg viewBox="0 0 '+w+' '+h+'" width="'+w+'" height="'+h+'" role="img" aria-label="点 '+state.pt+' 的拉曼指纹区光谱">'];
    [300,200,100,0,-100,-200,-300].forEach(function(g){
      s.push('<line x1="'+m.l+'" y1="'+fy(g)+'" x2="'+(m.l+iw)+'" y2="'+fy(g)+'" stroke="'+(g===0?'var(--line)':'var(--line2)')+'"/>');
      s.push('<text x="'+(m.l-8)+'" y="'+(fy(g)+4)+'" text-anchor="end" font-size="11" fill="var(--muted)">'+g+'</text>');
    });
    [0,10,20,30,40].forEach(function(i){
      s.push('<text x="'+fx(i)+'" y="'+(m.t+ih+18)+'" text-anchor="middle" font-size="11" fill="var(--muted)">'+AXIS[i]+'</text>');
    });
    s.push('<text x="'+fx(20)+'" y="'+(m.t+ih+38)+'" text-anchor="middle" font-size="12" fill="var(--ink)">波数 (cm-1)</text>');
    s.push('<text x="12" y="'+(m.t+ih/2)+'" font-size="12" fill="var(--ink)" transform="rotate(-90 12 '+(m.t+ih/2)+')" text-anchor="middle">强度 (计数)</text>');
    var d='';for(var i=0;i<n;i++) d+=(i?' L':'M')+fx(i).toFixed(1)+' '+fy(y[i]).toFixed(1);
    s.push('<path d="'+d+'" fill="none" stroke="var(--accent)" stroke-width="1.6"/>');
    var wch=WID[state.pt];
    if(wch>=1){
      var pi=0,j;
      for(j=0;j<n;j++) if(Math.abs(AXIS[j]-POS[state.pt])<Math.abs(AXIS[pi]-POS[state.pt])) pi=j;
      var col='var(--'+(verdict(wch)==='band'?'accent':verdict(wch)==='edge'?'edge':'spike')+')';
      s.push('<rect x="'+(fx(pi)-wch/2*pxc).toFixed(1)+'" y="'+(m.t+4)+'" width="'+(wch*pxc).toFixed(1)+'" height="'+(ih-8)+'" fill="'+col+'" opacity="0.14"/>');
      s.push('<circle cx="'+fx(pi).toFixed(1)+'" cy="'+fy(y[pi]).toFixed(1)+'" r="3.5" fill="'+col+'" stroke="var(--panel)" stroke-width="1.5"/>');
    }
    var lw=MINCH*pxc, lx=m.l+12, ly=m.t+ih-14;
    s.push('<line x1="'+lx+'" y1="'+ly+'" x2="'+(lx+lw)+'" y2="'+ly+'" stroke="var(--ink)" stroke-width="1.5"/>');
    s.push('<line x1="'+lx+'" y1="'+(ly-5)+'" x2="'+lx+'" y2="'+(ly+5)+'" stroke="var(--ink)" stroke-width="1.5"/>');
    s.push('<line x1="'+(lx+lw)+'" y1="'+(ly-5)+'" x2="'+(lx+lw)+'" y2="'+(ly+5)+'" stroke="var(--ink)" stroke-width="1.5"/>');
    s.push('<text x="'+(lx+lw+8)+'" y="'+(ly+4)+'" font-size="11" fill="var(--ink)">仪器可分辨下限</text>');
    s.push('</svg>');
    plot.innerHTML=s.join('');
  }
  function fill(){
    var w=WID[state.pt], v=verdict(w);
    document.getElementById('ptlabel').textContent=state.pt;
    if(w===0){ rows.innerHTML='<tr><td colspan="5" style="color:var(--muted)">该点在 800-1200 cm-1 内没有任何超过 3 sigma 的局部极大值。</td></tr>'; }
    else {
      rows.innerHTML='<tr><td class="num">'+POS[state.pt]+'</td><td class="num">'+HGT[state.pt]+
        '</td><td class="num">'+w+' / '+MINCH.toFixed(2)+'</td><td class="num">'+SNR[state.pt].toFixed(2)+
        '</td><td><span class="tag '+v+'">'+TAG[v]+'</span></td></tr>';
    }
  }
  function select(i){
    state.pt=i;
    var all=document.querySelectorAll('.pt');
    for(var k=0;k<all.length;k++) all[k].classList.toggle('on',+all[k].getAttribute('data-i')===i);
    draw(); fill();
  }
  build(); select(state.pt);
  if(window.ResizeObserver){ var last=0; new ResizeObserver(function(){
    if(plot.clientWidth!==last){ last=plot.clientWidth; draw(); } }).observe(plot); }
})();
"""

PAGE = """<!doctype html>
<html lang="zh"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>SHERLOC 质控台 · Sol __SOL__ __SEQ__</title>
<style>__CSS__</style></head>
<body><main>
<h1>SHERLOC 质控台</h1>
<p class="sub">Sol __SOL__ · 序列 <code>__SEQ__</code> · SCLK <code>__SCLK__</code> · __TIME__ · __SEGMENT__</p>

<div class="verdict __VCLASS__"><div>
<b>__HEADLINE__</b>
<div class="sub" style="margin-top:4px">判据：峰宽须达仪器分辨下限 <b>__MINCH__ 个通道</b>（__RESCM__ cm-1），且信噪比 &gt; 3。窄于此的候选峰按宇宙射线残留处理。</div>
</div></div>

<h2>质控指标</h2>
<div class="kpi">
<div><div class="k">光谱点数</div><div class="v">__NPTS__</div></div>
<div><div class="k">通道间距</div><div class="v">__DCHAN__ <span class="sub">cm-1</span></div></div>
<div><div class="k">仪器分辨本领</div><div class="v">__RESCM__ <span class="sub">cm-1</span></div></div>
<div><div class="k">单点噪声 sigma</div><div class="v">__SIGMA__</div></div>
<div><div class="k">拉曼中位强度</div><div class="v">__RLEVEL__</div></div>
<div><div class="k">荧光峰 / 拉曼中位</div><div class="v">__FRATIO__</div></div>
</div>
<p class="note">判定分布：宽峰候选 __NB__ · 临界 __NE__ · 窄尖刺 __NS__ · 无候选峰 __NN__　
中位谱最强特征 __MEDPEAK__ 计数（约 __MEDSIG__ sigma）</p>

<h2>矿物类别判定</h2>
__MINERALS__

<h2>逐点矿物归属</h2>
__POINTMAP__

<h2>质控台</h2>
<div class="grid2">
  <div class="panel">
    <p class="sub" style="margin:0 0 6px">已选点 <b id="ptlabel">-</b> · Region 1 指纹区</p>
    <div id="plot"></div>
    <table><thead><tr><th>峰位 cm-1</th><th>峰高</th><th>半高跨通道</th><th>SNR</th><th>判定</th></tr></thead>
    <tbody id="rows"></tbody></table>
  </div>
  <div class="panel">
    <p class="sub" style="margin:0 0 6px">__GRIDTITLE__</p>
    <div class="imgwrap" id="imgwrap"></div>
    <div class="legend">
      <span><i style="background:var(--accent)"></i>宽峰候选 __NB__</span>
      <span><i style="background:var(--edge)"></i>临界 __NE__</span>
      <span><i style="background:var(--spike)"></i>窄尖刺 __NS__</span>
      <span><i style="box-shadow:0 0 0 1px var(--none)"></i>无候选峰 __NN__</span>
    </div>
  </div>
</div>

<h2>文件解码</h2>
<p class="note" style="margin-top:0">以 PDS4 label 为准：sol、序列号、SCLK、光谱条数、每张表的字节偏移都取自 label；
处理标记（波长校正 / 去宇宙射线 / 增益 / 归一化）label 未提供，取自文件名，并以 label 的 alternate_id 交叉核对。</p>
<table><thead><tr><th>文件</th><th>产品</th><th>处理标记</th><th>label 条数</th><th>label 表定义</th></tr></thead>
<tbody>__FILEROWS__</tbody></table>

<p class="note">由 sherloc_qc 生成 · 观测含 __NFILES__ 个文件 · 数据源目录 <code>__FOLDER__</code><br>
生成环境 <code>__ENV__</code> · 判据：__JUDGE__</p>
</main><script>__JS__</script></body></html>
"""


def _grid_from_rls(pts, n_points):
    """把 RLS 的 samp/line 归一化成百分比坐标。

    一次扫描的各段可能对应不同的 ACI 图像，此时把面板按图像切成竖条，
    每段用它自己那张图的边界归一化，避免把两张图混在一个坐标系里。
    """
    if not pts:
        out = []
        for i in range(n_points):
            r, c = divmod(i, 10)
            if r % 2:
                c = 9 - c
            out.append([12 + c * 76 / 9, 12 + r * 76 / 9])
        return out, "无 RLS 坐标，按 10x10 点阵示意"

    us = [p for p in pts if p["number"] < n_points]
    images = []
    for p in us:
        if p["image"] not in images:
            images.append(p["image"])
    band_w = 92.0 / len(images)
    out = [None] * n_points
    for k, img in enumerate(images):
        grp = [p for p in us if p["image"] == img]
        xs = [p["samp"] for p in grp]
        ys = [p["line"] for p in grp]
        x0, x1 = min(xs), max(xs)
        y0, y1 = min(ys), max(ys)
        dx, dy = max(x1 - x0, 1e-6), max(y1 - y0, 1e-6)
        left = 4 + k * band_w
        for p in grp:
            i = p["number"]
            if out[i] is None:
                out[i] = [round(left + (p["samp"] - x0) / dx * (band_w - 4), 2),
                          round(8 + (p["line"] - y0) / dy * 84, 2)]
    for i in range(n_points):
        if out[i] is None:
            out[i] = [50, 50]
    if len(images) == 1:
        title = "RLS 真实坐标 · 来自 1 张 ACI 图像"
    else:
        title = "RLS 真实坐标 · 来自 %d 张 ACI 图像（各占一条竖条）" % len(images)
    return out, title


def _judge_text(res):
    """页脚那一行判据说明。把"用的是哪份谱、哪个窗口"写清楚。

    这行是给未来的自己看的：同一份数据在不同配置下点数能差十倍以上，
    报告里不写清楚配置，过两天就没人知道这张表是怎么来的。
    """
    from . import minerals as _min

    nw = res.get("noise_window_used") or (2500.0, 3000.0)
    parts = [
        "拉曼分析只用 Region 1" if res.get("raman_region_only", True)
        else "拉曼分析用三表相加（R1+R2+R3）",
        "噪声取 %.0f-%.0f cm-1 逐点 RMS" % (nw[0], nw[1]),
        "σ 稳健替代%s" % ("开" if res.get("robust_sigma", True) else "关"),
        "峰宽 ≥ %d cm-1" % int(res.get("min_fwhm_used", res.get("min_fwhm_cm", 30))),
        "强度 ≥ %.0f 倍噪声" % float(res.get("snr_min_used", 3.0)),
        "峰心/峰宽用伪 Voigt 拟合",
        "水带拟合 %.0f-%.0f、峰位限 %.0f-%.0f、FWHM %.0f-%.0f"
        % (_min.OH_FIT_WINDOW + _min.OH_PEAK_WINDOW
           + (_min.OH_MIN_FWHM, _min.OH_MAX_FWHM)),
    ]
    if res.get("n_split"):
        parts.append("有 %d 个簇做了双峰拆分" % res["n_split"])
    n_rep = len(res.get("sigma_replaced") or [])
    if n_rep:
        parts.append("本次有 %d 个点的 σ 被中位值替换" % n_rep)
    return " · ".join(parts)


def _median_spacing_pct(grid):
    """百分比坐标下，最近邻距离的中位数，单位统一成「容器宽度的百分比」。

    容器是 4:3，所以 1% 的高度 = 0.75% 的宽度；纵向乘 0.75 才能和横向比。
    返回值给 _point_map() 定圆点大小用。
    """
    pts = [(g[0], g[1] * 0.75) for g in (grid or []) if g]
    if len(pts) < 2:
        return None
    a = np.asarray(pts, dtype=float)
    d = np.sqrt((a[:, None, 0] - a[None, :, 0]) ** 2
                + (a[:, None, 1] - a[None, :, 1]) ** 2)
    np.fill_diagonal(d, np.inf)          # 自己不算自己的邻居
    return float(np.median(d.min(axis=1)))


def _point_map(res, points_csv="", grid=None, grid_title="",
               pixl_evidence=None):
    """逐点矿物归属图：按 RLS 真实坐标摆放，颜色和字母都表示矿物类别。

    grid 来自 _grid_from_rls()，是把 RLS 的 samp/line 归一化后的百分比坐标。
    没有坐标时退回示意点阵，至少还能看出点数分布。
    """
    from . import minerals as _min

    det = res.get("point_detail") or []
    if not det:
        return ("<p class='note'>这次观测没有跑矿物判定，逐点表为空。"
                "（常见原因：label 里没有波长表，或 Region 1 表缺失。）</p>")
    if not grid:
        grid = []
        for i in range(len(det)):
            r, c = divmod(i, 10)
            grid.append([12 + (c if r % 2 == 0 else 9 - c) * 76 / 9,
                         12 + r * 76 / 9])
        grid_title = "无 RLS 坐标，按 10x10 示意点阵摆放"

    # 圆点大小跟着点的疏密走（2026-09-17 修）——原来固定 20px，对 1296 点的
    # Survey（点间距只有零点几个百分点）会重叠成一坨。dot 的单位是
    # 「容器宽度的百分比」，写进 .mmap 的 --dot 变量里给 CSS 用。
    sp = _median_spacing_pct(grid)
    dot = max(0.22, min(2.0, 0.40 * sp)) if sp else 1.6
    show_letters = dot >= 1.0        # 点太小的时候字母根本读不出来，只留色点

    cells, used = [], []
    n_lab, n_unlab, n_none = 0, 0, 0
    for d in det:
        i = int(d["point"])
        x, y = grid[i] if i < len(grid) else (50.0, 50.0)
        keys = [b["key"] for b in d["classes"]]
        if keys:
            n_lab += 1
            short = "+".join(_min.SHORT.get(k, "?") for k in keys)
            names = "/".join(b["class"] for b in d["classes"])
            tip = "点 %d：%s" % (i, names)
            if d["pos"] is not None:
                tip += "；最佳峰位 %.0f cm-1，FWHM %.0f，SNR %.1f" % (
                    d["pos"], d["fwhm"], d["snr"])
            if pixl_evidence:
                ev = pixl_evidence.get(i)
                if ev:
                    tip += "；PIXL：%s" % ev
            cls = "mdot " + " ".join("pm-" + k for k in keys)
            if d["snr"] is not None and d["snr"] < 3.0:
                cls += " pm-faint"          # 未达 3σ 的，画淡一点
            for k in keys:
                if k not in used:
                    used.append(k)
        elif d["n_bands"]:
            n_unlab += 1
            if d["labeled"]:
                cls, short = "mdot pm-empty", "?"
                tip = "点 %d：有合格波段（%.0f cm-1），但峰位不落在任何已知矿物类别区间" % (
                    i, d["pos"])
            else:
                cls, short = "mdot pm-empty", "·"
                tip = ("点 %d：有合格波段（%.0f cm-1），但峰位相同的点太少、未成组，"
                       "按当前规则不计入任何类别" % (i, d["pos"]))
        else:
            n_none += 1
            cls, short = "mdot pm-empty", ""
            tip = "点 %d：800-1200 cm-1 内没有合格波段" % i
        cells.append("<div class='%s' style='left:%.2f%%;top:%.2f%%' title='%s'>%s</div>"
                     % (cls, x, y, html.escape(tip, quote=True),
                        html.escape(short if show_letters else "")))

    legend = []
    for k in used:
        n = sum(g["n_points"] for g in (res.get("mineral_groups") or [])
                if g["classes"] and g["classes"][0]["key"] == k)
        legend.append("<span><i class='pm-%s'></i>%s %s%s</span>"
                      % (k, html.escape(_min.SHORT.get(k, "?")),
                         html.escape(_min.LABELS.get(k, k)),
                         ("（%d 点）" % n) if n else ""))

    if points_csv:
        csv_html = ("<p class='note'>逐点明细已导出为 "
                    "<a href=\"%s\" download><code>%s</code></a>，可用 Excel 直接打开"
                    "（UTF-8 BOM，中文不乱码）。论文备注及一致性来自数据目录的用户标注，"
                    "未标注时留空，不代表自动真值核验。</p>"
                    % (html.escape(points_csv), html.escape(points_csv)))
    else:
        csv_html = ""

    return (
        "<p class='note' style='margin-top:0'>每个点按 <b>RLS 坐标表里的真实位置</b>摆放"
        "（samp/line 就是 ACI 图像的像素坐标）。颜色和字母都表示矿物类别；"
        "空心虚线小圈 = 有合格波段但没归入任何类别，虚线更小的 = 该点没有合格波段。"
        "未达 3σ 的点画成半透明。</p>"
        "<div class='mmap' style='--dot:%.2f'>%s</div>"
        "<div class='plegend'>%s</div>"
        "<p class='note'>本次共 %d 点：判出类别的 <b>%d</b> 点，有合格波段但未归入任何类别的 %d 点，"
        "无合格波段的 %d 点。空缺格 = 该点没有 FWHM ≥ 30 cm⁻¹ 且 ≥3σ 的波段。"
        "一个点可以同时归入多类（论文也是这个口径），所以各类点数之和会大于判出类别的点数。</p>"
        "<p class='note'>坐标来源：%s</p>%s"
        % (dot, "\n".join(cells), "".join(legend), len(det), n_lab, n_unlab, n_none,
           html.escape(grid_title or "-"), csv_html))


def paper_html(res):
    from instruments.notes import label
    observation=res.get('paper_observation') or {};entries=res.get('paper_points') or {}
    if not observation and not any(entries.values()):return ''
    rows=''.join('<tr><td>%s</td><td>%s</td><td>%s</td></tr>'%(i,html.escape(label(r)),html.escape(r['note'])) for i,r in entries.items() if r)
    return '<section class="panel"><h2>论文对照标注（用户记录，非自动真值）</h2><p>'+html.escape(label(observation)+'；'+observation.get('note',''))+'</p><table><tr><th>原始点号</th><th>对照状态</th><th>备注</th></tr>'+rows+'</table></section>'

def build_report(obs, res, rls_pts=None, folder="", points_csv="",
                 pixl_evidence=None):
    """生成一次观测的质控台 HTML。"""
    spec = obs.spectra
    pts = res["points"]
    axis = [int(round(v)) for v in res["x"]]
    spec_lines = []
    for p in pts:
        i = p["point"]
        spec_lines.append(",".join(str(int(round(v))) for v in res["sub"][i]))
    pos = ",".join(str(int(round(p["peak"]["pos"]))) if p["peak"] else "0" for p in pts)
    wid = ",".join(str(p["peak"]["channels"]) if p["peak"] else "0" for p in pts)
    snr = ",".join("%.2f" % p["peak"]["snr"] if p["peak"] else "0" for p in pts)
    hgt = ",".join(str(int(round(p["peak"]["height"]))) if p["peak"] else "0" for p in pts)
    grid, grid_title = _grid_from_rls(rls_pts, len(pts))
    det = res.get("point_detail") or []
    mgrid, mtitle = (_grid_from_rls(rls_pts, len(det)) if det else ([], ""))
    first = next((i for i, p in enumerate(pts) if p["verdict"] == "band"),
                 next((i for i, p in enumerate(pts) if p["verdict"] == "edge"), 0))

    js = (JS.replace("/*__AXIS__*/", "[" + ",".join(str(v) for v in axis) + "]")
            .replace("/*__SPEC__*/", "[" + ",".join("[" + s + "]" for s in spec_lines) + "]")
            .replace("/*__POS__*/", "[" + pos + "]")
            .replace("/*__WID__*/", "[" + wid + "]")
            .replace("/*__SNR__*/", "[" + snr + "]")
            .replace("/*__HGT__*/", "[" + hgt + "]")
            .replace("/*__GRID__*/", "[" + ",".join("[%s,%s]" % (g[0], g[1]) for g in grid) + "]")
            .replace("/*__MINCH__*/", "%.2f" % res["min_width_channels"])
            .replace("/*__FIRST__*/", str(first)))

    file_rows = []
    for p in sorted(obs.products, key=lambda x: (x.product_type, x.name)):
        tables = "<br>".join("%s <span class='sub'>(records=%s)</span>" % (html.escape(t.name), t.records)
                             for t in p.tables()[:4]) or "<span class='sub'>label 无表定义</span>"
        file_rows.append(
            "<tr><td><code>%s</code></td><td>%s<br><span class='sub'>%s</span></td>"
            "<td class='sub'>%s</td><td class='num'>%s</td><td class='sub'>%s</td></tr>" % (
                html.escape(p.name), html.escape(p.product_type),
                html.escape(p.pid.product_name), html.escape(p.pid.flag_text),
                p.n_spectra if p.n_spectra else "-", tables))

    c = res["counts"]
    vclass = "good" if c["band"] >= max(3, 0.05 * len(pts)) else "bad"
    fr = res["fluor_raman_ratio"]

    # 矿物类别判定表
    groups = res.get("mineral_groups") or []
    if groups:
        rows = []
        for g in groups:
            b = g["band"]
            cls = g["classes"][0] if g["classes"] else None
            name = cls["class"] if cls else "未归类"
            key = cls["key"] if cls else "unassigned"
            ev = "<br>".join(html.escape(e) for e in (cls["evidence"] if cls else []))
            # 谱库比对结果（Razzell Hollis et al. 2021 的 62 种矿物实测峰位库）
            lib = g.get("library") or []
            if lib:
                v, top = library.verdict(lib)
                tag = "band" if v == "可信" else ("edge" if v == "有可能" else "spike")
                obs_txt = "、".join(
                    ("%.0f" % p) if not isinstance(p, (tuple, list))
                    else "%.0f(%.0fσ)" % (p[0], p[1])
                    for p in (g.get("observed_peaks") or []))
                lib_html = ("<span class='tag %s'>%s</span> %s<br>"
                            "<span class='sub'>%.2f ｜ 峰 %s</span>"
                            % (tag, html.escape(v), html.escape(top["mineral"]),
                               top["score"], html.escape(obs_txt or "-")))
            else:
                lib_html = "<span class='sub'>无峰可对</span>"
            rows.append(
                "<tr><td class='num'>%.0f</td><td class='num'>%.0f</td>"
                "<td class='num'>%d</td><td class='num'>%.0f</td><td class='num'>%.0f</td>"
                "<td class='num'>%.1f</td><td><span class='tag %s'>%s</span></td>"
                "<td class='sub'>%s</td><td>%s</td></tr>"
                % (g["pos"], b["center"], g["n_points"], b["fwhm"], b["snr"],
                   b["snr"],
                   "band" if key not in ("unassigned",) else "edge",
                   html.escape(name), ev, lib_html))
        min_html = (
            "<p class='note' style='margin-top:0'>按 Corpolongo et al. (2023) 的口径："
            "逐点找 FWHM ≥ %d cm⁻¹ 且 ≥3σ 的波段，再把峰位相近的点分组，"
            "<b>只用组内有信号的点</b>求平均谱，用平均谱的峰位定类别。"
            "全扫描共 <b>%s</b> 个点进入峰组（不含未成组候选）。</p>"
            "<p class='note'><b>当前为自动候选，可能同时漏检与误报。</b>点数不能作为"
            "真实检出率或统一上下界。峰位受背景、拟合和采样影响，部分开发样例与论文接近"
            "不代表普适精度；应结合组内复核标记、原始谱与独立证据判断。</p>"
            "<table><thead><tr><th>聚类峰位</th><th>平均谱质心</th><th>点数</th>"
            "<th>FWHM</th><th>SNR</th><th></th><th>类别</th><th>依据</th>"
            "<th>谱库比对</th></tr></thead>"
            "<tbody>%s</tbody></table>"
            % (int(res.get("min_fwhm_cm", 30)),
               res.get("n_points_with_band", "-"), "\n".join(rows)))
    elif res.get("scan_level") and res["scan_level"].get("bands"):
        # 逐点聚类一个类别都没形成，但扫描级兜底有候选：别说成"没有信号"
        min_html = ("<p class='note'>逐点聚类没有形成任何矿物类别——过门的点太少，"
                    "聚类没有统计意义。<b>但全扫描平均谱上有候选带</b>，见下。</p>")
    else:
        min_html = ("<p class='note'>没有检测到任何峰位能对上已知矿物类别的波段。"
                    "这不是处理参数问题：逐点判据和扫描级兜底都没有在 "
                    "800–1800 cm⁻¹ 内找到达标的带。</p>")

    # 扫描级兜底的结果：逐点路径数据不足时，全扫描平均谱上的候选带。
    # 单独一块、单独措辞——它**不是矿物类别**，绝不能混进上面的类别表，
    # 否则"点数"这一栏就跟论文没法比了。
    sl = res.get("scan_level")
    if sl and sl.get("bands"):
        srows = []
        for b in sl["bands"]:
            cls = "、".join(c["class"] for c in b["classes"]) or "未归类"
            ev = "<br>".join(html.escape(e) for c in b["classes"] for e in c["evidence"])
            srows.append(
                "<tr><td class='num'>%.0f</td><td class='num'>%.0f</td>"
                "<td class='num'>%.1f</td><td class='num'>%.2f</td>"
                "<td><span class='tag edge'>%s</span></td><td class='sub'>%s</td></tr>"
                % (b["center"], b["fwhm"], b["snr"], b["ratio_to_point_sigma"],
                   html.escape(cls), ev))
        min_html += (
            "<h3 style='margin-top:18px'>扫描级弱信号候选</h3>"
            "<p class='note' style='margin-top:0'>%s。上面那条判据要求"
            "<b>单个点</b>就达到 ≥3σ 且 FWHM ≥ %d cm⁻¹；这次只有 <b>%s</b> 个点过门，"
            "逐点聚类已经没有统计意义。于是退回来看<b>全部 %d 个点</b>的平均谱——"
            "它在下面这些位置上有带。</p>"
            "<p class='note'><b>这些候选带不构成矿物判定</b>，也不参与与论文的逐点对照："
            "它们是「这里可能有东西，去看一眼」的提示。两条门槛同时满足才列出："
            "相对平均谱自身噪声 ≥3σ，且幅度 ≥ %.2f× 单点噪声中位数。"
            "第二条是必须的——平均谱的噪声按 1/√N 缩小，上千点的扫描里几个计数也能凑出 4σ。</p>"
            "<p class='note'>上面质控指标里的「宽峰候选 N」是<b>另一处判据</b>数的："
            "它在整条谱 800–4000 cm⁻¹ 上量半高宽、且没有上限，荧光的宽包也会被算进去；"
            "矿物判定这条要求 FWHM 落在 30–200 cm⁻¹ 之间，所以两个数字对不上是正常的。"
            "以<b>矿物判定这一栏</b>为准。</p>"
            "<table><thead><tr><th>峰位</th><th>FWHM</th><th>相对平均谱噪声</th>"
            "<th>幅度/单点噪声</th><th>类别提示</th><th>依据</th></tr></thead>"
            "<tbody>%s</tbody></table>"
            % (html.escape(sl.get("why", "")), int(res.get("min_fwhm_cm", 30)),
               res.get("n_points_with_band", "-"), int(sl.get("n_points", 0)),
               minerals.SCAN_LEVEL_MIN_RATIO,
               "\n".join(srows)))
    elif not groups:
        min_html += ("<p class='note'>扫描级兜底也检查了：全扫描平均谱上同样没有"
                     "≥3σ 且幅度达标的带。</p>")

    # 有机质 G 带通道（~1600 cm⁻¹）。单独一块：它**不是矿物类别**，
    # 判据和出处都和上面不同（Hurowitz et al. 2025 用它判有机质）。
    gb_scan = res.get("g_band_scan")
    grows = []
    for g in groups:
        gb = g.get("g_band")
        if gb:
            grows.append((gb, g["n_points"], g["band"]["center"]))
    if gb_scan or grows:
        body = []
        if gb_scan:
            gb = gb_scan
            body.append(
                "<tr><td>全扫描平均谱</td><td class='num'>%.0f</td><td class='num'>%.0f</td>"
                "<td class='num'>%.1f</td><td class='num'>%.2f</td>"
                "<td><span class='tag %s'>%s</span></td></tr>"
                % (gb["center"], gb["fwhm"], gb["snr"],
                   gb["ratio_to_point_sigma"] if np.isfinite(gb["ratio_to_point_sigma"]) else 0.0,
                   "band" if gb.get("passed") else "edge",
                   "过 3σ 门槛" if gb.get("passed") else "未过 3σ 门槛"))
        for gb, npts, gpos in sorted(grows, key=lambda r: -r[0]["snr"]):
            body.append(
                "<tr><td>簇心 %.0f（%d 个点）</td><td class='num'>%.0f</td><td class='num'>%.0f</td>"
                "<td class='num'>%.1f</td><td class='num'>%.2f</td>"
                "<td><span class='tag %s'>%s</span></td></tr>"
                % (gpos or 0, npts, gb["center"], gb["fwhm"], gb["snr"],
                   gb["ratio_to_point_sigma"],
                   "band" if gb.get("passed") else "edge",
                   "过 3σ 门槛" if gb.get("passed") else "未过 3σ 门槛"))
        min_html += (
            "<h3 style='margin-top:18px'>有机质 G 带（~1600 cm⁻¹）</h3>"
            "<p class='note' style='margin-top:0'>SHERLOC 判定<b>有机质</b>靠的就是这条带"
            "（Hurowitz et al. 2025 在 Bright Angel 组就是这么认出有机碳的），"
            "它和矿物类别是两条独立通道。判据：相对平均谱自身噪声 ≥3σ，"
            "且幅度 ≥ %.2f× 单点噪声中位数（只看显著性不够——上千点的平均谱里"
            "噪声也够 3σ）。</p>"
            "<p class='note'><b>必须提醒的边界</b>：同一位置有<b>熔融石英 Si–O 伸缩泛音</b>"
            "背景（仪器光学窗口，Hurowitz et al. 2025）。论文是用仪器背景谱去比、"
            "才把 G 带分出来的；这里只报候选和实测数字，不下有没有有机质的结论。</p>"
            "<table><thead><tr><th>来源</th><th>峰位</th><th>FWHM</th>"
            "<th>相对平均谱噪声</th><th>幅度/单点噪声</th><th>判读</th></tr></thead>"
            "<tbody>%s</tbody></table>"
            % (minerals.G_BAND_MIN_RATIO, "\n".join(body)))

    html_out = (PAGE
                .replace("__CSS__", CSS)
                .replace("__JS__", js)
                .replace("__SOL__", obs.sol)
                .replace("__SEQ__", obs.sequence_id)
                .replace("__SCLK__", obs.sclk)
                .replace("__TIME__", (obs.label_time or "")[:19].replace("T", " "))
                .replace("__SEGMENT__", html.escape(
                    ("扫描由 %d 段拼接（expID %s，各 %s 条）"
                     % (len(obs.parts), "+".join(obs.expids),
                        "/".join(str(p.n_spectra) for p in obs.parts)))
                    if obs.is_split else "单段扫描"))
                .replace("__VCLASS__", vclass)
                .replace("__HEADLINE__", html.escape(qc.headline(res)))
                .replace("__MINCH__", "%.2f" % res["min_width_channels"])
                .replace("__RESCM__", "%.1f" % res["resolution_cm"])
                .replace("__NPTS__", str(res["n_points_total"]))
                .replace("__DCHAN__", "%.2f" % res["d_channel_cm"])
                .replace("__SIGMA__", "%.1f" % res["noise_sigma_median"])
                .replace("__RLEVEL__", "%.1f" % res["raman_level"])
                .replace("__FRATIO__", ("%.0f 倍" % fr) if fr else "-")
                .replace("__NB__", str(c["band"])).replace("__NE__", str(c["edge"]))
                .replace("__NS__", str(c["spike"])).replace("__NN__", str(c["none"]))
                .replace("__MEDPEAK__", "%.1f" % res["median_peak_height"])
                .replace("__MEDSIG__", "%.1f" % res["median_peak_sigma"])
                .replace("__MINERALS__", min_html)
                .replace("__POINTMAP__", _point_map(res, points_csv, mgrid, mtitle,
                                    pixl_evidence))
                .replace("__GRIDTITLE__", html.escape(grid_title))
                .replace("__FILEROWS__", "\n".join(file_rows))
                .replace("__NFILES__", str(len(obs.products)))
                .replace("__ENV__", html.escape(RUN_TAG))
                .replace("__JUDGE__", _judge_text(res))
                .replace("__FOLDER__", html.escape(folder)))
    import appenv
    return appenv.clean_report(html_out.replace('</main>',paper_html(res)+'</main>'))
