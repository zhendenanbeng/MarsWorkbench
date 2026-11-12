#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""火星仪器操作台 · 启动器。

先开一个「选仪器」小窗口；新仪器（SuperCam / ChemCam / APXS）的工作台
服务由本进程内嵌（console.start），SHERLOC 工作台是独立程序，
按卡片以子进程拉起（它自带单实例逻辑，重复点也不会开一堆）。

    python launch.py [--port 8790] [--console-port 8766] [--browser] [--no-window]

看门狗口径：选择器窗口和工作台窗口都关了才退出——只关选择器不影响
正在用的工作台。
"""

from __future__ import annotations

import source_paths  # source roots; no-op for bundled imports

import argparse
import json
import os
import subprocess
import sys
import threading
import time
import urllib.parse
import urllib.request
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from product_theme import validation_label                           # noqa: E402

import appenv
import runtime_http

from sherloc.desktop import open_app_window                  # noqa: E402

PING_TIMEOUT = 180.0
VIEW_SCRIPT = os.path.join(HERE, "view.py")


def configured_sherloc_port():
    try:
        with open(appenv.settings('sherloc'), encoding='utf-8') as f:
            port = int(json.load(f).get("port", 8765))
        if 1 <= port <= 65535:
            return port
    except (OSError, ValueError, TypeError):
        pass
    return 8765


_sherloc_services=[]
_launcher_context=None

def ensure_sherloc(port,timeout=30):
    import view
    url='http://127.0.0.1:%d/'%port
    if port:
        try:
            with urllib.request.urlopen(url+'api/status',timeout=.5) as r:
                if json.load(r).get('build')==view.BUILD:return url
        except (OSError,ValueError):pass
    server,state=view.start(port)
    if _launcher_context:state['launcher_url']=_launcher_context.get('launcher_url')
    threading.Thread(target=server.serve_forever,daemon=True).start()
    _sherloc_services.append((server,state))
    return 'http://127.0.0.1:%d/'%server.server_port


def alert(title, message):
    try:
        import ctypes
        ctypes.windll.user32.MessageBoxW(0, message, title, 0x10)
    except Exception:                                       # noqa: BLE001
        print("%s\n%s" % (title, message))


PAGE = """<!doctype html><html lang="zh"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>MarsWorkbench · 火星光谱工作台</title>
<style>
:root{--bg:#F2F3F4;--panel:#FFFFFF;--panel2:#EDEEEF;--ink:#161616;--muted:#5C5C5C;
--line:#D7D9DA;--line2:#E8EAEA;--accent:#7DE6FF;--accent-ink:#0C3844;
--ok:#2E6B52;--radius:0px}
@media (prefers-color-scheme:dark){:root{--bg:#191B1C;--panel:#222526;--panel2:#2B2E30;
--ink:#EDEDED;--muted:#A3A7AA;--line:#3A3E40;--line2:#2E3133;--accent:#7DE6FF;--accent-ink:#0C3844}}
html[data-theme="dark"]{--bg:#191B1C;--panel:#222526;--panel2:#2B2E30;
--ink:#EDEDED;--muted:#A3A7AA;--line:#3A3E40;--line2:#2E3133;--accent:#7DE6FF;--accent-ink:#0C3844}
html[data-theme="light"]{--bg:#F2F3F4;--panel:#FFFFFF;--panel2:#EDEEEF;
--ink:#161616;--muted:#5C5C5C;--line:#D7D9DA;--line2:#E8EAEA;--accent:#7DE6FF;--accent-ink:#0C3844}
*{box-sizing:border-box;border-radius:0!important}
html,body{margin:0;height:100%}
body{background:var(--bg);color:var(--ink);display:flex;flex-direction:column;
align-items:center;justify-content:center;min-height:100vh;
font:14px/1.6 'Saira','Microsoft YaHei','Noto Sans SC','PingFang SC',sans-serif;
font-variant-numeric:tabular-nums}
body::before{content:'';position:fixed;top:0;left:0;right:0;height:3px;z-index:100;
background:linear-gradient(90deg,#FF00F0 0 33.3%,#7DE6FF 33.3% 66.6%,#00FFA2 66.6% 100%)}
.cmdbar{position:fixed;top:3px;left:0;right:0;height:46px;background:#161616;color:#EDEDED;
display:flex;align-items:stretch;border-bottom:1px solid #000;z-index:90}
.cmdbar .brand{display:flex;align-items:center;gap:10px;padding:0 22px;font-weight:800;
font-stretch:125%;font-size:14px;letter-spacing:.08em}
.cmdbar .brand .seal{display:inline-flex;align-items:center;justify-content:center;width:24px;height:24px;
background:var(--accent);color:#0C3844;font-size:13px;font-weight:900;
clip-path:polygon(0 0,100% 0,100% calc(100% - 6px),calc(100% - 6px) 100%,0 100%)}
.cmdbar .ver{margin-left:auto;display:flex;align-items:center;padding-right:20px;
font-size:11px;letter-spacing:.14em;font-stretch:125%;opacity:.4;text-transform:uppercase}
main{width:100%;max-width:74rem;margin:0 auto;padding:calc(49px + 2.5rem) 24px 4rem}
.hero{margin-bottom:2.75rem;border-bottom:1px solid var(--line);padding-bottom:2rem}
.eyebrow{font-size:10px;font-weight:600;letter-spacing:.14em;text-transform:uppercase;
color:var(--muted);font-stretch:125%}
.wordmark{font-stretch:125%;font-weight:800;font-size:clamp(2.4rem,6.5vw,4.8rem);
line-height:.98;letter-spacing:-.01em;margin-top:8px}
.wordmark .hollow{color:transparent;-webkit-text-stroke:2px var(--ink)}
.hero .sub{margin-top:1rem;font-size:15px;color:var(--muted);max-width:38em}
.hero .meta-row{display:flex;gap:2rem;flex-wrap:wrap;margin-top:1.4rem;
font-size:12px;font-stretch:112%;letter-spacing:.06em;color:var(--muted)}
.hero .meta-row b{color:var(--ink);font-weight:700}
.pair-head{display:flex;align-items:baseline;gap:16px;margin:0 0 14px}
.pair-head .cn{font-size:19px;font-weight:800}
.pair-head .cnt{font-size:12px;font-stretch:112%;color:var(--muted)}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,300px),1fr));
gap:1px;background:var(--line);border:1px solid var(--line)}
@media(min-width:860px){.grid{grid-template-columns:repeat(3,1fr)}
a.card[href$="app=apxs"]{grid-column:span 2}}
.rover{grid-column:1/-1;background:var(--bg);padding:10px 16px 6px;font-size:10px;
font-weight:600;letter-spacing:.14em;text-transform:uppercase;color:var(--muted);font-stretch:125%}
a.card{display:flex;flex-direction:column;gap:4px;text-decoration:none;position:relative;
background:var(--panel);padding:18px 18px 16px;color:var(--ink);min-height:118px;
transition:background-color .15s cubic-bezier(.7,0,.2,1)}
a.card .no{position:absolute;top:14px;right:16px;font-stretch:125%;font-weight:800;
font-size:2rem;line-height:1;color:transparent;-webkit-text-stroke:1.4px var(--line);
font-variant-numeric:tabular-nums;transition:-webkit-text-stroke-color .15s}
a.card:hover{background:color-mix(in srgb,var(--accent) 8%,var(--panel))}
a.card:hover .no{-webkit-text-stroke-color:var(--accent)}
a.card .nm{font-size:16px;font-weight:800;line-height:1.35;padding-right:56px}
a.card .nm .tag{display:inline-block;font-size:10px;font-weight:700;color:var(--accent-ink);
background:var(--accent);padding:0 6px;margin-left:8px;vertical-align:2px;
clip-path:polygon(0 0,100% 0,100% calc(100% - 5px),calc(100% - 5px) 100%,0 100%)}
a.card .ds{font-size:12px;color:var(--muted)}
.foot{margin-top:26px;font-size:11px;color:var(--muted);display:flex;gap:2rem;flex-wrap:wrap}
.foot b{color:var(--ink)}
</style>
__SKIN__
</head><body>
<header class="cmdbar"><a class="brand" href="#"><span class="seal">M</span>MARS&nbsp;WORKBENCH</a>
<span class="ver">INSTRUMENTS │ v__VERSION__</span><button id="themeBtn">◐ 深浅色</button></header>
<main>
<section class="hero">
<div class="eyebrow">// DATA WORKBENCH</div>
<div class="wordmark">MARS<br><span class="hollow">WORKBENCH</span></div>
<p class="sub">两辆火星车、五台仪器的数据分析工作台。下载 PDS 公开档案，逐点逐发看谱、
元素线识别、矿物候选判读、质控报告、空间分布图、跨仪器对照。全本地处理，不上传。</p>
<div class="meta-row">
<span>ROVERS │ <b>2</b></span>
<span>INSTRUMENTS │ <b>5</b></span>
<span>发行前工程检查 │ <b>__CHECKS__</b></span>
<span>记录日期 │ <b>__CHECKS_DATE__</b></span>
<span>DATA │ <b>NASA PDS</b></span>
</div>
</section>
<div class="pair-head">
<div><div class="eyebrow">// INSTRUMENTS</div><div class="cn">选择仪器</div></div>
<span class="cnt">5 个仪器入口</span>
</div>
<div class="grid">
  <div class="rover">毅力号 MARS 2020 · Perseverance</div>
  <a class="card" id="sherlocCard" href="__SHERLOC_URL__">
    <span class="no">01</span>
    <span class="nm">SHERLOC 紫外拉曼光谱</span>
    <span class="ds">下载 · 预处理 · 候选判读 · 质控台 · 分布图</span>
  </a>
  <a class="card" href="__CONSOLE_URL__?app=supercam" target="_blank">
    <span class="no">02</span>
    <span class="nm">SuperCam LIBS 激光击穿光谱</span>
    <span class="ds">下载 · 看谱 · 元素线识别 · PIXL 对照 · CI 图像</span>
  </a>
  <a class="card" href="/pixl" id="pixlCard">
    <span class="no">03</span>
    <span class="nm">PIXL 微区X射线成分</span>
    <span class="ds">本地观测 · MCC扫描点 · 成分热图 · 官方误差 · 报告导出</span>
  </a>
  <div class="rover">好奇号 MSL · Curiosity</div>
  <a class="card" href="__CONSOLE_URL__?app=chemcam" target="_blank">
    <span class="no">04</span>
    <span class="nm">ChemCam LIBS 激光击穿光谱</span>
    <span class="ds">下载 · 看谱 · 元素线识别 · RMI 元素分布图 · APXS 对照</span>
  </a>
  <a class="card" href="__CONSOLE_URL__?app=apxs" target="_blank">
    <span class="no">05</span>
    <span class="nm">APXS α 粒子 X 射线光谱</span>
    <span class="ds">下载 · 看谱 · 元素峰识别 · 氧化物成分表</span>
  </a>
</div>
<div class="foot">
<span>DATA SOURCE │ <b>NASA PDS</b> 公开档案</span>
<span>PROCESSING │ <b>LOCAL</b> 不上传</span>
</div>
<p id="startStatus" role="status" style="margin:8px 20px;color:var(--muted)"></p>
<script>
/* SHERLOC 是另一个程序：点卡片先让启动器把它拉起来，再等它就绪 */
var starting=false, card=document.getElementById('sherlocCard');
card.addEventListener('click',function(e){
  e.preventDefault();
  if(starting) return;
  starting=true; card.setAttribute('aria-disabled','true');
  var status=document.getElementById('startStatus');
  status.textContent='正在启动 SHERLOC，请稍候…';
  fetch('/api/open/sherloc').then(function(r){return r.json();}).then(function(r){
    if(r.error) throw new Error(r.error);
    location.href=r.url;
  }).catch(function(e){
    starting=false; card.removeAttribute('aria-disabled');
    status.textContent='启动失败：'+e.message+'。可以重试。';
  });
});
/* 记住上次用的仪器：顶部给出「继续」直达条 */
fetch('/api/last').then(function(r){return r.json();}).then(function(d){
  if(!d.last) return;
  var map={sherloc:'__SHERLOC_URL__',supercam:'__CONSOLE_URL__?app=supercam',chemcam:'__CONSOLE_URL__?app=chemcam',
           apxs:'__CONSOLE_URL__?app=apxs',pixl:'/pixl'};
  if(!Object.prototype.hasOwnProperty.call(map,d.last)) return;
  var bar=document.createElement('a');
  bar.href=map[d.last]; bar.target='_blank';
  if(d.last==='sherloc'){
    bar.removeAttribute('target');
    bar.addEventListener('click',function(e){e.preventDefault();card.click();});
  }
  bar.className='card'; bar.style.cssText='margin:0 0 14px;border-left:4px solid var(--ef-trace)';
  bar.innerHTML='<span class="nm">继续上次：'+d.last.toUpperCase()
    +' 工作台 <span style="float:right">→</span></span>';
  var grid=document.querySelector('.grid');
  grid.parentNode.insertBefore(bar, grid);
});
function ping(){ fetch('/api/ping').catch(function(){}); }
ping(); setInterval(ping, 10000);
</script></body></html>"""



def build_page(launcher_port, console_port, sherloc_port=None):
    return (PAGE
            .replace("__SKIN__", "")
            .replace("__VERSION__", appenv.VERSION)
            .replace("__CHECKS__", validation_label()[0])
            .replace("__CHECKS_DATE__", validation_label()[1])
            .replace("__SHERLOC_URL__", "http://127.0.0.1:%d/" % (
                sherloc_port if sherloc_port is not None else configured_sherloc_port()))
            .replace("__CONSOLE_URL__", "http://127.0.0.1:%d/" % console_port))


def make_handler(state):
    page = build_page(state["launcher_port"], state["console_port"],
                      state.get("sherloc_port")).encode("utf-8")
    start_lock = threading.Lock()

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
            except Exception:                               # noqa: BLE001
                pass

        def do_GET(self):
            path = urllib.parse.urlparse(self.path).path
            if path.startswith('/open/'):
                name=path.removeprefix('/open/')
                if name not in ('sherloc','supercam','pixl','chemcam','apxs'):
                    return self._send(404,'text/plain; charset=utf-8','未知仪器'.encode())
                try:
                    if name=='sherloc':
                        with start_lock:
                            target=ensure_sherloc(state.get('sherloc_port',configured_sherloc_port()))
                            state['sherloc_started']=time.time()
                    elif name=='pixl':target='/pixl'
                    else:target='http://127.0.0.1:%d/?app=%s'%(state['console_port'],name)
                    self.send_response(303)
                    self.send_header('Location',target)
                    self.send_header('Content-Length','0')
                    self.send_header('Cache-Control','no-store')
                    self.end_headers()
                    return
                except Exception as exc:
                    import html
                    message='<!doctype html><html lang="zh"><meta charset="utf-8"><body><h1>仪器暂未打开</h1><p>'+html.escape(str(exc))+'</p><p><a href="/">回到仪器选择页</a></p></body></html>'
                    return self._send(503,'text/html; charset=utf-8',message.encode())
            if path in ("/", "/index.html"):
                return self._send(200, "text/html; charset=utf-8", page)
            if path == "/api/last":
                try:
                    from console import load_settings
                    last = load_settings().get("last_app", "")
                except Exception:                               # noqa: BLE001
                    last = ""
                return self._send(200, "application/json",
                                  json.dumps({"last": last}).encode("utf-8"))
            if path == "/api/ping":
                state["last_ping"] = time.time()
                return self._send(200, "application/json", b'{"ok":true}')
            if path == "/favicon.ico":
                return self._send(204, "image/x-icon", b"")
            if path == "/api/open/sherloc":
                try:
                    with start_lock:
                        url = ensure_sherloc(state.get("sherloc_port", configured_sherloc_port()))
                        state["sherloc_started"] = time.time()
                    return self._send(200, "application/json",
                                      json.dumps({"ok": True, "url": url}).encode())
                except Exception as exc:                     # noqa: BLE001
                    return self._send(503, "application/json",
                                      json.dumps({"error": str(exc)}).encode())
            return self._send(404, "text/plain; charset=utf-8", b"not found")

    return runtime_http.decorate(Handler,state)


def _watchdog(state, console_state):
    """选择器与工作台窗口都关掉才退出（SHERLOC 是独立进程，不受影响）。"""
    while True:
        time.sleep(5)
        pings = [p for p in (state.get("last_ping"),
                             (console_state or {}).get("last_ping"))
                 if p]+[st['last_ping'] for _,st in _sherloc_services if st.get('last_ping')]
        if pings and time.time() - max(pings) > PING_TIMEOUT:
            print("\n所有窗口都已关闭，自动退出。")
            try:
                state["shutdown"]()
            except Exception:                               # noqa: BLE001
                pass
            return


def main(argv=None):
    global _launcher_context
    ap=argparse.ArgumentParser(description='MarsWorkbench · 火星光谱工作台')
    ap.add_argument('--version',action='version',version=appenv.build_id())
    ap.add_argument('--port',type=int,default=8790)
    ap.add_argument('--console-port',type=int,default=8766)
    ap.add_argument('--sherloc-port',type=int,default=8765)
    ap.add_argument('--browser',action='store_true')
    ap.add_argument('--no-window',action='store_true')
    ap.add_argument('--no-watchdog',action='store_true')
    args=ap.parse_args(argv)
    if any(not 0<=p<=65535 for p in (args.port,args.console_port,args.sherloc_port)):ap.error('端口须在0到65535之间')
    runtime_http.initialize()
    import console
    common,common_state=console.start(args.console_port)
    threading.Thread(target=common.serve_forever,daemon=True).start()
    state={'last_ping':None,'launcher_port':args.port,'console_port':common.server_port,
           'sherloc_port':args.sherloc_port,'shutdown':lambda:None,'ensure_sherloc':ensure_sherloc,'console_state':common_state}
    server=runtime_http.bind(make_handler(state),args.port)
    state['launcher_port']=server.server_port;state['shutdown']=server.shutdown
    # Rebuild the page using the actual port if the requested launcher port was occupied.
    server.RequestHandlerClass=make_handler(state)
    url='http://127.0.0.1:%d/'%server.server_port
    state['launcher_url']=url;common_state['launcher_url']=url;_launcher_context=state
    print('MarsWorkbench '+appenv.build_id()+' · 选择器：'+url,flush=True)
    print('通用工作台：http://127.0.0.1:%d/'%common.server_port,flush=True)
    if not args.no_watchdog:threading.Thread(target=_watchdog,args=(state,common_state),daemon=True).start()
    if not args.no_window:open_app_window(url,not args.browser,size=(1440,900))
    try:server.serve_forever()
    except KeyboardInterrupt:pass
    finally:
        common.shutdown();common.server_close()
        for child,_ in _sherloc_services:child.shutdown();child.server_close()
        server.server_close()
    return 0

if __name__=='__main__':
    try:raise SystemExit(main())
    except Exception as exc:
        print(appenv.crash(exc));raise SystemExit(1)
