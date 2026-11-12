"""Theme integration contracts; data colors, DOM controls and numerical kernels are preserved."""

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout
import base64
import json
import logging
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import urllib.error
import urllib.request
from html.parser import HTMLParser
from http.server import BaseHTTPRequestHandler
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import appenv,console,launch,pixl_view,product_theme as theme,product_ui,view,runtime_http
from instruments import REGISTRY,pixl_browser as pixl

N=BAD=0
def check(ok,name):
    global N,BAD
    N+=1;BAD+=not ok;print(('OK ' if ok else 'FAIL ')+name,flush=True)

def chrome_checks(pages,names):
    class Surface(HTMLParser):
        def __init__(self):super().__init__();self.ids=[];self.links=[]
        def handle_starttag(self,tag,attrs):
            values=dict(attrs)
            if 'id' in values:self.ids.append(values['id'])
            if tag=='a' and 'data-instrument' in values:self.links.append(values)
    parsed=[]
    for page in pages:
        parser=Surface();parser.feed(page);parsed.append(parser)
    check(all(len(p.links)==5 and {a['data-instrument'] for a in p.links}=={i[0] for i in theme.INSTRUMENTS} and [a['data-instrument'] for a in p.links if a.get('aria-current')=='page']==([] if name=='launcher' else [name]) for p,name in zip(parsed,names)),'六页真实仪器导航完整，当前项唯一，启动器不伪装成已选仪器')
    check(all(len(p.ids)==len(set(p.ids)) and p.ids.count('mwRail')==1 and p.ids.count('themeBtn')==1 for p in parsed) and all(id in parsed[-1].ids for id in ('frameComposition','framePreviousPage','frameNextPage','scanRegion')),'导航和既有控件无重复ID，PIXL图内成分表与全图/扫描区入口有效')
    report=theme.report_style('<html><body><h1>报告</h1><svg><text>41.6797</text></svg></body></html>')
    check(product_ui.inject(report.encode(),'sherloc','/').decode()==report and 'id="mwRail"' not in report and '41.6797' in report,'报告自带阅读主题，不附加工作台导航或依赖本机接口')
    class NoRedirect(urllib.request.HTTPRedirectHandler):
        def redirect_request(self,*args,**kwargs):return None
    class Bare(BaseHTTPRequestHandler):
        def log_message(self,*args):pass
        def _send(self,code,ctype,body):self.send_response(code);self.send_header('Content-Type',ctype);self.end_headers();self.wfile.write(body)
        def do_GET(self):self._send(404,'text/plain',b'not found')
    with tempfile.TemporaryDirectory() as temporary,patch.dict(os.environ,{'MARSWORKBENCH_DATA_HOME':temporary}):
        from instruments import history
        state={'launcher_port':8970,'console_port':8971,'sherloc_port':0}
        servers=[runtime_http.bind(launch.make_handler(state),0),runtime_http.bind(runtime_http.decorate(Bare,{}),0)]
        for server in servers:threading.Thread(target=server.serve_forever,daemon=True).start()
        opener=urllib.request.build_opener(NoRedirect)
        def request(path,body=None,server=0):
            req=urllib.request.Request('http://127.0.0.1:%d'%servers[server].server_port+path,json.dumps(body).encode() if body is not None else None,{'Content-Type':'application/json'})
            try:
                with opener.open(req,timeout=10) as response:return response.status,dict(response.headers),response.read()
            except urllib.error.HTTPError as error:return error.code,dict(error.headers),error.read()
        try:
            with patch.object(launch,'ensure_sherloc',return_value='http://127.0.0.1:8982/') as ensure:
                expected={'sherloc':'http://127.0.0.1:8982/','supercam':'http://127.0.0.1:8971/?app=supercam','pixl':'/pixl','chemcam':'http://127.0.0.1:8971/?app=chemcam','apxs':'http://127.0.0.1:8971/?app=apxs'}
                answers={key:request('/open/'+key) for key in expected}
                check(all(answers[k][0]==303 and answers[k][1].get('Location')==target for k,target in expected.items()) and ensure.call_count==1,'真实HTTP五仪器导航指向各自工作台，SHERLOC采用实际分配端口')
                invalid=[request('/open/'+name) for name in ('unknown','../outside','https://example.invalid')]
                check(all(code==404 and 'Location' not in headers for code,headers,_ in invalid) and ensure.call_count==1,'导航仅接受已知仪器，不接受路径或外部地址作为跳转目标')
            path=appenv.settings('appearance');history.update(path,lambda c:dict(c,extra='preserve'))
            request('/api/appearance',{'theme':'dark','tone':'soft'})
            first=json.loads(request('/api/appearance',server=1)[2])
            request('/api/appearance',{'tone':'bright'},server=1)
            invalid=request('/api/appearance',{'theme':'light','tone':'unknown'})[0]
            final=json.loads(request('/api/appearance')[2])
            check(first=={'theme':'dark','tone':'soft'} and final=={'theme':'dark','tone':'bright'} and invalid==400 and history.read(path)['extra']=='preserve','深浅色及黄色强度跨端口同步，坏值整批拒绝，已有配置保留')
        finally:
            for server in servers:server.shutdown();server.server_close()
            log=logging.getLogger('marsworkbench')
            for handler in list(log.handlers):
                target=getattr(handler,'baseFilename',None)
                if target and Path(target).is_relative_to(Path(temporary)):
                    handler.close();log.removeHandler(handler)
def main():
    original=[launch.build_page(8970,8971),view.build_page([]),*[console.build_page(b) for b in REGISTRY.values()],pixl_view.PAGE]
    names=['launcher','sherloc',*REGISTRY,'pixl'];pages=[product_ui.inject(p.encode(),name,'/').decode() for p,name in zip(original,names)]
    check(all(p.count('data-mw-theme="endfield-yellow-paper"')==1 and p.count('data:font/woff2;base64,')==1 for p in pages),'启动器和五台各注入一份皮肤与离线字体，无重复CSS')
    check(all(product_ui.inject(p.encode(),name,'/').decode()==p for p,name in zip(pages,names)),'重复路由装饰保持幂等，不重绑控件')
    check(all(all('data-ef-'+k+'="'+v+'"' in p for k,v in theme.VARIANT.items()) for p in pages),'六页共享单一产品变体配置，全部标记覆盖')
    controls=('browseBtn','dlOpenBtn','scan','pt','qcBtn')
    check(all(all('id="'+i+'"' in p for i in controls) for p in pages[1:]) and all('id="'+i+'"' in pages[0] for i in ('sherlocCard','pixlCard')),'原工作台DOM操作入口及仪器卡全部保留')
    check(all('mw-kicker' in p and '// '+theme.LABELS[name] in p for p,name in zip(pages[1:],names[1:])),'五台标题使用真实仪器标签，无占位技术文字')
    node=shutil.which('node');syntax=[]
    for p in pages:
        scripts='\n'.join(re.findall(r'<script[^>]*>(.*?)</script>',p,re.S));syntax.append(subprocess.run([node,'--check'],input=scripts,text=True,capture_output=True).returncode==0)
    check(all(syntax),'六页实际操作脚本保持有效，皮肤不改动作语义')
    raw=base64.b64decode(__import__('ui_font').FONT_B64);license=appenv.resources()/'licenses/Saira-OFL.txt'
    check(raw.startswith(b'wOF2') and license.is_file() and 'SIL OPEN FONT LICENSE' in license.read_text(encoding='utf-8'),'字体可离线解码且完整OFL随资源分发')
    svg='<svg><path stroke="#c23135"/><path stroke="#287cc0"/><rect fill="#31ad62"/><text>2809</text></svg>'
    html='<html><head></head><body><h1>数据</h1>'+svg+'</body></html>'
    report=theme.report_style(html)
    implicit=theme.report_style('<!doctype html><html lang="zh"><meta charset="utf-8"><h1>特征与分型</h1>'+svg+'</html>')
    check(svg in report and svg in theme.decorate(html,'sherloc') and theme.report_style(report)==report and implicit.count('data-mw-theme=')==1 and svg in implicit and theme.report_style(implicit)==implicit,'科学曲线/色标原样保留，质控及隐式body分析报告均自含主题且幂等')
    with tempfile.TemporaryDirectory() as name,patch.dict(os.environ,{'MARSWORKBENCH_DATA_HOME':str(Path(name)/'state')}):
        root=Path(name);resources=root/'resources';resources.mkdir()
        snapshot={'scope':'developer','failures':0,'assertions':1234,'finished_utc':'2026-11-11T00:00:00Z','results':[{'ok':True}]}
        (resources/'_last_run.json').write_text(json.dumps(snapshot),encoding='utf-8')
        with patch.object(appenv,'resources',return_value=resources):rendered=launch.build_page(8970,8971)
        check('1234' in rendered and '2026-11-11' in rendered and 'v'+appenv.VERSION in rendered and 'v11<' not in rendered and '@font-face' not in re.sub(r'<style.*?</style>','',rendered,flags=re.S),'首页版本/记录读取真实来源，CSS不会露成正文')
        with patch.object(appenv,'resources',return_value=root/'missing'):missing=theme.validation_label()
        check(missing==('待运行','无发行记录'),'没有发行记录时不宣称检查通过')
        folder,obs,_=pixl.scan(appenv.resources()/'pixl-demo');data=pixl.load(folder,obs[0]);out=root/'report';pixl.export(data,out,'SiO2')
        body=(out/'index.html').read_text(encoding='utf-8')
        check('data-mw-surface="report"' in body and body.count('data-mw-theme=')==1 and '41.6797' in body and 'composition.csv' in body and not re.search(r'(?:src|href)="https?://',body),'可携带报告自含主题/字体、数值和本地链接，无网络资源依赖')
    check('prefers-reduced-motion:reduce' in theme.SKIN_CSS and 'animation:none!important' in theme.SKIN_CSS and 'focus-visible' in theme.SKIN_CSS and all('@media(max-width:'+str(w)+'px)' in theme.SKIN_CSS for w in (600,900)),'低动效、键盘焦点和窄窗口布局规则覆盖')
    chrome_checks(pages,names)
    print('共 %d 项断言，%d 项失败'%(N,BAD));return int(BAD>0)
if __name__=='__main__':raise SystemExit(main())
