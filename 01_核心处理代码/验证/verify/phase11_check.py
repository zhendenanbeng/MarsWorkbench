"""Productization contract guards; science result guards remain in the original suites."""

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout
import builtins
import json
import logging
from pathlib import Path
import os
import re
import subprocess
import sys
import tempfile
import threading
import urllib.error
import urllib.request
from http.server import BaseHTTPRequestHandler
from unittest.mock import patch

APP=_layout.project_layout.APP
PROJECT=_layout.project_layout.ROOT
sys.path.insert(0,str(APP))
import appenv
import runtime_http
import app_network

TOTAL=FAIL=0
def check(ok,message):
    global TOTAL,FAIL
    TOTAL+=1;FAIL+=not ok
    print(('OK ' if ok else 'FAIL ')+message,flush=True)

def get(server,path):
    with urllib.request.urlopen('http://127.0.0.1:%d'%server.server_port+path,timeout=20) as response:
        return json.load(response)

def serve(server):
    thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start();return thread

def stop(server,thread):
    server.shutdown();server.server_close();thread.join()


def checks_a(temp):
    from instruments import history
    path=Path(appenv.settings('console'));history.update(path,lambda _:dict(theme='dark'))
    check(path.is_relative_to(Path(os.environ['MARSWORKBENCH_DATA_HOME'])) and history.read(path)['theme']=='dark' and not path.is_relative_to(APP),'A1 设置在用户数据目录读写，软件目录外置')
    legacy=temp/'legacy';legacy.mkdir();old=legacy/'.console_settings.json';old.write_text('{"last_app":"apxs"}',encoding='utf-8')
    with patch.object(appenv,'app_dir',return_value=legacy),patch.dict(os.environ,{'MARSWORKBENCH_DATA_HOME':str(temp/'migrated')}):
        target=Path(appenv.settings('console'));first=target.read_bytes()
        with patch.object(appenv,'app_dir',return_value=temp/'moved'):
            moved=Path(appenv.settings('console')).read_bytes()
    check(first==moved and json.loads(first)['last_app']=='apxs' and old.read_text(encoding='utf-8')=='{"last_app":"apxs"}','A2 旧配置复制迁移且原件保留，软件搬家不丢设置')
    portable=temp/'portable';portable.mkdir();(portable/'portable.flag').touch()
    with patch.dict(os.environ,{},clear=False),patch.object(appenv,'app_dir',return_value=portable):
        saved=os.environ.pop('MARSWORKBENCH_DATA_HOME')
        try:home=appenv.data_home()
        finally:os.environ['MARSWORKBENCH_DATA_HOME']=saved
    check(home==portable/'userdata' and all((home/n).is_dir() for n in ['settings','cache','logs','reports_selfcheck']),'A3 便携标志切换及四类状态目录齐全')
    from instruments import report,REGISTRY
    fixture=_layout.project_layout.VERIFY/'fixtures/apxs';out=temp/'report'
    report.build_report(REGISTRY['apxs'],REGISTRY['apxs'].scan_folder(str(fixture)),str(out),data_dir=str(fixture))
    page=(out/'index.html').read_text(encoding='utf-8')
    synthetic=appenv.clean_report('<a href="../../04_使用文档/a.md">说明</a> C:\\Users\\author\\Desktop\\tool')
    check(not re.search(r'[A-Za-z]:[\\/]|href=[\"\'][^\"\']*04_使用文档',page+synthetic) and str(Path.home()) not in page,'A4 真实报告与开发路径正例均不残留开发机绝对地址/文档外链')
    import tomllib
    meta=tomllib.loads((PROJECT/'pyproject.toml').read_text(encoding='utf-8'))
    result=subprocess.run([sys.executable,str(APP/'launch.py'),'--version'],capture_output=True,text=True)
    check(meta['project']['dynamic']==['version'] and result.returncode==0 and result.stdout.strip().startswith(appenv.VERSION+'+'),'A5 pyproject、CLI和构建号使用单一版本源')
    class Broken(BaseHTTPRequestHandler):
        def log_message(self,*a):pass
        def _send(self,code,ctype,body):
            self.send_response(code);self.send_header('Content-Type',ctype);self.end_headers();self.wfile.write(body)
        def do_GET(self):raise RuntimeError('test stack provenance')
    server=runtime_http.bind(runtime_http.decorate(Broken,{}),0);thread=serve(server)
    try:
        try:get(server,'/fault')
        except urllib.error.HTTPError as exc:error=json.load(exc)
    finally:stop(server,thread)
    try:raise ValueError('startup proof')
    except ValueError as exc:human=appenv.crash(exc)
    log=(appenv.location('logs','app.log')).read_text(encoding='utf-8');crash=appenv.location('logs','crash.txt').read_text(encoding='utf-8')
    check(error['request_id'] in log and 'Traceback' in log and 'test stack provenance' not in error['error'] and 'startup proof' in crash and '日志' in human,'A6 请求编号对应完整栈，启动失败有crash文件与友好说明')
    occupied=runtime_http.bind(Broken,0);fallback=runtime_http.bind(Broken,occupied.server_port)
    check(occupied.server_port<fallback.server_port<occupied.server_port+20,'A7 占用端口自动递增，不关闭旧服务')
    check(occupied.server_address[0]==fallback.server_address[0]=='127.0.0.1','A8 所有绑定均只开放本机回环')
    occupied.server_close();fallback.server_close()
    import console,view
    source_writes=[];original=builtins.open
    def read_only(file,mode='r',*args,**kwargs):
        if not isinstance(file,int) and any(c in mode for c in 'wax+'):
            path=Path(file).resolve()
            if path.is_relative_to(APP):source_writes.append(str(path));raise PermissionError('software directory read only')
        return original(file,mode,*args,**kwargs)
    with patch('builtins.open',read_only):
        common,state=console.start(0)
        for name,b in REGISTRY.items():
            folder=_layout.project_layout.VERIFY/'fixtures'/name;state['app'][name]={'folder':str(folder),'observations':b.scan_folder(str(folder))}
        common_thread=serve(common)
        graph=[]
        try:
            for name in REGISTRY:
                payload=get(common,'/api/plot?app='+name+'&scan=0&point=0&shot=-1')
                graph.append(bool(payload.get('series')))
        finally:stop(common,common_thread)
        sh,shstate=view.start(0,str(appenv.resources()/'sample-data/sherloc'))
        shthread=serve(sh)
        try:graph.append(bool(get(sh,'/api/plot?scan=0&point=0&mode=all').get('series')))
        finally:stop(sh,shthread)
    check(len(graph)==4 and all(graph) and not source_writes,'A9 四仪器真实浏览API在禁止软件目录写入时可用')
    oldroot=temp/'old-cache';(oldroot/'.cache/edr_labels').mkdir(parents=True)
    label=oldroot/'.cache/edr_labels/a.lbl';label.write_bytes(b'original cache')
    with patch.object(appenv,'SOURCE',oldroot),patch.dict(os.environ,{'MARSWORKBENCH_DATA_HOME':str(temp/'cache-user')}):
        read=appenv.cache_read('edr_labels/a.lbl');migrated=appenv.copy_cached('edr_labels')/'a.lbl'
        compatible=read==label and migrated.read_bytes()==label.read_bytes()==b'original cache'
    with patch.dict(os.environ,{'MARSWORKBENCH_OFFLINE':'1'}):offline=app_network.probe()
    check(compatible and not offline['ok'] and '健康检查' in offline['message'],'A10 旧缓存只读兼容，离线有统一健康检查提示')


def checks_b(temp):
    import hashlib
    import shutil
    import samples,health,product_ui
    from sherloc import watson
    data=samples.manifest();rows=data['files']
    intact=all((appenv.resources()/r['path']).stat().st_size==r['size'] and hashlib.sha256((appenv.resources()/r['path']).read_bytes()).hexdigest()==r['sha256'] for r in rows)
    check(intact and set(data['instruments'])=={'sherloc','chemcam','supercam','apxs'} and data['total_bytes']==sum(r['size'] for r in rows),'B1 四仪器示例与缓存种子清单、字节数、SHA完整')
    with patch.object(watson,'CACHE',str(temp/'new-watson')):
        installed=samples.install(temp/'examples')
        user=temp/'examples/chemcam/user-note.txt';user.write_text('keep',encoding='utf-8')
        again=samples.install(temp/'examples')
        confirmed=samples.install(temp/'examples',True)
    check(not installed['requires_confirmation'] and again['requires_confirmation'] and not confirmed['requires_confirmation'] and user.read_text(encoding='utf-8')=='keep','B2 重复加载必须确认，只替换同名示例并保留用户文件')
    pending=samples.welcome('pending');skipped=samples.welcome('skipped');complete=samples.welcome('complete')
    check(pending['show'] and not skipped['show'] and complete['state']=='complete' and not samples.welcome()['show'],'B3 首跑向导可跳过、可完成、状态持久化')
    import console,view
    from instruments import REGISTRY
    node=shutil.which('node');assert node,'前端检查需要Node.js'
    pages=[console.build_page(b) for b in REGISTRY.values()]+[view.build_page([])]
    syntax=[]
    for page in pages:
        injected=product_ui.inject(page.encode(),'sherloc').decode()
        scripts='\n'.join(re.findall(r'<script[^>]*>(.*?)</script>',injected,re.S))
        syntax.append(subprocess.run([node,'--check'],input=scripts,text=True,capture_output=True).returncode==0 and 'id="mwHelp"' in injected)
    with patch.object(urllib.request,'urlopen',side_effect=AssertionError('帮助不得联网')):
        import manuals
        help_parts=manuals.help_sections()
        sections=[help_parts[n] for n in ['quickstart','faq','scope','about']]
    check(all(syntax) and all('<h2>' in s for s in sections) and '快速浏览' in sections[2] and ('不等于定量' in sections[2] or '≠定量' in sections[2]),'B4 四台离线帮助来自保留手册、四节齐全，实际页面脚本语法正确')
    with patch.object(health.app_network,'probe',return_value={'ok':True,'message':'fixture online'}):
        good=health.diagnose(8766)
        with patch.object(appenv,'data_home',side_effect=PermissionError('fixture readonly')):bad=health.diagnose(8766)
    with patch.object(health.app_network,'probe',return_value={'ok':False,'message':'fixture offline'}):offline=health.diagnose(8766)
    check(good['writable'] and good['bind']=='127.0.0.1' and 'numpy' in good['components'] and bad['status']=='error' and offline['pds']['ok'] is False and 'fixture offline' in offline['text'],'B5 健康检查区分必需环境错误、离线警告并生成可复制文本')
    sources=(appenv.resources()/'sample-data/SOURCES.md').read_text(encoding='utf-8')
    check(all(r['path'] in sources and r['url'] in sources and r['bundle'] in sources and r['sha256'] in sources for r in rows),'B6 SOURCES逐文件覆盖来源、bundle及SHA，不遗漏缓存种子')


def checks_c_smoke(temp):
    sys.path.insert(0,str(_layout.project_layout.INSTALL/'构建工具'))
    import smoke_fresh
    result=smoke_fresh.run()
    for index,name in enumerate(('冷启动、四仪器和PNG','SHERLOC报告','ChemCam地图及报告链接')):
        rows=result.get('checks',[])
        check(len(rows)>index and rows[index]['ok'] and (index<2 or result['ok']),'C%d 源码隔离端到端：%s；%s'%(index+1,name,result.get('error','通过')))
    return result


def checks_c_legal(temp):
    import hashlib
    import health,product_api
    base=appenv.resources()/'licenses'
    manifest=json.loads((base/'manifest.json').read_text(encoding='utf-8'))
    check(all((base/name).is_file() for name in ('LICENSE','THIRD_PARTY_NOTICES.md','DATA_SOURCES.md','CPython-LICENSE.txt')) and all(hashlib.sha256((base/r['file']).read_bytes()).hexdigest()==r['sha256'] for r in manifest['files']) and json.loads((_layout.project_layout.INSTALL/'版本信息.json').read_text(encoding='utf-8'))['version']==appenv.VERSION,'C4 许可原文与清单SHA一致，版本与安装元数据对应')
    with patch.dict(sys.modules,{'cv2':None,'astropy':None}):
        degraded=health.diagnose(network=False)
    spec=(_layout.project_layout.INSTALL/'构建工具/win.spec').read_text(encoding='utf-8')
    check(degraded['status']=='warning' and degraded['components']['numpy']['available'] and not degraded['components']['astropy']['available'] and '全功能版' in degraded['text'] and 'collect_submodules' in spec and 'instruments' in spec and 'sherloc' in spec,'C5 标准版缺可选库仅警告，打包显式覆盖惰性模块')
    class Handler:
        path='/api/version/check';command='GET'
        def _send(self,code,ctype,body):self.response=(code,json.loads(body))
    handler=Handler()
    with patch.object(urllib.request,'urlopen',side_effect=AssertionError('更新未配置不许联网')):
        handled=product_api.route(handler,{})
    check(handled and handler.response[0]==200 and handler.response[1]['enabled'] is False and handler.response[1]['configured'] is False,'C6 更新源未配置时默认关闭且零联网')


def main():
    with tempfile.TemporaryDirectory(prefix='mars-product-') as name:
        temp=Path(name)
        with patch.dict(os.environ,{'MARSWORKBENCH_DATA_HOME':str(temp/'runtime')}):
            try:
                checks_a(temp)
                checks_b(temp)
                smoke=checks_c_smoke(temp)
                checks_c_legal(temp)
                reports=list((Path(smoke['root'])/'profile/Documents/MarsWorkbench示例').rglob('*.html'))
                check(len(reports)>=4 and all(not re.search(r'(?i)(?:(?<![A-Za-z0-9])[A-Za-z]:[\\/]|file://|href=[\"\'][^\"\']*04_使用文档)',p.read_text(encoding='utf-8')) for p in reports),'E1 实际端到端全部报告不含盘符、file链接或开发文档外链')
            finally:
                logger=logging.getLogger('marsworkbench')
                for handler in list(logger.handlers):handler.close();logger.removeHandler(handler)
    print('共 %d 项断言，%d 项失败'%(TOTAL,FAIL))
    return int(FAIL>0)

if __name__=='__main__':raise SystemExit(main())
