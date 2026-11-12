"""Workflow contracts: user-visible export scope, report coverage, history and shared preferences."""

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout
import csv
import json
import logging
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time
import urllib.error
import urllib.request
from http.server import BaseHTTPRequestHandler
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import appenv,pixl_view,product_ui,runtime_http
from instruments import history,pixl_browser as b,pixl_workflow as w

N=BAD=0
def check(ok,name):
    global N,BAD
    N+=1;BAD+=not ok;print(('OK ' if ok else 'FAIL ')+name,flush=True)
def fails(fn):
    try:fn();return False
    except (ValueError,PermissionError,FileNotFoundError):return True

def main():
    with tempfile.TemporaryDirectory() as name,patch.dict(os.environ,{'MARSWORKBENCH_DATA_HOME':str(Path(name)/'state'),'MARSWORKBENCH_HEADLESS':'1'}):
        temp=Path(name);folder=temp/'data';shutil.copytree(appenv.resources()/'pixl-demo',folder)
        st=pixl_view.state_for({});root,obs,_=b.scan(folder);data=b.load(root,obs[0])
        history.update(appenv.settings('pixl'),lambda c:dict(c,download_dest='keep'))
        for i in range(6):
            new=temp/('data'+str(i));shutil.copytree(folder,new);pixl_view.load_folder(st,str(new))
        config=history.read(appenv.settings('pixl'))
        check(len(history.rows(config))==5 and history.rows(config)[0]['path']==str(new) and config['download_dest']=='keep','换目录复用共享最近五次记录且保留其他偏好')
        before=history.read(appenv.settings('pixl'));bad=fails(lambda:pixl_view.load_folder(st,str(temp/'missing')))
        check(bad and history.read(appenv.settings('pixl'))==before,'失败目录不污染最近记录或覆盖已打开数据')
        out,csvname=w.export_single(data,'csv',str(temp/'csv'),pmc=6)
        with (out/csvname).open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f))
        check(len(rows)==1 and rows[0]['PMC']=='6' and rows[0]['SiO2_wt%']=='41.6797' and rows[0]['SiO2_wt%_err']=='2.1058' and not (out/'index.html').exists(),'当前点CSV仅一行，保持数值/误差且不捆绑报告')
        out,csvname=w.export_single(data,'csv',str(temp/'csv'),scope='all')
        with (out/csvname).open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f))
        check(len(rows)==2809 and (out/csvname).read_bytes().startswith(b'\xef\xbb\xbf'),'全部点CSV仅当前观测的全部RQC点，BOM保留')
        out,png=w.export_single(data,'png',str(temp/'png'),'SO3')
        check((out/png).read_bytes().startswith(b'\x89PNG') and sorted(p.name for p in out.iterdir())==['map_SO3.png','sources.json'],'PNG独立导出当前指标，不捆绑CSV/报告')
        check(fails(lambda:w.export_single(data,'csv',str(temp/'csv'),pmc=5)) and fails(lambda:w.export_single(data,'png',str(temp/'png'),'SO3',2)) and fails(lambda:w.export_single(data,'csv',str(temp/'csv'),scope='invalid')),'无成分PMC、缺图PNG与无效范围均明确拒绝')
        with obs[0]['path'].open(encoding='utf-8-sig',newline='') as f:direct=list(csv.DictReader(f))
        valid=sum(math.isfinite(float(row[e+'_wt%'])) and float(row[e+'_wt%'])>=0 for row in direct for e in data['elements'])
        summary=w.check_summary(data)
        check(summary['valid_values']==valid and summary['chemical_points']==len(direct) and summary['matched_points']==2809,'归档数据检查统计与独立CSV逐格复算相符')
        second=dict(obs[0],id='second');report,entries=w.build_report(root,[obs[0],second],str(temp/'qc'),obs[0]['id'],None,0)
        nested=folder/'下载';shutil.copytree(appenv.resources()/'pixl-demo',nested);rqc=next(nested.glob('*rqc*.csv'));rqc.write_bytes(rqc.read_bytes()+b'\n')
        check(len(entries)==2 and all(r.get('report') for r in entries) and all((report/r['report']).exists() for r in entries) and '归档数据检查' in (report/entries[0]['report']).read_text(encoding='utf-8') and len(b.scan(folder)[1])==1 and len(b.scan(nested)[1])==1,'全目录报告不受筛选影响；父目录不混入下载副本，下载目录可单独打开')
        broken=dict(second,path=temp/'absent.csv');report,entries=w.build_report(root,[broken,obs[0]],str(temp/'qc'))
        check('error' in entries[0] and entries[1].get('report') and str(temp) not in (report/'index.html').read_text(encoding='utf-8'),'报告单观测失败保留错误并继续，分享页无本机路径')
        stop=[False];report,entries=w.build_report(root,[obs[0],second],str(temp/'qc'),cancel=lambda:stop[0],progress=lambda *a:stop.__setitem__(0,True))
        check(len(entries)==1 and '已取消' in (report/'index.html').read_text(encoding='utf-8'),'报告取消保留已完成的逐观测产物并标明部分结果')
        node=shutil.which('node');scripts=history.NAV_JS+"\nconsole.log(JSON.stringify([navigationCommand('PageDown','',0,0,3,9,0),navigationCommand('End','',0,0,3,9,0),navigationCommand('1','100',20,0,1,2816,0,4),navigationCommand('1','100',20,0,1,2816,0)]));"
        result=json.loads(subprocess.check_output([node,'-e',scripts],text=True))
        check(result[0]['index']==1 and result[1]['index']==9 and result[2]['index']==1001 and result[3]['index']==1,'共享快捷键保留原台三位语义，PIXL支持四位实际PMC')
        scripts=history.NAV_JS+"\nconsole.log(JSON.stringify([retainedPointId(1001,[5,8,1001],6),retainedPointId(9,[5,8,1001],8),retainedPointId(9,[],8)]));"
        check(json.loads(subprocess.check_output([node,'-e',scripts],text=True))==[1001,8,None],'切帧保留实际PMC，缺失时只落到真实存在的点')
        page=product_ui.inject(pixl_view.PAGE.encode(),'pixl','/').decode();scripts='\n'.join(re.findall(r'<script[^>]*>(.*?)</script>',page,re.S))
        checked=subprocess.run([node,'--check'],input=scripts,text=True,capture_output=True)
        check(checked.returncode==0 and all('id="'+i+'"' in page for i in ('folderPath','browseBtn','dlOpenBtn','scan','pt','ptPrev','ptNext','exportBtn','csvBtn','qcBtn','themeBtn','browser','recentFolders','dlCard','dlSol')),'实际页面使用通用操作入口与双向翻点，全部脚本语法有效')
        class Handler(BaseHTTPRequestHandler):
            def log_message(self,*a):pass
            def _send(self,code,ctype,body):self.send_response(code);self.send_header('Content-Type',ctype);self.end_headers();self.wfile.write(body)
            def do_GET(self):self._send(404,'text/plain',b'not found')
        states=[{},{}];servers=[runtime_http.bind(runtime_http.decorate(Handler,s),0) for s in states]
        threads=[threading.Thread(target=s.serve_forever,daemon=True) for s in servers]
        for t in threads:t.start()
        def request(path,body=None,server=0):
            url='http://127.0.0.1:%d'%servers[server].server_port+path
            req=urllib.request.Request(url,json.dumps(body).encode() if body is not None else None,{'Content-Type':'application/json'})
            with urllib.request.urlopen(req,timeout=30) as r:
                raw=r.read();return json.loads(raw) if 'json' in r.headers.get('Content-Type','') else raw
        try:
            history.update(appenv.settings('appearance'),lambda c:dict(c,extra='preserve'))
            request('/api/appearance',{'theme':'dark'});dark=request('/api/appearance',server=1)
            request('/api/appearance',{'theme':'light'},server=1);light=request('/api/appearance')
            try:request('/api/appearance',{'theme':'bad'});invalid=False
            except urllib.error.HTTPError as e:invalid=e.code==400
            check(dark['theme']=='dark' and light['theme']=='light' and invalid and history.read(appenv.settings('appearance'))['extra']=='preserve','两工作台共享主题偏好，坏值拒绝且不破坏其他设置')
            info=request('/api/samples');installed=request('/api/samples/install',{'instrument':'pixl','destination':str(temp/'examples')});example=Path(installed['folders']['pixl']);(example/'user-note.txt').write_text('keep')
            again=request('/api/samples/install',{'instrument':'pixl','destination':str(temp/'examples')});confirmed=request('/api/samples/install',{'instrument':'pixl','destination':str(temp/'examples'),'confirm':True})
            check('pixl' in info['instruments'] and again['requires_confirmation'] and not confirmed['requires_confirmation'] and (example/'user-note.txt').read_text()=='keep','PIXL复用统一示例选择/目录/覆盖询问，保留用户其他文件')
            loaded=request('/api/pixl/load',{'path':str(folder)});oid=loaded['observations'][0]['id']
            exported=request('/api/pixl/export-selected',{'kind':'csv','obs':oid,'pmc':6,'scope':'current'});assert b'41.6797' in request(exported['url'])
            job=request('/api/pixl/report/start',{'obs':oid});end=time.monotonic()+20
            while time.monotonic()<end:
                status=request('/api/pixl/report/status')
                if not status['running']:break
                time.sleep(.05)
            body=request(status['url']);sub=next(r['report'] for r in status['rows'] if r.get('report'));sub_body=request(status['url'].rsplit('/',1)[0]+'/'+sub)
            try:request('/api/pixl/openfolder',{'kind':'export','token':'unknown'});denied=False
            except urllib.error.HTTPError as e:denied=e.code==400
            from sherloc import desktop
            with patch.dict(os.environ,{'MARSWORKBENCH_HEADLESS':'0'}),patch.object(desktop,'open_app_window') as opened:
                request('/api/pixl/openreport',{})
            check('数据检查与成分报告' in body.decode() and b'41.6797' in sub_body and status['done']==1 and denied and opened.call_count==1 and opened.call_args[0][0].endswith('index.html'),'实际HTTP独立CSV和后台报告完成，自动打开已生成报告且只允许已生成目录')
        finally:
            for s in servers:s.shutdown();s.server_close()
            for t in threads:t.join()
            logger=logging.getLogger('marsworkbench')
            for h in list(logger.handlers):h.close();logger.removeHandler(h)
    print('共 %d 项断言，%d 项失败'%(N,BAD));return int(BAD>0)

if __name__=='__main__':raise SystemExit(main())
