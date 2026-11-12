"""Independent PIXL local-view contracts, real archive sample and negative cases."""

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout
import csv
import hashlib
import io
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
import urllib.request
import urllib.error
from http.server import BaseHTTPRequestHandler
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import appenv,launch,pixl_view,runtime_http
from instruments import pixl_browser as p
from PIL import Image
import numpy as np

N=BAD=0
def check(ok,name):
    global N,BAD
    N+=1;BAD+=not ok;print(('OK ' if ok else 'FAIL ')+name,flush=True)
def fails(fn):
    try:fn();return False
    except (ValueError,FileNotFoundError,PermissionError):return True

def main():
    demo=appenv.resources()/'pixl-demo'
    manifest=json.loads((demo/'manifest.json').read_text(encoding='utf-8'))
    check(len(manifest['files'])==6 and all(hashlib.sha256((demo/r['file']).read_bytes()).hexdigest()==r['sha256'] for r in manifest['files']),'示例六个原件SHA锚定')
    with tempfile.TemporaryDirectory() as name,patch.dict(os.environ,{'MARSWORKBENCH_DATA_HOME':str(Path(name)/'state')}):
        temp=Path(name);folder=temp/'local';shutil.copytree(demo,folder)
        root,obs,w=p.scan(folder);data=p.load(root,obs[0]);payload=p.payload(data)
        check(len(obs)==1 and obs[0]['sol']==187 and len(data['chem'])==2809 and len(data['raw']['pts'])==2811 and payload['matched_points']==2809,'真实处理sol与PMC覆盖数不混淆')
        with obs[0]['path'].open(encoding='utf-8-sig',newline='') as stream:direct=next(row for row in csv.DictReader(stream) if row['PMC']=='6')
        check(data['chem'][6]['SiO2']=={'value':float(direct['SiO2_wt%']),'error':float(direct['SiO2_wt%_err'])} and data['chem'][6]['SiO2']['value']==41.6797 and data['chem'][6]['SiO2']['error']==2.1058 and '1σ' in payload['error_label'],'独立CSV复算PMC6成分与1σ绝对误差')
        point=next(x for x in payload['points'] if x['pmc']==6);raw=data['raw']['pts'][6]['PMC_0005_MCC']
        check(point['x']==raw[1] and point['y']==raw[0] and point['x']==229.584213 and point['y']==419.844055,'SIS/XML约定：sample=j、line=i，禁止沿用相反轴')
        pcd='pcd_0186_0683484439_000rcm_n007000006992743100050luj03.tif';pcd_csv=temp/'pcd.csv';pcd_csv.write_bytes(obs[0]['rxl'][0].read_bytes());pcd_csv.with_suffix('.xml').write_text('<root><file_name>'+pcd+'</file_name></root>',encoding='utf-8')
        pcd_raw=p.read_coordinates(pcd_csv);pcd_frames=p.frames(root,pcd_csv,pcd_raw)
        check(len(payload['frames'])==4 and sum(not f['missing'] for f in payload['frames'])==2 and [f['key'] for f in payload['frames'][:2]]==['PMC_0005_MCC','PMC_0006_MCC'] and pcd_raw['mcc_tifs']==[pcd] and len(pcd_frames)==1 and pcd_frames[0]['key']=='PMC_0005_MCC','pcw/pcd及数字/lu压缩字段均按XML和PMC关联，缺图不按时间猜测')
        check(all(p.finite(v) is None for v in ('nan','inf','-999','',None)) and p.finite('0')==0 and p.finite('2.4')==2.4,'缺值/非有限/负值不伪装为零；真实零保留')
        synthetic=temp/'synthetic.csv';synthetic.write_text('PMC,SiO2_wt%,SiO2_wt%_err\n1,0,nan\n2,-999,1\n',encoding='utf-8')
        chem,_=p.chemistry(synthetic)
        check(chem[1]['SiO2']=={'value':0.,'error':None} and chem[2]['SiO2']=={'value':None,'error':1.},'缺成分与缺误差分别保留')
        synthetic.write_text('PMC,SiO2_wt%\n1,4\n1,6\n',encoding='utf-8');check(fails(lambda:p.chemistry(synthetic)),'重复成分PMC拒绝静默覆盖')
        synthetic.write_text('PMC,PMC_0005_MCC_i,PMC_0005_MCC_j\n1,2,3\n1,4,5\n',encoding='utf-8');check(fails(lambda:p.read_coordinates(synthetic)),'重复坐标PMC拒绝静默覆盖')
        conflict=temp/'conflict';conflict.mkdir();stem='ps__0100_0674838381_000rqc__00417120483005510000___j0%d.csv'
        (conflict/(stem%1)).write_text('PMC,SiO2_wt%\n1,1\n',encoding='utf-8');(conflict/(stem%2)).write_text('PMC,SiO2_wt%\n1,2\n',encoding='utf-8');(conflict/'a').mkdir();(conflict/'a'/(stem%1)).write_text('PMC,SiO2_wt%\n1,999\n',encoding='utf-8')
        _,picked,warned=p.scan(conflict)
        check(picked and picked[0]['file'].endswith('j02.csv') and any('异字节' in w for w in warned),'低版本异字节副本留痕警告并按最高版本选择，不再依赖目录顺序')
        (conflict/'b').mkdir();(conflict/'b'/(stem%2)).write_text('PMC,SiO2_wt%\n1,5\n',encoding='utf-8')
        check(fails(lambda:p.scan(conflict)),'最高版本异字节冲突仍拒绝')
        stats=p.composition_stats(data,salt_free={'SO3':1.0,'Cl':1.0});s=stats['SiO2']
        check(s['n']==2809 and abs(s['mean']-43.8548)<5e-4 and abs(s['median']-45.2226)<5e-4 and abs(s['std']-10.6188)<5e-4 and s['salt_free_n']==670 and abs(s['salt_free_mean']-50.6479)<5e-4,'成分统计均值/中位/标准差/盐剔除均值锚定真实件')
        check(abs(stats['Na2O']['mean']-3.8478)<5e-4,'Na2O逐点均值锚定（bulk sum口径差的对照基准）')
        check(fails(lambda:p.composition_stats(data,salt_free={'XX':1.0})),'盐剔除阈值元素无效拒绝')
        labelled=temp/'labelled';p.export(data,labelled,'SiO2');page=(labelled/'index.html').read_text(encoding='utf-8')
        check('档案逐点算术均值' in page and 'bulk sum' in page and '中位数' in page and '45.2226' in page,'报告统计表携带估计量口径与Na2O注记')
        empty=temp/'empty';empty.mkdir()
        check(fails(lambda:p.scan('')) and p.scan(empty)[1]==[],'空路径不扫描开发目录，空数据明确返回')
        with patch.object(p.pixl,'HAS_CV2',False):png=p.render(data,'SiO2')
        image=Image.open(io.BytesIO(png));check(image.size==(752,676) and png.startswith(b'\x89PNG') and p.color(0,0,1)!=p.color(1,0,1),'标准版无OpenCV可导出原分辨率实测散点图')
        out=temp/'export';p.export(data,out,'SiO2')
        with (out/'composition.csv').open(encoding='utf-8-sig',newline='') as stream:exported={int(r['PMC']):r for r in csv.DictReader(stream)}
        check(len(exported)==2809 and float(exported[6]['SiO2_wt%'])==41.6797 and float(exported[6]['SiO2_wt%_err'])==2.1058 and (out/'composition.csv').read_bytes().startswith(b'\xef\xbb\xbf'),'CSV完整精度、误差、UTF8-BOM及全部点')
        report=(out/'index.html').read_text(encoding='utf-8');sources=json.loads((out/'sources.json').read_text(encoding='utf-8'))
        check('composition.csv' in report and 'map_SiO2.png' in report and '41.6797' in report and 'sample' not in str(sources['axis']) and sources['axis']=='label' and len(sources['files'])>=5 and str(temp) not in report,'报告相对链接、来源SHA和坐标约定完整且无本机路径')
        _,empty,_=p.canvas_data(data,2);noimage=p.export(data,temp/'noimage','SiO2',2)
        with patch.object(p.pixl,'mcc_gray',side_effect=ValueError('损坏图像')):broken=p.payload(data)
        check(not empty and noimage['map'] is None and (temp/'noimage/index.html').exists() and '缺' in ''.join(noimage['warnings']) and not broken['has_image'] and broken['chemical_points']==2809 and '读取失败' in ''.join(broken['warnings']),'缺/坏MCC仍可导出成分和报告，禁止生成虚构底图')
        extra=root/(obs[0]['rxl'][0].name.replace('10005','10006'));shutil.copy2(obs[0]['rxl'][0],extra)
        _,ambiguous,_=p.scan(root);amb=p.load(root,ambiguous[0]);check(len(ambiguous[0]['rxl'])==2 and amb['rxl_path'] is None and amb['chem'],'同RTT多坐标表不自动选最近或混拼')
        extra.unlink()
        other=root/(obs[0]['rxl'][0].name.replace('69927431','69927439'));shutil.copy2(obs[0]['rxl'][0],other)
        _,different,_=p.scan(root);check(len(different[0]['rxl'])==1,'相同定位、相近sol但不同RTT不配对')
        other.unlink()
        scripts='\n'.join(re.findall(r'<script[^>]*>(.*?)</script>',pixl_view.PAGE,re.S))
        syntax=subprocess.run([shutil.which('node'),'--check'],input=scripts,text=True,capture_output=True)
        check(syntax.returncode==0 and 'id="pixlCard"' in launch.PAGE and 'PIXL' in pixl_view.PAGE,'启动器入口和实际前端脚本语法')
        class Handler(BaseHTTPRequestHandler):
            def log_message(self,*a):pass
            def _send(self,status,ctype,body):self.send_response(status);self.send_header('Content-Type',ctype);self.end_headers();self.wfile.write(body)
            def do_GET(self):self._send(404,'text/plain',b'not found')
        state={};server=runtime_http.bind(runtime_http.decorate(Handler,state),0);thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start();url='http://127.0.0.1:%d'%server.server_port
        def request(path,body=None,origin=None):
            headers={'Content-Type':'application/json'}
            if origin:headers['Origin']=origin
            req=urllib.request.Request(url+path,data=json.dumps(body).encode() if body is not None else None,headers=headers)
            with urllib.request.urlopen(req,timeout=30) as r:
                raw=r.read();return json.loads(raw) if 'json' in r.headers.get('Content-Type','') else raw
        def bad(path,origin=None):
            try:request(path,origin=origin);return 0
            except urllib.error.HTTPError as e:return e.code
        try:
            with patch.object(urllib.request,'urlopen',wraps=urllib.request.urlopen):
                loaded=request('/api/pixl/load',{'path':str(folder)});oid=loaded['observations'][0]['id'];q='?obs='+oid
                point=request('/api/pixl/point'+q+'&pmc=6');export=request('/api/pixl/export',{'obs':oid,'element':'SiO2'})
                body=request(export['links']['报告']);pic=request(export['links']['热图PNG'])
            check(point['values']['SiO2']['value']==41.6797 and b'41.6797' in body and pic.startswith(b'\x89PNG') and bad('/api/pixl/point'+q+'&pmc=999999')==400,'真实HTTP看点、导出及报告/PNG链接闭环，越界PMC拒绝')
            check(bad('/api/pixl/init',origin='https://attacker.example')==403 and bad('/data/pixl/unknown/../../secret')==400,'跨源拒绝与导出路径白名单')
        finally:
            server.shutdown();server.server_close();thread.join()
            logger=logging.getLogger('marsworkbench')
            for handler in list(logger.handlers):handler.close();logger.removeHandler(handler)
        target=Path(pixl_view.demo());first=target/manifest['files'][0]['file'];first.write_bytes(b'local changed')
        check(fails(pixl_view.demo) and first.read_bytes()==b'local changed','再次加载示例不覆盖已修改副本')
    print('共 %d 项断言，%d 项失败'%(N,BAD));return int(BAD>0)

if __name__=='__main__':raise SystemExit(main())
