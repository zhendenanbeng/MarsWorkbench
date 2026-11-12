# -*- coding: utf-8 -*-
"""多产品报告、导出、成分缓存和会话状态的离线回归检查。"""
from __future__ import annotations

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout
import csv
import json
from pathlib import Path
import shutil
import sys
import tempfile
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from http.server import ThreadingHTTPServer
from unittest.mock import patch

ROOT=_layout.project_layout.APP
sys.path.insert(0,str(ROOT))
import console
from instruments import REGISTRY, report, supercam
from sherloc.desktop import remember_instrument

FIX=_layout.project_layout.VERIFY/'fixtures'
TOTAL=FAIL=0


def check(value,message):
    global TOTAL,FAIL
    TOTAL+=1
    if not value:
        FAIL+=1;print(' FAIL',message)


def mixed_products(temp):
    b=REGISTRY['supercam'];obs=b.scan_folder(str(FIX/'supercam'))
    passive=[o for o in obs if o.product in ('CP2','CP3')]
    active=next(o for o in obs if o.product=='CL1')
    with patch.object(b,'annotate',side_effect=AssertionError('被动谱调用了发射线判定')),patch.object(b,'moc_payload',side_effect=AssertionError('被动谱关联了LIBS成分')):
        for o in passive:
            p=console.plot_for({'observations':[o]},b,0,0,-1,{})
            check(p['passive'] and not p['annotate']['marks'],'被动谱API不得返回发射线')
            check('moc' not in p,'被动谱API不得返回LIBS成分')
            console.render_point_png(b,o,0,None,{},str(temp/(o.product+'.png')))
            check((temp/(o.product+'.png')).is_file(),'被动谱PNG必须生成')
        res=report.build_report(b,passive,str(temp/'passive'))
    check(len(set(res['csv']))==2,'CP2/CP3同序列CSV必须分开')
    check(not res['lines_csv'],'被动谱不应生成发射线SNR表')
    for name in res['csv']:
        with (temp/'passive'/name).open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f))
        check(len(rows)==1 and not rows[0]['检出元素'],'被动谱不能假报元素')
        check(not rows[0]['发数'] and '被动' in rows[0]['判定口径'],'被动谱不得虚构激光发数')
        check(all(not v for k,v in rows[0].items() if k.startswith('官方预测')),'被动谱成分列必须为空')
    page=Path(res['index']).read_text(encoding='utf-8')
    check('波长 (µm)' in page and 'Spectrum DN' in page,'CP2报告必须保留实际单位和谱列')
    check('不适用（被动光谱）' in page,'被动报告必须明确成分判定不适用')
    mixed=report.build_report(b,obs,str(temp/'mixed'))
    check(len(mixed['csv'])==len(set(mixed['csv']))==len(obs),'混合产品报告不覆盖CSV')
    check(len(mixed['lines_csv'])==1,'只有CL1应有发射线明细')
    with (temp/'mixed'/mixed['lines_csv'][0]).open(encoding='utf-8-sig',newline='') as f:line_rows=list(csv.DictReader(f))
    check(len(line_rows)>3 and any(r['达到门槛']=='False' for r in line_rows),'谱线明细应包含完整判定而非只列前三条')
    state={'app':{'supercam':{'observations':obs}},'export':{'supercam':{}}}
    console._export_worker(state,b,str(temp/'png'),{'dust':True,'despike':True})
    check(not state['export']['supercam'].get('error'),'混合产品PNG导出不得失败')
    check(len(list((temp/'png').glob('*.png')))==sum(o.n_points for o in obs),'导出完成数必须等于实际PNG文件数')
    check(len({console.export_name(o,i) for o in obs for i in range(o.n_points)})==sum(o.n_points for o in obs),'混合产品导出名称不得冲突')
    missing={'matched':True,'rows':[{'oxide':'SiO2','wt':10.27,'err':None}]}
    captured=[]
    with patch.object(b,'moc_payload',return_value=missing),patch.object(console,'render_png',side_effect=lambda *a,**kw:captured.extend(kw['notes'])):
        console.render_point_png(b,active,0,None,{},str(temp/'missing_sigma.png'))
    check(any('未提供' in line for line in captured),'缺失标准差不得伪装成零误差')
    with patch.object(b,'moc_payload',return_value=missing):
        numeric=report.build_report(b,[active],str(temp/'numeric'))
    with (temp/'numeric'/numeric['csv'][0]).open(encoding='utf-8-sig',newline='') as f: row=next(csv.DictReader(f))
    check(float(row['官方预测 SiO2 wt%'])==10.27 and row['官方预测 SiO2 σ wt%']=='','成分与标准差须分列为数值，缺失误差留空')
    apxs=REGISTRY['apxs'];aobs=apxs.scan_folder(str(FIX/'apxs'))
    numeric=report.build_report(apxs,aobs,str(temp/'numeric_apxs'))
    with (temp/'numeric_apxs'/numeric['csv'][0]).open(encoding='utf-8-sig',newline='') as f: row=next(csv.DictReader(f))
    check(float(row['APXS FeO wt%'])>5 and float(row['APXS FeO σ wt%'])>0,'APXS官方成分及误差须同步导出CSV')


def moc_cache(temp):
    original=(FIX/'supercam/_moc/supercam_moc.csv').read_bytes()
    target=temp/'moc.csv';target.write_bytes(original)
    b=supercam.SuperCamBackend()
    with patch.object(supercam,'moc_cache_path',return_value=str(target)):
        with patch('instruments.pds_listing.http_get',return_value=b'<html>temporary server error</html>'):
            try:
                b.ensure_moc(force=True);check(False,'错误页必须拒绝缓存')
            except ValueError:
                check(target.read_bytes()==original,'坏下载不能覆盖原有成分表')
        check(not Path(str(target)+'.part').exists(),'无效下载应清理临时文件')
        with patch('instruments.pds_listing.http_get',return_value=original):b.ensure_moc(force=True)
        check(len(b.moc_table()[0])==3,'有效成分表可以更新缓存')
        stem='SCAM_0998_0755522849_085_CL1_SCAM02998_AEGIS_0998A__________06P01.FITS'
        payload=b.moc_payload(stem)
        check(payload.get('matched'),'成分键应容忍产品文件大小写及扩展名')
    check(supercam._to_f('NaN') is None and supercam._to_f('Infinity') is None,'成分缺值不可进入JSON为NaN/Infinity')


def settings(temp):
    path=temp/'settings.json';value={'folders':{'apxs':'original'},'other':42}
    path.write_text(json.dumps(value),encoding='utf-8')
    for bid in ('chemcam','supercam','apxs','sherloc'):
        remember_instrument(bid,str(path));d=json.loads(path.read_text(encoding='utf-8'))
        check(d['last_app']==bid and d['folders']==value['folders'] and d['other']==42,'记录仪器应保留其他设置：'+bid)
    before=path.read_bytes();remember_instrument('invalid',str(path))
    check(path.read_bytes()==before,'未知仪器不能修改设置')
    path.write_text('{bad json',encoding='utf-8');remember_instrument('apxs',str(path))
    check(path.read_text(encoding='utf-8')=='{bad json','最近仪器记录不能吞掉损坏的配置')
    check(not list(temp.glob('.instrument-*.tmp')),'原子设置写入不得留下临时文件')


def report_api(temp):
    b=REGISTRY['apxs'];observations=b.scan_folder(str(FIX/'apxs'))
    data=temp/'data';report_dir=data/'_qc/consoles';report_dir.mkdir(parents=True)
    index=report_dir/'index.html';index.write_text('<html>report</html>')
    wrong=data/'outside.html';wrong.write_text('outside')
    csvfile=report_dir/'table.csv';csvfile.write_text('x,y')
    state={'app':{'apxs':{'folder':str(data),'observations':observations}},'export':{'apxs':{}},'dl':{'apxs':{}},'qc':{'apxs':{}}}
    server=ThreadingHTTPServer(('127.0.0.1',0),console.make_handler(state))
    threading.Thread(target=server.serve_forever,daemon=True).start()
    base='http://127.0.0.1:%d'%server.server_port
    def get(endpoint,**params):
        try:
            with urllib.request.urlopen(base+endpoint+'?'+urllib.parse.urlencode({'app':'apxs',**params}),timeout=5) as r:return r.status,json.load(r)
        except urllib.error.HTTPError as exc:return exc.code,json.load(exc)
    try:
        with patch.object(console,'open_app_window') as opened:
            check(get('/api/openreport',path=str(index))[0]==200 and opened.call_count==1,'真实质控HTML应可打开')
            check(get('/api/openreport',path=str(wrong))[0]==400,'数据目录里非报告HTML也应拒绝')
            check(get('/api/openreport',path=str(csvfile))[0]==400,'打开报告接口只允许HTML')
            check(get('/api/openreport',path=str(report_dir/'../outside.html'))[0]==400,'报告路径不能越过产物目录')
        entered,release=threading.Event(),threading.Event();captured={}
        def build(backend,obs,out_dir,progress,data_dir=None):
            captured.update(obs=obs,out=out_dir);entered.set();release.wait(3)
            return {'out':out_dir,'index':str(Path(out_dir)/'index.html'),'count':len(obs),'points':1}
        with patch.object(report,'build_report',side_effect=build):
            check(get('/api/qc')[1].get('started'),'后台报告应能启动')
            entered.wait(2)
            state['app']['apxs']={'folder':str(temp/'other'),'observations':[]}
            release.set();deadline=time.monotonic()+5
            while state['qc']['apxs'].get('running') and time.monotonic()<deadline:time.sleep(.01)
            check(captured['out']==str(report_dir) and len(captured['obs'])==len(observations),'后台报告应绑定启动时目录与观测快照')
    finally:
        server.shutdown();server.server_close()


def main():
    with tempfile.TemporaryDirectory(prefix='sherloc-multi-guards-') as folder:
        temp=Path(folder)
        mixed_products(temp);moc_cache(temp);settings(temp);report_api(temp)
    print('共 %d 项断言，%d 项失败'%(TOTAL,FAIL))
    return 1 if FAIL else 0


if __name__=='__main__':raise SystemExit(main())
