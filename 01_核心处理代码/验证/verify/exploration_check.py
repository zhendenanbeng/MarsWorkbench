# -*- coding: utf-8 -*-
"""特征积分、分型与任务接口的行为验证；不联网、不覆盖既有报告。"""
from __future__ import annotations

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout
import csv
import json
import os
from pathlib import Path
import sys
import tempfile
import threading
import time
import urllib.request
import urllib.error
from http.server import ThreadingHTTPServer
from unittest.mock import patch

import numpy as np

ROOT = _layout.project_layout.APP
sys.path.insert(0, str(ROOT))
from sherloc import explore
from sherloc.organize import scan_folder
import view

TOTAL = FAIL = 0


def check(value, message):
    global TOTAL, FAIL
    TOTAL += 1
    if not value:
        FAIL += 1
        print(' FAIL', message)


def numeric_checks():
    x = np.array([0., .4, 1.2, 2., 3.])
    w = explore.integration_weights(x, .3, 2.7)
    exact = (2.7**2+2.7) - (.3**2+.3)
    check(np.isclose((2*x+1)@w, exact), '不等间距与插值端点的线性谱积分必须精确')
    check(np.isclose(w.sum(), 2.4), '常数谱积分必须等于窗口宽度')
    spectra = np.vstack((2*x+1, -2*x-1, np.zeros(len(x))))
    a = explore.band_features(x, spectra, np.ones(3)*.01, (.3,2.7))
    check(a['area'][0] > 0 and a['area'][1] < 0, '面积必须保留正负号')
    check(not a['valid'][2], '零信号不能通过积分SNR门槛')
    ratio = explore.ratio_features(a,a,np.ones(3)*.01)
    check(ratio['value'][0] == 1, '相同窗口面积比为1')
    check(ratio['sigma'][0] < 1e-12, '相同窗口比值误差必须考虑共享通道协方差')
    check(not ratio['valid'][1] and not ratio['valid'][2], '负面积和零分母不得报告比值')
    low = explore.band_features(x, spectra, np.ones(3)*1e9, (.3,2.7))
    check(not low['valid'].any(), '低SNR积分保留数值但必须标记无效')
    for xx, lo, hi in [(x, -1,2),(x[::-1],.3,2.7),(x,2,1)]:
        try:
            explore.integration_weights(xx,lo,hi)
            check(False,'无效积分窗口/波数轴应拒绝')
        except ValueError:
            check(True,'无效积分输入已拒绝')
    rng = np.random.default_rng(123)
    axis = np.linspace(800,1800,101)
    prototypes = [np.exp(-((axis-c)/30)**2) for c in (900,1150,1650)]
    y = np.vstack([p + rng.normal(0,.003,(20,len(axis))) for p in prototypes])
    y *= np.linspace(1,20,len(y))[:,None]
    original = y.copy()
    model = explore.spectral_types(axis,y,np.ones(len(y),bool),3)
    check(model['status']=='ok' and model['effective_k']==3,'三个明确谱形应识别为三组')
    check(all(len(set(model['labels'][i:i+20]))==1 for i in (0,20,40)), '同形不同强度的光谱应归为同组')
    check(len(set(model['labels'][[0,20,40]]))==3,'三个不同谱形不可混为一组')
    again = explore.spectral_types(axis,y,np.ones(len(y),bool),3)
    check(np.array_equal(model['labels'],again['labels']),'重复运行分型可复现')
    check(np.array_equal(original,y),'分型不得修改输入矩阵')
    check(0 < model['retained_variance'] <= 1+1e-12,'PCA方差比例必须合法')
    same = explore.spectral_types(axis,np.vstack([prototypes[0]*i for i in (1,2,3)]),np.ones(3,bool),3)
    check(same['status']=='single_pattern' and same['effective_k']==1,'同形谱不得被强制拆分')
    check(len(same['groups'])==1,'同形谱也应有平均谱和组摘要')
    tiny = explore.spectral_types(axis,y[:1],np.ones(1,bool),3)
    check(tiny['status']=='insufficient','一个候选点须报告样本不足')
    bad = np.zeros((3,len(axis)));bad[0,0]=np.nan
    check(explore.spectral_types(axis,bad,np.ones(3,bool),3)['status']=='insufficient','非有限值和零谱不参加分型')
    encoded = json.dumps(explore.clean_json({'x':np.array([np.nan,np.inf,1])}),allow_nan=False)
    check(json.loads(encoded)['x']==[None,None,1.0],'JSON不得输出NaN或Infinity')
    for bands,snr,k,scope in [(((900,800),(990,1030)),3,3,'all'),(explore.DEFAULT_BANDS,0,3,'all'),(explore.DEFAULT_BANDS,3,1,'all'),(explore.DEFAULT_BANDS,3,2.5,'all'),(explore.DEFAULT_BANDS,3,3,'bad')]:
        try:
            explore.validate_parameters(bands,snr,k,scope)
            check(False,'非法参数应拒绝')
        except ValueError:
            check(True,'非法参数已拒绝')


def real_data_checks(folder):
    data = ROOT.parent/'03_下载数据/SHERLOC/练习数据/processed/sol_00370'
    _,scans = scan_folder(str(data))
    scan = next(s for s in scans if s.spectra is not None and (s.n_points or 0)>10)
    before = {str(p.data_path):explore.file_hash(p.data_path) for p in scan.products if Path(p.data_path).is_file()}
    result = explore.run(scan,folder)
    out = Path(result['out'])
    with (out/'points.csv').open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f))
    manifest=json.loads((out/'analysis.json').read_text(encoding='utf-8'),parse_constant=lambda v: (_ for _ in ()).throw(ValueError(v)))
    check(len(rows)==scan.n_points,'所有点都必须写入特征表，不能只保留参与分型的点')
    check(result['candidate_points']>=2 and result['type_status']=='ok','真实未归类数据应产生分型结果')
    check(all(not r['known_classes'] for r in rows if r['type']),'未归类模式不得混入已知矿物类别点')
    check(all(r['excluded_reason'] for r in rows if not r['type']),'未参与分型的点必须有排除原因')
    check((out/'index.html').is_file() and (out/'type_spectra.csv').is_file(),'须生成报告与类型平均谱')
    check(result['maps']==4,'单ACI坐标系须生成A/B/比值/类型四张图')
    check(any(p['path'].lower().endswith('.xml') for p in manifest['input_files']),'元数据XML也必须记录校验值')
    check(all(explore.file_hash(p)==digest for p,digest in before.items()),'真实输入文件不得被更改')
    class NoPositions:
        def __getattr__(self,name):return getattr(scan,name)
        def read_shot_positions(self):return []
    no_positions = explore.run(NoPositions(),folder)
    check(no_positions['maps']==0,'无RLS时不得编造地图坐标')
    check(Path(no_positions['index']).is_file(),'无RLS时数据分析报告仍应可用')
    try:
        explore.run(scan,folder,cancel=lambda:True)
        check(False,'取消任务应停止')
    except explore.CancelledError as exc:
        check(Path(exc.output_path,'job.json').is_file(),'取消后须保留任务状态和产物位置')
        check(not Path(exc.output_path,'index.html').exists(),'取消不能伪装为完整报告')


def api_checks(folder):
    state={'folder':folder,'listing':[],'scans':[object()], 'export':{},'dl':{},'aci':{},'watson':{}}
    server=ThreadingHTTPServer(('127.0.0.1',0),view.make_handler(state))
    threading.Thread(target=server.serve_forever,daemon=True).start()
    base='http://127.0.0.1:%d'%server.server_port
    def request(path,payload=None):
        req=urllib.request.Request(base+path,data=json.dumps(payload or {}).encode() if payload is not None else None,headers={'Content-Type':'application/json'})
        try:
            with urllib.request.urlopen(req,timeout=5) as r:return r.status,r.read()
        except urllib.error.HTTPError as e:return e.code,e.read()
    def wait_done():
        until=time.monotonic()+5
        while state.get('analysis',{}).get('running') and time.monotonic()<until:time.sleep(.01)
        check(not state.get('analysis',{}).get('running'),'任务应结束')
    entered,release=threading.Event(),threading.Event()
    def blocked(*args,**kwargs):
        entered.set();release.wait(3);raise explore.CancelledError('测试取消')
    try:
        bad_request=urllib.request.Request(base+'/api/analysis/start',data=b'{',headers={'Content-Type':'application/json'})
        try:
            urllib.request.urlopen(bad_request,timeout=5)
            check(False,'无效JSON不能启动默认任务')
        except urllib.error.HTTPError as exc:
            check(exc.code==400 and not state.get('analysis'),'无效JSON必须拒绝且不创建任务')
        code,_=request('/api/analysis/start',{'scan':9})
        check(code==400,'无效扫描须在启动前拒绝')
        code,_=request('/api/analysis/start',{'scan':0,'k':1})
        check(code==400,'无效参数须在启动前拒绝')
        with patch.object(explore,'run',side_effect=blocked):
            code,body=request('/api/analysis/start',{'scan':0,'out':folder})
            check(code==200 and json.loads(body)['started'],'任务接口应启动后台任务')
            entered.wait(2)
            code,_=request('/api/analysis/start',{'scan':0,'out':folder})
            check(code==409,'同时只允许一个分析任务')
            code,body=request('/api/analysis/cancel',{})
            check(code==200 and json.loads(body)['cancelling'],'取消请求应传到任务')
            release.set();wait_done()
            check(state['analysis'].get('cancelled'),'取消须显示取消状态')
        counter=[0]
        def completed(*args,**kwargs):
            counter[0]+=1;out=Path(folder)/('report%d'%counter[0]);out.mkdir()
            (out/'index.html').write_text('<h1>analysis report</h1>',encoding='utf-8')
            return {'out':str(out),'index':str(out/'index.html'),'maps':0,'points':1,'types':0,'candidate_points':1,'type_status':'insufficient'}
        with patch.object(explore,'run',side_effect=completed):
            _,body=request('/api/analysis/start',{'scan':0,'out':folder});first=json.loads(body)['id'];wait_done()
            code,body=request('/analysis-result/'+first+'/index.html')
            check(code==200 and b'analysis report' in body,'报告通过同源地址可读取')
            code,_=request('/analysis-result/'+first+'/%2e%2e/secret.json')
            check(code==403,'报告接口不可越过产物目录')
            request('/api/analysis/start',{'scan':0,'out':folder});wait_done()
            check(request('/analysis-result/'+first+'/index.html')[0]==200,'新任务不能令旧报告链接失效')
    finally:
        release.set();server.shutdown();server.server_close()


def main():
    if not os.path.isdir(os.path.join(str(ROOT), "..", "03_下载数据", "SHERLOC", "练习数据", "processed", "sol_00370")):
        print("SKIP 数据不在位：PDS 16.0 未发布 sol_00370 的 processed RRS（exploration 需它在位）")
        print("共 0 项断言，0 项失败")
        return 0

    numeric_checks()
    with tempfile.TemporaryDirectory(prefix='sherloc-explore-guards-') as tmp:
        real_data_checks(tmp)
        api_checks(tmp)
    print('共 %d 项断言，%d 项失败'%(TOTAL,FAIL))
    return 1 if FAIL else 0


if __name__=='__main__':raise SystemExit(main())
