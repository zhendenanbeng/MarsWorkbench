"""第十期：图像浏览、报告闭环与可信度产物的离线契约。"""

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout
import io
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import urllib.request
import urllib.error
from http.server import ThreadingHTTPServer
from unittest.mock import patch

import numpy as np
from PIL import Image

ROOT=_layout.project_layout.APP
sys.path.insert(0,str(ROOT))
from instruments import supercam_images as images
from instruments import REGISTRY
import console
sys.path.insert(0,str(_layout.project_layout.VERIFY))

TOTAL=FAIL=0
def check(ok, text):
    global TOTAL,FAIL
    TOTAL+=1;FAIL+=not ok
    print(('OK ' if ok else 'FAIL ')+text)


def image_checks(temp):
    name='scam_0187_0683540642_630_ci__scam01187_bellegarde_187_scam__01p01.fits'
    r=images.parse_ci(name)
    check(r['sol']==187 and r['point']==1 and r['sclk']=='0683540642.630' and r['seq']=='scam01187' and not images.parse_ci('../'+name), 'CI名称须解析Sol/SCLK/点/序列并拒绝路径')
    stems=['sif_0186_0683466183_675eby_n0070000srlc00621_0000lmj','sif_0186_0683466184_675eby_n0070000srlc00621_0000lmj','sif_0188_0683643528_171eby_n0070000srlc00678_0000lmj']
    inventory=temp/'inventory';inventory.touch()
    with patch.object(images.watson,'_inventory_path',return_value=str(inventory)),patch.object(images.watson,'load_inventory',return_value=stems):
        rows=images.watson_groups(186)
        check(len(rows)==1 and rows[0]['count']==2 and rows[0]['placement']=='n0070000' and not images.watson_groups(187),'WATSON严格按Sol/序列/placement归组，不借邻Sol')
    # Native big-endian image, physical scaling and a masked integer BLANK.
    cards=[('SIMPLE','T'),('BITPIX','16'),('NAXIS','2'),('NAXIS1','3'),('NAXIS2','2'),('BSCALE','2'),('BZERO','100'),('BLANK','-1')]
    header=''.join((k.ljust(8)+'= '+v).ljust(80) for k,v in cards)+'END'.ljust(80)
    stored=np.array([[-1,0,1],[2,3,4]],dtype='>i2')
    fits=temp/'small.fits';fits.write_bytes(header.ljust(2880).encode()+stored.tobytes())
    a=images.fits_image(fits)
    label='ODL_VERSION_ID = ODL3\n^IMAGE = 1025 <BYTES>\nOBJECT = IMAGE\nLINES = 2\nLINE_SAMPLES = 3\nSAMPLE_BITS = 32\nSAMPLE_TYPE = IEEE_REAL\nBANDS = 2\nBAND_STORAGE_TYPE = BAND_SEQUENTIAL\nEND_OBJECT = IMAGE\nEND\n'
    vicar=temp/'small.IMG';vicar.write_bytes(label.ljust(1024).encode()+np.arange(12,dtype='>f4').tobytes())
    v=images.vicar_image(vicar)
    png=np.asarray(Image.open(io.BytesIO(images.stretch_png(v))));
    check(np.isnan(a[0,0]) and a[1,2]==108 and np.array_equal(v,np.arange(6).reshape(2,3)) and png[0,0]==0 and png[-1,-1]==255,'FITS物理缩放/BLANK及VICAR浮点大端BSQ波段0、显示拉伸须正确')
    with patch.object(images,'CACHE',temp/'no_network'),patch.object(images.pds_listing,'list_files',side_effect=OSError('offline')),patch.object(images.watson,'_inventory_path',return_value=str(temp/'missing')),patch.object(images.watson,'load_inventory',side_effect=AssertionError('不应抓全任务清单')):
        d=images.listing(187)
        check('CI目录列表失败' in d['ci_error'] and '本地清单未缓存' in d['watson_error'] and d['ci']==[],'断网与缺清单分别降级，不猜图、不触发全量下载')
    node=shutil.which('node');assert node,'前端护栏需要Node.js'
    page=console.build_page(REGISTRY['supercam'])
    syntax=subprocess.run([node,'--check'],input='\n'.join(re.findall(r'<script[^>]*>(.*?)</script>',page,re.S)),text=True,capture_output=True)
    check(syntax.returncode==0 and all('id="'+key+'"' in page for key in ['scamImagePanel','scamCI','scamWatsonGroup','scamWatson','scamImageViewport']) and 'id="scamImagePanel"' not in console.build_page(REGISTRY['chemcam']), 'SuperCam面板包含两类选择/缩放，其他台无面板，实际页JS语法正确')
    state={'app':{name:{'folder':str(temp),'observations':[]} for name in REGISTRY}}
    server=ThreadingHTTPServer(('127.0.0.1',0),console.make_handler(state));thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start()
    codes=[]
    try:
        for path in ['/api/scamimg?app=chemcam&sol=187','/api/scamimg/image?app=apxs&kind=ci']:
            try:urllib.request.urlopen('http://127.0.0.1:%d'%server.server_port+path)
            except urllib.error.HTTPError as e:codes.append(e.code)
        check(codes==[400,400],'两个图像API均须拒绝非SuperCam仪器')
    finally:server.shutdown();server.server_close();thread.join()


def map_checks(temp):
    from instruments import report
    b=REGISTRY['chemcam'];obs=b.scan_folder(str(_layout.project_layout.VERIFY/'fixtures/chemcam'))
    data=temp/'mapdata';maps=data/'_maps';maps.mkdir(parents=True)
    for clock in ['486263086','486264264']:
        stem='map_sol00013_ccam02013_'+clock
        Image.fromarray(np.zeros((2,2),dtype=np.uint8)).save(maps/(stem+'.png'))
        (maps/(stem+'_points.csv')).write_text('point,line,sample\n',encoding='utf-8')
        (maps/(stem+'_points.json')).write_text(json.dumps({'basis':'EDR <原片>'}),encoding='utf-8')
    (maps/'random.png').touch();(maps/'map_sol00013_ccam02013_111111111.png').touch()
    rows=report._map_rows(data,temp/'report')
    report.build_report(b,obs,str(temp/'report'),data_dir=str(data))
    page=(temp/'report/index.html').read_text(encoding='utf-8')
    check(len(rows)==2 and '元素分布图 2 帧' in page and 'EDR &lt;原片&gt;' in page and page.count('map_sol00013_ccam02013_486263086.png')>=3 and all((temp/'report'/r['url']).exists() for r in rows),'地图卡与对应观测链接、相对路径和底图口径须出现')
    check('random.png' not in page and '111111111' not in page and '_points.csv' in page,'乱名与缺坐标CSV的地图不得进入报告')
    none=temp/'emptydata';none.mkdir()
    report.build_report(b,obs,str(temp/'emptyreport'),data_dir=str(none))
    absent=(temp/'emptyreport/index.html').read_text(encoding='utf-8')
    with patch.object(report,'_map_rows',side_effect=AssertionError('非ChemCam不扫描地图')):
        report.build_report(REGISTRY['apxs'],[],str(temp/'apxs'),data_dir=str(data))
    check('elementMaps' not in absent and 'elementMapLinks' not in absent,'无地图无卡片/占位；非ChemCam完全不扫描_maps')


def trust_checks(temp):
    import run_all
    import make_trust_report as trust
    from datetime import datetime,timezone,timedelta
    now=datetime.now(timezone.utc);path=temp/'last.json'
    results=[{'script':name,'what':what,'code':0,'n':n,'expect':n} for name,what,n in run_all.CHECKS]
    with patch.object(run_all,'LAST_RUN',str(path)):
        snap=run_all.write_snapshot(results,now.isoformat())
    text=trust.generate(path,now)
    check(all('## '+title in text for title in ['护栏总览','金标准与差异台账','已知边界','各期验证文档索引','口径申明']) and '上限60点' in text and '是否当前可用：是' in text,'可信度报告须包含活台账、边界、文档索引及完整护栏记录')
    old=dict(snap,finished_utc=(now-timedelta(days=8)).isoformat());path.write_text(json.dumps(old),encoding='utf-8')
    stale=trust.generate(path,now)
    snap['source_sha256']={};path.write_text(json.dumps(snap),encoding='utf-8');changed=trust.generate(path,now)
    check('数据过期，请重跑 run_all' in stale and 'color:#c62828' in stale and '源码已变化' in changed,'超过七天或源码不对应不得显示当前通过')
    missing=trust.generate(temp/'absent.json',now)
    partial=dict(snap,results=[]);path.write_text(json.dumps(partial),encoding='utf-8')
    check('没有可核验的总检记录' in missing and '数据过期' in trust.generate(path,now),'缺失或不完整总检记录须保守标记且不崩溃')


def library_checks():
    import hashlib
    from collections import Counter
    from sherloc import library
    from phase10_library_diagnosis import score_auc,METHODS
    path=_layout.project_layout.VERIFY/'fixtures/历史数值基准/C_交付诊断.json'
    data=json.loads(path.read_text(encoding='utf-8'));rows=data['rows']
    # Product paths/HTML may change; freeze the actual numerical/selection/read kernels.
    science={'library.py','minerals.py','qc.py','preprocess.py','peaks.py','fitting.py','reader.py','label.py','organize.py','repeatability.py'}
    sources=all(hashlib.sha256(_layout.project_layout.runtime_path(name).read_bytes()).hexdigest()==sha for name,sha in data['production_source_sha256'].items() if Path(name).name in science)
    inputs=all(hashlib.sha256(_layout.project_layout.resolve_legacy(ROOT.parent/name).read_bytes()).hexdigest()==sha for name,sha in data['input_sha256'].items())
    check(sources and inputs and len(rows)==233 and sum(r['named'] for r in rows)==67 and data['unnamed']==166,'谱库诊断须绑定冻结来源和233个规则组（2026-11-19 当前树重算），不能套到新内核/新数据')
    visible=[len([v for v in library.sample_peaks(modes) if 800<=v<=1800]) for cls in library.REFERENCE.values() for _,_,modes in cls['samples']]
    check(data['global_snr_scale_invariant'] and len(visible)==103 and visible.count(0)==38 and visible.count(1)==19,'锁定统一SNR缩放不改变当前分数、参考窗口信息不足的诊断事实')
    cross=Counter((r['rule_key'],r['library_class']) for r in rows)
    check(cross[('unassigned','oxychlorine')]==51 and cross[('unassigned','borate')]==50 and cross[('unassigned','sulfate')]==18 and all(abs(score_auc(rows,m)-data['full_data_auc'][m])<1e-12 for m in METHODS),'交叉表与原始代理AUC须可从冻结峰组独立复算，不套用预期误分对')
    folds=data['folds'];tests=[s for f in folds for s in f['test_sols']]
    check(len(tests)==len(set(tests)) and all(not(set(f['test_sols'])&set(f['train_sols'])) for f in folds) and data['decision']=='close_no_heldout_path' and max(data['sol_blocked_macro_auc'][m] for m in METHODS[1:])<.70 and not data['default_changed'],'Sol留出不能泄漏；三项均未达继续线时必须止损，不接默认')


def main():
    import os as _os
    if not _os.path.isdir(_os.path.join(str(ROOT), "..", "03_下载数据", "SHERLOC", "练习数据", "processed")):
        print("SKIP 数据不在位：练习数据 processed 树不在本仓")
        print("共 0 项断言，0 项失败")
        return 0

    with tempfile.TemporaryDirectory(prefix='sherloc-phase10-') as folder:
        image_checks(Path(folder))
        map_checks(Path(folder))
        trust_checks(Path(folder))
        # library_checks 与冻结 C_交付诊断.json（2026-11-19 当前练习树重算，
        # 旧 0.16.0 时代基线存档于 C_交付诊断_0160旧树存档.json）对账；
        # 练习树或 sherloc 内核再变更时指纹比对 FAIL——按纪律重新诊断冻结，
        # 不静默跳过。
        library_checks()
    print('共 %d 项断言，%d 项失败'%(TOTAL,FAIL))
    return int(FAIL>0)

if __name__=='__main__':raise SystemExit(main())
