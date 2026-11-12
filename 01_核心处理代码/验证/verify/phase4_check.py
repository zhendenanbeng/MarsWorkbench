"""Fourth-cycle feature contracts, expanded one independently delivered task at a time."""

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout
from pathlib import Path
import csv
import sys
import tempfile
import json
import threading
import subprocess
import shutil
import urllib.request
from http.server import ThreadingHTTPServer
from unittest.mock import patch
import numpy as np
from PIL import Image
ROOT=_layout.project_layout.APP;sys.path.insert(0,str(ROOT))
from instruments import REGISTRY,lines,evidence,report,notes,spectrum_export,absorption,fitsio,rmi,filtering,history,observation_summary,download_queue
from sherloc import qc,workbench,reference
import console
import view
NODE=shutil.which('node')
if NODE is None:raise RuntimeError('第四期前端逻辑护栏需要Node.js（不属于软件运行依赖）；请将node加入PATH')
TOTAL=FAIL=0
def check(ok,name):
    global TOTAL,FAIL
    TOTAL+=1
    if not ok:FAIL+=1;print(' FAIL '+name)

def task_a(temp):
    for bid in ['chemcam','supercam']:
        b=REGISTRY[bid];obs=next(o for o in b.scan_folder(str(_layout.project_layout.VERIFY/'fixtures'/bid)) if not getattr(o,'passive',False) and o.product not in ['CP2','CP3'])
        pt=obs.points[0];d=b.load_point(pt,None,{'dust':True});marks,elems=b.annotate_data(d)
        check(d['shot_matrix'].shape==(len(d['x']),25),bid+'保留25发矩阵与波长轴匹配')
        raw=b.load_point(pt,None,{})
        check(raw['shot_matrix'].shape[1]==30,bid+'不除尘时分母为30')
        base_marks,base_elements=b.annotate(d['x'],d['y'])
        check(elems==base_elements and [m['detected'] for m in marks]==[m['detected'] for m in base_marks],bid+'经验统计不得更改当前谱分类')
        individual=[b.annotate(d['x'],d['shot_matrix'][:,i]) for i in range(25)]
        check(all(m['shot_detected']==sum(mm[j]['detected'] for mm,_ in individual) for j,m in enumerate(marks)),bid+'逐线分子可逐发独立复算')
        check(all(m['element_shot_detected']==sum(m['elem'] in es for _,es in individual) for m in marks),bid+'元素支持须逐发满足两条独立线')
        check(all(m['shot_total']==25 and 0<=m['shot_detected']<=25 for m in marks),bid+'分母统一，分子不能越界')
        check(all(m['fwhm'] is None or 0<m['fwhm']<=1.4 for m in marks if m['detected']),bid+'有效检出宽度遵守原形状门槛')
        again=evidence.shot_statistics(d['x'],d['shot_matrix']);again[0]['shot_detected']=-7
        check(evidence.shot_statistics(d['x'],d['shot_matrix'])[0]['shot_detected']>=0,bid+'缓存返回值不能污染后续结果')
        n=console.plot_for({'observations':[obs]},b,0,0,-1,{'dust':True})
        check(bool(n['annotate']['element_evidence']),bid+'API带元素逐发证据')
        result=report.build_report(b,[obs],str(temp/bid))
        rows=list(csv.DictReader((temp/bid/result['lines_csv'][0]).open(encoding='utf-8-sig')))
        check(all(int(r['保留发数'])==25 for r in rows),bid+'报告谱线表使用默认预处理分母')
        point=next(csv.DictReader((temp/bid/result['csv'][0]).open(encoding='utf-8-sig')))
        check(bool(point['逐发元素支持']) and 'SNR中位' in point['逐线经验统计'],bid+'逐点CSV保留检出率、SNR与宽度')
        console.render_point_png(b,obs,0,None,{'dust':True},str(temp/(bid+'.png')))
        with Image.open(temp/(bid+'.png')) as im:check(im.height>860,bid+'PNG必须带实际逐线证据注记')
    x=np.linspace(270,290,1200);y=100*np.exp(-.5*((x-280)/.1)**2)+np.sin(x*30)*.02
    def truncated(*args,**kw):
        return {'fwhm':.4,'left':279.8,'right':280.2,'censored':True} if kw.get('allow_censored') else None
    with patch('sherloc.peaks.half_maximum',side_effect=truncated):
        mark=lines.match_lines(x,y,[('X',280.,'test')],tol=.35)[0][0]
    check(mark['fwhm_censored'] and mark['fwhm_observed']>0 and not mark['detected'],'截断诊断须向上传递，不能借†取得元素支持')
    passive=next(o for o in REGISTRY['chemcam'].scan_folder(str(_layout.project_layout.VERIFY/'fixtures/chemcam')) if getattr(o,'passive',False))
    with patch.object(evidence,'enrich',side_effect=AssertionError('被动谱不许逐发打分')):
        p=console.plot_for({'observations':[passive]},REGISTRY['chemcam'],0,0,-1,{})
    check(not p['annotate']['element_evidence'],'被动谱不能生成LIBS逐发统计')

def task_c(temp):
    folder=temp/'notes';folder.mkdir();b=REGISTRY['chemcam'];obs=next(o for o in b.scan_folder(str(_layout.project_layout.VERIFY/'fixtures/chemcam')) if not getattr(o,'passive',False));pt=obs.points[0]
    check(notes.read(folder,b.id,obs,pt)['point'] is None and not (folder/'_notes').exists(),'未标注读取不能创建文件或伪造状态')
    notes.save(folder,b.id,obs,pt,'point','consistent','<script>论文图3</script>')
    notes.save(folder,b.id,obs,pt,'observation','pending','观测待查')
    saved=notes.read(folder,b.id,obs,pt)
    check(saved['point']['status']=='consistent' and saved['observation']['status']=='pending','观测和点标注不得互相覆盖')
    worker="from instruments import notes;from types import SimpleNamespace;import sys;obs=SimpleNamespace(sol='13',seq='ccam02013',product='RDR');[notes.save(sys.argv[1],'chemcam',obs,int(sys.argv[2])*10+i,'point','pending','并发') for i in range(6)]"
    processes=[subprocess.Popen([sys.executable,'-c',worker,str(folder),str(i)],cwd=str(ROOT),stdout=subprocess.PIPE,stderr=subprocess.PIPE) for i in [10,11]]
    results=[p.communicate(timeout=15) for p in processes]
    for process,output in zip(processes,results):
        if process.returncode:print(output[1].decode('utf-8','replace')[-1800:])
    keys=notes.load(folder)['annotations']
    check(all(p.returncode==0 for p in processes) and sum('并发'==r['note'] for r in keys.values())==12,'两个工作台进程并发保存不能丢失彼此标注')
    check(console.obs_listing([obs],str(folder),b.id)[0]['label'].startswith('● '),'有点或观测标注时下拉必须标记')
    check(notes.point_key(b.id,obs,pt)==notes.point_key(b.id,obs,obs.points[0]),'标注身份不能依赖下拉排序索引')
    before=(folder/'_notes/notes.json').read_bytes()
    try:notes.save(folder,b.id,obs,pt,'point','unknown','x');check(False,'无效标注须拒绝')
    except ValueError:check((folder/'_notes/notes.json').read_bytes()==before,'无效标注须拒绝并保持原文件')
    r=report.build_report(b,[obs],str(temp/'annotated'),data_dir=str(folder))
    row=next(csv.DictReader((temp/'annotated'/r['csv'][0]).open(encoding='utf-8-sig')))
    check(row['论文对照']=='与论文一致' and row['观测论文对照']=='待查','报告CSV必须回填两级标注')
    page=Path(r['index']).read_text(encoding='utf-8')
    check('&lt;script&gt;' in page and '<script>论文' not in page,'用户备注在HTML中必须转义')
    class Scan:
        sol='293';sequence_id='srlc11360';sclk='123';n_points=1
    scan=Scan();notes.save(folder,'sherloc',scan,0,'point','inconsistent','人工对照备注')
    result={'point_detail':[{'point':0,'classes':[],'pos':None,'labeled':False,'n_bands':0,'fwhm':None,'snr':None,'sigma':1.,'groups':[]}]}
    notes.attach_sherloc(result,folder,scan);header,rows=qc.point_table(result)
    check(rows[0][header.index('是否一致')]=='与论文不一致' and rows[0][header.index('论文备注')]=='人工对照备注','SHERLOC逐点表须使用同一标注存储')
    check('人工对照备注' in workbench.paper_html(result) and workbench.paper_html({})=='','SHERLOC报告展示人工对照，空标注不增加伪内容')
    state={'app':{'chemcam':{'folder':str(folder),'observations':[obs]}}}
    for handler,url in [(console.make_handler(state),'/api/paper?app=chemcam&scan=0&point=0'),(view.make_handler({'folder':str(folder),'scans':[scan],'listing':[{'label':'test'}]}),'/api/paper?scan=0&point=0')]:
        server=ThreadingHTTPServer(('127.0.0.1',0),handler);thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start()
        try:
            endpoint='http://127.0.0.1:%d'%server.server_port+url
            body=json.dumps({'scope':'point','status':'pending','note':'HTTP往返中文'}).encode('utf-8')
            req=urllib.request.Request(endpoint,data=body,headers={'Content-Type':'application/json'},method='POST')
            posted=json.load(urllib.request.urlopen(req));got=json.load(urllib.request.urlopen(endpoint))
            check(posted['point']['note']==got['point']['note']=='HTTP往返中文','两工作台POST/GET需按真实HTTP往返中文标注')
        finally:server.shutdown();thread.join();server.server_close()
    file=folder/'_notes/notes.json';file.write_text('{broken',encoding='utf-8')
    try:notes.save(folder,b.id,obs,pt,'point','pending','x');check(False,'损坏存储不得被覆盖')
    except ValueError:check(file.read_text()=='{broken','损坏存储须报错并保留原件')

def csv_table(path):
    rows=list(csv.reader(Path(path).open(encoding='utf-8-sig',newline='')))
    metadata={r[0][2:]:json.loads(r[1]) for r in rows if r[0].startswith('# ')}
    table=[r for r in rows if not r[0].startswith('# ')]
    return metadata,table[0],np.array([[float(v) if v else np.nan for v in r] for r in table[1:]])

def task_d(temp):
    b=REGISTRY['chemcam'];obs=next(o for o in b.scan_folder(str(_layout.project_layout.VERIFY/'fixtures/chemcam')) if not getattr(o,'passive',False))
    exported=spectrum_export.from_backend(b,obs,0,2,{'dust':True},temp/'csv');meta,header,table=csv_table(exported['path'])
    raw=b.load_point(obs.points[0],2,{'dust':True})
    check(Path(exported['path']).read_bytes().startswith(b'\xef\xbb\xbf'),'CSV必须UTF8-BOM供Excel读取')
    check(meta['shot']==8 and header[2]=='原始第6发' and header[-1]=='原始第30发','当前发号与全部逐发列均使用原始编号')
    check(np.array_equal(table[:,1],raw['y']) and np.array_equal(table[:,2:],raw['shot_matrix']),'当前谱和各发必须无显示舍入地精确导出')
    check(meta['preprocessing']=={'dust':True} and len(meta['sources'])==2,'头部需保留处理开关与原件SHA')
    second=spectrum_export.from_backend(b,obs,0,None,{},temp/'csv')
    check(second['path']!=exported['path'] and Path(exported['path']).exists(),'连续导出不能覆盖先前文件')
    passive=next(o for o in b.scan_folder(str(_layout.project_layout.VERIFY/'fixtures/chemcam')) if getattr(o,'passive',False))
    result=spectrum_export.from_backend(b,passive,0,None,{},temp/'csv');meta,header,_=csv_table(result['path'])
    check(len(header)==2 and meta['shot']=='not_applicable','被动谱导出不能伪造逐发列')
    class Scan:sol='293';sequence_id='srlc11360';sclk='123';n_points=1;parts=[]
    y=np.arange(7,dtype=float)[None,:]+.123456789123
    pipeline={'r1':y,'r2':None,'r3':None,'shift':np.linspace(200,4200,7),'log':[],'product':'RRS','wl_source':'test'}
    with patch.object(view.preprocess,'run_pipeline',return_value=pipeline):
        result=view.export_current_csv(Scan(),{},0,'r1',temp/'csv')
    _,header,table=csv_table(result['path'])
    check(header[1]=='r1' and np.array_equal(table[:,1],y[0]),'SHERLOC导出必须直接取管线全精度，不复用两位小数显示JSON')
    try:spectrum_export.write(temp/'invalid','bad',[1,2],[('y',[1])],{});check(False,'错位列应拒绝')
    except ValueError:check(not (temp/'invalid').exists(),'错位列须在写文件前拒绝')

def task_b():
    script=reference.JS+'''
var original={x:[1,2,3],series:[{key:'r1',label:'R1',y:[10,20,30]}],xlabel:'nm',product:'RRS'};
var ref=freezeReference(original,'first');original.series[0].y[0]=999;
var kept=ref.series[0].y[0]===10;
var matched=referenceSelection(ref,{x:[1,2,3],series:[{key:'r1'}],xlabel:'nm',product:'RRS'});
var missing=referenceSelection(ref,{x:[1,2,3],series:[{key:'r2'}],xlabel:'nm',product:'RRS'});
var units=referenceSelection(ref,{x:[1,2,3],series:[{key:'r1'}],xlabel:'channel',product:'RRS'});
var read=referenceRead(ref,ref.series[0],2.1);
console.log(JSON.stringify([kept,matched.series.length===1,missing.series.length===0&&missing.message.indexOf('r2')>=0&&units.series.length===0,read.x===2&&read.y===20&&referenceRead(ref,ref.series[0],7)===null]));
'''
    p=subprocess.run([NODE,'-'],input=script,text=True,capture_output=True,check=True)
    for ok,name in zip(json.loads(p.stdout),['参考捕获须深拷贝防后续数据改变','同键同轴参考可叠加','缺键或异轴参考须隐藏并说明','悬停按参考自身轴读数且不越界外推']):check(ok,name)

def task_e():
    x=np.linspace(1.30,2.6,500);cont=1+.2*x;y=cont*(1-.2*np.exp(-.5*((x-1.93)/.025)**2))
    rows=absorption.payload(x,y,'synthetic')['references'];water=next(r for r in rows if r['key']=='water19')
    check(water['detected'] and abs(water['depth']-.2)<.005,'已知20%倒峰须恢复连续谱相对深度')
    check(not any(r['detected'] for r in absorption.payload(x,np.ones_like(x))['references']),'平谱不能报吸收')
    check(not any(r['detected'] for r in absorption.payload(x,cont)['references']),'纯斜率不能报吸收')
    emission=cont*(1+.2*np.exp(-.5*((x-1.93)/.025)**2))
    check(not next(r for r in absorption.payload(x,emission)['references'] if r['key']=='water19')['detected'],'发射峰不能被当吸收倒峰')
    check(not rows[-1]['available'] and '窗口外' in rows[-1]['reason'],'Fe²⁺约1µm位于CP2窗口外须明确标记')
    scaled=next(r for r in absorption.payload(x,3*y)['references'] if r['key']=='water19')
    check(abs(scaled['depth']-water['depth'])<1e-12,'相对吸收深度须对整体强度尺度不变')
    b=REGISTRY['supercam'];obs=next(o for o in b.scan_folder(str(_layout.project_layout.VERIFY/'fixtures/supercam')) if o.product=='CP2');pt=obs.points[0]
    default=b.load_point(pt,None,{})
    check(not default['absorption']['enabled'] and not any(r['detected'] for r in default['absorption']['references']),'默认DN只能画参考，不能冒充反射评分')
    d=b.load_point(pt,None,{'reflectance':True,'baseline':True});plain=b.load_point(pt,None,{'reflectance':True})
    check(np.array_equal(d['y'],plain['y']) and d['meta']['spectral_column'].lower()=='i_f_atm','反射模式必须保留原产品I/F而非全局扣背景结果')
    raw=fitsio.read_bintable(next(iter(pt.files.values())),'SPECTRA')['I_F_atm'];axis=fitsio.read_bintable(next(iter(pt.files.values())),'WAVELENGTH')['Wavelength (um)']
    check(np.array_equal(d['y'],raw[np.argsort(axis)]),'反射强度与波长必须按相同顺序配对')
    short=next(r for r in absorption.payload(x[x<1.98],y[x<1.98])['references'] if r['key']=='water19')
    check(not short['available'],'缺右连续谱不能给出水带评分')

def task_j(temp):
    folder=temp/'rmi/sol01000';folder.mkdir(parents=True);name='cr0_486263086prc_f0481570ccam02000l1'
    original=np.linspace(1000,50000,128*128).astype(np.uint16).reshape(128,128)
    Image.fromarray(original).save(folder/(name+'.tif'))
    (folder/(name+'.lbl')).write_bytes((_layout.project_layout.VERIFY/'fixtures/chemcam/sol01000'/ (name+'.lbl')).read_bytes())
    b=REGISTRY['chemcam'];images=rmi.inventory(b,str(folder))
    check(len(images)==1 and images[0]['sol']=='1000' and images[0]['seq']=='ccam02000','RMI清单须从本图标签识别Sol/序列')
    d=rmi.details(images[0],'486263230')
    check(d['marker'] is None and not d['special_raw'] and d['delta_seconds']==-144,'本图无SPECIAL时不得借用LIBS坐标画十字')
    import io
    with Image.open(io.BytesIO(rmi.preview(str(folder/(name+'.tif'))))) as image:
        check(image.size==(128,128) and np.ptp(np.asarray(image))>200,'16位TIF预览须显示拉伸，不能convertL全白')
    with Image.open(folder/(name+'.tif')) as image:check(np.array_equal(np.asarray(image),original),'浏览不得更改原TIF像素')

def task_f():
    b=REGISTRY['supercam'];observations=b.scan_folder(str(_layout.project_layout.VERIFY/'fixtures/supercam'))
    with patch.object(b,'load_point',side_effect=AssertionError('元数据过滤不能读谱')):
        all_rows=filtering.select(b,observations,'',{})
        by_target=filtering.select(b,observations,'sharp',{})
        by_sol=filtering.select(b,observations,'Sol 700',{})
    check(len(all_rows['indices'])==3,'空搜索不计算元素且保留全部观测')
    check(len(by_target['indices'])==2 and len(by_sol['indices'])==2,'目标名和Sol过滤须只读元数据')
    filtering._CACHE.clear();original=b.load_point
    with patch.object(b,'load_point',side_effect=original) as loaded:
        first=filtering.select(b,observations,'Fe',{'dust':True});second=filtering.select(b,observations,'Fe',{'dust':True})
        check(first['indices']==second['indices']==[2] and loaded.call_count==1,'元素过滤惰性计算并复用同内容/处理配置缓存')
    check(all(observations[i].product=='CL1' for i in first['indices']),'被动谱不得因参考带文字匹配元素')
    with patch.object(b,'load_point',side_effect=ValueError('bad product')):
        failed=filtering.select(b,observations,'Fe',{'baseline':True})
    check(not failed['indices'] and failed['unresolved']==[2] and failed['errors'],'读取失败必须保留未判定状态，不当阴性')
    check(filtering.select(b,observations,'Fe',{},lambda:True)['cancelled'],'过期查询须可停止，不能回填新搜索')

def task_h(temp):
    original={'mode':'r1','custom':42};config=original
    for i in range(6):config=history.remember(config,str(temp/('folder%d'%i)))
    check(len(history.rows(config))==5 and original=={'mode':'r1','custom':42},'最近目录上限5且不污染调用方设置对象')
    config=history.remember(config,str(temp/'folder3'))
    check(config['recent_folders'][0]==str(temp/'folder3') and len(config['recent_folders'])==5,'重复目录应去重并移到首位')
    profile=temp/'settings.json';profile.write_text(json.dumps({'custom':42}),encoding='utf-8')
    history.store_folder(profile,str(temp),'chemcam');history.store_folder(profile,str(temp/'other'),'supercam')
    saved=json.loads(profile.read_text());check(saved['custom']==42 and len(saved['recent_folders'])==2,'各仪器最近目录分开并保留不相关设置')
    test=history.NAV_JS+'''console.log(JSON.stringify([navigationCommand('PageDown','',900,1,3,99,0),navigationCommand('Home','',900,1,3,99,0),navigationCommand('End','',900,1,3,99,0),navigationCommand('2','1',20,0,1,99,1),navigationCommand('2','1',20,0,1,99,0),navigationCommand('9','99',20,0,1,99,1)]));'''
    process=subprocess.run([NODE,'-'],input=test,text=True,capture_output=True,check=True)
    result=json.loads(process.stdout)
    check([r['index'] for r in result[:3]]==[2,0,99],'翻观测和Home/End必须使用可见选项及首末点')
    check(result[3]['index']==11 and result[4]['index']==12 and result[5]['index'] is None,'数字连续输入须区分界面点号基准且拒绝越界')

def task_i(temp):
    rows=[{'elements':['Fe'],'n_shots':25,'all_shot_support':{'Fe':10},'moc':{'rows':[{'oxide':'SiO2','wt':50.},{'oxide':'FeOT','wt':20.}]}},
          {'elements':[],'n_shots':5,'all_shot_support':{'Fe':0},'moc':None}]
    s=observation_summary.summarize(rows,'supercam')
    check(s['top'][0]['numerator']==10 and s['top'][0]['denominator']==30,'观测检出率须合计发数，不能平均各点百分比')
    check(s['official']['SiO2']['mean']==50 and s['official']['SiO2']['available_points']==1,'官方成分缺值不补零且保留覆盖数')
    apxs=[{'elements':['Fe'],'apxs':{'rows':[{'oxide':'FeO','wt':4}]}},{'elements':[],'apxs':{'rows':[{'oxide':'FeO','wt':12}]}}]
    a=observation_summary.summarize(apxs,'apxs')
    check(a['official']['FeO']['mean']==8 and a['top'][0]['unit']=='点','APXS使用各自配对成分均值与点出现率')
    b=REGISTRY['chemcam'];mean=observation_summary.MeanSpectrum();grid=np.array([1.,1.5,2.,10.,10.5,11.])
    mean.add({'x':grid,'y':np.ones(6)},b);mean.add({'x':[1.1,1.9,10.1,10.9],'y':[3.,3.,3.,3.]},b)
    result,note=mean.result()
    check(np.array_equal(result[1],[1,2,1,1,2,1]) and '不跨空档' in note,'平均谱须分段对齐且不外推覆盖边缘')
    mean.add({'x':grid,'y':np.ones(6),'ylabel':'另一种单位'},b)
    check(mean.result()[0] is None,'不同强度单位禁止混合成平均谱')
    cp3=observation_summary.MeanSpectrum();native={'x':[1.,2.,2.,3.],'y':[1.,2.,4.,3.]};cp3.add(native,b)
    check(np.array_equal(cp3.result()[0][1],native['y']),'同一原生重复波长网格须保留不同通道，不能抹掉CP3接缝')
    obs=next(o for o in b.scan_folder(str(_layout.project_layout.VERIFY/'fixtures/chemcam')) if not getattr(o,'passive',False))
    r=report.build_report(b,[obs],str(temp/'overview'));page=Path(r['index']).read_text(encoding='utf-8')
    check('observationSummary' in page and 'data-sort="elements"' in page and '观测平均谱' in page and len(report._spectrum_svg(b,np.arange(6144),np.arange(6144),w=260,h=85).encode())<30000,'报告须有可排序汇总、平均缩略图且SVG增量不超预算')

def task_g(temp):
    check(download_queue.normalize(['sol_00001','sol_00001','sol00002'])==['sol_00001','sol00002'],'多Sol队列保持选择顺序并去重')
    try:download_queue.normalize('../escape');check(False,'Sol路径越界须拒绝')
    except ValueError:check(True,'Sol路径越界须在下载前拒绝')
    job={'cancel':False};seen=[]
    def one(sol,out):
        seen.append(sol)
        if sol=='sol00002':raise IOError('mock failure')
        Path(out).mkdir(parents=True);(Path(out)/'done.bin').write_bytes(b'ok');job.update(done=1,total=1)
    download_queue.run(job,['sol00001','sol00002','sol00003'],temp/'queue',one)
    check(len(seen)==3 and job['failed_sols']==['sol00002'] and len(job['completed_sols'])==2 and not job['running'],'单Sol失败须继续后面的Sol并正确汇总')
    cancelled={'cancel':False};seen=[]
    def stopping(sol,out):
        seen.append(sol);Path(out).mkdir(parents=True);(Path(out)/'kept.bin').write_bytes(b'ok');cancelled.update(done=1,total=1 if len(seen)==1 else 2)
        if len(seen)==2:cancelled['cancel']=True
    download_queue.run(cancelled,['sol00001','sol00002','sol00003'],temp/'cancelled',stopping)
    check(len(seen)==2 and cancelled['queue_status']=='cancelled' and cancelled['files_done']==2 and (temp/'cancelled/sol00001/kept.bin').exists(),'取消须停后续队列并保留已完成与当前段完成文件')
    b=REGISTRY['chemcam'];state={'dl':{'chemcam':{'cancel':False}}};name='cl5_398645626ccs_f0030004ccam02013p3.csv'
    entries=[{'name':name,'size':4,'url':'https://example.invalid/mock.csv'}]
    def mock_files(files,dest,job):
        Path(dest).mkdir(parents=True,exist_ok=True)
        for f in files:(Path(dest)/f['name']).write_bytes(b'mock')
        job.update(done=len(files),total=len(files),bytes_done=4*len(files),bytes_total=4*len(files),log=['mock transport'])
    with patch.object(console.pds_listing,'list_files',return_value=entries),patch.object(console.ibase,'download_files',side_effect=mock_files):
        console._dl_worker(state,b,['sol00010','sol00013'],{'CL5_CCS'},str(temp/'common'))
    check(len(state['dl']['chemcam']['completed_sols'])==2 and (temp/'common/sol00010'/name).exists() and (temp/'common/sol00013'/name).exists(),'通用台实际worker须分目录执行多Sol，避免同名覆盖')
    import io
    class Response(io.BytesIO):
        def __init__(self,declared=4):super().__init__(b'mock');self.headers={'Content-Length':str(declared)}
    files=[{'name':'mock.csv','size':4,'ext':'csv','url':'https://example.invalid/mock.csv'},{'name':'mock.xml','size':4,'ext':'xml','url':'https://example.invalid/mock.xml'}]
    info={'groups':[{'key':'RRS','files':files}]};state={'dl':{'cancel':False}}
    with patch.object(view,'list_sol_files',return_value=info),patch.object(urllib.request,'urlopen',side_effect=lambda *a,**k:Response()):
        view._dl_worker(state,'processed',['sol_00001','sol_00002'],{'RRS'},str(temp/'sherloc'))
    check(len(state['dl']['completed_sols'])==2 and (temp/'sherloc/sol_00002/mock.xml').exists(),'SHERLOC实际worker须逐Sol下载并保留配套XML')
    data={'cancel':False}
    from instruments import base
    with patch.object(urllib.request,'urlopen',side_effect=lambda *a,**k:Response(9)),patch('time.sleep'):
        base.download_files(entries,str(temp/'short'),data)
    check(data['done']==0 and data.get('error') and not (temp/'short'/name).exists() and any('短读' in s for s in data['log']),'短读重试后仍失败时不得发布截断文件或记为完成')

def main():
    with tempfile.TemporaryDirectory(prefix='sherlock-phase4-') as tmp:task_a(Path(tmp));task_c(Path(tmp));task_d(Path(tmp));task_b();task_e();task_j(Path(tmp));task_f();task_h(Path(tmp));task_i(Path(tmp));task_g(Path(tmp))
    print('共 %d 项断言，%d 项失败'%(TOTAL,FAIL));return int(FAIL>0)
if __name__=='__main__':raise SystemExit(main())
