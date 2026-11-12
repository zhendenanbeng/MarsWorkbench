"""Full precision, provenance-carrying current-spectrum CSV export."""
import csv
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import re
import numpy as np

def sources(paths):
    return [{'path':str(Path(p).resolve()),'sha256':hashlib.sha256(Path(p).read_bytes()).hexdigest()} for p in sorted(set(str(p) for p in paths))]

def write(out,stem,x,columns,metadata,shots=None,shot_numbers=None):
    x=np.asarray(x,float)
    if x.ndim!=1 or not len(x) or not np.isfinite(x).all():raise ValueError('导出坐标必须为非空有限一维数组')
    arrays=[np.asarray(y,float) for _,y in columns]
    if not arrays or any(y.shape!=x.shape for y in arrays):raise ValueError('导出强度列须与当前横轴等长')
    labels=[name for name,_ in columns]
    if shots is not None:
        shots=np.asarray(shots,float)
        if shots.ndim!=2 or shots.shape[0]!=len(x) or len(shot_numbers or [])!=shots.shape[1]:raise ValueError('逐发矩阵与原始发号不一致')
        arrays.extend(shots.T);labels.extend('原始第%d发'%n for n in shot_numbers)
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    stem=re.sub(r'[^\w.-]','_',stem)[:160]
    path=out/(stem+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%f')+'.csv')
    meta=dict(metadata,schema_version=1,generated_utc=datetime.now(timezone.utc).isoformat(),nonfinite_intensity='空单元格；保留行位置，禁止补零',
              processing_note='当前平均谱先平均后扣背景；各发列单独扣背景，非线性预处理后其平均未必等于当前平均谱')
    with path.open('x',encoding='utf-8-sig',newline='') as stream:
        writer=csv.writer(stream)
        for key,value in meta.items():writer.writerow(['# '+key,json.dumps(value,ensure_ascii=False,allow_nan=False)])
        writer.writerow([meta.get('xlabel','坐标'),*labels])
        for i,v in enumerate(x):writer.writerow([float(v),*(float(a[i]) if np.isfinite(a[i]) else '' for a in arrays)])
    return {'path':str(path.resolve()),'rows':len(x),'columns':len(arrays)+1,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}

def from_backend(backend,obs,index,shot,steps,out):
    if not 0<=index<len(obs.points):raise ValueError('点号无效')
    pt=obs.points[index];d=backend.load_point(pt,shot,steps)
    numbers=d.get('shot_numbers') or []
    original_shot=(numbers[shot] if shot is not None and 0<=shot<len(numbers) else 'mean') if not d.get('passive') else 'not_applicable'
    metadata={'instrument':backend.id,'sol':obs.sol,'sequence':obs.seq,'product':obs.product,'point_index':index,
              'sclk':pt.sclk,'shot':original_shot,'preprocessing':steps,'xlabel':d.get('xlabel') or backend.xlabel,
              'ylabel':d.get('ylabel') or backend.ylabel,'sources':sources(pt.files.values()),'product_metadata':d.get('meta',{}),'processing_log':d.get('log',[])}
    return write(out,'%s_sol%s_%s_%s_%s_pt%d'%(backend.id,obs.sol,obs.seq,obs.product,pt.sclk,index),d['x'],[('当前谱',d['y'])],metadata,d.get('shot_matrix'),numbers)

HTML='<button id="csvBtn">导出当前谱 CSV</button><span id="csvMessage" class="sub"></span>'
JS='''
$('csvBtn').addEventListener('click',function(){var b=this;b.disabled=true;$('csvMessage').textContent='正在导出…';fetch('/api/csv?'+qs()+'&out='+encodeURIComponent($('outDir').value)).then(function(r){return r.json();}).then(function(d){b.disabled=false;$('csvMessage').textContent=d.error?d.error:('已导出 '+d.rows+' 行：'+d.path);}).catch(function(e){b.disabled=false;$('csvMessage').textContent='导出失败：'+e;});});
'''
