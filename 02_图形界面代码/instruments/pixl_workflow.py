"""PIXL UI workflows: separate exports and an all-observation data-check report."""
import copy
import csv
import hashlib
import html
import json
from pathlib import Path
import threading
import uuid
import appenv
from . import pixl_browser as backend

def writable(path):
    if not str(path).strip():raise ValueError('请选择可写输出目录')
    root=Path(path).expanduser().resolve()
    if root.is_relative_to(appenv.app_dir().resolve()) and not root.is_relative_to(appenv.data_home().resolve()):raise PermissionError('请选择软件目录外的可写输出目录')
    return root

def settings(folder):
    root=Path(folder)
    return {'png_out':str(root/'_png'),'csv_out':str(root/'_csv'),'report_out':str(root/'_qc/pixl')}

def source_manifest(data,folder,frame):
    paths=[data['observation']['path']]
    if data['rxl_path']:paths.append(data['rxl_path'])
    if data['frames'] and data['frames'][frame]['path']:paths.append(data['frames'][frame]['path'])
    paths+=[p.with_suffix('.xml') for p in list(paths) if p.with_suffix('.xml').exists()]
    (folder/'sources.json').write_text(json.dumps({'files':[{'file':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in paths],
        'error_label':data['error_label'],'axis':'sample=j,line=i','frame':frame,'note':backend.NOTE},ensure_ascii=False,indent=2),encoding='utf-8')

def export_single(data,kind,out,element='points',frame=0,pmc=None,scope='current'):
    if kind not in ('png','csv'):raise ValueError('导出类型无效')
    if element not in ('points',*data['elements']):raise ValueError('成分指标无效')
    if scope not in ('current','all'):raise ValueError('CSV范围无效')
    if frame<0 or (data['frames'] and frame>=len(data['frames'])) or (not data['frames'] and frame!=0):raise ValueError('图像索引无效')
    if kind=='png':raw=backend.render(data,element,frame)
    else:
        if scope=='current' and pmc not in data['chem']:raise ValueError('当前PMC没有RQC记录，不能导出伪零；可选择全部点CSV')
        rows={pmc:data['chem'][pmc]} if scope=='current' else data['chem']
    folder=writable(out)/(data['observation']['id']+'_'+uuid.uuid4().hex[:12]);folder.mkdir(parents=True)
    if kind=='png':name='map_'+element+'.png';(folder/name).write_bytes(raw)
    else:
        name='current_point.csv' if scope=='current' else 'composition.csv'
        with (folder/name).open('w',encoding='utf-8-sig',newline='') as stream:
            writer=csv.writer(stream);writer.writerow(['PMC']+[e+s for e in data['elements'] for s in ('_wt%','_wt%_err')])
            for key,values in sorted(rows.items()):writer.writerow([key]+[values[e][field] if values[e][field] is not None else '' for e in data['elements'] for field in ('value','error')])
    source_manifest(data,folder,frame)
    return folder,name

def check_summary(data):
    valid=sum(row['value'] is not None for values in data['chem'].values() for row in values.values())
    errors_missing=sum(row['value'] is not None and row['error'] is None for values in data['chem'].values() for row in values.values())
    return {'chemical_points':len(data['chem']),'coordinate_points':len(data['raw']['pts']),
        'matched_points':len(set(data['chem'])&set(data['raw']['pts'])),'valid_values':valid,
        'missing_values':len(data['chem'])*len(data['elements'])-valid,'missing_errors':errors_missing,
        'note':'仅统计归档字段与配套覆盖；检查通过不代表矿物鉴定、浓度准确率或定位精度。'}

def build_report(root,observations,out,current=None,rxl=None,frame=0,cancel=lambda:False,progress=lambda *args:None):
    folder=writable(out)/('report_'+uuid.uuid4().hex[:12]);folder.mkdir(parents=True)
    rows=[]
    for index,ob in enumerate(observations):
        if cancel():break
        entry={'file':ob['file'],'sol':ob['sol'],'id':ob['id']}
        try:
            data=backend.load(root,ob,rxl if ob['id']==current else None)
            selected=frame if ob['id']==current else 0
            sub=folder/ob['id'];backend.export(data,sub,'points',selected)
            summary=check_summary(data);entry.update(summary,report=ob['id']+'/index.html')
            subpage=sub/'index.html';text=subpage.read_text(encoding='utf-8')
            title='<h2>归档数据检查</h2><p>成分点 {chemical_points}；坐标点 {coordinate_points}；PMC交集 {matched_points}；有效成分值 {valid_values}；缺值 {missing_values}；有成分但缺err {missing_errors}。</p><p>{note}</p>'.format(**summary)
            subpage.write_text(text.replace('<table>',title+'<table>',1),encoding='utf-8')
        except Exception as exc:entry['error']=appenv.clean_report(str(exc))
        rows.append(entry);progress(index+1,len(observations),entry)
    esc=lambda s:html.escape(str(s))
    body=''.join('<tr><td>'+str(row['sol'])+'</td><td>'+esc(row['file'])+'</td><td>'+str(row.get('chemical_points','—'))+'</td><td>'+str(row.get('matched_points','—'))+'</td><td>'+('<a href="'+row['report']+'">打开</a>' if row.get('report') else esc(row.get('error','')))+'</td></tr>' for row in rows)
    title='已取消，以下仅为已完成部分' if cancel() else '本目录全部观测的归档数据检查'
    page='<!doctype html><html lang="zh"><meta charset="utf-8"><title>PIXL 数据质控报告</title><style>body{font:14px/1.6 sans-serif;margin:30px;color:#173147}td,th{border:1px solid #ccd;padding:8px}table{border-collapse:collapse}a{color:#126c93}</style><body><h1>PIXL 数据检查与成分报告</h1><p>'+title+'</p><p>'+esc(backend.NOTE)+'</p><p>多坐标候选的其他观测仅列成分，不自动猜坐标；当前观测使用明确选择的坐标表及MCC帧。字段覆盖检查不是准确率证明。</p><table><tr><th>Sol</th><th>产品</th><th>成分点</th><th>PMC交集</th><th>报告/读取错误</th></tr>'+body+'</table></body></html>'
    from product_theme import report_style
    (folder/'index.html').write_text(report_style(page),encoding='utf-8');(folder/'checks.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2),encoding='utf-8')
    return folder,rows

class ReportJob:
    def __init__(self):self.lock=threading.RLock();self.state={'running':False,'status':'idle'}
    def snapshot(self):
        with self.lock:return copy.deepcopy(self.state)
    def cancel(self):
        with self.lock:self.state['cancel']=True
        return self.snapshot()
    def start(self,root,observations,out,current,rxl,frame):
        with self.lock:
            if self.state['running']:raise ValueError('质控报告正在生成')
            if not observations:raise ValueError('当前目录没有RQC观测')
            writable(out)
            self.state={'running':True,'status':'running','cancel':False,'done':0,'total':len(observations),'source_folder':str(root)}
            self.thread=threading.Thread(target=self.run,args=(Path(root),copy.deepcopy(observations),out,current,rxl,frame),daemon=True);self.thread.start()
        return self.snapshot()
    def run(self,*args):
        try:
            def progress(done,total,row):
                with self.lock:self.state.update(done=done,total=total,current=row['file'])
            folder,rows=build_report(*args,cancel=lambda:self.snapshot().get('cancel',False),progress=progress)
            with self.lock:self.state.update(folder=str(folder),rows=rows,status='cancelled' if self.state['cancel'] else 'completed',errors=sum('error' in row for row in rows))
        except Exception as exc:
            with self.lock:self.state.update(status='failed',error=str(exc))
        finally:
            with self.lock:self.state['running']=False
