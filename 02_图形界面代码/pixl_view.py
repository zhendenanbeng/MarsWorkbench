"""PIXL adapter for the common workbench workflow, hosted by the launcher."""

import source_paths  # source roots; no-op for bundled imports
import hashlib
import json
import os
from pathlib import Path
import shutil
import threading
import urllib.parse
import uuid
import xml.etree.ElementTree as ET
import appenv
from instruments import history
from instruments import pixl_browser as backend
from instruments import pixl_download,pixl_workflow
from instruments.pixl_workflow_ui import PAGE

# /api/pixl/sfs 采样谱缓存：key=(扫描目录, 文件名/大小/mtime 元组, 轴参数)。
# 全扫描 CSV 重读约 4s 且页面每次加载都调用——缓存后重复请求毫秒级。
_SFS_SAMPLE_CACHE = {'key': None, 'spectrum': None}

def state_for(state):
    if 'pixl' not in state:
        try:folder=json.loads(Path(appenv.settings('pixl')).read_text(encoding='utf-8')).get('folder','')
        except (OSError,ValueError):folder=''
        import project_layout
        folder=str(project_layout.resolve_legacy(folder)) if folder else ''
        state['pixl']={'folder':folder,'observations':[],'loaded':{},'exports':{},'report':pixl_workflow.ReportJob(),'piquant':pixl_workflow.PiquantJob(),'lock':threading.RLock()}
    return state['pixl']

def load_folder(state,path):
    root,observations,warnings=backend.scan(path)
    with state['lock']:state.update(folder=str(root),observations=observations,loaded={})
    history.store_folder(appenv.settings('pixl'),str(root))
    from sherloc.desktop import remember_instrument
    remember_instrument('pixl',appenv.settings('console'))
    for ob in observations:
        try:ob['title']=next((n.text or '' for n in ET.parse(ob['path'].with_suffix('.xml')).getroot().iter() if n.tag.rsplit('}',1)[-1]=='title'),'')
        except (OSError,ET.ParseError):ob['title']=''
    return {'folder':str(root),'observations':[{k:([p.name for p in v] if k=='rxl' else v) for k,v in ob.items() if k!='path'} for ob in observations],'warnings':warnings}

def dataset(state,p):
    ob=next((o for o in state['observations'] if o['id']==p.get('obs')),None)
    if ob is None:raise ValueError('观测已变更或不存在，请重新打开数据目录')
    rxl=int(p['rxl']) if p.get('rxl') not in ('',None) else None
    key=(ob['id'],rxl)
    paths=[ob['path'],*ob['rxl']];signature=tuple((p.stat().st_mtime_ns,p.stat().st_size) for p in paths)
    if key not in state['loaded'] or state['loaded'][key][0]!=signature:state['loaded'][key]=(signature,backend.load(Path(state['folder']),ob,rxl))
    return state['loaded'][key][1]

def demo():
    source=appenv.resources()/'pixl-demo';manifest=json.loads((source/'manifest.json').read_text(encoding='utf-8'))
    import project_layout
    dest=project_layout.instrument_data('pixl','工作示例')
    for row in manifest['files']:
        path=(source/row['file']).resolve();target=(dest/row['file']).resolve()
        if not path.is_relative_to(source.resolve()) or not target.is_relative_to(dest.resolve()):raise ValueError('示例清单路径无效')
        raw=path.read_bytes()
        if hashlib.sha256(raw).hexdigest()!=row['sha256']:raise ValueError('PIXL示例SHA校验失败')
        if target.exists() and target.read_bytes()!=raw:raise ValueError('PIXL示例副本已修改，请保留文件并改用自己的数据目录：'+target.name)
        if not target.exists():shutil.copy2(path,target)
    return str(dest)

def install_demo(destination,confirm=False):
    import samples
    root=Path(destination or samples.default_destination()).expanduser().resolve()
    dest=root/samples._relative_destination(root,'pixl')
    if dest.is_relative_to(appenv.app_dir().resolve()) and not dest.is_relative_to(appenv.data_home().resolve()):raise ValueError('请选择软件目录之外的示例位置')
    source=appenv.resources()/'pixl-demo';rows=json.loads((source/'manifest.json').read_text(encoding='utf-8'))['files']
    existing=[r for r in rows if (dest/r['file']).exists()]
    if existing and not confirm:return {'requires_confirmation':True,'existing_files':len(existing),'destination':str(root)}
    for row in rows:
        original=(source/row['file']).resolve();target=(dest/row['file']).resolve()
        if not original.is_relative_to(source.resolve()) or not target.is_relative_to(dest.resolve()):raise ValueError('示例文件路径无效')
        samples._copy(original,target,row['sha256'])
    return {'requires_confirmation':False,'files':len(rows),'destination':str(root),'folders':{'pixl':str(dest)}}

def route(handler,state):
    url=urllib.parse.urlparse(handler.path)
    if url.path!='/pixl' and not url.path.startswith(('/api/pixl/','/data/pixl/')):return False
    st=state_for(state)
    def send(value,ctype='application/json; charset=utf-8'):
        body=json.dumps(value,ensure_ascii=False,allow_nan=False).encode() if isinstance(value,dict) else value
        handler._send(200,ctype,body);return True
    if url.path=='/pixl':
        import product_ui
        return send(product_ui.inject(PAGE.encode(),'pixl',state.get('launcher_url','/')),'text/html; charset=utf-8')
    if url.path == '/api/pixl/sfs':
        import project_layout
        axis_p = project_layout.data_root() / 'PIXL' / '缓存' / 'sfs_energy_axis.json'
        if not axis_p.is_file():
            return send({'available': False,
                         'reason': 'SFS 能量轴未生成：先下载 sol 的 SFS 数据并运行谱层自证'})
        doc = json.loads(axis_p.read_text(encoding='utf-8'))
        # 谱采样（能量轴自证后的 mean spectrum，~400 点供 SVG 渲染）；
        # 扫描目录按轴 doc 记录优先（Bellegarde 全扫描集 sfs_sol0186），
        # 旧安装回退 sfs_sol187 尾部。
        # 重读全扫描 CSV 约 4s/次且页面每次加载都调本接口——按文件指纹
        # 缓存采样结果（轴参数进 key，轴重冻结或扫描变更自动失效）。
        spectrum = None
        scan_used = None
        try:
            import numpy as np
            from instruments import sfs_spectra
            sfs_dir = None
            for cand in (doc.get('scan_dir'), 'sfs_sol0186', 'sfs_sol187'):
                if not cand:
                    continue
                d = project_layout.data_root() / 'PIXL' / '下载' / cand
                if d.is_dir():
                    sfs_dir = d
                    scan_used = cand
                    break
            if sfs_dir is None:
                sfs_dir = project_layout.data_root() / 'PIXL' / '下载' / 'sfs_sol187'
            files = sorted(sfs_dir.glob('*.csv'))
            axis = doc.get('axis', {})
            key = (str(sfs_dir),
                   tuple((f.name, f.stat().st_size, f.stat().st_mtime_ns)
                         for f in files),
                   axis.get('keV_per_ch'), axis.get('zero_keV'))
            if _SFS_SAMPLE_CACHE.get('key') == key:
                spectrum = _SFS_SAMPLE_CACHE.get('spectrum')
            else:
                specs = []
                for f in files:
                    sp, _ = sfs_spectra.read_sfs(f)
                    if sp is not None:
                        specs.append(sp)
                if specs:
                    mean = np.vstack(specs).mean(axis=0)
                    kev = (axis['keV_per_ch'] * np.arange(len(mean))
                           + axis['zero_keV'])
                    step = max(1, len(mean) // 400)
                    sel = slice(0, len(mean), step)
                    spectrum = {'kev': [round(float(v), 3) for v in kev[sel]],
                                'counts': [float(v) for v in mean[sel]],
                                'ka_lines': {k: v for k, v in
                                             sfs_spectra.KA_KEV.items()
                                             if axis['zero_keV'] <= v <= kev[-1]}}
                _SFS_SAMPLE_CACHE['key'] = key
                _SFS_SAMPLE_CACHE['spectrum'] = spectrum
        except Exception:
            spectrum = None
        ratios = []
        try:
            mat_p = (project_layout.data_root() / 'PIXL' / '缓存'
                     / 'sfs_ratio_matrix.json')
            if mat_p.is_file():
                mat = json.loads(mat_p.read_text(encoding='utf-8'))
                ratios = [{'pair': r['pair'],
                           'corr': round(float(r['loglog_corr']), 3)}
                          for r in (mat.get('ratios') or [])[:6]]
        except Exception:
            ratios = []
        # PIQUANT 官方开源定量层（第三条对账线）：只读缓存状态，不在此
        # 起 subprocess——批量对账走 数据准备脚本/02_piquant_quantify.py
        piquant = {'available': False}
        try:
            from instruments import piquant_bridge
            pq_status = piquant_bridge.status()
            if pq_status.get('available'):
                piquant = {'available': True, 'bin': pq_status.get('bin', ''),
                           'n_quantified': None, 'majors': [], 'scope': None}
                pq_p = (project_layout.data_root() / 'PIXL' / '缓存'
                        / 'piquant_vs_rqc.json')
                if pq_p.is_file():
                    pq = json.loads(pq_p.read_text(encoding='utf-8'))
                    majors = ('FeO-T', 'SiO2', 'CaO', 'Al2O3')
                    stats = pq.get('stats', {})
                    piquant.update({
                        'n_quantified': pq.get('n_quantified'),
                        'scan_dir': pq.get('scan_dir'),
                        'majors': [{'oxide': o,
                                    'median_dev_pct': round(
                                        stats[o]['median_rel_dev_pct'], 1),
                                    'corr': round(stats[o]['corr'], 3)}
                                   for o in majors if o in stats],
                        'scope': pq.get('scope'),
                    })
            else:
                piquant = {'available': False, 'reason': pq_status.get('reason')}
        except Exception:
            pass
        return send({'available': True, 'axis': doc.get('axis'),
                     'mode': doc.get('mode'), 'scan_dir': scan_used or doc.get('scan_dir'),
                     'ca_corr': doc.get('ca_corr'), 'fe_corr': doc.get('fe_corr'),
                     'n_pmc': doc.get('n_pmc'), 'scope': doc.get('scope'),
                     'ratios': ratios, 'spectrum': spectrum, 'piquant': piquant})
    if url.path.startswith('/data/pixl/'):
        parts=urllib.parse.unquote(url.path).split('/')
        folder=st['exports'].get(parts[3] if len(parts)>=5 else '')
        if not folder:raise FileNotFoundError('导出文件链接已失效，请从数据目录打开或重新导出')
        path=(folder/Path(*parts[4:])).resolve()
        if not path.is_relative_to(folder.resolve()) or not path.is_file() or path.suffix.lower() not in ('.csv','.png','.html','.json'):raise FileNotFoundError('没有该导出文件')
        import mimetypes
        return send(path.read_bytes(),mimetypes.guess_type(path)[0] or 'application/octet-stream')
    action=url.path.rsplit('/',1)[-1];p={k:v[0] for k,v in urllib.parse.parse_qs(url.query).items()}
    if action in ('load','demo','export','export-selected','openfolder','openreport','start','cancel'):
        if handler.command!='POST':raise ValueError('此操作需POST')
        size=int(handler.headers.get('Content-Length','0'))
        if not 0<size<=65536:raise ValueError('请求大小无效')
        p=json.loads(handler.rfile.read(size))
    if url.path.startswith('/api/pixl/download/'):
        with st['lock']:
            if 'downloader' not in st:st['downloader']=pixl_download.Controller()
            controller=st['downloader']
        if action=='catalog':return send(controller.query(p.get('refresh')=='1'))
        if action=='status':return send(dict(controller.snapshot(),destination=pixl_download.destination()))
        if action=='start':
            result=controller.start(p.get('catalog_id'),p.get('scan_ids'),p.get('destination'),bool(p.get('images',True)))
            history.update(appenv.settings('pixl'),lambda config:dict(config,download_dest=p['destination']))
            return send(result)
        if action=='cancel':return send(controller.cancel())
        raise FileNotFoundError('PIXL下载接口不存在')
    if action=='init':return send({'folder':st['folder'],'cache':backend.pixl.CACHE,'download_root':appenv.empty_data('pixl'),'preferences':history.read(appenv.settings('pixl'))})
    if action=='recent':return send({'folders':history.rows(history.read(appenv.settings('pixl')))})
    if action=='load':return send(load_folder(st,p.get('path','')))
    if action=='demo':return send({'folder':demo()})
    if url.path.startswith('/api/pixl/report/'):
        job=st['report']
        if action=='start':
            with st['lock']:
                if p.get('obs') and not any(o['id']==p['obs'] for o in st['observations']):raise ValueError('观测已变化，请重新选择后生成报告')
                return send(job.start(st['folder'],list(st['observations']),pixl_workflow.settings(st['folder'])['report_out'],p.get('obs'),int(p['rxl']) if p.get('rxl') not in ('',None) else None,int(p.get('frame',0))))
        if action=='cancel':return send(job.cancel())
        if action=='status':
            result=job.snapshot()
            if result.get('folder'):
                with st['lock']:
                    folder=Path(result['folder'])
                    if st.get('report_folder')!=folder:
                        token=uuid.uuid4().hex;st['exports'][token]=folder;st.update(report_folder=folder,report_token=token)
                    result['url']='/data/pixl/'+st['report_token']+'/index.html'
            return send(result)
        raise ValueError('报告操作无效')
    if url.path.startswith('/api/pixl/piquant/'):
        job=st['piquant']
        if action=='start':
            return send(job.start(p.get('scan_dir'),p.get('rqc') or None,int(p.get('stride',10))))
        if action=='cancel':return send(job.cancel())
        if action=='status':return send(job.snapshot())
        raise ValueError('PIQUANT 对账操作无效')
    if action=='openreport':
        result=st['report'].snapshot();folder=result.get('folder')
        if result.get('running') or not folder or not (Path(folder)/'index.html').is_file():raise ValueError('报告尚未生成')
        if os.environ.get('MARSWORKBENCH_HEADLESS')!='1':
            from sherloc.desktop import open_app_window
            open_app_window((Path(folder)/'index.html').as_uri(),True)
        return send({'ok':True})
    if action=='openfolder':
        kind=p.get('kind')
        if kind=='download':
            folder=st['downloader'].snapshot().get('folder') if st.get('downloader') else None
            target=Path(folder) if folder else None
        elif kind=='report':target=Path(st['report'].snapshot().get('folder','')) if st['report'].snapshot().get('folder') else None
        elif kind=='export':target=st['exports'].get(p.get('token'))
        else:target=None
        if target is None or not target.is_dir():raise ValueError('尚无已生成的输出目录')
        if os.environ.get('MARSWORKBENCH_HEADLESS')!='1' and hasattr(os,'startfile'):os.startfile(str(target))
        return send({'ok':True})
    if action=='browse':
        root=Path(p.get('path') or Path.home()).expanduser().resolve()
        if not root.is_dir():raise ValueError('目录不存在')
        return send({'path':str(root),'parent':str(root.parent),'dirs':[{'name':v.name,'path':str(v)} for v in sorted(root.iterdir()) if v.is_dir() and not v.name.startswith(('.', '$'))][:300]})
    with st['lock']:
        data=dataset(st,p);frame=int(p.get('frame',0))
        if action=='observation':return send(dict(backend.payload(data,frame),frame=frame))
        if action=='point':
            pmc=int(p.get('pmc',-1))
            if pmc not in data['chem'] and pmc not in data['raw']['pts']:raise ValueError('该PMC不存在')
            return send({'pmc':pmc,'present':pmc in data['chem'],'values':data['chem'].get(pmc,{})})
        if action=='image':
            image,_,_=backend.canvas_data(data,frame)
            if image is None:raise FileNotFoundError('缺MCC原图')
            import io
            stream=io.BytesIO();backend.Image.fromarray(image).save(stream,format='PNG');return send(stream.getvalue(),'image/png')
        if action=='export-selected':
            defaults=pixl_workflow.settings(st['folder']);kind=p.get('kind')
            out=p.get('out') or defaults['png_out' if kind=='png' else 'csv_out']
            folder,name=pixl_workflow.export_single(data,kind,out,p.get('element','points'),frame,int(p['pmc']) if p.get('pmc') is not None else None,p.get('scope','current'))
            token=uuid.uuid4().hex;st['exports'][token]=folder
            return send({'folder':str(folder),'token':token,'url':'/data/pixl/'+token+'/'+name})
        if action=='export':
            root=Path(st['folder']).resolve();out=(root/'_pixl_exports').resolve()
            if not out.is_relative_to(root) or (out.is_relative_to(appenv.app_dir().resolve()) and not out.is_relative_to(appenv.data_home().resolve())):raise PermissionError('导出需选择软件目录外的可写数据目录')
            token=uuid.uuid4().hex;folder=out/(data['observation']['id']+'_'+token[:8]);result=backend.export(data,folder,p.get('element','points'),frame)
            st['exports'][token]=folder
            names={'成分CSV':result['csv'],'坐标CSV':'coordinates.csv','来源SHA':'sources.json','报告':result['report']}
            if result['map']:names['热图PNG']=result['map']
            return send({'folder':str(folder),'links':{name:'/data/pixl/'+token+'/'+file for name,file in names.items()},'warnings':result['warnings']})
    raise FileNotFoundError('PIXL接口不存在')
