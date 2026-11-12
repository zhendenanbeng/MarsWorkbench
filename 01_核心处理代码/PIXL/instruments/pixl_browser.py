"""Local PIXL RQC/RXL browsing. No cross-instrument matching or new quantification."""
import csv
import hashlib
import html
import io
import math
from pathlib import Path
import re
import xml.etree.ElementTree as ET
import numpy as np
from PIL import Image, ImageDraw
from sherloc import pixl

ELEMENTS=pixl.ELEMENTS
NOTE='官方RQC成分直接展示，不重新定量、不判矿物；err保留原列定义。热图只画实测点，不插值、不补零、不跨仪器配准，无物理比例尺。'

def finite(value):
    try:
        number=float(value)
        return number if math.isfinite(number) and number>=0 else None
    except (TypeError,ValueError):return None

def info(path):
    return pixl.parse_lid('urn:nasa:pds:mars2020_pixl:local:'+path.stem.lower())

def scan(folder):
    if not str(folder).strip():raise ValueError('请先选择本地数据目录或加载PIXL示例')
    root=Path(folder).expanduser().resolve()
    if not root.is_dir():raise ValueError('请选择包含PIXL RQC/RXL CSV及XML的本地目录')
    groups={};warnings=[]
    for path in sorted(root.rglob('*')):
        if not path.is_file() or path.suffix.lower()!='.csv' or not path.resolve().is_relative_to(root):continue
        if any(p.startswith(('.', '_')) or p in ('下载','userdata') for p in path.relative_to(root).parts[:-1]):continue
        meta=info(path)
        if not meta or meta['type'] not in ('rqc','rxl'):continue
        logical=path.stem.lower()[:-2]
        version=int(path.stem[-2:]) if path.stem[-2:].isdigit() else 0
        groups.setdefault(logical,{}).setdefault(version,[]).append((path,meta))
    products={}
    for logical,by_version in groups.items():
        top=max(by_version)
        for version,entries in by_version.items():
            if len(entries)>1 and len({entry[0].read_bytes() for entry in entries})>1:
                # 同版本异字节：若恰为被选用的最高版本则无法确定取舍，拒绝；低版本冲突按版本策略忽略但必须留痕
                if version==top:raise ValueError('同名同版本产品字节冲突：'+entries[0][0].name)
                warnings.append('发现同名同版本异字节副本，已按最高版本选择，冲突件：'+entries[0][0].name)
        chosen=sorted(by_version[top])[0]
        products[logical]=(top,chosen[0],chosen[1])
    rows=list(products.values());observations=[]
    for _,path,meta in rows:
        if meta['type']!='rqc':continue
        matches=[p for _,p,m in rows if m['type']=='rxl' and m['pos']==meta['pos'] and m['seq'][:9]==meta['seq'][:9]]
        observations.append({'id':hashlib.sha256(str(path.relative_to(root)).encode()).hexdigest()[:16],
            'file':path.name,'path':path,'sol':int(meta['sol']),'rtt':meta['seq'][:9],'placement':meta['pos'],
            'rxl':sorted(matches),'warning':'' if len(matches)==1 else '缺RXL坐标表，成分仍可浏览' if not matches else '多份RXL候选，请明确选择坐标表'})
    if not observations:warnings.append('未找到RQC成分表；只有RXL不能显示成分。采集sol与成分处理sol可能不同，请打开包含两者的上级目录。')
    return root,sorted(observations,key=lambda o:(o['sol'],o['file'])),warnings

def chemistry(path):
    # Keep missing/error values separate; duplicate PMCs must never be silently overwritten.
    result={}
    with path.open(encoding='utf-8-sig',newline='') as stream:
        reader=csv.DictReader(stream)
        if not reader.fieldnames or 'PMC' not in reader.fieldnames:raise ValueError('成分表缺PMC列')
        available=[e for e in ELEMENTS if e+'_wt%' in reader.fieldnames]
        if not available:raise ValueError('不支持的成分表：缺官方wt%列')
        for row in reader:
            try:pmc=int(row['PMC'])
            except (TypeError,ValueError):raise ValueError('成分表有无效PMC')
            if pmc in result:raise ValueError('成分表PMC重复：'+str(pmc))
            result[pmc]={e:{'value':finite(row.get(e+'_wt%')),'error':finite(row.get(e+'_wt%_err'))} for e in available}
    return result,available

def metadata(path):
    xml=path.with_suffix('.xml')
    if not xml.exists():return {},'缺XML，err定义尚未核实'
    try:
        root=ET.parse(xml).getroot();fields={}
        for node in root.iter():
            if node.tag.rsplit('}',1)[-1]!='Field_Delimited':continue
            row={c.tag.rsplit('}',1)[-1]:(c.text or '') for c in node}
            fields[row.get('name','')]=row.get('description','')
        definitions=[v for k,v in fields.items() if k.endswith('_err')]
        error='官方1σ绝对误差（wt%）' if definitions and all('1-sigma' in s.lower() and 'absolute' in s.lower() for s in definitions) else '官方err列（wt%，σ定义未核实）'
        return fields,error
    except ET.ParseError:return {},'XML无法解析，err定义尚未核实'

def read_coordinates(path):
    raw=pixl.read_rxl(str(path));seen=set()
    # MCC references include both white-light pcw and dark/other-illumination pcd products.
    xml=path.with_suffix('.xml')
    if xml.exists():
        tree=ET.parse(xml).getroot()
        raw['mcc_tifs']=sorted({(n.text or '').strip() for n in tree.iter() if n.tag.rsplit('}',1)[-1]=='file_name'
            and (n.text or '').strip().lower().startswith(('pcw_','pcd_','pc__'))
            and (n.text or '').strip().lower().endswith(('.tif','.tiff','.img'))})
    with path.open(encoding='utf-8-sig',newline='') as stream:
        for row in csv.DictReader(stream):
            try:pmc=int(row['PMC'])
            except (TypeError,ValueError):raise ValueError('坐标表有无效PMC')
            if pmc in seen:raise ValueError('坐标表PMC重复：'+str(pmc))
            seen.add(pmc)
    return raw

def frames(root,rxl,raw):
    output=[]
    for name in raw['mcc_tifs']:
        match=re.search(r'n\d{7}(\d{13})[0-9_][a-z0-9_]{2}j',name.lower())
        if not match:continue
        key='PMC_%04d_MCC'%int(match[1][-4:])
        if key not in raw['images']:continue
        # Only exact referenced stems, including the cache's IMG alias. No nearest-time guess.
        candidates=[p for p in root.rglob('*') if p.is_file() and p.suffix.lower() in ('.img','.tif','.tiff') and p.stem.lower()==Path(name).stem.lower() and p.resolve().is_relative_to(root)]
        if not candidates:output.append({'name':name,'key':key,'path':None,'missing':True});continue
        candidate=next((p for p in candidates if p.parent==rxl.parent),sorted(candidates)[0])
        if len(candidates)>1 and len({hashlib.sha256(p.read_bytes()).hexdigest() for p in candidates})>1:
            output.append({'name':name,'key':key,'path':None,'missing':True,'reason':'同名MCC字节冲突'});continue
        output.append({'name':candidate.name,'key':key,'path':candidate,'missing':False})
    return output

def load(root,observation,rxl_index=None):
    chem,elements=chemistry(observation['path']);_,error_label=metadata(observation['path'])
    available=observation['rxl']
    if rxl_index is None and len(available)==1:rxl_index=0
    if rxl_index is not None and not 0<=rxl_index<len(available):raise ValueError('坐标表索引无效')
    rxl=available[rxl_index] if rxl_index is not None else None
    raw=read_coordinates(rxl) if rxl else {'pts':{},'mcc_tifs':[],'images':[]}
    images=frames(root,rxl,raw) if rxl else []
    return {'observation':observation,'chem':chem,'elements':elements,'error_label':error_label,
        'rxl_path':rxl,'raw':raw,'frames':images}

def canvas_data(data,frame=0,axis='label',origin=0):
    if axis not in ('label','ij') or origin not in (0,1):raise ValueError('坐标约定无效')
    images=data['frames'];image=None;key=None;warning=[]
    if images:
        if not 0<=frame<len(images):raise ValueError('图像索引无效')
        item=images[frame];key=item['key']
        if item['path']:
            try:image=pixl.mcc_gray(str(item['path']))
            except (OSError,ValueError,IndexError) as exc:warning.append('MCC读取失败，成分表仍可使用：'+str(exc))
        else:warning.append(item.get('reason','缺该帧MCC本地原图'))
    if image is None:warning.append('无MCC底图：保留成分表，不生成虚构点阵或猜图')
    points=[];outside=0
    if image is not None:
        h,w=image.shape
        for pmc,coordinates in data['raw']['pts'].items():
            pair=coordinates.get(key)
            if not pair or not all(math.isfinite(v) for v in pair):continue
            # RXL XML explicitly declares i=line/y, j=sample/x. Alternate is manual, disclosed.
            x,y=(pair[1],pair[0]) if axis=='label' else pair
            x-=origin;y-=origin
            if not (0<=x<w and 0<=y<h):outside+=1;continue
            points.append({'pmc':pmc,'x':x,'y':y,'values':data['chem'].get(pmc,{})})
    if outside:warning.append('%d个坐标在图像外，已排除绘制（原值保留）'%outside)
    return image,points,warning

def composition_stats(data,salt_free=None):
    """逐元素统计（全部来自归档逐点定量值）。

    均值是**逐点算术平均**，不是任务方报告的 bulk sum（对全部点光谱求和后拟合一次）；
    贫钠镁铁质岩上两种估计量的 Na2O 可差 ~0.8–1.0 wt%（2026-11-06 论文对账：
    Alfalfa 富钠目标零差，Dourbes/Guillaumes/Quartier/Máaz 系统偏低）。
    salt_free 如 {'SO3':1.0,'Cl':1.0}：按阈值逐点剔除后另给均值与点数；**仅剔除**，
    不做文献配方（Farley 2022 Methods）的扣 Ca/Na 与归一化。
    """
    if salt_free is not None:
        if not isinstance(salt_free,dict) or not salt_free or not all(e in data['elements'] and isinstance(v,(int,float)) and v>=0 for e,v in salt_free.items()):raise ValueError('盐剔除阈值无效')
    keys=list(data['chem'])
    keep=None if salt_free is None else [all((data['chem'][p].get(e,{}).get('value') is None) or data['chem'][p][e]['value']<=t for e,t in salt_free.items()) for p in keys]
    out={}
    for e in data['elements']:
        values=[data['chem'][p][e]['value'] for p in keys]
        v=np.array([x for x in values if x is not None])
        entry={'n':int(v.size),'mean':float(v.mean()),'median':float(np.median(v)),'std':float(v.std(ddof=1)) if v.size>1 else 0.0}
        if keep is not None:
            w=np.array([x for x,k in zip(values,keep) if x is not None and k])
            entry['salt_free_n']=int(w.size)
            entry['salt_free_mean']=float(w.mean()) if w.size else None
        out[e]=entry
    return out

STATS_NOTE='统计为档案逐点算术均值，非任务方 bulk sum（求和谱拟合）；Na2O 在贫钠镁铁质岩上两者可有 ~0.8–1.0 wt% 系统差（估计量差异，非数据缺陷）。'

def payload(data,frame=0,axis='label',origin=0):
    image,points,warnings=canvas_data(data,frame,axis,origin)
    return {'elements':data['elements'],'initial_pmc':min(data['chem']) if data['chem'] else None,'pmcs':sorted(set(data['chem'])|set(data['raw']['pts'])),
        'frames':[{'name':f['name'],'key':f['key'],'missing':f['missing']} for f in data['frames']],
        'points':points,'width':image.shape[1] if image is not None else 752,'height':image.shape[0] if image is not None else 580,
        'has_image':image is not None,'warnings':warnings,'chemical_points':len(data['chem']),
        'coordinate_points':len(data['raw']['pts']),'matched_points':len(set(data['chem'])&set(data['raw']['pts'])),
        'error_label':data['error_label'],'rxl_file':data['rxl_path'].name if data['rxl_path'] else None,'note':NOTE}

# Compatibility entry points; actual graphical output is implemented in the GUI layer.
def color(value,lo,hi):
    from .pixl_output import color as implementation
    return implementation(value,lo,hi)

def render(data,element,frame=0,axis='label',origin=0):
    from .pixl_output import render as implementation
    return implementation(data,element,frame,axis,origin)

def export(data,out,element,frame=0,axis='label',origin=0):
    from .pixl_output import export as implementation
    return implementation(data,out,element,frame,axis,origin)
