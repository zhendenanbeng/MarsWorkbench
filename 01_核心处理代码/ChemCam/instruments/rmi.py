"""Single-image ChemCam RMI context browsing. No cross-image laser projection."""
import io
from pathlib import Path
import re
import numpy as np
from PIL import Image
from . import pds3

RX=re.compile(r'^cr[01]_(\d{9})prc_.*?(ccam\d{5})[a-z]\d\.tif$',re.I)

def inventory(backend,root):
    rows=[]
    for path,name,dirname in backend.walk_files(root,r'.*\.tif$'):
        match=RX.match(name)
        if not match:continue
        label=Path(path).with_suffix('.lbl');data=pds3.parse_lbl(label)
        sol=data.get_int('PLANET_DAY_NUMBER')
        if sol is None:sol=backend.sol_from_dirname(dirname)
        if sol is None:continue
        rows.append({'path':path,'label':str(label),'sol':str(int(sol)),
                     'seq':str(data.get('SEQUENCE_ID') or match.group(2)).lower(),'sclk':match.group(1),'name':name})
    return sorted(rows,key=lambda r:(r['sclk'],r['name'],r['path']))

def details(image,clock):
    path=Path(image['label']);text=path.read_text(encoding='ascii') if path.exists() else ''
    state=re.search(r'GROUP\s*=\s*INSTRUMENT_STATE_PARMS(.*?)END_GROUP\s*=\s*INSTRUMENT_STATE_PARMS',text,re.S)
    raw={}
    for key in ['SPECIAL_LINE','SPECIAL_SAMPLE']:
        match=re.search(r'MSL:'+key+r'\s*=\s*([-+0-9.]+)',state.group(1) if state else '')
        if match:raw[key]=float(match.group(1))
    delta=int(image['sclk'])-int(clock)
    return {'name':image['name'],'sclk':image['sclk'],'delta_seconds':delta,
            'relationship':'同Sol/序列上下文；较当前点约%s%d秒，不代表实际LIBS光斑共帧'%('晚' if delta>=0 else '早',abs(delta)),
            'marker':None,'special_raw':raw,
            'marker_reason':'本图标签没有完整SPECIAL光斑坐标，未画十字' if len(raw)<2 else '本图给出坐标原值，但像素起算约定尚未核验，未画十字'}

def payload(obs,index):
    if not 0<=index<len(obs.points):raise ValueError('点号无效')
    return {'images':[dict(details(row,obs.points[index].sclk),index=i) for i,row in enumerate(getattr(obs,'rmi_images',[]))],
            'scope':'仅浏览同Sol/序列图像，不作跨图投影；不能把候选时间关联当光斑映射。'}

def preview(path):
    with Image.open(path) as source:
        # PIL.convert(L) clips most 16-bit RMI values to white. Stretch for display only.
        if len(source.getbands())==1:array=np.asarray(source,dtype=float)
        else:array=np.asarray(source.convert('L'),dtype=float)
    valid=np.isfinite(array)
    if not valid.any():raise ValueError('RMI没有可显示的有限像素')
    low,high=np.percentile(array[valid],[1,99])
    display=np.zeros(array.shape,dtype=np.uint8)
    if high>low:display[valid]=np.clip((array[valid]-low)/(high-low)*255,0,255).astype(np.uint8)
    else:display[valid]=128
    buffer=io.BytesIO();Image.fromarray(display).save(buffer,format='PNG')
    return buffer.getvalue()

