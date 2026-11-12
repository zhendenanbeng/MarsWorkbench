"""ChemCam CL1/CL9 PSV single-column spectra and the official channel wavelength map."""
from __future__ import annotations
import hashlib
import io
import os
from pathlib import Path
import re
import tempfile
import numpy as np
from . import pds_listing

WAVE_URL=('https://pds-geosciences.wustl.edu/msl/msl-m-chemcam-libs-4_5-rdr-v1/'
          'mslccm_1xxx/document/ccam_default_wave.tab')
WAVE_SHA256='d161471c3fc8f35dbcd19579242177502b1e0bf0c79791df55d9736ec4ef6e90'
CHANNELS=6444

def cache_path():
    import appenv
    return appenv.copy_cached('ccam_default_wave.tab')

def wavelengths(raw):
    if hashlib.sha256(raw).hexdigest()!=WAVE_SHA256:
        raise ValueError('被动谱波长表与已验证的PDS版本不一致；不能猜测波长轴')
    x=np.loadtxt(io.BytesIO(raw),ndmin=2)
    if x.shape!=(CHANNELS,1) or not np.isfinite(x).all():
        raise ValueError('被动谱波长表须为6444行、1列且数值有限')
    valid=x[:,0]>0
    if valid.sum()!=6144:raise ValueError('波长表有效通道数量不符')
    return x[:,0],valid

def ensure_wavelengths(path=None):
    path=Path(path) if path is not None else cache_path()
    if path.exists():
        wavelengths(path.read_bytes())
        return path
    if path!=cache_path():raise ValueError('指定的被动谱波长表不存在')
    raw=pds_listing.http_get(WAVE_URL,timeout=45)
    wavelengths(raw)  # 校验完成后才进入缓存；失败不覆盖任何已有表。
    path.parent.mkdir(parents=True,exist_ok=True)
    fd,name=tempfile.mkstemp(prefix='chemcam-wave-',suffix='.part',dir=path.parent)
    try:
        with os.fdopen(fd,'wb') as stream:stream.write(raw)
        os.replace(name,path)
    finally:
        if os.path.exists(name):os.remove(name)
    return path

def read_spectrum(path,label_path,wavelength_path=None):
    path,label_path=Path(path),Path(label_path)
    if not label_path.is_file():raise ValueError('被动谱需要配套LBL以核验表偏移、行数和单位')
    text=label_path.read_text(encoding='ascii')
    def integer(key,block=text):
        match=re.search(r'^\s*'+re.escape(key)+r'\s*=\s*(\d+)\s*$',block,re.M)
        if not match:raise ValueError('被动谱LBL缺少或不支持字段 '+key)
        return int(match.group(1))
    table=re.search(r'OBJECT\s*=\s*TABLE\s*\n(.*?)END_OBJECT\s*=\s*TABLE',text,re.S)
    pointer=re.search(r'\^TABLE\s*=\s*\("([^"]+)"\s*,\s*(\d+)\s*\)',text)
    if not table or not pointer or Path(pointer.group(1)).name.lower()!=path.name.lower():
        raise ValueError('被动谱TABLE指针与数据文件不匹配或格式未支持')
    block=table.group(1)
    if integer('ROWS',block)!=CHANNELS or integer('COLUMNS',block)!=1:
        raise ValueError('本版仅支持6444行、单列平均谱PSV')
    if not re.search(r'PRODUCT_TYPE\s*=\s*"CHEMCAM-PSV"',text):
        raise ValueError('不是CHEMCAM-PSV产品')
    unit=re.search(r'UNIT\s*=\s*"?([A-Z_]+)',block)
    if not unit or unit.group(1)!='DIGITAL_NUMBER':raise ValueError('尚未支持此PSV强度单位')
    record_bytes=integer('RECORD_BYTES');start=int(pointer.group(2));records=integer('FILE_RECORDS')
    raw=path.read_bytes()
    if start<1 or len(raw)!=record_bytes*records or records-start+1!=CHANNELS:
        raise ValueError('PSV文件长度、记录偏移或行数与LBL不一致')
    data=np.loadtxt(io.BytesIO(raw[(start-1)*record_bytes:]),ndmin=2)
    if data.shape!=(CHANNELS,1) or not np.isfinite(data).all():
        raise ValueError('PSV数值表不是有限的6444行单列数据')
    wave_path=ensure_wavelengths(wavelength_path)
    x,valid=wavelengths(wave_path.read_bytes())
    # 原始通道顺序是VNIR/VIS/UV，每段2148含两端各50个无效通道。
    # 必须同一掩码筛波长和强度，再用同一顺序排序；绝不能截取前6144行。
    order=np.argsort(x[valid],kind='stable');x=x[valid][order];y=data[:,0][valid][order]
    return x,y,{'wavelength_table':str(wave_path),'wavelength_sha256':WAVE_SHA256,
                'input_channels':CHANNELS,'valid_channels':len(x),'discarded_channels':int((~valid).sum()),
                'unit':'DN','axis_source':WAVE_URL}
