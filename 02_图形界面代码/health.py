"""User-copyable diagnostics; optional features and network outages are warnings."""

import source_paths  # source roots; no-op for bundled imports
import importlib
import importlib.metadata
import platform
import tempfile
from pathlib import Path
import sys
import appenv
import app_network

def components():
    result={}
    for name,module,required in [('numpy','numpy',True),('Pillow','PIL',True),('opencv-python','cv2',False),('astropy','astropy',False)]:
        try:
            loaded=importlib.import_module(module)
            try:version=importlib.metadata.version(name)
            except importlib.metadata.PackageNotFoundError:version=str(getattr(loaded,'__version__','unknown'))
            result[name]={'available':True,'version':version,'required':required}
        except Exception:result[name]={'available':False,'version':None,'required':required}
    return result

def diagnose(port=None,network=True):
    warnings=[];home=None;writable=False;size=0
    try:
        home=appenv.data_home()
        with tempfile.TemporaryFile(dir=home) as probe:probe.write(b'write probe')
        writable=True
        for path in (home/'cache').rglob('*'):
            try:
                if path.is_file():size+=path.stat().st_size
            except OSError:pass
    except OSError:warnings.append('用户数据目录不可写；检查权限、磁盘空间或便携模式位置。')
    deps=components()
    for name,row in deps.items():
        if not row['available']:warnings.append(name+('缺失（必需）' if row['required'] else '未安装（可选；相关功能需全功能版）'))
    connection=app_network.probe(10) if network else {'ok':None,'message':'未检查PDS连接'}
    if connection['ok'] is False:warnings.append(connection['message'])
    result={'version':appenv.VERSION,'build':appenv.build_id(),'python':platform.python_version(),
            'platform':platform.platform(),'frozen':bool(getattr(sys,'frozen',False)),'port':port,'bind':'127.0.0.1',
            'data_home':str(home) if home else '不可用','writable':writable,'cache_bytes':size,
            'components':deps,'pds':connection,'warnings':warnings,
            'status':'error' if not writable or any(not r['available'] and r['required'] for r in deps.values()) else 'warning' if warnings else 'ok'}
    lines=['MarsWorkbench '+result['build'],'Python '+result['python']+' · '+result['platform'],
           '本机端口：'+str(port)+' / 127.0.0.1','用户数据：'+result['data_home'],
           '可写：'+str(writable)+'；新缓存：%.2f MB'%(size/1048576),
           '旧PIXL/WATSON缓存继续保留在用户目录，未计入新缓存大小。']
    for name,row in deps.items():lines.append(name+'：'+(row['version'] if row['available'] else '未安装')+('（可选）' if not row['required'] else ''))
    lines.append('PDS：'+connection['message']);lines.extend('提示：'+v for v in warnings)
    result['text']='\n'.join(lines)
    return result
