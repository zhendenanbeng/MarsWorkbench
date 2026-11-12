"""Source roots, runtime assets and instrument data. No scientific computation."""
from pathlib import Path
import hashlib
import json
import os
import shutil
import sys

_SELF=Path(__file__).resolve()
ROOT=None
CONFIG={}
for _parent in _SELF.parents:
    _marker=_parent/'.project-layout.json'
    if _marker.is_file() and (_parent/'01_核心处理代码/共享/project_layout.py').resolve()==_SELF:
        ROOT=_parent
        CONFIG=json.loads(_marker.read_text(encoding='utf-8'))
        break
APP=ROOT/CONFIG['gui'] if ROOT else _SELF.parent
CORE=ROOT/CONFIG['core'] if ROOT else APP
VERIFY=ROOT/CONFIG['verification'] if ROOT else APP/'verify'
INSTALL=ROOT/CONFIG['install'] if ROOT else APP
MANUALS=ROOT/CONFIG['manuals'] if ROOT else APP/'resources/manuals'
GROUPS={'sherloc':'SHERLOC','supercam':'SuperCam','chemcam':'ChemCam','apxs':'APXS','pixl':'PIXL'}

def activate():
    if ROOT and not getattr(sys,'frozen',False):
        paths=[ROOT/p for p in CONFIG['namespace_roots']]+[VERIFY.parent,VERIFY]
        for path in reversed(paths):
            if str(path) not in sys.path:sys.path.insert(0,str(path))

def runtime_sources():
    """(logical import-compatible name, canonical file); no copied source tree."""
    roots=[ROOT/p for p in CONFIG['namespace_roots']] if ROOT else [APP]
    result={}
    for root in roots:
        for directory in ('','instruments','sherloc'):
            for path in sorted((root/directory).glob('*.py')):
                key=(Path(directory)/path.name).as_posix()
                if key in result and result[key]!=path:raise ValueError('重复运行模块：'+key)
                result[key]=path
    return sorted(result.items())

def source_hashes(include_verification=True):
    paths=[p for _,p in runtime_sources()]
    if include_verification and VERIFY.is_dir():paths.extend(VERIFY.glob('*.py'))
    if ROOT:paths.extend(ROOT/name for name in ('.project-layout.json','setup.py','pyproject.toml'))
    base=ROOT or APP
    return {p.relative_to(base).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(set(paths))}

def source_path(relative):
    return (ROOT or APP)/relative

def runtime_path(logical):
    value=dict(runtime_sources()).get(Path(logical).as_posix())
    if value is None:raise FileNotFoundError('运行模块不存在：'+str(logical))
    return value

def migration_context():
    if ROOT:return ROOT,CONFIG
    override=os.environ.get('MARSWORKBENCH_DATA_ROOT')
    if override:
        data=Path(override).expanduser().resolve();candidate=data.parent;marker=candidate/'.project-layout.json'
        try:
            config=json.loads(marker.read_text(encoding='utf-8'))
            if (candidate/config['data']).resolve()==data:return candidate,config
        except (OSError,ValueError,KeyError):pass
    return None,{}

def resolve_legacy(value):
    path=Path(value).expanduser()
    root,config=migration_context()
    if not root:return path
    if path.is_absolute():
        try:relative=path.resolve().relative_to(root).as_posix()
        except ValueError:return path
    else:relative=path.as_posix()
    aliases=dict(config.get('sources',{}),**config.get('legacy_data',{}))
    for old,new in sorted(aliases.items(),key=lambda row:-len(row[0])):
        if relative==old or relative.startswith(old+'/'):
            candidate=(root/(new+relative[len(old):])).resolve()
            if candidate.is_relative_to(root):return candidate
    return root/relative if not path.is_absolute() else path

def remap_preferences(value):
    if isinstance(value,dict):return {key:remap_preferences(item) for key,item in value.items()}
    if isinstance(value,list):return [remap_preferences(item) for item in value]
    root,_=migration_context()
    if root and isinstance(value,str) and Path(value).is_absolute():
        candidate=resolve_legacy(value)
        if candidate.exists():return str(candidate)
    return value

def data_root():
    override=os.environ.get('MARSWORKBENCH_DATA_ROOT')
    if override:return Path(override).expanduser().resolve()
    state=os.environ.get('MARSWORKBENCH_DATA_HOME')
    if state:return Path(state).expanduser().resolve()/'data'
    if ROOT and not getattr(sys,'frozen',False):return ROOT/CONFIG['data']
    # Installed applications keep user data writable and separate from binaries.
    from appenv import data_home
    return data_home()/'data'

def instrument_data(name,section='下载'):
    if name not in GROUPS:raise ValueError('未知仪器：'+name)
    path=data_root()/GROUPS[name]/section
    path.mkdir(parents=True,exist_ok=True)
    return path

def legacy_cache(relative):
    if ROOT:
        value=CONFIG.get('legacy_cache',{}).get(Path(relative).as_posix())
        if value:return ROOT/value
    return None

def cache_entries(relative):
    if not ROOT:return []
    prefix=Path(relative).as_posix().rstrip('/')
    return [(key,ROOT/value) for key,value in CONFIG.get('legacy_cache',{}).items() if key==prefix or key.startswith(prefix+'/')]

def legacy_cache_root():
    """Read-only test input view of formerly mixed caches; canonical data stays grouped."""
    if not ROOT:
        from appenv import data_home
        return data_home()/'cache'
    target=INSTALL/'构建缓存/回归缓存视图'
    inputs=[(key,ROOT/value) for key,value in CONFIG.get('legacy_cache',{}).items() if (ROOT/value).is_file()]
    signature=hashlib.sha256('\n'.join(key+'|'+str(path.stat().st_size)+'|'+str(path.stat().st_mtime_ns) for key,path in inputs).encode()).hexdigest()
    stamp=target/'.input-index.json'
    try:ready=json.loads(stamp.read_text(encoding='utf-8')).get('signature')==signature
    except (OSError,ValueError):ready=False
    if not ready:
        target.mkdir(parents=True,exist_ok=True)
        for relative,path in inputs:
            dest=target/relative;dest.parent.mkdir(parents=True,exist_ok=True)
            if not dest.exists() or dest.read_bytes()!=path.read_bytes():shutil.copy2(path,dest)
        stamp.write_text(json.dumps({'signature':signature}),encoding='utf-8')
    return target

def resources():
    if not ROOT or getattr(sys,'frozen',False):return APP/'resources'
    template=INSTALL/'资源模板'
    target=INSTALL/'构建缓存/运行资源'
    # Fingerprint canonical inputs. Runtime copies are disposable installation assets.
    files=[p for p in template.rglob('*') if p.is_file() and 'docs' not in p.relative_to(template).parts and p.name!='trust-report.txt']
    pairs=[(p,p.relative_to(template)) for p in files]
    for logical,relative in CONFIG.get('resource_sources',{}).items():
        p=ROOT/relative
        if p.is_file():pairs.append((p,Path(logical)))
    if MANUALS.is_dir():
        pairs.extend((p,Path('manuals')/p.relative_to(MANUALS)) for p in MANUALS.rglob('*') if p.is_file())
    signature=hashlib.sha256('\n'.join(str(p)+'|'+str(p.stat().st_size)+'|'+str(p.stat().st_mtime_ns)+'|'+str(rel) for p,rel in sorted(pairs)).encode()).hexdigest()
    stamp=target/'.asset-index.json'
    try:ready=json.loads(stamp.read_text(encoding='utf-8')).get('signature')==signature
    except (OSError,ValueError):ready=False
    if not ready:
        target.mkdir(parents=True,exist_ok=True)
        for path,relative in pairs:
            dest=target/relative;dest.parent.mkdir(parents=True,exist_ok=True)
            if not dest.exists() or dest.read_bytes()!=path.read_bytes():shutil.copy2(path,dest)
        manifest=target/'sample-data/manifest.json'
        if manifest.exists():
            rows=json.loads(manifest.read_text(encoding='utf-8'))['files']
            text='\n\n'.join('\n\n'.join(str(row.get(k,'')) for k in ('path','url','bundle','sha256')) for row in rows)
            (target/'sample-data/SOURCES.md').write_text(text+'\n',encoding='utf-8')
        stamp.write_text(json.dumps({'signature':signature},ensure_ascii=False),encoding='utf-8')
    return target

activate()
