"""Verified examples copied out of installed resources; never run analyses in resources."""

import source_paths  # source roots; no-op for bundled imports
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile
import appenv

NAMES={'sherloc':'SHERLOC拉曼','chemcam':'ChemCam九点光栅','apxs':'APXS成分','supercam':'SuperCam Bellegarde'}

def manifest():
    return json.loads((appenv.resources()/'sample-data/manifest.json').read_text(encoding='utf-8'))

def documents():
    if os.name=='nt':
        try:
            import ctypes
            buffer=ctypes.create_unicode_buffer(32768)
            if ctypes.windll.shell32.SHGetFolderPathW(None,5,None,0,buffer)==0:return Path(buffer.value)
        except Exception:pass
    return Path.home()/'Documents'

def default_destination():
    import project_layout
    return project_layout.data_root()

def _relative_destination(root,relative):
    import project_layout
    parts=Path(relative).parts
    if root==project_layout.data_root().resolve() and parts and parts[0] in project_layout.GROUPS:
        return Path(project_layout.GROUPS[parts[0]])/'工作示例'/Path(*parts[1:])
    return Path(relative)

def _path(root,relative):
    root=Path(root).resolve();path=(root/relative).resolve()
    if not path.is_relative_to(root):raise ValueError('示例路径越出选定目录')
    return path

def _copy(source,dest,expected):
    raw=source.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=expected:raise ValueError('随包示例校验失败，请重新解压完整发行包：'+source.name)
    dest.parent.mkdir(parents=True,exist_ok=True)
    fd,temp=tempfile.mkstemp(dir=dest.parent,suffix='.part')
    try:
        with os.fdopen(fd,'wb') as f:f.write(raw)
        os.replace(temp,dest)
    finally:
        if os.path.exists(temp):os.unlink(temp)

def seed_cache():
    for row in manifest()['files']:
        if not row['path'].startswith('cache-seed/'):continue
        relative=row['path'][len('cache-seed/'):]
        if relative=='watson_inventory.csv':
            from sherloc import watson
            target=Path(watson._inventory_path())
        else:target=appenv.cache_write(relative)
        if not target.exists():_copy(_path(appenv.resources(),row['path']),target,row['sha256'])

def install(destination=None,confirm=False):
    root=Path(destination or default_destination()).expanduser().resolve()
    if root.is_relative_to(appenv.app_dir()) and not root.is_relative_to(appenv.data_home()):raise ValueError('请选择软件目录之外的示例位置，例如文档目录')
    rows=[r for r in manifest()['files'] if r['path'].startswith('sample-data/')]
    existing=[r for r in rows if _path(root,_relative_destination(root,r['path'][len('sample-data/'):])).exists()]
    if existing and not confirm:
        return {'requires_confirmation':True,'existing_files':len(existing),'destination':str(root)}
    for row in rows:
        relative=row['path'][len('sample-data/'):]
        _copy(_path(appenv.resources(),row['path']),_path(root,_relative_destination(root,relative)),row['sha256'])
    seed_cache()
    return {'requires_confirmation':False,'destination':str(root),'files':len(rows),
            'folders':{name:str(root/_relative_destination(root,name)) for name in NAMES}}

def welcome(mark=None):
    from instruments import history
    path=appenv.settings('launcher')
    if mark is not None:
        if mark not in ('skipped','complete','pending'):raise ValueError('向导状态无效')
        history.update(path,lambda data:dict(data,onboarding=mark))
    value=history.read(path).get('onboarding','pending')
    return {'state':value,'show':value=='pending'}
