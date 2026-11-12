"""Single runtime path/version authority. No writes to the installed application."""
from __future__ import annotations
import hashlib
import json
import logging
from logging.handlers import RotatingFileHandler
import os
from pathlib import Path
import shutil
import sys
import tempfile
import traceback

VERSION = '0.17.0'
PRODUCT = 'MarsWorkbench'
import project_layout
SOURCE = project_layout.APP
_DEFAULT_SOURCE = SOURCE


def app_dir():
    return Path(sys.executable).resolve().parent if getattr(sys,'frozen',False) else SOURCE


def resources():
    return project_layout.resources() if SOURCE==_DEFAULT_SOURCE else SOURCE/'resources'


def data_home():
    override = os.environ.get('MARSWORKBENCH_DATA_HOME')
    if override:
        root = Path(override).expanduser().resolve()
    elif (app_dir()/'portable.flag').is_file():
        root = app_dir()/'userdata'
    elif os.name=='nt':
        root = Path(os.environ.get('LOCALAPPDATA',str(Path.home()/'AppData/Local')))/PRODUCT
    else:
        root = Path.home()/'.marsworkbench'
    root.mkdir(parents=True,exist_ok=True)
    for name in ('settings','cache','logs','reports_selfcheck'):(root/name).mkdir(exist_ok=True)
    return root


def location(group, *parts):
    root=data_home()/group
    root.mkdir(parents=True,exist_ok=True)
    return root.joinpath(*parts)


def settings(name):
    path=location('settings',name+'.json')
    legacy=app_dir()/('.console_settings.json' if name=='console' else '.gui_settings.json' if name=='sherloc' else '.'+name+'_settings.json')
    if not path.exists() and legacy.is_file():
        # Never remove or rewrite the original. Malformed legacy config remains recoverable.
        try:
            value=json.loads(legacy.read_text(encoding='utf-8'))
            if isinstance(value,dict):
                atomic_json(path,value)
                logger().info('migrated configuration copy %s -> %s',legacy,path)
        except (OSError,ValueError):
            logger().warning('legacy configuration unreadable; retained %s',legacy)
    return str(path)


def atomic_json(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    fd,temp=tempfile.mkstemp(dir=path.parent,suffix='.part')
    try:
        with os.fdopen(fd,'w',encoding='utf-8') as f:json.dump(value,f,ensure_ascii=False,indent=2)
        os.replace(temp,path)
    finally:
        if os.path.exists(temp):os.unlink(temp)


def cache_write(relative):
    path=location('cache',relative)
    path.parent.mkdir(parents=True,exist_ok=True)
    return path


def cache_read(relative):
    new=location('cache',relative)
    if new.exists():return new
    old=SOURCE/'.cache'/relative
    indexed=project_layout.legacy_cache(relative)
    return old if old.exists() else indexed if indexed and indexed.exists() else new


def copy_cached(relative):
    """Seed a mutable cache from legacy and bundled originals without overwriting it."""
    source=SOURCE/'.cache'/relative;dest=cache_write(relative)
    if source.is_dir():
        for old in source.rglob('*'):
            if old.is_file():
                new=dest/old.relative_to(source)
                if not new.exists():new.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(old,new)
    elif source.is_file() and not dest.exists():shutil.copy2(source,dest)
    elif not source.exists():
        for logical,canonical in project_layout.cache_entries(relative):
            suffix=logical[len(str(relative).replace('\\\\','/').rstrip('/')):].lstrip('/')
            target=dest/suffix if suffix else dest
            if not target.exists() and canonical.is_file():target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(canonical,target)
    bundled=resources()/'cache-seed'/relative
    if bundled.is_dir():
        for original in bundled.rglob('*'):
            if original.is_file():
                target=dest/original.relative_to(bundled)
                if not target.exists():target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(original,target)
    elif bundled.is_file() and not dest.exists():shutil.copy2(bundled,dest)
    return dest


def empty_data(instrument):
    return str(project_layout.instrument_data(instrument))


def build_id():
    info=resources()/'build-info.json'
    if getattr(sys,'frozen',False) and info.exists():
        try:return VERSION+'+'+json.loads(info.read_text(encoding='utf-8'))['build']
        except (OSError,ValueError,KeyError):pass
    h=hashlib.sha256(VERSION.encode())
    for logical,path in project_layout.runtime_sources():
        h.update(logical.encode());h.update(path.read_bytes())
    return VERSION+'+'+h.hexdigest()[:12]


def logger():
    log=logging.getLogger('marsworkbench')
    target=location('logs','app.log')
    if not any(getattr(h,'baseFilename',None)==str(target) for h in log.handlers):
        handler=RotatingFileHandler(target,maxBytes=1024*1024,backupCount=4,encoding='utf-8')
        handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s %(message)s'))
        log.addHandler(handler);log.setLevel(logging.INFO)
    return log


def crash(exc):
    text=''.join(traceback.format_exception(type(exc),exc,exc.__traceback__))
    try:
        path=location('logs','crash.txt');path.write_text(text,encoding='utf-8')
        logger().error('startup failed',exc_info=(type(exc),exc,exc.__traceback__))
        return '启动未完成，请查看日志：'+str(path)
    except OSError:return '启动未完成，用户数据目录不可写；请改用可写目录或取消便携模式。'


def clean_report(text):
    """Reports travel with data: no developer-home paths or development-doc links."""
    import re
    for path in (str(SOURCE),str(SOURCE.parent),str(Path.home())):
        for variant in {path,path.replace('\\','/')}:
            text=text.replace(variant,'[本地数据]')
    text=re.sub(r'(?i)\bhref=([\"\'])([^\"\']*(?:04_使用文档|file://)[^\"\']*)\1','',text)
    text=re.sub(r'(?<![A-Za-z0-9])[A-Za-z]:[\\/]', '[数据路径]/',text)
    if '</body>' in text.lower():
        from product_theme import report_style
        text=report_style(text)
    return text
