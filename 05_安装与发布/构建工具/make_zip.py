"""Build pinned onedir distributions. Publishing is deliberately not performed."""

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[2]/'02_图形界面代码'))
import source_paths as _layout
import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile

PROJECT=_layout.project_layout.ROOT
APP=_layout.project_layout.APP
sys.path.insert(0,str(APP))
import appenv

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def prepare(variant):
    if sys.version_info[:2]!=(3,12):raise RuntimeError('请使用Python3.12构建')
    resources=appenv.resources();resources.mkdir(exist_ok=True)
    snapshot=json.loads((resources/'_last_run.json').read_text(encoding='utf-8'))
    actual=_layout.project_layout.source_hashes()
    if snapshot.get('scope')!='developer' or snapshot.get('failures') or snapshot.get('assertions',0)<992 or snapshot.get('source_sha256')!=actual:
        raise RuntimeError('发行快照与当前源码不一致，请全量总检并重新冻结；不允许拿子集或旧记录充当992项通过')
    from PIL import Image,ImageDraw
    icon=Image.new('RGBA',(256,256),(0,0,0,0));draw=ImageDraw.Draw(icon)
    draw.polygon([(8,8),(248,8),(248,212),(212,248),(8,248)],fill='#FFFA00')
    draw.line([(8,60),(8,8),(60,8)],fill='#161616',width=10)
    draw.line([(52,180),(52,82),(128,146),(204,82),(204,180)],fill='#161616',width=18,joint='curve')
    draw.line((52,211,168,211),fill='#161616',width=9)
    icon.save(resources/'brand/mark.ico',sizes=[(16,16),(24,24),(32,32),(48,48),(64,64),(128,128),(256,256)])
    build=appenv.build_id().split('+',1)[-1]+'.'+variant
    versions={name:importlib.metadata.version(name) for name in ('numpy','pillow')+(('opencv-python','astropy') if variant=='full' else ())}
    (resources/'build-info.json').write_text(json.dumps({'version':appenv.VERSION,'build':build,'variant':variant,'python':sys.version.split()[0],
        'dependencies':versions,'sample_manifest_sha256':digest(resources/'sample-data/manifest.json')},ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    return build

def build(variant):
    identity=prepare(variant)
    dist=_layout.project_layout.INSTALL/'安装包'/variant;work=_layout.project_layout.INSTALL/'构建缓存'/('work-'+variant)
    assert dist.resolve().is_relative_to((_layout.project_layout.INSTALL/'安装包').resolve()) and work.resolve().is_relative_to((_layout.project_layout.INSTALL/'构建缓存').resolve())
    env=dict(os.environ,MW_VARIANT=variant,MARSWORKBENCH_DATA_HOME=str(_layout.project_layout.INSTALL/'构建缓存/构建用户'),PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1')
    evidence=_layout.project_layout.INSTALL/'机器验收'/('v'+appenv.VERSION);evidence.mkdir(parents=True,exist_ok=True)
    log=evidence/('build-'+variant+'.txt')
    command=[sys.executable,'-m','PyInstaller','--noconfirm','--clean','--distpath',str(dist),'--workpath',str(work),str(_layout.project_layout.INSTALL/'构建工具/win.spec')]
    with log.open('w',encoding='utf-8') as stream:
        result=subprocess.run(command,cwd=APP,env=env,stdout=stream,stderr=subprocess.STDOUT)
    if result.returncode:raise RuntimeError('构建失败，见'+str(log))
    folder=dist/'MarsWorkbench'
    (folder/'启动.bat').write_bytes('@echo off\r\nchcp 65001 >nul\r\ncd /d "%~dp0"\r\nstart "" /min "%~dp0MarsWorkbench.exe" %*\r\n'.encode('utf-8'))
    for name in ('LICENSE','THIRD_PARTY_NOTICES.md','DATA_SOURCES.md'):
        shutil.copy2(appenv.resources()/'licenses'/name,folder/name)
    manuals=_layout.project_layout.MANUALS
    for name in ('操作手册.md','代码作用手册.md'):shutil.copy2(manuals/name,folder/name)
    if (manuals/'插图').exists():shutil.copytree(manuals/'插图',folder/'插图',dirs_exist_ok=True)
    (folder/'版本信息.json').write_text(json.dumps({'version':appenv.VERSION,'build':identity,'variant':variant,'algorithm_protection':'未作逆向保护加固'},ensure_ascii=False,indent=2),encoding='utf-8')
    checks=[]
    for p in sorted(folder.rglob('*')):
        if p.is_file() and p.name!='校验.txt':checks.append(digest(p)+'  '+p.relative_to(folder).as_posix())
    (folder/'校验.txt').write_text('\n'.join(checks)+'\n',encoding='utf-8')
    archive=_layout.project_layout.INSTALL/'安装包'/('MarsWorkbench-'+appenv.VERSION+'-windows-x64-'+variant+'.zip')
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(folder.rglob('*')):
            if p.is_file():z.write(p,p.relative_to(folder).as_posix())
    result={'variant':variant,'version':appenv.VERSION,'build':identity,'zip':str(archive),'zip_bytes':archive.stat().st_size,
            'sha256':digest(archive),'directory':str(folder),'files':len(checks)+1}
    (evidence/('artifact-'+variant+'.json')).write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    archives=sorted((_layout.project_layout.INSTALL/'安装包').glob('*.zip'))
    (_layout.project_layout.INSTALL/'安装包/SHA256SUMS.txt').write_text('\n'.join(digest(p)+'  '+p.name for p in archives)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False),flush=True)

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--variant',choices=['standard','full','both'],default='both');args=parser.parse_args()
    for variant in (['standard','full'] if args.variant=='both' else [args.variant]):build(variant)

if __name__=='__main__':main()
