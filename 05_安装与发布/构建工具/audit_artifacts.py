"""Verify every archived byte against the shipped checksums and frozen run."""

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[2]/'02_图形界面代码'))
import source_paths as _layout
import hashlib
import json
from pathlib import Path
import zipfile
import sys

PROJECT=_layout.project_layout.ROOT
sys.path.insert(0,str(_layout.project_layout.APP))
import appenv

def main():
    rows=[];outer=[]
    for variant in ('standard','full'):
        archive=_layout.project_layout.INSTALL/'安装包'/('MarsWorkbench-'+appenv.VERSION+'-windows-x64-'+variant+'.zip')
        with zipfile.ZipFile(archive) as z:
            names=set(z.namelist());checks={line.split('  ',1)[1]:line.split('  ',1)[0] for line in z.read('校验.txt').decode('utf-8').splitlines()}
            assert names==set(checks)|{'校验.txt'}
            assert all(hashlib.sha256(z.read(name)).hexdigest()==sha for name,sha in checks.items())
            snap=json.loads(z.read('_internal/resources/_last_run.json'))
            assert snap['assertions']>=992 and snap['failures']==0 and snap['scope']=='developer'
            assert not any('/verify/' in name or '/03_下载数据/SHERLOC/练习数据/' in name or '.console_settings' in name for name in names)
            assert all(name in names for name in ('启动.bat','操作手册.md','代码作用手册.md','LICENSE','DATA_SOURCES.md','版本信息.json'))
            info=json.loads(z.read('_internal/resources/build-info.json'))
            assert info['variant']==variant and info['version']==appenv.VERSION
        sha=hashlib.sha256(archive.read_bytes()).hexdigest();outer.append(sha+'  '+archive.name)
        rows.append({'variant':variant,'bytes':archive.stat().st_size,'sha256':sha,'files':len(names),'all_inner_sha256_ok':True,'developer_data_excluded':True,'frozen_assertions':snap['assertions'],'build':info['build']})
    wheel=_layout.project_layout.INSTALL/'安装包/wheels'/('mars_workbench-'+appenv.VERSION+'-py3-none-any.whl');wheel_sha=hashlib.sha256(wheel.read_bytes()).hexdigest()
    outer.append(wheel_sha+'  wheels/'+wheel.name)
    (_layout.project_layout.INSTALL/'安装包/SHA256SUMS.txt').write_text('\n'.join(outer)+'\n',encoding='utf-8')
    data={'zip':rows,'wheel':{'name':wheel.name,'bytes':wheel.stat().st_size,'sha256':wheel_sha}}
    evidence=_layout.project_layout.INSTALL/'机器验收'/('v'+appenv.VERSION);evidence.mkdir(exist_ok=True)
    (evidence/'发行物逐文件校验.json').write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(data,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
