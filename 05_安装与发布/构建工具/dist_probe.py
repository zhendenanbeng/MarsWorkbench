"""Run distribution checks from a minimal, reorganized source tree."""
from pathlib import Path
import sys,json,os,shutil,subprocess,tempfile
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'02_图形界面代码'))
import source_paths,appenv
L=source_paths.project_layout

def main():
    base=L.INSTALL/'构建缓存/发行检查';base.mkdir(parents=True,exist_ok=True)
    root=Path(tempfile.mkdtemp(prefix='dist-layout-',dir=base))
    shutil.copy2(L.ROOT/'.project-layout.json',root/'.project-layout.json')
    shutil.copy2(L.ROOT/'pyproject.toml',root/'pyproject.toml')
    shutil.copy2(L.ROOT/'setup.py',root/'setup.py')
    shutil.copy2(L.ROOT/'启动.bat',root/'启动.bat')
    for relative in L.CONFIG['namespace_roots']+[str(L.VERIFY.relative_to(L.ROOT)),L.CONFIG['manuals']]:
        src=L.ROOT/relative
        shutil.copytree(src,root/relative,ignore=shutil.ignore_patterns('__pycache__','.cache','build','.gui_settings.json','.console_settings.json'))
    for relative in [str((L.INSTALL/'构建工具').relative_to(L.ROOT)),str((L.INSTALL/'资源模板').relative_to(L.ROOT)),str((L.INSTALL/'依赖').relative_to(L.ROOT))]:
        shutil.copytree(L.ROOT/relative,root/relative,ignore=shutil.ignore_patterns('work*','lib','bdist*','__pycache__'))
    shutil.copy2(L.INSTALL/'版本信息.json',root/L.CONFIG['install']/'版本信息.json')
    for relative in L.CONFIG['resource_sources'].values():
        src=L.ROOT/relative;dest=root/relative
        if src.is_file():dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dest)
    assert not (root/'03_下载数据/SHERLOC/练习数据').exists()
    assert not (root/'03_下载数据/SHERLOC/E2E测试数据').exists()
    env=dict(os.environ,MARSWORKBENCH_DATA_HOME=str(root/'test-user'),PYTHONUTF8='1',PYTHONUNBUFFERED='1',PYTHONDONTWRITEBYTECODE='1')
    env.pop('MARSWORKBENCH_DATA_ROOT',None)
    evidence=L.INSTALL/'机器验收'/('v'+appenv.VERSION);evidence.mkdir(parents=True,exist_ok=True)
    with (evidence/'发行独立检查.txt').open('w',encoding='utf-8') as stream:
        result=subprocess.run([sys.executable,str(root/L.CONFIG['verification']/'run_all.py'),'--dist'],cwd=root,env=env,stdout=stream,stderr=subprocess.STDOUT)
    snapshot=root/'test-user/reports_selfcheck/dist_last_run.json'
    data=json.loads(snapshot.read_text(encoding='utf-8')) if snapshot.exists() else {}
    proof={'minimal_root':str(root),'developer_data_excluded':True,'code':result.returncode,'assertions':data.get('assertions'),'failures':data.get('failures'),'scope':data.get('scope')}
    (evidence/'发行独立检查.json').write_text(json.dumps(proof,ensure_ascii=False,indent=2),encoding='utf-8')
    if snapshot.exists():shutil.copy2(snapshot,evidence/'dist_last_run.json')
    print(json.dumps(proof,ensure_ascii=False))
    return int(result.returncode!=0 or data.get('failures')!=0)
if __name__=='__main__':raise SystemExit(main())
