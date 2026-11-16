"""Reorganization contracts: real source ownership, preserved inputs and usable handbooks."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths
import project_layout as L
import appenv
import hashlib,json,os,subprocess,tempfile,threading,urllib.request,urllib.error
from http.server import BaseHTTPRequestHandler
from unittest.mock import patch

N=BAD=0
def check(ok,text):
    global N,BAD
    N+=1;BAD+=not ok;print(('OK ' if ok else 'FAIL ')+text,flush=True)

def main():
    distribution='--dist' in sys.argv
    from instruments import REGISTRY,pixl_browser,absorption,filtering,rmi,supercam_images
    from sherloc import qc
    import console,pixl_view,runtime_http
    core=L.CORE
    check(all((core/name).is_dir() for name in ('SHERLOC','SuperCam','ChemCam','APXS','PIXL')) and qc.__file__.startswith(str(core/'SHERLOC')) and pixl_browser.__file__.startswith(str(core/'PIXL')),'五仪器拥有实际源码目录，导入的是对应核心文件')
    modules=[('supercam',REGISTRY['supercam']),('chemcam',REGISTRY['chemcam']),('apxs',REGISTRY['apxs'])]
    check(all(Path(sys.modules[type(b).__module__].__file__).is_relative_to(core/L.GROUPS[name]) for name,b in modules) and console.__file__.startswith(str(L.APP)),'通用工作台与三个后端的物理归属分离，未从旧目录或临时副本导入')
    check(all(not hasattr(module,'HTML') and not hasattr(module,'JS') for module in (absorption,filtering,rmi,supercam_images)) and (L.APP/'instruments/pixl_output.py').is_file(),'界面HTML/JS及PIXL出图实现移出核心数据读取模块')
    baseline=json.loads((L.VERIFY/'fixtures/历史数值基准/重构前科学源码.json').read_text(encoding='utf-8'))
    check(len(baseline['files'])==10 and all(hashlib.sha256((L.ROOT/row['new']).read_bytes()).hexdigest()==row['sha256'] for row in baseline['files']),'10个冻结科学核心文件与基线指纹逐字节一致（演进须经 revisions 记录）')
    if distribution:
        rows=json.loads((appenv.resources()/'sample-data/manifest.json').read_text(encoding='utf-8'))['files']
        valid=all(hashlib.sha256((appenv.resources()/row['path']).read_bytes()).hexdigest()==row['sha256'] for row in rows)
        check(valid,'发行最小树内随包PDS原件SHA保持，未依赖开发练习集')
    elif not (L.INSTALL/'机器验收/数据迁移校验.json').is_file():
        check(True,'SKIP 数据不在位：机器验收迁移校验为构建产物（源码树不含；发行树全量校验）')
    else:
        audit=json.loads((L.INSTALL/'机器验收/数据迁移校验.json').read_text(encoding='utf-8'))
        removed=json.loads((L.INSTALL/'机器验收/文档移除清单.json').read_text(encoding='utf-8'))
        retired={row['path'] for row in removed['files']}
        rows=[row for row in audit['files'] if row['new'] not in retired]
        check(len(rows)>1000 and all((L.ROOT/row['new']).is_file() and hashlib.sha256((L.ROOT/row['new']).read_bytes()).hexdigest()==row['sha256'] for row in rows),'全部保留数据逐文件SHA一致；只排除已列明移入回收站的旧文档')
    with tempfile.TemporaryDirectory() as folder,patch.dict(os.environ,{},clear=False):
        os.environ.pop('MARSWORKBENCH_DATA_ROOT',None);os.environ.pop('MARSWORKBENCH_DATA_HOME',None)
        expected=L.ROOT/L.CONFIG['data']
        normal=L.data_root()==expected and all(L.instrument_data(name).is_relative_to(expected/L.GROUPS[name]) for name in L.GROUPS)
        os.environ['MARSWORKBENCH_DATA_HOME']=folder
        isolated=L.instrument_data('pixl').is_relative_to(Path(folder))
        seeded=appenv.copy_cached('edr_labels')
        originals=list((appenv.resources()/'cache-seed/edr_labels').rglob('*.img'))
        complete=len(originals)==2 and all((seeded/p.relative_to(appenv.resources()/'cache-seed/edr_labels')).read_bytes()==p.read_bytes() for p in originals)
        try:L.instrument_data('../outside');rejected=False
        except ValueError:rejected=True
        check(normal and isolated and complete and rejected,'五仪器下载根统一；干净缓存含原始图且状态隔离，未知仪器拒绝')
    old=L.ROOT/'02_练习数据与验证基准/processed'
    expected=L.ROOT/'03_下载数据/SHERLOC/练习数据/processed'
    source_alias=L.resolve_legacy(old)==expected
    with patch.object(L,'ROOT',None),patch.dict(os.environ,{'MARSWORKBENCH_DATA_ROOT':str(expected.parents[2])}):
        packaged_alias=L.resolve_legacy(old)==expected
    check(source_alias and packaged_alias and L.remap_preferences({'tone':'soft','folder':str(old)})['tone']=='soft','源码及工程内exe均可映射旧路径，不借源目录导入，非路径偏好不改')
    manuals={p.name for p in L.MANUALS.glob('*.md')}
    check(manuals=={'操作手册.md','代码作用手册.md'} and all(not (L.ROOT/name).exists() for name in ('04_使用文档','05_汇报材料','07_发布验收','docs','01_软件工作台')),'只保留两本手册，旧资料目录及旧代码根已经移除')
    check((L.INSTALL/'构建工具/win.spec').is_file() and len(list((L.INSTALL/'依赖').glob('*-lock.txt')))==3 and '05_安装与发布' in (L.ROOT/'启动.bat').read_text(encoding='utf-8') and not list((L.ROOT/L.CONFIG['data']).rglob('*.py')),'安装/依赖独立归位，启动入口指向新位置，数据目录不混入源码')
    class Bare(BaseHTTPRequestHandler):
        def log_message(self,*args):pass
        def _send(self,code,ctype,body):self.send_response(code);self.send_header('Content-Type',ctype);self.end_headers();self.wfile.write(body)
        def do_GET(self):self._send(404,'text/plain',b'not found')
    server=runtime_http.bind(runtime_http.decorate(Bare,{}),0);threading.Thread(target=server.serve_forever,daemon=True).start()
    try:
        base='http://127.0.0.1:%d'%server.server_port
        with urllib.request.urlopen(base+'/api/help',timeout=20) as response:help_data=json.load(response)
        with urllib.request.urlopen(base+'/manuals/%E6%93%8D%E4%BD%9C%E6%89%8B%E5%86%8C.md',timeout=20) as response:page=response.read().decode()
        quickstart=help_data['sections']['quickstart']
        check(set(help_data['sections'])=={'quickstart','faq','scope','about'} and '各工作台统一的操作方式' in quickstart and '文件夹与启动入口' not in quickstart and '操作手册' in page and 'endfield-yellow-paper' in page and 'id="mwRail"' not in page,'实际HTTP帮助按手册章节读取；完整阅读页不依赖工作台导航')
    finally:server.shutdown();server.server_close()
    entries=L.runtime_sources();book=(L.MANUALS/'代码作用手册.md').read_text(encoding='utf-8')
    check(len({name for name,_ in entries})==len(entries) and all(path.relative_to(L.ROOT).as_posix() in book for _,path in entries) and not (appenv.resources()/'docs').exists(),'运行模块唯一且都有作用索引；不维护额外旧帮助文档树')
    result=subprocess.run([sys.executable,str(L.ROOT/'setup.py'),'--version'],cwd=L.ROOT,capture_output=True,text=True)
    check(result.returncode==0 and result.stdout.strip().endswith(appenv.VERSION),'实际源码wheel构建元数据读取统一版本，重组后仍可打包')
    print('共 %d 项断言，%d 项失败'%(N,BAD));return int(BAD>0)
if __name__=='__main__':raise SystemExit(main())
