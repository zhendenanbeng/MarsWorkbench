"""Explicit opt-in real PDS download test for a built Windows executable."""

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[2]/'02_图形界面代码'))
import source_paths as _layout
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile
import time
import urllib.request

PROJECT=_layout.project_layout.ROOT

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--exe',type=Path,required=True);parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    root=Path(tempfile.mkdtemp(prefix='pixl-network-',dir=PROJECT/'.work'))
    env=dict(os.environ,MARSWORKBENCH_DATA_HOME=str(root/'state'),MARSWORKBENCH_HEADLESS='1',PYTHONUTF8='1',PATH=os.environ.get('SystemRoot','C:\\Windows')+'\\System32')
    env.pop('MARSWORKBENCH_OFFLINE',None);env.pop('PYTHONPATH',None);env.pop('PYTHONHOME',None)
    def api(path,data=None):
        request=urllib.request.Request(url+path,json.dumps(data).encode() if data is not None else None,{'Content-Type':'application/json'})
        with urllib.request.urlopen(request,timeout=90) as response:
            raw=response.read()
            if 'json' not in response.headers.get('Content-Type',''):return raw
            result=json.loads(raw)
            if result.get('error'):raise RuntimeError(str(result))
            return result
    result={'ok':False,'root':str(root)};url=None
    with (root/'server.txt').open('wb') as stream:process=subprocess.Popen([str(args.exe.resolve()),'--no-window','--no-watchdog','--port','0','--console-port','0','--sherloc-port','0'],env=env,stdout=stream,stderr=subprocess.STDOUT,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
    try:
        end=time.monotonic()+30
        while time.monotonic()<end:
            matches=re.findall(r'http://127\.0\.0\.1:\d+/',(root/'server.txt').read_text(encoding='utf-8',errors='replace'))
            if matches:url=matches[0];break
            if process.poll() is not None:raise RuntimeError('构建程序启动失败')
            time.sleep(.1)
        if url is None:raise TimeoutError('没有取得启动器地址')
        health=api('api/health?network=0');assert health['frozen']
        catalog=api('api/pixl/download/catalog?refresh=1');result['inventories']=catalog['inventories']
        assert all(row['source']=='在线PDS清单' for row in catalog['inventories'])
        target=next(row for row in catalog['scans'] if row['sol']==187 and row['rtt']=='069927431')
        api('api/pixl/download/start',{'catalog_id':catalog['catalog_id'],'scan_ids':[target['id']],'destination':str(root/'download'),'images':True})
        end=time.monotonic()+180
        while time.monotonic()<end:
            status=api('api/pixl/download/status')
            if not status['running']:break
            time.sleep(.3)
        result['job']=status
        assert status['status']=='completed' and status['ready'] and status['files_done']>=8,status
        local=api('api/pixl/load',{'path':status['folder']});oid=local['observations'][0]['id']
        point=api('api/pixl/point?obs='+oid+'&pmc=6');result['point6']=point
        assert point['values']['SiO2']=={'value':41.6797,'error':2.1058}
        observation=api('api/pixl/observation?obs='+oid)
        assert len(observation['frames'])>=4 and not any(f['missing'] for f in observation['frames'])
        assert api('api/pixl/image?obs='+oid+'&frame=3').startswith(b'\x89PNG')
        result.update(ok=True,build=health['build'],frames=len(observation['frames']))
    except Exception as exc:result['error']=repr(exc)
    finally:
        if url:
            try:api('api/quit')
            except Exception:pass
        try:process.wait(timeout=12)
        except subprocess.TimeoutExpired:process.terminate();process.wait(timeout=10)
        args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2),flush=True)
    return int(not result['ok'])

if __name__=='__main__':raise SystemExit(main())
