"""Exercise real HTTP pipelines in a fresh user profile, offline and without PATH Python."""

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[2]/'02_图形界面代码'))
import source_paths as _layout
import argparse
import csv
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.request
import zipfile

PROJECT=_layout.project_layout.ROOT
APP=_layout.project_layout.APP

def run(archive=None,proof=None,command_override=None,readonly=False,bat=False):
    owned=_layout.project_layout.INSTALL/'构建缓存/work-smoke';owned.mkdir(exist_ok=True)
    for stale in owned.glob('fresh-*'):          # 上次失败留下的调试目录
        shutil.rmtree(stale,ignore_errors=True)
    root=Path(tempfile.mkdtemp(prefix='fresh-',dir=owned))
    profile=root/'profile';profile.mkdir()
    package=root/'package';package.mkdir()
    if archive:
        with zipfile.ZipFile(archive) as z:
            for member in z.infolist():
                if not (package/member.filename).resolve().is_relative_to(package.resolve()):raise ValueError('Unsafe ZIP member')
            z.extractall(package)
        command=['cmd.exe','/d','/c',str(package/'启动.bat')] if bat else [str(package/'MarsWorkbench.exe')]
    else:command=command_override or [sys.executable,str(APP/'launch.py')]
    env=dict(os.environ,USERPROFILE=str(profile),HOME=str(profile),LOCALAPPDATA=str(profile/'AppData/Local'),
             MARSWORKBENCH_DATA_HOME=str(profile/'runtime'),MARSWORKBENCH_HEADLESS='1',MARSWORKBENCH_OFFLINE='1',
             PYTHONUTF8='1',PYTHONUNBUFFERED='1',PYTHONDONTWRITEBYTECODE='1',
             PATH=os.environ.get('SystemRoot','C:\\Windows')+'\\System32;'+os.environ.get('SystemRoot','C:\\Windows'))
    for key in ('PYTHONHOME','PYTHONPATH'):env.pop(key,None)
    denied_sid=None
    if readonly:
        if not archive or os.name!='nt':raise ValueError('只读ACL验收仅用于Windows ZIP')
        assert package.resolve().is_relative_to(owned.resolve()) and package.name=='package'
        # The only ACL target is this newly extracted, owned package directory.
        identity=subprocess.check_output(['whoami','/user','/fo','csv','/nh'],text=True).strip()
        denied_sid='*'+next(csv.reader([identity]))[1]
        subprocess.run(['icacls',str(package),'/deny',denied_sid+':(OI)(CI)(WD,AD,WEA,WA)'],check=True,capture_output=True)
        try:
            (package/'write-probe.tmp').write_text('must fail')
        except PermissionError:pass
        else:
            subprocess.run(['icacls',str(package),'/remove:d',denied_sid],check=True,capture_output=True)
            raise AssertionError('ACL未阻止写入')
    output=root/'server.txt';launcher=None;process=None;started=time.monotonic()
    shell=Path(os.environ.get('SystemRoot','C:\\Windows'))/'System32/WindowsPowerShell/v1.0/powershell.exe'
    exact_exe=str(package/'MarsWorkbench.exe').replace("'","''")
    own_filter="Get-CimInstance Win32_Process -Filter \"Name='MarsWorkbench.exe'\" | Where-Object { $_.ExecutablePath -eq '"+exact_exe+"' }"
    def bat_endpoints():
        script='$owned=@('+own_filter+'); $ports=@(Get-NetTCPConnection -State Listen | Where-Object {$_.OwningProcess -in $owned.ProcessId} | Select-Object -ExpandProperty LocalPort); @{pids=@($owned.ProcessId);ports=$ports} | ConvertTo-Json -Compress'
        answer=subprocess.run([str(shell),'-NoProfile','-Command',script],capture_output=True,text=True,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0),timeout=12)
        return json.loads(answer.stdout) if answer.returncode==0 and answer.stdout.strip() else {'pids':[],'ports':[]}
    def request(url,data=None):
        req=urllib.request.Request(url,json.dumps(data).encode() if data is not None else None,{'Content-Type':'application/json'})
        with urllib.request.urlopen(req,timeout=180) as r:
            body=r.read()
            if 'json' in r.headers.get('Content-Type',''):
                value=json.loads(body)
                if value.get('error'):raise RuntimeError(str(value))
                return value
            return body
    result={'mode':'zip' if archive else 'source','archive':str(archive) if archive else None,'root':str(root),'readonly_acl':readonly,'entry':'bat' if bat else 'exe' if archive else 'source','checks':[]}
    try:
        with output.open('wb') as stream:
            process=subprocess.Popen(command+['--no-window','--no-watchdog','--port','0','--console-port','0','--sherloc-port','0'],
                cwd=package,env=env,stdout=stream,stderr=subprocess.STDOUT,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
        while time.monotonic()-started<30:
            text=output.read_text(encoding='utf-8',errors='replace')
            urls=re.findall(r'http://127\.0\.0\.1:\d+/',text)
            if bat and len(urls)<2:
                endpoints=bat_endpoints();result['owned_pids']=endpoints['pids']
                for port in endpoints['ports']:
                    candidate='http://127.0.0.1:%d/'%port
                    try:page=request(candidate).decode('utf-8')
                    except Exception:continue
                    if 'id="sherlocCard"' in page:
                        common_urls=re.findall(r'href="(http://127\.0\.0\.1:\d+/)\?app=chemcam',page)
                        if common_urls:urls=[candidate,common_urls[0]];break
            if len(urls)>=2:
                launcher,common=urls[:2]
                request(launcher+'api/ping');break
            if process.poll() is not None and (not bat or process.returncode):raise RuntimeError(text)
            time.sleep(.1)
        if launcher is None:raise TimeoutError('启动超过30秒\n'+output.read_text(encoding='utf-8',errors='replace'))
        result['startup_seconds']=round(time.monotonic()-started,3)
        result['launcher']=launcher;result['console']=common
        install=request(launcher+'api/samples/install',{'instrument':'sherloc','destination':str(profile/'Documents/MarsWorkbench示例')})
        view=install['url'];folders=install['folders'];result['sherloc']=view
        for url in (launcher,view,common+'?app=supercam',common+'?app=chemcam',common+'?app=apxs',launcher+'pixl'):
            page=request(url).decode('utf-8')
            assert page.count('data-mw-theme="endfield-yellow-paper"')==1 and 'data:font/woff2;base64,' in page,url
        result['checks'].append({'name':'启动器与五仪器单一主题、字体离线自包含','ok':True})
        spectrum=request(view+'api/plot?scan=0&point=0&mode=r1')
        png=request(view+'api/png?scan=0&point=0&mode=r1')
        assert spectrum.get('series') and png.startswith(b'\x89PNG\r\n\x1a\n')
        for instrument in ('apxs','supercam','chemcam'):
            assert request(common+'api/plot?app='+instrument+'&scan=0&point=0&shot=-1').get('series'),instrument
        health=request(launcher+'api/health?network=0');result['health']=health
        assert health['frozen']==bool(archive) and health['writable'] and health['status']!='error'
        if archive:
            full='full' in Path(archive).name
            assert health['components']['astropy']['available']==full
            assert health['components']['opencv-python']['available']==full
        result['checks'].append({'name':'冷启动<30秒、四仪器看谱、PNG、组件隔离','ok':True})
        qc=request(view+'api/qc');assert Path(qc['index']).is_file(),qc
        result['checks'].append({'name':'SHERLOC真实质控报告生成','ok':True})
        maps=request(common+'api/rmi/map?app=chemcam&scan=0');result['maps']=maps
        map_files=list((Path(folders['chemcam'])/'_maps').glob('map_*.png'))
        assert len(map_files)>=2,maps
        request(common+'api/qc?app=chemcam')
        deadline=time.monotonic()+180
        while time.monotonic()<deadline:
            status=request(common+'api/status?app=chemcam')['qc']
            if not status.get('running'):break
            time.sleep(.2)
        assert status.get('index') and not status.get('error'),status
        html=request(common+'data/chemcam/_qc/consoles/index.html').decode('utf-8')
        assert 'elementMaps' in html and '../../_maps/' in html
        assert html.count('data-mw-theme="endfield-yellow-paper"')==1
        assert request(common+'data/chemcam/_maps/'+map_files[0].name).startswith(b'\x89PNG')
        result['checks'].append({'name':'ChemCam离线共帧地图及报告链接闭环','ok':True})
        demo=request(launcher+'api/pixl/demo',{})
        pixl=request(launcher+'api/pixl/load',{'path':demo['folder']})
        oid=pixl['observations'][0]['id']
        point=request(launcher+'api/pixl/point?obs='+oid+'&pmc=6')
        assert point['values']['SiO2']=={'value':41.6797,'error':2.1058}
        exported=request(launcher+'api/pixl/export',{'obs':oid,'element':'SiO2'})
        assert request(launcher.rstrip('/')+exported['links']['热图PNG']).startswith(b'\x89PNG')
        assert b'41.6797' in request(launcher.rstrip('/')+exported['links']['报告'])
        result['pixl_export']=exported
        # Unified workflow endpoints remain separate from the backwards-compatible bundle API.
        request(launcher+'api/appearance',{'theme':'dark'})
        assert request(common+'api/appearance')['theme']=='dark'
        request(common+'api/appearance',{'theme':'light'})
        assert request(launcher+'api/appearance')['theme']=='light'
        standalone=request(launcher+'api/pixl/export-selected',{'kind':'csv','obs':oid,'pmc':6,'scope':'current'})
        assert b'41.6797' in request(launcher.rstrip('/')+standalone['url'])
        request(launcher+'api/pixl/report/start',{'obs':oid})
        deadline=time.monotonic()+90
        while time.monotonic()<deadline:
            report=request(launcher+'api/pixl/report/status')
            if not report['running']:break
            time.sleep(.1)
        assert report['done']==report['total']==1 and not report.get('error') and report['url']
        report_html=request(launcher.rstrip('/')+report['url']).decode('utf-8')
        assert '数据检查与成分报告' in report_html and report_html.count('data-mw-theme="endfield-yellow-paper"')==1
        assert request(launcher+'api/pixl/openreport',{})['ok']
        result['workflow_report']=report
        result['checks'].append({'name':'PIXL离线示例、官方成分/误差、热图与报告','ok':True})
        if proof:
            proof=Path(proof);proof.mkdir(parents=True,exist_ok=True)
            (proof/'spectrum.png').write_bytes(png)
            for name in ('sherloc','chemcam'):
                for derived in ('_qc','_maps'):
                    folder=Path(folders[name])/derived
                    if folder.exists():shutil.copytree(folder,proof/name/derived,dirs_exist_ok=True)
        result['ok']=True
    except Exception as exc:
        result.update(ok=False,error=repr(exc))
    finally:
        # 成功运行删临时目录（每次 run 约 76 MB，曾累积到 6.2 GB）；
        # 失败时保留 root 供调试，由下次 run 前的清扫兜底。
        if result.get('ok'):
            shutil.rmtree(root,ignore_errors=True)
        if launcher:
            try:request(launcher+'api/quit')
            except Exception:pass
        if process:
            try:process.wait(timeout=12)
            except subprocess.TimeoutExpired:process.terminate();process.wait(timeout=10)
        if denied_sid:
            subprocess.run(['icacls',str(package),'/remove:d',denied_sid],check=True,capture_output=True)
        if bat:
            # BAT starts a detached console. Scope cleanup to this extracted exe only.
            assert package.resolve().is_relative_to(owned.resolve())
            subprocess.run([str(shell),'-NoProfile','-Command',own_filter+' | ForEach-Object { Stop-Process -Id $_.ProcessId -Force -ErrorAction SilentlyContinue }'],capture_output=True,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0),timeout=12)
        result['server_log']=output.read_text(encoding='utf-8',errors='replace')
        result['elapsed_seconds']=round(time.monotonic()-started,3)
        # Keep the isolated profile as reproducible evidence; it is ignored by Git.
    return result

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--zip',type=Path);parser.add_argument('--source',action='store_true');parser.add_argument('--out',type=Path,required=True);parser.add_argument('--proof',type=Path);parser.add_argument('--readonly',action='store_true');parser.add_argument('--bat',action='store_true')
    args=parser.parse_args()
    if not args.zip and not args.source:parser.error('提供--zip或--source')
    result=run(args.zip,args.proof,readonly=args.readonly,bat=args.bat);args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k not in ('server_log','health','maps')},ensure_ascii=False),flush=True)
    return int(not result['ok'])

if __name__=='__main__':raise SystemExit(main())
