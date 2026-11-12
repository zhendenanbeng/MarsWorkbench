"""Offline acquisition contracts: identity, file integrity, preservation and cancellation."""

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout
import copy
import hashlib
import io
import json
import logging
import os
from pathlib import Path
import sys
import tempfile
import threading
import time
import urllib.request
import urllib.error
from http.server import BaseHTTPRequestHandler
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import appenv,pixl_view,runtime_http
from instruments import pixl_download as d
from instruments import pixl_browser as browser

N=BAD=0
def check(ok,message):
    global N,BAD
    N+=1;BAD+=not ok;print(('OK ' if ok else 'FAIL ')+message,flush=True)
def fails(fn,typ=Exception):
    try:fn();return False
    except typ:return True
class Response(io.BytesIO):
    def __init__(self,raw,url,length=None):super().__init__(raw);self.url=url;self.headers={'Content-Length':str(len(raw) if length is None else length)}
    def geturl(self):return self.url

def main():
    demo=appenv.resources()/'pixl-demo';rqc=next(demo.glob('*rqc*.csv'));rxl=next(demo.glob('*rxl*.csv'))
    clid='urn:nasa:pds:mars2020_pixl:data_oxides_pmc:'+rqc.stem;plid='urn:nasa:pds:mars2020_pixl:data_processed:'+rxl.stem
    cmeta=d.pixl.parse_lid(clid);pmeta=d.pixl.parse_lid(plid)
    inv=('P,'+clid+'::1.0\nS,'+clid+'::1.0\nP,'+clid[:-2]+'02::1.0\n').encode()
    check(len(d.parse_inventory(inv,'data_oxides_pmc','rqc'))==1 and d.parse_inventory(inv,'data_oxides_pmc','rqc')[0]['stem'].endswith('02'),'清单仅主成员、仅RQC、同逻辑产品取最新版')
    check(not d.parse_inventory(('P,'+plid).encode(),'data_oxides_pmc','rqc'),'坐标产品不可冒充成分清单')
    entries={};raws={};refs={}
    for lid,path,kind in ((clid,rqc,'rqc'),(plid,rxl,'rxl')):
        mapping={}
        for src in (path,path.with_suffix('.xml')):
            raw=src.read_bytes();url='https://pds-geosciences.wustl.edu/test/'+src.name
            mapping[src.name.lower()]={'name':src.name,'url':url,'lid':lid,'md5':hashlib.md5(raw).hexdigest()};raws[url]=raw
        if kind=='rxl':
            for name in browser.read_coordinates(path)['mcc_tifs']:
                raw=next(demo.glob('*.IMG')).read_bytes();url='https://pds-geosciences.wustl.edu/test/'+name
                mapping[name.lower()]={'name':name,'url':url,'lid':lid,'md5':hashlib.md5(raw).hexdigest()};raws[url]=raw
        refs[lid]=mapping
    entry=refs[clid][rqc.name.lower()]
    record={'summary':{'hits':1},'data':[{'properties':{'lid':[clid],
        'ops:Data_File_Info.ops:file_name':[rqc.name],'ops:Data_File_Info.ops:file_ref':[entry['url']],
        'ops:Data_File_Info.ops:md5_checksum':[entry['md5']],
        'ops:Label_File_Info.ops:file_name':[rqc.with_suffix('.xml').name],
        'ops:Label_File_Info.ops:file_ref':[refs[clid][rqc.with_suffix('.xml').name.lower()]['url']]}}]}
    mapping=d.refs_from_record(record,clid)
    check(mapping[rqc.name.lower()]['md5']==entry['md5'] and any(k.endswith('.xml') for k in mapping),'官方API同时解析数据与XML引用及公开MD5')
    check(fails(lambda:d.refs_from_record(record,plid),ValueError),'注册表LID不匹配拒绝')
    bad=copy.deepcopy(record);bad['data'][0]['properties']['ops:Data_File_Info.ops:file_name']=['wrong.csv']
    check(fails(lambda:d.refs_from_record(bad,clid),ValueError),'文件名与直链不一致拒绝')
    check(all(fails(lambda u=u:d.safe_url(u),ValueError) for u in ('http://pds.nasa.gov/a','https://pds.nasa.gov.evil.example/a','https://user:pass@pds.nasa.gov/a','file:///etc/a')) and d.safe_url('https://pds.mcp.nasa.gov/api/search/1/products'),'仅核验过的PDS HTTPS主机，支持实测NASA跳转主机')
    with tempfile.TemporaryDirectory() as name,patch.dict(os.environ,{'MARSWORKBENCH_DATA_HOME':str(Path(name)/'state')}):
        root=Path(name);folder=root/'out';folder.mkdir();manifest={}
        def urlopen(req,**kwargs):
            url=req if isinstance(req,str) else req.full_url
            return Response(raws[url],url)
        with patch.object(urllib.request,'urlopen',side_effect=urlopen):
            row=d.get_file(folder,entry,'rqc',manifest,lambda:False,lambda n,t:None)
        check((folder/rqc.name).read_bytes()==rqc.read_bytes() and row['sha256']==hashlib.sha256(rqc.read_bytes()).hexdigest() and not list(folder.glob('*.part')),'下载原字节、公开MD5和SHA通过后才成为正式文件')
        with patch.object(urllib.request,'urlopen',side_effect=AssertionError('不得重下')):
            reused=d.get_file(folder,entry,'rqc',manifest,lambda:False,lambda n,t:None)
        check(reused['reused'],'已有文件官方哈希一致时核验复用，不重复读取数据网络')
        (folder/rqc.name).write_bytes(b'user file')
        check(fails(lambda:d.get_file(folder,entry,'rqc',manifest,lambda:False,lambda n,t:None)) and (folder/rqc.name).read_bytes()==b'user file','已有不同字节文件不覆盖')
        (folder/rqc.name).unlink()
        with patch.object(urllib.request,'urlopen',return_value=Response(rqc.read_bytes(),entry['url'],len(rqc.read_bytes())+5)):
            short=fails(lambda:d.get_file(folder,entry,'rqc',manifest,lambda:False,lambda n,t:None))
        check(short and not (folder/rqc.name).exists() and not list(folder.glob('*.part')),'短读失败清理临时体，不留下伪完成文件')
        with patch.object(urllib.request,'urlopen',return_value=Response(b'<html>oops</html>',entry['url'])):
            errorpage=fails(lambda:d.get_file(folder,entry,'rqc',manifest,lambda:False,lambda n,t:None))
        check(errorpage and not (folder/rqc.name).exists(),'HTML错误页不成为CSV')
        cancel=[False]
        with patch.object(urllib.request,'urlopen',side_effect=urlopen):
            stopped=fails(lambda:d.get_file(folder,entry,'rqc',manifest,lambda:cancel[0],lambda n,t:cancel.__setitem__(0,True)),d.Cancelled)
        check(stopped and not (folder/rqc.name).exists() and not list(folder.glob('*.part')),'读取中取消清理part，不破坏其他已完成文件')
        escape=dict(entry,name='../escape.csv')
        check(fails(lambda:d.get_file(folder,escape,'rqc',manifest,lambda:False,lambda n,t:None),ValueError),'远端文件名不得越出保存目录')
        xml=refs[clid][rqc.with_suffix('.xml').name.lower()];x=root/'wrong.xml';x.write_bytes(rqc.with_suffix('.xml').read_bytes())
        check(fails(lambda:d.validate(x,dict(xml,lid=plid),'rqc'),ValueError),'XML标签必须对应选定LID')
        with patch.object(d.pixl,'CACHE',str(root/'no-legacy')),patch.object(d,'read_url',return_value=('P,'+clid+'::1.1\n').encode()):
            values,meta=d.inventory('data_oxides_pmc','rqc',True)
        with patch.object(d.pixl,'CACHE',str(root/'no-legacy')),patch.object(d,'read_url',side_effect=OSError('offline')):
            values,cached=d.inventory('data_oxides_pmc','rqc',True)
        check(len(values)==1 and cached['warning'] and cached['source']=='本地缓存' and meta['source']=='在线PDS清单','离线清单降级如实标注来源，不宣称联网成功')
        scan={'id':'a','sol':187,'lid':clid,'rxl_lids':[plid],'rtt':cmeta['seq'][:9],'file':rqc.name,'placement':cmeta['pos'],'coordinate_sols':[186]}
        c=d.Controller();c.catalog={'scans':[scan]};c.token='token'
        check(fails(lambda:c.start('stale',['a'],str(root/'job')),ValueError) and fails(lambda:c.start('token',['unknown'],str(root/'job')),ValueError),'过期清单和不存在扫描不启动后台任务')
        with patch.object(d,'resolve',side_effect=lambda lid:refs[lid]),patch.object(urllib.request,'urlopen',side_effect=urlopen):
            c.start('token',['a'],str(root/'job'),True);c.thread.join(20)
        first=c.snapshot();files=list((root/'job').rglob('*.csv'))
        local,observations,w=browser.scan(root/'job');data=browser.load(local,observations[0])
        check(first['status']=='completed' and first['ready'] and first['files_done']==8 and len(files)==2 and len(data['frames'])==4 and all(f['path'] for f in data['frames']),'整套配套下载：成分sol187、坐标sol186、四帧MCC均可独立浏览')
        with patch.object(d,'resolve',side_effect=lambda lid:refs[lid]),patch.object(urllib.request,'urlopen',side_effect=AssertionError('已完成文件不应再下载')):
            c.start('token',['a'],str(root/'job'),True);c.thread.join(20)
        check(c.snapshot()['status']=='completed' and c.snapshot()['reused']==8,'重复任务复用全部核验过的文件，清单SHA持续保留')
        def broken(req,**kwargs):
            url=req if isinstance(req,str) else req.full_url
            if url.endswith('.tif'):raise OSError('MCC unavailable')
            return Response(raws[url],url)
        with patch.object(d,'resolve',side_effect=lambda lid:refs[lid]),patch.object(urllib.request,'urlopen',side_effect=broken):
            c.start('token',['a'],str(root/'partial'),True);c.thread.join(20)
        check(c.snapshot()['status']=='completed_with_errors' and c.snapshot()['ready'] and c.snapshot()['files_done']==4,'图像失败不丢已完成表，状态明确部分失败且允许打开数据')
        class Handler(BaseHTTPRequestHandler):
            def log_message(self,*args):pass
            def _send(self,code,ctype,body):self.send_response(code);self.send_header('Content-Type',ctype);self.end_headers();self.wfile.write(body)
            def do_GET(self):self._send(404,'text/plain',b'not found')
        state={};server=runtime_http.bind(runtime_http.decorate(Handler,state),0);thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start()
        url='http://127.0.0.1:%d/api/pixl/download/'%server.server_port
        try:
            with urllib.request.urlopen(url+'status') as r:status=json.load(r)
            try:urllib.request.urlopen(url+'cancel');wrong_method=False
            except urllib.error.HTTPError as exc:wrong_method=exc.code==400
            check(status['status']=='idle' and Path(status['destination']).is_dir() and wrong_method,'真实HTTP进度接口及POST动作限制，默认目录可选择')
        finally:
            server.shutdown();server.server_close();thread.join()
            logger=logging.getLogger('marsworkbench')
            for handler in list(logger.handlers):handler.close();logger.removeHandler(handler)
    print('共 %d 项断言，%d 项失败'%(N,BAD));return int(BAD>0)

if __name__=='__main__':raise SystemExit(main())
