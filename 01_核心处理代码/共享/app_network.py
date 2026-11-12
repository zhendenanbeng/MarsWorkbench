"""Shared timeout/offline semantics for outbound PDS requests; loopback is untouched."""
import os
import socket
import threading
import urllib.error
import urllib.parse
import urllib.request

OFFLINE_MESSAGE='网络暂不可用或服务未响应；已有本地数据仍可使用。请打开“帮助 → 关于 → 健康检查”查看连接状态。'
_local=threading.local()
_original=urllib.request.urlopen
_installed=False

def canonical_url(url):
    old='pds-imaging.jpl.nasa.gov/archive/m20/cumulative/mars2020_imgops/'
    if old in url:return url.replace(old,'planetarydata.jpl.nasa.gov/img/data/mars2020/mars2020_imgops/').replace('http://','https://',1)
    return url

class NetworkUnavailable(OSError):pass

class Response:
    def __init__(self,raw):self.raw=raw
    def __getattr__(self,name):return getattr(self.raw,name)
    def __enter__(self):self.raw.__enter__();return self
    def __exit__(self,*args):return self.raw.__exit__(*args)
    def read(self,*args,**kwargs):
        try:return self.raw.read(*args,**kwargs)
        except (TimeoutError,socket.timeout,ConnectionError,OSError) as exc:
            _local.failed=True;raise NetworkUnavailable(OFFLINE_MESSAGE) from exc

def begin_request():_local.failed=False
def failed():return bool(getattr(_local,'failed',False))

def install():
    global _installed
    if _installed:return
    def request(url,*args,**kwargs):
        address=url if isinstance(url,str) else url.full_url
        mapped=canonical_url(address)
        if mapped!=address:
            if isinstance(url,str):url=mapped
            else:url.full_url=mapped
            address=mapped
        host=urllib.parse.urlparse(address).hostname
        if host in ('localhost','127.0.0.1','::1'):return _original(url,*args,**kwargs)
        if os.environ.get('MARSWORKBENCH_OFFLINE')=='1':
            _local.failed=True;raise NetworkUnavailable(OFFLINE_MESSAGE)
        # urllib's default socket timeout is otherwise unbounded.
        if len(args)>=2:
            args=list(args);args[1]=min(float(args[1] or 30),30)
        else:kwargs['timeout']=min(float(kwargs.get('timeout') or 30),30)
        try:return Response(_original(url,*args,**kwargs))
        except urllib.error.HTTPError as exc:
            if exc.code in (408,429,500,502,503,504):
                _local.failed=True;raise NetworkUnavailable(OFFLINE_MESSAGE) from exc
            raise
        except (urllib.error.URLError,TimeoutError,socket.timeout,ConnectionError) as exc:
            _local.failed=True;raise NetworkUnavailable(OFFLINE_MESSAGE) from exc
    urllib.request.urlopen=request;_installed=True


def probe(timeout=10):
    if os.environ.get('MARSWORKBENCH_OFFLINE')=='1':return {'ok':False,'message':OFFLINE_MESSAGE,'offline_mode':True}
    try:
        req=urllib.request.Request('https://pds-geosciences.wustl.edu/',method='HEAD',headers={'User-Agent':'MarsWorkbench/0.11'})
        with urllib.request.urlopen(req,timeout=min(timeout,10)) as r:
            return {'ok':200<=r.status<400,'message':'PDS连接可用'}
    except Exception:return {'ok':False,'message':OFFLINE_MESSAGE}
