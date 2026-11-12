"""Transport-only concerns: correlation IDs, logs, local binding and version footer."""

import source_paths  # source roots; no-op for bundled imports
import errno
import json
import sys
import uuid
import socket
from http.server import ThreadingHTTPServer
import appenv
import app_network

class LocalServer(ThreadingHTTPServer):
    allow_reuse_address=False
    def server_bind(self):
        if hasattr(socket,'SO_EXCLUSIVEADDRUSE'):
            self.socket.setsockopt(socket.SOL_SOCKET,socket.SO_EXCLUSIVEADDRUSE,1)
        return super().server_bind()


def bind(handler,port,attempts=20):
    last=None
    for candidate in ([0] if port==0 else range(port,min(port+attempts,65536))):
        try:return LocalServer(('127.0.0.1',candidate),handler)
        except OSError as exc:
            last=exc
            if exc.errno not in (errno.EADDRINUSE,errno.EACCES) and getattr(exc,'winerror',None) not in (10048,10013):raise
    raise OSError('连续%d个本机端口不可用，请关闭占用程序或用--port指定其他端口：%s'%(attempts,last))


def initialize():
    if getattr(sys,'frozen',False):
        for stream in (sys.stdout,sys.stderr):
            if stream is not None and hasattr(stream,'reconfigure'):stream.reconfigure(encoding='utf-8',errors='replace')
    appenv.logger();app_network.install()
    import os
    if os.environ.get('MARSWORKBENCH_HEADLESS')=='1' and hasattr(os,'startfile'):
        os.startfile=lambda *args,**kwargs:None
    if (appenv.resources()/'sample-data/manifest.json').exists():
        import samples
        samples.seed_cache()


def decorate(Original,state):
    class Handler(Original):
        def do_GET(self):
            return self._dispatch('do_GET')
        def do_POST(self):
            return self._dispatch('do_POST')
        def _dispatch(self,method):
            self.request_id=uuid.uuid4().hex[:12]
            app_network.begin_request()
            try:
                origin=self.headers.get('Origin')
                expected={'http://127.0.0.1:%d'%self.server.server_port,'http://localhost:%d'%self.server.server_port}
                if (origin and origin not in expected) or self.headers.get('Sec-Fetch-Site')=='cross-site':
                    return self._send(403,'application/json',json.dumps({'error':'请从本机工作台页面操作。'}).encode())
                if self.path.split('?')[0]=='/api/quit' and 'launcher_port' in state:
                    import threading
                    self._send(200,'application/json',b'{"ok":true}')
                    threading.Thread(target=state['shutdown'],daemon=True).start()
                    return
                try:import product_api
                except ImportError:product_api=None
                if product_api is not None and product_api.route(self,state):return
                return getattr(super(),method)()
            except (ValueError,FileNotFoundError,PermissionError) as exc:
                return self._send(400,'application/json; charset=utf-8',json.dumps({'error':str(exc)},ensure_ascii=False).encode())
            except Exception:
                appenv.logger().exception('request_id=%s path=%s',self.request_id,self.path)
                return self._send(500,'application/json; charset=utf-8',json.dumps({'error':'操作未完成，请查看健康检查及日志，并提供下方请求编号。','request_id':self.request_id},ensure_ascii=False).encode())

        def _send(self,code,ctype,body):
            if 'application/json' in ctype:
                try:
                    payload=json.loads(body)
                    if isinstance(payload,dict):
                        error=payload.get('error')
                        if error:
                            rid=getattr(self,'request_id',uuid.uuid4().hex[:12])
                            appenv.logger().error('request_id=%s path=%s error=%s',rid,self.path,error,exc_info=sys.exc_info()[0] is not None)
                            payload['request_id']=rid
                            if code==500:payload['error']='操作未完成，请查看健康检查及日志，并提供请求编号。'
                            payload['error']=str(payload['error'])+'（请求编号：'+rid+'）'
                        if app_network.failed():
                            for key in list(payload):
                                if key=='error' or key.endswith('_error'):
                                    if payload[key] and app_network.OFFLINE_MESSAGE not in str(payload[key]):payload[key]=str(payload[key])+' '+app_network.OFFLINE_MESSAGE
                        body=json.dumps(payload,ensure_ascii=False).encode('utf-8')
                except (ValueError,TypeError):pass
            elif 'text/html' in ctype and b'</body>' in body and not self.path.startswith('/data/'):
                try:
                    import product_ui
                    import urllib.parse
                    params=urllib.parse.parse_qs(urllib.parse.urlparse(self.path).query)
                    kind='launcher' if 'launcher_port' in state else 'sherloc' if 'scans' in state else params.get('app',['chemcam'])[0]
                    body=product_ui.inject(body,kind,state.get('launcher_url'))
                except ImportError:
                    footer=('<footer style="text-align:center;padding:12px;color:#677">MarsWorkbench '+appenv.build_id()+'</footer>').encode()
                    body=body.replace(b'</body>',footer+b'</body>')
            return super()._send(code,ctype,body)
    return Handler
