"""First-run/help/diagnostics routes shared by the launcher and all workbenches."""

import source_paths  # source roots; no-op for bundled imports
import json
from pathlib import Path
import urllib.parse
import urllib.request
import mimetypes
import html
import appenv
import health
import samples

ROUTES={'/api/help','/api/health','/api/welcome','/api/samples','/api/samples/install','/api/trust','/api/version/check','/api/appearance','/api/paths'}

def _activate(state,selected,folders):
    from instruments import REGISTRY,history
    common=state.get('console_state') or (state if 'app' in state else None)
    if common is not None:
        for name,b in REGISTRY.items():
            current=common['app'][name]
            if name==selected or not current.get('observations'):
                current.update(folder=folders[name],observations=b.scan_folder(folders[name]),filter_token=None)
                history.store_folder(appenv.settings('console'),folders[name],name)
    if selected=='sherloc':
        history.store_folder(appenv.settings('sherloc'),folders['sherloc'])
        if 'scans' in state:
            import view
            _,observations,listing=view.scan_listing(folders['sherloc'])
            state.update(folder=folders['sherloc'],scans=observations,listing=listing,mineral_cache={})
            return None
        ensure=state.get('ensure_sherloc')
        if ensure is not None:
            url=ensure(state.get('sherloc_port',8765))
            with urllib.request.urlopen(url+'api/load?path='+urllib.parse.quote(folders['sherloc']),timeout=30) as response:
                loaded=json.load(response)
            if loaded.get('error'):raise ValueError(loaded['error'])
            return url
        # A common workbench can also open SHERLOC without a system Python.
        import launch
        return launch.ensure_sherloc(8765)
    port=state.get('console_port')
    return ('http://127.0.0.1:%d/?app=%s'%(port,selected)) if port else '?app='+selected

def route(handler,state):
    url=urllib.parse.urlparse(handler.path)
    if url.path.startswith('/manuals/'):
        import manuals
        import product_theme
        relative=urllib.parse.unquote(url.path[len('/manuals/'):])
        if relative in manuals.NAMES:
            page='<!doctype html><html lang="zh"><meta charset="utf-8"><body>'+manuals.render(manuals.read(relative))+'</body></html>'
            handler._send(200,'text/html; charset=utf-8',product_theme.report_style(page).encode());return True
        root=appenv.resources()/'manuals';path=(root/relative).resolve()
        if path.is_relative_to(root.resolve()) and path.is_file() and relative.startswith('插图/') and path.suffix.lower() in ('.png','.jpg','.svg'):
            handler._send(200,mimetypes.guess_type(path)[0] or 'application/octet-stream',path.read_bytes());return True
        handler._send(404,'text/plain',b'manual not found');return True
    if url.path=='/pixl' or url.path.startswith(('/api/pixl/','/data/pixl/')):
        import pixl_view
        return pixl_view.route(handler,state)
    if url.path.startswith('/data/'):
        parts=urllib.parse.unquote(url.path).split('/')
        common=state.get('console_state') or (state if 'app' in state else None)
        instrument=parts[2] if len(parts)>3 else ''
        folder=(state.get('folder') if instrument=='sherloc' and 'scans' in state else
                (common['app'].get(instrument,{}).get('folder') if common else None))
        relative=Path(*parts[3:])
        if folder and parts[3] in ('_qc','_maps','_png'):
            root=Path(folder).resolve();path=(root/relative).resolve();allowed=(root/parts[3]).resolve()
            if path.is_relative_to(allowed) and allowed.is_relative_to(root) and path.is_file() and path.suffix.lower() in ('.html','.png','.svg','.csv','.json'):
                handler._send(200,mimetypes.guess_type(path)[0] or 'application/octet-stream',path.read_bytes());return True
        handler._send(404,'application/json',b'{"error":"Report file not available"}');return True
    if url.path not in ROUTES:return False
    query={k:v[0] for k,v in urllib.parse.parse_qs(url.query).items()}
    data={}
    if handler.command=='POST':
        size=int(handler.headers.get('Content-Length','0'))
        if not 0<size<=65536:raise ValueError('请求大小无效')
        data=json.loads(handler.rfile.read(size))
    def send(value,code=200):
        handler._send(code,'application/json; charset=utf-8',json.dumps(value,ensure_ascii=False).encode('utf-8'))
        return True
    if url.path=='/api/appearance':
        from instruments import history
        path=appenv.settings('appearance')
        if handler.command=='POST':
            changes={k:data[k] for k in ('theme','tone') if k in data}
            if not changes:raise ValueError('请选择主题或黄色强度')
            if 'theme' in changes and changes['theme'] not in ('light','dark','system'):raise ValueError('主题需为light、dark或system')
            if 'tone' in changes and changes['tone'] not in ('bright','soft'):raise ValueError('黄色强度需为bright或soft')
            history.update(path,lambda config:dict(config,**changes))
        config=history.read(path)
        return send({'theme':config.get('theme') or 'light','tone':config.get('tone','bright')})
    if url.path=='/api/help':
        import manuals
        sections=manuals.help_sections()
        for filename,title in [('LICENSE','代码许可'),('THIRD_PARTY_NOTICES.md','第三方许可'),('DATA_SOURCES.md','数据来源与引用')]:
            notice=appenv.resources()/'licenses'/filename
            if notice.exists():sections['about']+='<details><summary>'+title+'</summary><pre>'+html.escape(notice.read_text(encoding='utf-8'))+'</pre></details>'
        return send({'sections':sections,'version':appenv.VERSION,'build':appenv.build_id()})
    if url.path=='/api/health':return send(health.diagnose(handler.server.server_port,query.get('network','1')!='0'))
    if url.path=='/api/welcome':return send(samples.welcome(data.get('mark') if handler.command=='POST' else None))
    if url.path=='/api/samples':
        pixl_meta=json.loads((appenv.resources()/'pixl-demo/manifest.json').read_text(encoding='utf-8'))
        return send({'instruments':dict(samples.NAMES,pixl='PIXL成分扫描'),'destination':str(samples.default_destination()),'bytes':samples.manifest()['total_bytes'],'pixl_bytes':sum(row['bytes'] for row in pixl_meta['files'])})
    if url.path=='/api/samples/install':
        if handler.command!='POST':return send({'error':'请使用加载示例按钮'},405)
        instrument=data.get('instrument','sherloc')
        if instrument=='pixl':
            import pixl_view
            result=pixl_view.install_demo(data.get('destination'),bool(data.get('confirm')))
            if not result['requires_confirmation']:
                pixl_view.load_folder(pixl_view.state_for(state),result['folders']['pixl'])
                samples.welcome('complete')
                result['url']=(state.get('launcher_url') or '/')+'pixl'
            return send(result)
        if instrument not in samples.NAMES:return send({'error':'示例仪器无效'},400)
        result=samples.install(data.get('destination'),bool(data.get('confirm')))
        if not result['requires_confirmation']:
            result['url']=_activate(state,instrument,result['folders'])
            from sherloc.desktop import remember_instrument
            remember_instrument(instrument,appenv.settings('console'))
            samples.welcome('complete')
        return send(result)
    if url.path=='/api/trust':
        import manuals
        return send({'text':manuals.trust_text()})
    if url.path=='/api/paths':
        import project_layout
        instrument=query.get('instrument','sherloc')
        return send({'instrument':instrument,'data_root':str(project_layout.data_root()),'download':appenv.empty_data(instrument)})
    if url.path=='/api/version/check':
        # External publishing was explicitly deferred. Never contact an unconfigured source.
        return send({'enabled':False,'configured':False,'message':'本地预发行版：在线更新默认关闭，尚未配置发布地址。'})
    return False
