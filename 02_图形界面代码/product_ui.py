"""Shared non-blocking onboarding, examples and offline help."""

import source_paths  # source roots; no-op for bundled imports
import json
import appenv

HTML='''<style>
#mwFooter{position:fixed;bottom:8px;right:12px;z-index:30;background:var(--panel,#fff);color:var(--muted,#576170);border:1px solid var(--line,#ccd);border-radius:8px;padding:6px 10px;font:12px/1.5 sans-serif}
#mwFooter button{margin-left:8px;padding:3px 9px;cursor:pointer}#mwDialog{width:min(860px,92vw);max-height:85vh;overflow:auto;border:1px solid var(--line,#aab5c5);border-radius:12px;padding:20px;color:var(--ink,#162235);background:var(--panel,#fff);font:14px/1.7 "Segoe UI","Microsoft YaHei",sans-serif}
#mwDialog::backdrop{background:#10203088}#mwDialog button{padding:6px 12px;margin:3px;cursor:pointer}#mwDialog input,#mwDialog select{padding:6px;max-width:95%}#mwClose{float:right}#mwContent pre{white-space:pre-wrap;overflow-wrap:anywhere;font:12px/1.6 monospace}#mwContent h2{margin-top:0}#mwContent details{padding:7px 0;border-bottom:1px solid #ddd}#mwContent section{margin:14px 0}#mwWelcome{max-width:700px;width:88vw;background:var(--panel,#fff);border:1px solid var(--accent,#2465a8);border-radius:12px;padding:16px;margin:15px 0;text-align:left}#mwWelcome ol{padding-left:24px}#mwWelcome button{padding:6px 12px;margin-right:8px;cursor:pointer}#mwStatus{white-space:pre-wrap}
#folderPath.mwSampleAlias{font-size:0}#folderPath.mwSampleAlias::before{content:attr(data-mw-label);font-size:12px}
</style><div id="mwFooter"><span id="mwBuildLabel" title="MarsWorkbench __BUILD__">MarsWorkbench __BUILD__</span><button id="mwSamples">加载示例数据</button><button id="mwHelp">帮助 / 关于</button></div>
<dialog id="mwDialog"><button id="mwClose">关闭</button><nav id="mwTabs"></nav><p id="mwStatus" role="status"></p><div id="mwContent"></div></dialog>'''

JS=r'''
(function(){
 const config=__CONFIG__,dialog=document.getElementById('mwDialog'),content=document.getElementById('mwContent'),status=document.getElementById('mwStatus'),tabs=document.getElementById('mwTabs');
 function element(tag,text){var e=document.createElement(tag);if(text!=null)e.textContent=text;return e;}
 async function request(url,data){var r=await fetch(url,data?{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(data)}:undefined);var d=await r.json();if(d.error)throw Error(d.error);return d;}
 function open(){if(!dialog.open)dialog.showModal();status.textContent='';content.replaceChildren();tabs.replaceChildren();}
 function fail(error){status.textContent='未完成：'+error.message;}
 document.getElementById('mwClose').onclick=()=>dialog.close();dialog.addEventListener('keydown',e=>e.stopPropagation());
 async function showSamples(){open();content.appendChild(element('h2','加载示例数据'));status.textContent='正在准备示例信息…';
  try{var d=await request('/api/samples');status.textContent='公开PDS示例会复制到可写数据目录；之后可以独立分析和标注。';
   var select=element('select');select.id='mwSampleInstrument';Object.keys(d.instruments).forEach(k=>select.add(new Option(d.instruments[k],k)));if(config.instrument!=='launcher'){select.value=config.instrument;select.disabled=true;}
   var label=element('label','打开仪器 ');label.appendChild(select);content.appendChild(label);content.appendChild(element('br'));
   var target=element('input');target.id='mwSampleDestination';target.value=d.destination;target.style.width='90%';content.appendChild(element('p','示例保存位置（可修改）：'));content.appendChild(target);
   var note=element('p');function sampleNote(){note.textContent=(select.value==='pixl'?'PIXL示例约 '+(d.pixl_bytes/1048576).toFixed(1):'完整光谱示例约 '+(d.bytes/1048576).toFixed(1))+' MB。首次无需网络；额外图像按需联网。';}select.onchange=sampleNote;sampleNote();content.appendChild(note);
   var button=element('button','复制示例并打开');button.id='mwInstallSamples';content.appendChild(button);
   button.onclick=async function(){button.disabled=true;status.textContent='正在校验并复制示例，请稍候…';try{
    var body={instrument:select.value,destination:target.value,confirm:false};var result=await request('/api/samples/install',body);
    if(result.requires_confirmation){if(!confirm('该位置已有'+result.existing_files+'个同名示例文件。只覆盖这些示例文件，保留其他文件及标注，是否继续？')){status.textContent='已取消覆盖';button.disabled=false;return;}body.confirm=true;result=await request('/api/samples/install',body);}
    status.textContent='示例已就绪，正在打开工作台…';if(result.url)location.href=result.url;else location.reload();
   }catch(e){fail(e);button.disabled=false;}};
  }catch(e){fail(e);}}
 async function showHelp(){open();try{var d=await request('/api/help');var names={quickstart:'快速开始',faq:'常见问题',scope:'科学口径',about:'关于'};
   function section(key){content.innerHTML=d.sections[key];status.textContent=key==='about'?'版本 '+d.version+' · 构建 '+d.build:'';
    if(key==='about'){
     var health=element('button','健康检查（最多约10秒联网）');health.id='mwHealth';content.appendChild(health);
     var trust=element('button','可信度与发行自检');trust.id='mwTrust';content.appendChild(trust);
     var welcome=element('button','下次显示首跑向导');content.appendChild(welcome);welcome.onclick=()=>request('/api/welcome',{mark:'pending'}).then(()=>status.textContent='已设置，下次打开启动器显示向导').catch(fail);
     health.onclick=async function(){health.disabled=true;status.textContent='正在检查本机环境与PDS连接…';try{var result=await request('/api/health');var pre=element('pre',result.text);pre.id='mwHealthText';content.appendChild(pre);var copy=element('button','复制诊断文本');content.appendChild(copy);copy.onclick=()=>navigator.clipboard.writeText(result.text).then(()=>status.textContent='已复制').catch(()=>status.textContent='可选中文本后按Ctrl+C复制');status.textContent=result.status==='ok'?'检查完成':result.status==='error'?'检查完成：必需环境异常，请先处理诊断中的问题':'检查完成：有可选组件或联网提示';}catch(e){fail(e);}finally{health.disabled=false;}};
     trust.onclick=async function(){try{var result=await request('/api/trust');var pre=element('pre',result.text);pre.id='mwTrustText';content.appendChild(pre);}catch(e){fail(e);}};
    }}
   Object.keys(names).forEach(k=>{var b=element('button',names[k]);b.id='mwTab_'+k;b.onclick=()=>section(k);tabs.appendChild(b);});section('quickstart');
  }catch(e){fail(e);}}
 document.getElementById('mwSamples').onclick=showSamples;document.getElementById('mwHelp').onclick=showHelp;
 function applyAppearance(theme,tone){if(theme==='system')document.documentElement.removeAttribute('data-theme');else if(theme)document.documentElement.setAttribute('data-theme',theme);if(tone)document.documentElement.setAttribute('data-tone',tone);document.querySelectorAll('[data-mw-theme-choice]').forEach(b=>b.setAttribute('aria-pressed',String(b.dataset.mwThemeChoice===theme)));document.querySelectorAll('[data-mw-tone-choice]').forEach(b=>b.setAttribute('aria-pressed',String(b.dataset.mwToneChoice===tone)));document.dispatchEvent(new Event('workbench-theme'));}
 var appearanceVersion=0,appearanceBusy=false,appearanceWrites=Promise.resolve();
 async function refreshAppearance(){if(appearanceBusy)return;var tick=appearanceVersion;try{var d=await request('/api/appearance');if(tick===appearanceVersion&&!appearanceBusy)applyAppearance(d.theme,d.tone);}catch(e){}}
 async function saveAppearance(change){var tick=++appearanceVersion;applyAppearance(change.theme,change.tone);appearanceBusy=true;appearanceWrites=appearanceWrites.catch(()=>{}).then(()=>request('/api/appearance',change));try{var d=await appearanceWrites;if(tick===appearanceVersion)applyAppearance(d.theme,d.tone);}catch(error){status.textContent='配色暂未保存：'+error.message;}finally{if(tick===appearanceVersion)appearanceBusy=false;}}
 var themeButton=document.getElementById('themeBtn');if(themeButton)themeButton.addEventListener('click',function(e){e.preventDefault();e.stopImmediatePropagation();var active=document.documentElement.getAttribute('data-theme')||(matchMedia('(prefers-color-scheme:dark)').matches?'dark':'light');saveAppearance({theme:active==='dark'?'light':'dark'});},{capture:true});
 async function showAppearance(){open();content.appendChild(element('h2','配色设置'));content.appendChild(element('p','主题与黄色强度在各工作台之间共用。'));var options=element('div');options.id='mwAppearanceOptions';[['light','浅色'],['dark','深色'],['system','跟随系统']].forEach(([value,label])=>{var b=element('button',label);b.dataset.mwThemeChoice=value;b.setAttribute('aria-pressed','false');b.onclick=()=>saveAppearance({theme:value});options.append(b);});content.append(options);content.appendChild(element('h3','黄色强度'));var tones=element('div');tones.className='bar';[['bright','亮黄'],['soft','柔黄']].forEach(([value,label])=>{var b=element('button');b.dataset.mwToneChoice=value;b.setAttribute('aria-pressed','false');var swatch=element('i');swatch.className='mwToneSwatch';swatch.dataset.value=value;b.append(swatch,document.createTextNode(label));b.onclick=()=>saveAppearance({tone:value});tones.append(b);});content.append(tones);await refreshAppearance();}
 refreshAppearance();window.addEventListener('focus',refreshAppearance);document.addEventListener('visibilitychange',()=>{if(!document.hidden)refreshAppearance();});
 if(config.instrument!=='launcher'&&config.launcherURL&&!document.getElementById('backBtn')){var back=element('button','选择仪器');var top=document.querySelector('.top');back.id=top?'backBtn':'mwBack';back.className='ghost';back.onclick=()=>location.href=config.launcherURL;if(top)top.prepend(back);else document.getElementById('mwFooter').appendChild(back);}
 var folder=document.getElementById('folderPath');if(folder){function labelFolder(){var path=folder.textContent,example=path.includes('MarsWorkbench示例');folder.classList.toggle('mwSampleAlias',example);folder.dataset.mwLabel='示例数据 · '+config.instrument.toUpperCase();folder.title=path;}labelFolder();new MutationObserver(labelFolder).observe(folder,{childList:true,characterData:true,subtree:true});}
 if(config.launcherURL&&document.getElementById('backBtn'))document.getElementById('backBtn').addEventListener('click',function(e){e.preventDefault();e.stopImmediatePropagation();location.href=config.launcherURL;},{capture:true});
 if(config.instrument==='launcher')request('/api/last').then(d=>{if(!d.last)return;document.querySelectorAll('a.card').forEach(card=>{if(card.textContent.toLowerCase().includes(d.last)){card.style.borderColor='var(--ef-trace,#161A1D)';card.setAttribute('aria-label',card.textContent+' · 最近打开');}});}).catch(()=>{});
 if(config.instrument==='launcher')request('/api/welcome').then(d=>{if(!d.show)return;var welcome=element('section');welcome.id='mwWelcome';welcome.innerHTML='<h2>欢迎使用火星光谱工作台</h2><ol><li>加载随包示例，或选择自己的数据</li><li>选择观测和点，查看光谱</li><li>导出谱图与质控报告，记录复核结论</li></ol>';
  var load=element('button','加载示例数据');load.id='mwWelcomeSamples';load.onclick=showSamples;welcome.appendChild(load);var skip=element('button','跳过，直接使用');skip.id='mwWelcomeSkip';welcome.appendChild(skip);skip.onclick=()=>request('/api/welcome',{mark:'skipped'}).then(()=>welcome.remove()).catch(fail);
  var grid=document.querySelector('.grid');if(grid)grid.parentNode.insertBefore(welcome,grid);else document.body.prepend(welcome);
 }).catch(()=>{});
 function setupChrome(){var rail=document.getElementById('mwRail');if(!rail)return;rail.addEventListener('keydown',e=>e.stopPropagation());function control(id,target,key,label){var b=document.getElementById(id),slot=document.getElementById(target);if(!b||!slot)return;b.className='mwRailControl';b.innerHTML=config.icons[key]+'<span>'+label+'</span>';slot.append(b);}control('themeBtn','mwRailTheme','moon','主题');control('mwHelp','mwRailHelp','help','帮助');document.getElementById('mwPalette').onclick=showAppearance;
  var top=document.querySelector('.top');if(top){var actions=element('div');actions.className='mwHeadActions';['browseBtn','dlOpenBtn'].forEach(id=>{var b=document.getElementById(id);if(b){if(id==='browseBtn')b.textContent='打开本地数据';actions.append(b);}});top.append(actions);}
  var row=document.querySelector('main>.row'),scan=document.getElementById('scan');if(!row||!scan)return;row.classList.add('mwWorkbenchRow');var point=row.firstElementChild,plot=row.lastElementChild;point.classList.add('mwPointPanel');plot.classList.add('mwPlotPanel');function heading(parent,chinese,english){var h=element('h2',chinese);h.className='mwPanelHeading';h.append(element('small',english));parent.prepend(h);}heading(point,'当前点位','SELECTED POINT');heading(plot,config.instrument==='pixl'?'MCC 图像与扫描点':'光谱与证据',config.instrument==='pixl'?'IMAGE & POINTS':'SPECTRUM & EVIDENCE');
  var observation=element('section');observation.id='mwObservation';observation.className='card mwObservation';var field=element('div');field.className='mwObservationField';var label=scan.previousElementSibling;if(label&&label.classList.contains('lbl'))field.append(label);scan.setAttribute('aria-label','观测');field.append(scan);observation.append(field);var filter=document.getElementById('scanFilter');if(filter){var f=element('div');f.className='mwObservationFilter';var l=filter.previousElementSibling;if(l&&l.tagName==='LABEL'&&l.htmlFor==='scanFilter'){l.textContent='筛选观测';f.append(l);}else{l=element('label','筛选观测');l.htmlFor='scanFilter';l.className='lbl';f.append(l);}f.append(filter);var info=document.getElementById('filterInfo');if(info)f.append(info);observation.append(f);}var ptlabel=document.getElementById('ptlabel');if(ptlabel){var meta=element('div');meta.className='mwObservationMeta';function mirror(){meta.textContent=ptlabel.textContent;}mirror();new MutationObserver(mirror).observe(ptlabel,{childList:true,characterData:true,subtree:true});observation.append(meta);}row.before(observation);
  if(config.instrument==='pixl'){var composition=document.getElementById('composition'),metric=document.getElementById('element');function mark(){if(!composition||!metric)return;composition.querySelectorAll('tr').forEach(tr=>tr.classList.toggle('mwMetricRow',tr.firstElementChild?.textContent===metric.value));}if(composition)new MutationObserver(mark).observe(composition,{childList:true,subtree:true});if(metric)metric.addEventListener('change',mark);mark();}
 }
 setupChrome();
})();
'''

def inject(body,instrument='launcher',launcher_url=None):
    if b'id="mwFooter"' in body or b'data-mw-surface="report"' in body:return body
    from product_theme import decorate,rail_html,icon
    html=HTML.replace('__BUILD__',appenv.build_id())
    script=JS.replace('__CONFIG__',json.dumps({'instrument':instrument,'launcherURL':launcher_url,'icons':{name:icon(name) for name in ('moon','help')}},ensure_ascii=False))
    import re
    body=re.sub(br'(<body\b[^>]*>)',lambda m:m[0]+rail_html(instrument,launcher_url).encode('utf-8'),body,count=1,flags=re.I)
    body=body.replace(b'</body>',(html+'<script>'+script+'</script></body>').encode('utf-8'))
    return decorate(body,instrument)
