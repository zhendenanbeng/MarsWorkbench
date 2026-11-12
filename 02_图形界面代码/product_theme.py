"""Approved yellow-paper shell. Scientific palettes, values and geometry are retained."""

import source_paths  # source roots; no-op for bundled imports
import base64
import html
import json
import math
import re
from ui_font import FONT_B64

ACCENT='#FFFA00'
VARIANT={'variant':'yellow-paper','accent':'yellow','paper':'neutral','cut':'tight','ground':'contours','density':'standard','motion':'calm'}
LABELS={'launcher':'DATA WORKBENCH','sherloc':'01 / SHERLOC','supercam':'02 / SUPERCAM','pixl':'03 / PIXL','chemcam':'04 / CHEMCAM','apxs':'05 / APXS','report':'DATA REPORT'}
INSTRUMENTS=[('sherloc','拉曼','SHERLOC','spectrum'),('supercam','激光谱','SUPERCAM','laser'),('pixl','微区成分','PIXL','grid'),('chemcam','激光谱','CHEMCAM','laser'),('apxs','X 射线','APXS','atoms')]
_ICONS={
'spectrum':'<path d="M3 18h18M4 15h3l2-10 3 11 3-6 2 5h3"/>',
'laser':'<path d="m4 20 7-7m-3 0 3 3m1-12v4m8 4h-4m2-6-3 3M4 7l3 2M17 20l-2-3m-4-7 3 3-3 3-3-3Z"/>',
'grid':'<path d="M4 4h6v6H4Zm10 0h6v6h-6ZM4 14h6v6H4Zm10 0h6v6h-6Z"/>',
'atoms':'<circle cx="12" cy="12" r="2"/><ellipse cx="12" cy="12" rx="10" ry="4" transform="rotate(55 12 12)"/><ellipse cx="12" cy="12" rx="10" ry="4" transform="rotate(-55 12 12)"/>',
'moon':'<path d="M20 15A8.5 8.5 0 0 1 9 4a8.5 8.5 0 1 0 11 11Z"/>',
'help':'<circle cx="12" cy="12" r="9"/><path d="M9 9a3 3 0 1 1 5 2c-1.5 1-2 1-2 3m0 3v.1"/>',
'palette':'<path d="M4 6h16M4 12h16M4 18h16M8 3v6m8 0v6m-5 0v6"/>'}

def icon(name):
    return '<svg class="mwIcon" viewBox="0 0 24 24" aria-hidden="true">'+_ICONS[name]+'</svg>'

def rail_html(instrument,launcher_url=None):
    base=(launcher_url or '/').rstrip('/')+'/'
    home='<a class="mwBrand" href="'+html.escape(base,quote=True)+'" aria-label="选择仪器"><span class="mwMark">M</span><span>MARS / WB</span></a>'
    links=[]
    for key,chinese,english,symbol in INSTRUMENTS:
        current=' aria-current="page"' if key==instrument else ''
        links.append('<a class="mwInstrument" data-instrument="'+key+'" href="'+html.escape(base+'open/'+key,quote=True)+'"'+current+' aria-label="'+english+' '+chinese+'">'+icon(symbol)+'<span>'+chinese+'</span><small>'+english+'</small></a>')
    return '<aside id="mwRail" aria-label="仪器导航">'+home+'<nav>'+''.join(links)+'</nav><div class="mwRailBottom"><div id="mwRailHelp"></div><button id="mwPalette" class="mwRailControl" aria-label="配色设置">'+icon('palette')+'<span>配色</span></button><div id="mwRailTheme"></div><small class="mwRailVersion">MARSWORKBENCH</small></div></aside>'

def _contours():
    paths=[]
    for i in range(14):
        points=['%s%.1f,%d'%('M' if y==-30 else 'L',190+i*24+38*math.sin(y/90+i*.19)+20*math.sin(y/39+i*.13),y) for y in range(-30,600,12)]
        paths.append('<path d="'+' '.join(points)+'"/>')
    svg='<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 640 540"><g fill="none" stroke="#777772" stroke-width="1">'+''.join(paths)+'</g></svg>'
    return base64.b64encode(svg.encode()).decode()

FONT_CSS="@font-face{font-family:Saira;src:url(data:font/woff2;base64,"+FONT_B64+") format('woff2');font-weight:100 900;font-stretch:50% 125%;font-display:swap}"
TOKEN_CSS=r'''
:root{--bg:#FEFEFD;--panel:#F1F1EF;--panel2:#E8E8E5;--field:#FFFFFF;--ink:#191919;--muted:#626260;--line:#DDDDDA;--line2:#E5E5E1;--accent:#FFFA00;--accent-ink:#161616;--signal-soft:#FFFEC9;--ok:#42423D;--warn:#42423D;--err:#42423D;--radius:0;--ef-signal:#FFFA00;--ef-on-signal:#161616;--ef-trace:#191919;--ef-cut:7px;--ef-motion:160ms;--ef-mech:cubic-bezier(.7,0,.2,1);--ef-font:Saira,'Microsoft YaHei','Noto Sans SC','PingFang SC',sans-serif;color-scheme:light}
@media(prefers-color-scheme:dark){:root{--bg:#1C1C1B;--panel:#272726;--panel2:#343432;--field:#222221;--ink:#EEEEEA;--muted:#B5B5AC;--line:#42423C;--line2:#373733;--signal-soft:#454424;--ok:#D8D8CF;--warn:#D8D8CF;--err:#D8D8CF;--ef-trace:#EEEEEA;color-scheme:dark}}
html[data-theme="light"]{--bg:#FEFEFD;--panel:#F1F1EF;--panel2:#E8E8E5;--field:#FFFFFF;--ink:#191919;--muted:#626260;--line:#DDDDDA;--line2:#E5E5E1;--signal-soft:#FFFEC9;--ok:#42423D;--warn:#42423D;--err:#42423D;--ef-trace:#191919;color-scheme:light}
html[data-theme="dark"]{--bg:#1C1C1B;--panel:#272726;--panel2:#343432;--field:#222221;--ink:#EEEEEA;--muted:#B5B5AC;--line:#42423C;--line2:#373733;--signal-soft:#454424;--ok:#D8D8CF;--warn:#D8D8CF;--err:#D8D8CF;--ef-trace:#EEEEEA;color-scheme:dark}
html[data-tone="soft"]{--ef-signal:#EBE45B;--signal-soft:#F5F3D5}html[data-tone="soft"][data-theme="dark"]{--signal-soft:#434128}
@media(prefers-color-scheme:dark){html[data-tone="soft"]:not([data-theme="light"]){--signal-soft:#434128}}
html:root{--accent:var(--ef-signal);--accent-ink:var(--ef-on-signal)}
'''
SHELL_CSS=r'''
html{scroll-padding-top:20px}body{font-family:var(--ef-font);font-variant-numeric:tabular-nums;background:var(--bg);color:var(--ink);position:relative;isolation:isolate;height:auto;min-height:100svh}html[data-mw-surface]:not([data-mw-surface="report"]) body{padding-left:86px}
*,*::before,*::after{box-sizing:border-box}button,input,select,textarea{font-family:var(--ef-font);font-variant-numeric:tabular-nums}button{font-weight:600;font-stretch:112%}body::before{display:none}body::after{content:'';position:fixed;right:-100px;top:-40px;width:640px;height:540px;background:var(--mw-contours) no-repeat center/contain;opacity:.09;z-index:-1;pointer-events:none}
.card,.panel,.top,.cmdbar,#mwDialog,#mwFooter,button,input,select,textarea,code.path,.prog,.dllog,.charttip{border-radius:0!important;box-shadow:none!important}
.mwIcon{width:19px;height:19px;fill:none;stroke:currentColor;stroke-width:1.6;stroke-linecap:square;stroke-linejoin:miter;display:block;flex:none}
#mwRail{position:fixed;inset:0 auto 0 0;width:86px;background:var(--panel);border-right:1px solid var(--line);display:flex;align-items:center;flex-direction:column;z-index:35;overflow-y:auto;scrollbar-width:thin}.mwBrand{margin:21px 0 24px;display:flex;flex-direction:column;align-items:center;text-decoration:none!important;font-size:8px;letter-spacing:.1em;gap:6px;color:var(--ink)}.mwMark{background:var(--accent);color:var(--accent-ink);display:grid;place-items:center;width:40px;height:42px;font-size:25px;font-weight:850;font-stretch:125%;clip-path:polygon(0 0,100% 0,100% 85%,85% 100%,0 100%)}
#mwRail nav{display:grid;gap:6px;width:100%;padding:0 7px}.mwInstrument{position:relative;display:flex;flex-direction:column;align-items:center;justify-content:center;gap:4px;min-height:65px;color:var(--muted);font-size:11px;text-decoration:none!important;border-left:3px solid transparent}.mwInstrument small{font-size:8px;font-stretch:125%;letter-spacing:.04em}.mwInstrument:hover{background:var(--panel2);color:var(--ink)}.mwInstrument[aria-current="page"]{background:#191919;color:white;border-left-color:var(--accent);font-weight:700}.mwInstrument[aria-current="page"]:after{content:'';position:absolute;width:6px;height:6px;right:8px;top:8px;background:var(--accent)}
.mwRailBottom{margin-top:auto;display:grid;gap:5px;width:100%;padding:14px 7px 12px}.mwRailBottom>div{display:contents}.mwRailControl,#mwRail button{display:flex;flex-direction:column;align-items:center;justify-content:center;gap:4px;min-height:47px;width:100%;padding:5px!important;background:transparent;border:0;clip-path:none;font-size:10px;color:var(--muted)}.mwRailControl:hover,#mwRail button:hover{background:var(--accent);color:var(--accent-ink)}.mwRailVersion{text-align:center;color:var(--muted);font-size:7px;letter-spacing:.05em;margin-top:5px}
.top{position:relative!important;top:0;z-index:10;background:transparent!important;backdrop-filter:none!important;border:0!important;padding:20px 30px 17px!important;display:grid!important;grid-template-columns:minmax(0,1fr) auto;gap:6px 20px;align-items:end;max-width:1600px;margin:auto}
.top h1{font-size:28px!important;font-weight:800;line-height:1.3;margin:0;grid-column:1;grid-row:1;border-left:4px solid var(--accent);padding-left:13px;min-width:0;overflow-wrap:anywhere}.top h1:before{display:none!important}.top .sub{grid-column:1;grid-row:2;color:var(--muted);font-size:11px;line-height:1.6;margin:0 0 0 17px}.mw-kicker{display:block;font-size:12px;font-weight:550;font-stretch:125%;letter-spacing:.12em;line-height:1.5;color:var(--muted);margin-bottom:8px;text-transform:uppercase}.top .right,.top>#backBtn,.top #themeBtn{display:none}
.mwHeadActions{grid-column:2;grid-row:1/3;display:flex;gap:7px;align-self:end}.mwHeadActions button{font-size:12px;padding:8px 14px;background:var(--panel2);border:0;white-space:nowrap}
main{max-width:1600px!important;margin:0 auto;padding:0 30px 70px!important}.card,.panel{border:0!important;background:var(--panel);padding:15px 17px;margin-bottom:16px}main>.card:first-child{background:var(--panel);padding:12px 16px}.bar{gap:9px}.row>*{min-width:0}.bar>*{max-width:100%}.bar .grow{min-width:0}.row{gap:20px}
.mwWorkbenchRow{grid-template-columns:minmax(0,1fr) 300px!important}.mwWorkbenchRow>.mwPointPanel{grid-column:2;grid-row:1}.mwWorkbenchRow>.mwPlotPanel{grid-column:1;grid-row:1;position:relative}.mwPlotPanel:before{content:'';height:3px;width:40px;position:absolute;left:17px;top:0;background:var(--ink)}
.mwPanelHeading{font-size:14px!important;font-weight:750;margin:0 0 13px!important;display:flex;align-items:baseline;justify-content:space-between;gap:8px}.mwPanelHeading small{font-size:9px;font-weight:500;letter-spacing:.12em;font-stretch:125%;color:var(--muted)}
.mwObservation{display:flex;align-items:end;gap:18px;flex-wrap:wrap;padding:11px 16px!important}.mwObservationField{min-width:0;flex:2 1 280px}.mwObservationFilter{min-width:0;flex:1 1 200px;max-width:430px}.mwObservation .lbl{font-size:10px;margin:0 0 4px}.mwObservation #scan,.mwObservation #scanFilter{width:100%!important;min-width:0;margin:0!important}.mwObservation #filterInfo{margin:3px 0 0;font-size:10px;line-height:1.5}.mwObservation #filterInfo:empty{display:none}.mwObservationMeta{margin:0 0 5px auto;border-left:1px solid var(--line);padding-left:17px;font-size:11px;color:var(--muted);min-width:140px;font-stretch:112%}
.mwPointPanel #ptlabel{font-stretch:112%;font-size:12px!important;margin:6px 0 12px!important;color:var(--ink)!important}.mwPointPanel h2{font-size:13px}.mwPointPanel .stepper button{padding:5px 0}.mwPointPanel .bar input[type=number]{min-width:0;flex:1;width:75px!important}.mwPointPanel .log{font-size:11px;max-height:170px}
button,.iconbtn,.buse,a.card{clip-path:polygon(0 0,100% 0,100% calc(100% - var(--ef-cut)),calc(100% - var(--ef-cut)) 100%,0 100%);padding-right:calc(12px + var(--ef-cut));transition:background-color var(--ef-motion) var(--ef-mech),color var(--ef-motion) var(--ef-mech),border-color var(--ef-motion) var(--ef-mech)}
button{min-height:35px;background:var(--panel2);color:var(--ink);border:1px solid transparent;font-size:12px}button:hover,.iconbtn:hover,.buse:hover{background:var(--accent);color:var(--accent-ink);border-color:transparent;filter:none!important}button.primary,.buse{background:var(--ink);color:var(--bg);border-color:transparent}button.primary:hover{background:var(--accent);color:var(--accent-ink)}#qcBtn,#mwWelcomeSamples{background:var(--accent);color:var(--accent-ink);font-weight:750;border-color:transparent}button.ghost{background:transparent;color:var(--muted)}button.ghost:hover{background:var(--accent);color:var(--accent-ink)}button:disabled{color:var(--muted);opacity:.55;cursor:default;background:repeating-linear-gradient(135deg,var(--panel2) 0 5px,var(--line2) 5px 6px)}
button:focus-visible,.iconbtn:focus-visible,a:focus-visible{outline:2px solid currentColor;outline-offset:-4px}input:focus-visible,select:focus-visible,textarea:focus-visible{outline:2px solid var(--ink);outline-offset:2px}input,select,textarea{background:var(--field);color:var(--ink);border:1px solid var(--line);font-size:12px}select,input[type=text],input[type=search]{min-height:34px}input[type=checkbox],input[type=range]{accent-color:var(--ink)}input[type=range]:focus-visible{outline-color:var(--ink)}
a{color:var(--ink);text-decoration-thickness:1px;text-underline-offset:3px}a:hover{text-decoration-thickness:2px}.lbl,.hint,.sub{color:var(--muted)}.lbl{font-size:11px}.hint{font-size:11px;line-height:1.7}.tage,.tagm{background:var(--signal-soft);color:var(--ink);border:0;border-radius:0!important}.dot{color:var(--ink)}.stepper button:hover{color:var(--accent-ink);background:var(--accent);border-color:transparent}
.msg.err,[role=alert]:not([hidden]),#warnings:not(:empty){background:var(--panel2);color:var(--ink);border-left:5px solid var(--accent);padding:8px 10px}.msg.ok{color:var(--ink);border-left:3px solid var(--ink);padding-left:8px}.prog,.prog>i{border-radius:0}.prog>i{background:var(--ink)}.dllog,.log{color:var(--muted)}summary{white-space:normal;color:var(--muted)}summary:hover{color:var(--ink)}
.msg{max-width:100%;overflow-wrap:anywhere;word-break:normal}
code.path{background:var(--field);border:0;font:11px/1.6 var(--ef-font);overflow-wrap:anywhere}.path:empty{min-height:30px}#folderPath{min-width:0}.charttip{font-family:var(--ef-font)!important;border:1px solid var(--line);color:var(--ink);background:var(--field)}.ev,.oxtable,.elebox{min-width:0;max-width:100%;overflow-wrap:anywhere}.ev{font-size:11px}.legend{font-family:var(--ef-font);font-size:11px}
/* Existing scientific category palettes and the quantitative heatmap remain unchanged. */
svg [stroke="var(--accent)"]{stroke:var(--ef-trace)}svg [fill="var(--accent)"]{fill:var(--ef-trace)}.legend i[style*="var(--accent)"]{background:var(--ef-trace)!important}#zoomrect{border-color:var(--ink)!important}
html[data-mw-instrument="sherloc"][data-theme="dark"] #plot svg,html[data-mw-instrument="sherloc"][data-theme="dark"] .legend i{filter:brightness(1.2)}
#compositionBox{max-height:285px!important;overflow:auto;scrollbar-width:thin}#compositionBox table{font-size:11px}#compositionBox th{position:sticky;top:0;background:#252525;color:white;padding:6px 7px}#compositionBox td{padding:6px 7px}.mwMetricRow td{font-weight:700}.mwMetricRow td:first-child{background:var(--accent);color:var(--accent-ink)}#viewport{height:300px!important;background:#20201f;border:0!important}#legend{font-family:var(--ef-font);font-size:10px}#pointNote{font-size:11px}
#mwFrameTable{margin-top:16px;border-top:1px solid var(--line);padding-top:13px}.mwTableScroll{overflow:auto;max-width:100%;scrollbar-width:thin}.mwDataTable{font-size:11px;min-width:620px;width:100%;border-collapse:collapse}.mwDataTable th{background:#252525;color:white;font-size:10px;font-weight:550;padding:7px 10px;white-space:nowrap}.mwDataTable td{padding:6px 10px;border-bottom:1px solid var(--line);white-space:nowrap;text-align:right;background:var(--field)}.mwDataTable td:first-child{text-align:left}.mwDataTable tbody tr{cursor:pointer}.mwDataTable tbody tr:hover td{background:var(--panel2)}.mwDataTable .mwMetricCell{background:var(--signal-soft)}.mwDataTable th.mwMetricCell{background:var(--accent);color:var(--accent-ink)}.mwDataTable tr[aria-selected="true"] td{background:var(--signal-soft)}.mwDataTable tr[aria-selected="true"] td:first-child{background:var(--accent);color:var(--accent-ink);font-weight:750;border-left:3px solid var(--ink)}.mwTableFooter{display:flex;justify-content:space-between;gap:10px;align-items:center;margin-top:8px;font-size:10px;color:var(--muted)}.mwTablePager{display:flex;align-items:center;gap:7px}.mwTablePager button{min-height:27px;min-width:27px;padding:3px 8px}
#mwFooter{position:fixed;bottom:0;left:86px;right:0;z-index:30;display:flex;align-items:center;gap:10px;max-width:none;border:0;border-top:1px solid var(--line);background:var(--panel);padding:4px 30px;font:10px/1.6 var(--ef-font)!important;min-height:32px}#mwBuildLabel{overflow:hidden;text-overflow:ellipsis;white-space:nowrap;min-width:0}#mwFooter #mwSamples{margin-left:auto;font-size:10px;min-height:26px;padding:3px 9px}#mwFooter #mwHelp{font-size:10px}
#mwDialog{border:1px solid var(--line);border-top:5px solid var(--accent);background:var(--field);color:var(--ink);padding:22px;font-family:var(--ef-font);width:min(860px,calc(100vw - 24px));max-height:85dvh;overflow:auto}#mwDialog::backdrop{background:#19191999}#mwContent pre{font-family:var(--ef-font);border-top:1px solid var(--line);padding-top:10px}#mwContent h2{font-weight:750}#mwTabs{display:flex;flex-wrap:wrap;gap:4px}#mwTabs button{margin:0}#mwStatus{overflow-wrap:anywhere}#mwWelcome{width:100%;max-width:100%;margin:16px 0;padding:16px;border:0;border-left:5px solid var(--accent);background:var(--panel);box-shadow:none}
#mwAppearanceOptions{display:flex;gap:7px;flex-wrap:wrap;margin-bottom:18px}#mwAppearanceOptions button[aria-pressed="true"]{background:var(--ink);color:var(--bg)}.mwToneSwatch{display:inline-block;width:11px;height:11px;background:var(--accent);margin-right:7px;outline:1px solid #16161644}.mwToneSwatch[data-value="bright"]{background:#FFFA00}.mwToneSwatch[data-value="soft"]{background:#EBE45B}
.cmdbar{position:relative!important;top:0!important;left:auto!important;right:auto!important;height:46px!important;background:transparent!important;color:var(--ink)!important;border:0!important;border-bottom:1px solid var(--line)!important;backdrop-filter:none!important}.cmdbar .brand{color:var(--ink)!important;padding-left:30px!important;text-decoration:none}.cmdbar .brand .seal{background:var(--accent)!important;color:var(--accent-ink)!important;border:0}.cmdbar .ver{color:var(--muted)!important;font-size:10px;letter-spacing:.08em}
html[data-mw-surface="launcher"] body{display:block;height:auto;min-height:100svh}html[data-mw-surface="launcher"] main{padding-top:26px!important;max-width:1360px!important}
.hero{display:grid;grid-template-columns:1.1fr 1fr;column-gap:40px;align-items:center;margin-bottom:24px!important;padding-bottom:22px!important;border-bottom:1px solid var(--line)}.hero>.eyebrow{grid-column:1/-1;color:var(--muted)!important;background:transparent!important}.hero .wordmark{font-size:clamp(36px,4.4vw,64px);max-width:100%;font-stretch:125%;line-height:1.05}.hero .wordmark .hollow{-webkit-text-stroke:0;color:var(--ink)}.hero .sub{margin:8px 0 0;font-size:13px;line-height:1.9}.hero .meta-row{grid-column:1/-1;column-gap:28px;row-gap:10px;margin-top:20px}.hero .meta-row b{font-stretch:125%;font-weight:650}
.grid{background:var(--line);border:0}.grid .rover{color:var(--muted);background:var(--panel2);font-size:11px}a.card{background:var(--panel);border:0!important;border-left:4px solid transparent!important;box-shadow:none;min-width:0;padding:20px 18px 22px}a.card:hover,a.card[aria-current="page"]{background:var(--signal-soft);border-left-color:var(--ink)!important}a.card .no{-webkit-text-stroke:0!important;color:var(--muted)!important;font-size:22px;opacity:1!important;right:18px!important;top:18px!important}a.card .nm{padding-right:48px;font-size:17px;font-weight:750}a.card .ds{line-height:1.7;font-size:11px}.pair-head{flex-wrap:wrap}.pair-head .cn{font-size:23px;font-weight:800}.pair-head .eyebrow{color:var(--muted)}
html[data-mw-surface="launcher"] main>a.card{min-height:0;padding:12px 18px;background:var(--signal-soft)}html[data-mw-surface="launcher"] #mwWelcome h2{font-size:17px;margin:0 0 8px}html[data-mw-surface="launcher"] #mwWelcome ol{display:grid;grid-template-columns:repeat(3,1fr);gap:12px;margin:8px 0 12px}
@media(min-width:1650px){.mwWorkbenchRow{grid-template-columns:minmax(0,1fr) 330px!important}#viewport{height:360px!important}}
@media(max-width:1100px){.top{padding-left:20px!important;padding-right:20px!important;gap:6px 12px}.top h1{font-size:25px!important}main{padding-left:20px!important;padding-right:20px!important}.mwWorkbenchRow{grid-template-columns:minmax(0,1fr) 265px!important;gap:14px}.card,.panel{padding:13px}.mwHeadActions button{padding:8px 10px;font-size:11px}.mwObservationMeta{display:none}.mwObservation{gap:13px}}
@media(max-width:900px){.mwWorkbenchRow{grid-template-columns:minmax(0,1fr)!important}.mwWorkbenchRow>.mwPointPanel{grid-column:1;grid-row:2}.mwWorkbenchRow>.mwPlotPanel{grid-column:1;grid-row:1}.mwPointPanel .log{max-height:200px}.top h1{font-size:24px!important}.mwHeadActions{flex-direction:column;align-self:center}.typelist{grid-template-columns:minmax(0,1fr)!important}#viewport{height:390px!important}.hero{grid-template-columns:1fr;gap:10px}.hero .wordmark{font-size:48px}.hero .sub{max-width:640px}.hero .meta-row{gap:10px 22px}}
@media(max-width:600px){html[data-mw-surface]:not([data-mw-surface="report"]) body{padding-left:57px}#mwRail{width:57px}.mwBrand{margin:16px 0 19px;font-size:6px;letter-spacing:.02em}.mwMark{width:33px;height:35px;font-size:22px}#mwRail nav{padding:0 3px;gap:4px}.mwInstrument{min-height:59px;font-size:9px;gap:4px}.mwInstrument small{font-size:6px}.mwIcon{width:17px;height:17px}.mwInstrument[aria-current="page"]:after{width:4px;height:4px;right:5px;top:6px}.mwRailBottom{padding:12px 3px 9px;gap:4px}.mwRailControl,#mwRail button{min-height:44px;font-size:9px}.mwRailVersion{font-size:5px;letter-spacing:0}
.top{padding:15px 12px!important;display:flex!important;gap:8px}.top h1{font-size:22px!important;flex-basis:100%;padding-left:10px;border-left-width:3px}.mw-kicker{font-size:9px;margin-bottom:6px}.top .sub{margin-left:13px;font-size:10px}.mwHeadActions{flex-direction:row;margin-top:4px;gap:6px;flex-wrap:wrap}.mwHeadActions button{min-height:39px;font-size:10px}
main{padding-left:12px!important;padding-right:12px!important;padding-bottom:95px!important}.card,.panel{padding:12px 10px;margin-bottom:14px}.bar{gap:7px}.bar .grow{flex-basis:100%}.bar>select{max-width:100%;width:auto;flex:1 1 140px}.bar>button{flex:0 1 auto;min-height:39px;font-size:11px}.mwPanelHeading{font-size:12px!important}.mwPanelHeading small{font-size:7px}.mwObservation{padding:11px!important;gap:10px}.mwObservationField,.mwObservationFilter{flex-basis:100%;max-width:none}.mwDataTable{font-size:10px;min-width:600px}.mwTableFooter{font-size:9px;align-items:flex-start}.mwTablePager{gap:3px}.mwPlotPanel .bar{gap:6px}#viewport{height:300px!important}
#mwFooter{left:57px;padding:5px 12px;min-height:39px;font-size:8px!important;gap:6px}#mwBuildLabel{max-width:calc(100% - 100px)}#mwFooter #mwSamples{font-size:9px;min-height:30px;padding:3px 8px}#mwDialog{padding:15px}.cmdbar .brand{padding-left:12px!important;font-size:11px!important}.cmdbar .ver{max-width:45%;font-size:8px!important;padding-right:10px!important}.hero{display:block}.hero .wordmark{font-size:clamp(28px,8vw,45px)!important}.hero .sub{font-size:12px}.hero .meta-row{gap:8px 13px}.pair-head{gap:10px}.pair-head .cn{font-size:20px}.grid{grid-template-columns:minmax(0,1fr)!important}a.card{min-height:100px;padding:17px 13px}a.card .nm{font-size:15px}.typelist{grid-template-columns:1fr!important}#scanFilter,#scan,#rxl,#frame,#dlSol,#dlScans{min-width:0!important;max-width:100%}#mwSampleDestination{width:100%}html[data-mw-surface="launcher"] #mwWelcome ol{display:block}}
@media(pointer:coarse){button,.iconbtn,summary{min-height:44px}.mwInstrument{min-height:60px}}
@media(prefers-reduced-motion:reduce){html:root{--ef-motion:0ms}*,*::before,*::after{animation:none!important;transition:none!important;scroll-behavior:auto!important}}
'''
SHELL_CSS=':root{--mw-contours:url(data:image/svg+xml;base64,'+_contours()+')}'+SHELL_CSS
READING_CSS=r'''
html[data-mw-surface="report"] body{font:15px/1.75 var(--ef-font);margin:0 auto;padding:32px 20px;max-width:1240px;background:var(--bg);color:var(--ink);overflow-wrap:anywhere}html[data-mw-surface="report"] body::after{display:none}html[data-mw-surface="report"] h1{font-size:28px;line-height:1.4;border-bottom:2px solid var(--ink);padding-bottom:14px}html[data-mw-surface="report"] h1::before{content:'// DATA REPORT';display:block;font-size:11px;font-stretch:125%;letter-spacing:.14em;color:var(--muted);font-weight:600}html[data-mw-surface="report"] h2{border-top:1px solid var(--line);padding-top:18px}html[data-mw-surface="report"] table{display:block;overflow:auto;max-width:100%;border-collapse:collapse;background:var(--field);font-size:12px}html[data-mw-surface="report"] td,html[data-mw-surface="report"] th{border:0;border-bottom:1px solid var(--line);padding:7px 9px}html[data-mw-surface="report"] th{background:#252525;color:white;font-weight:600}html[data-mw-surface="report"] img{max-width:100%;height:auto}html[data-mw-surface="report"] .card{border:0;border-top:1px solid var(--line);box-shadow:none}html[data-mw-surface="report"] a{color:var(--ink);text-decoration:underline;text-decoration-color:var(--ink);text-underline-offset:3px}html[data-mw-surface="report"] a:hover{background:var(--accent);color:var(--accent-ink)}
'''
SKIN_CSS=FONT_CSS+TOKEN_CSS+SHELL_CSS+READING_CSS

def style_tag():
    return '<style data-mw-theme="endfield-yellow-paper">'+SKIN_CSS+'</style>'

def decorate(body,instrument='launcher',surface=None):
    is_bytes=isinstance(body,bytes)
    text=body.decode('utf-8') if is_bytes else body
    closing=r'</body\s*>' if re.search(r'</body\s*>',text,re.I) else r'</html\s*>'
    if not re.search(closing,text,re.I):return body
    text=re.sub(r'<style\b[^>]*data-mw-theme=[^>]*>.*?</style>','',text,flags=re.S|re.I)
    mode=surface or ('launcher' if instrument=='launcher' else 'workbench')
    attributes={'data-mw-surface':mode,'data-mw-instrument':instrument,**{'data-ef-'+k:v for k,v in VARIANT.items()}}
    def tag(match):
        value=match[0]
        for key in attributes:value=re.sub(r'\s+'+re.escape(key)+r'=[\"\'][^\"\']*[\"\']','',value,flags=re.I)
        return value[:-1]+''.join(' '+k+'="'+html.escape(v,quote=True)+'"' for k,v in attributes.items())+'>'
    text=re.sub(r'<html\b[^>]*>',tag,text,count=1,flags=re.I)
    text=re.sub(closing,lambda m:style_tag()+m[0],text,count=1,flags=re.I)
    if mode=='workbench' and 'class="mw-kicker"' not in text:
        kicker='<span class="mw-kicker">// '+LABELS.get(instrument,'DATA WORKBENCH')+'</span>'
        text=re.sub(r'(<h1\b[^>]*>)',lambda m:m[0]+kicker,text,count=1,flags=re.I)
    return text.encode('utf-8') if is_bytes else text

def report_style(text):
    return decorate(text,'report','report')

def validation_label():
    import appenv
    try:
        data=json.loads((appenv.resources()/'_last_run.json').read_text(encoding='utf-8'))
        if data.get('failures')==0 and data.get('scope')=='developer' and data.get('results'):
            return str(int(data['assertions'])),str(data.get('finished_utc',''))[:10]
    except (OSError,ValueError,KeyError,TypeError):pass
    return '待运行','无发行记录'

