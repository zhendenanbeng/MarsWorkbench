"""PIXL原分辨率出图、CSV及HTML报告；数值读取来自核心后端。"""
import csv,hashlib,html,io,math
from pathlib import Path
import numpy as np
from PIL import Image,ImageDraw
from .pixl_browser import *

def color(value,lo,hi):
    t=max(0,min(1,(value-lo)/(hi-lo))) if hi>lo else .5
    stops=[(37,75,150),(31,164,186),(234,218,68),(188,49,42)]
    index=min(2,int(t*3));f=t*3-index
    return tuple(round(a*(1-f)+b*f) for a,b in zip(stops[index],stops[index+1]))

def render(data,element,frame=0,axis='label',origin=0):
    if element!='points' and element not in data['elements']:raise ValueError('成分指标无效')
    image,points,warnings=canvas_data(data,frame,axis,origin)
    if image is None:raise ValueError('缺MCC图，不能导出叠图；成分CSV和报告仍可导出')
    h,w=image.shape;out=Image.new('RGB',(w,h+96),'white');out.paste(Image.fromarray(image).convert('RGB'),(0,0));draw=ImageDraw.Draw(out)
    values=[p['values'].get(element,{}).get('value') for p in points];valid=[v for v in values if v is not None];lo=min(valid) if valid else 0;hi=max(valid) if valid else 0
    for p,v in zip(points,values):
        fill=(237,164,47) if element=='points' else color(v,lo,hi) if v is not None else (145,145,145)
        x,y=p['x'],p['y'];draw.ellipse((x-2.5,y-2.5,x+2.5,y+2.5),fill=fill)
    draw.text((8,h+6),'PIXL '+element+' | raw wt% range: '+('%.6g .. %.6g'%(lo,hi) if valid else 'N/A'),fill='black')
    draw.text((8,h+27),'Missing=gray; measured points only; no interpolation; no scale bar',fill='black')
    draw.text((8,h+48),'Axes: '+('XML: sample=j, line=i' if axis=='label' else 'MANUAL: sample=i, line=j')+'; pixel offset='+str(origin),fill='black')
    draw.text((8,h+69),'Official RQC values; no mineral identification / cross-instrument alignment',fill='black')
    buffer=io.BytesIO();out.save(buffer,format='PNG');return buffer.getvalue()

def export(data,out,element,frame=0,axis='label',origin=0):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    columns=['PMC']+[e+s for e in data['elements'] for s in ('_wt%','_wt%_err')]
    with (out/'composition.csv').open('w',encoding='utf-8-sig',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(columns)
        for pmc,values in sorted(data['chem'].items()):writer.writerow([pmc]+[values[e][k] if values[e][k] is not None else '' for e in data['elements'] for k in ('value','error')])
    png=None
    if element not in ('points',*data['elements']):raise ValueError('成分指标无效')
    image,points,warnings=canvas_data(data,frame,axis,origin)
    if image is not None:png='map_'+element+'.png';(out/png).write_bytes(render(data,element,frame,axis,origin))
    point_by_pmc={p['pmc']:p for p in points}
    with (out/'coordinates.csv').open('w',encoding='utf-8-sig',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(['PMC','MCC_frame','raw_i','raw_j','display_x','display_y','visible','axis','pixel_origin_offset'])
        key=data['frames'][frame]['key'] if data['frames'] else None
        for pmc,row in sorted(data['raw']['pts'].items()):
            pair=row.get(key,(None,None));p=point_by_pmc.get(pmc,{})
            writer.writerow([pmc,key,*pair,p.get('x',''),p.get('y',''),bool(p),axis,origin])
    esc=lambda v:html.escape(str(v));obs=data['observation'];rows=[]
    stats=composition_stats(data)
    stat_cells=''.join('<tr><td>'+esc(e)+'</td><td>'+('%.4f'%stats[e]['mean'])+'</td><td>'+('%.4f'%stats[e]['median'])+'</td><td>'+('%.4f'%stats[e]['std'])+'</td><td>'+str(stats[e]['n'])+'</td></tr>' for e in data['elements'])
    for pmc,values in sorted(data['chem'].items()):rows.append('<tr><td>'+str(pmc)+'</td>'+''.join('<td>'+('—' if v['value'] is None else esc(v['value']))+' ± '+('—' if v['error'] is None else esc(v['error']))+'</td>' for v in values.values())+'</tr>')
    report='<!doctype html><html lang="zh"><meta charset="utf-8"><title>PIXL 成分报告</title><style>body{font:14px/1.6 sans-serif;margin:30px;color:#173147}td,th{border:1px solid #ccd;padding:5px}table{border-collapse:collapse}img{max-width:100%}</style><body><h1>PIXL 官方成分浏览报告</h1><p>'+esc(NOTE)+'</p><p>'+esc(obs['file'])+'</p><p>坐标：'+esc(data['rxl_path'].name if data['rxl_path'] else '缺失')+'；按定位+RTT配对，PMC连接；误差：'+esc(data['error_label'])+'</p><p>坐标约定：'+esc(axis)+'；原点偏移：'+str(origin)+'；只按归档坐标浏览，未独立核验像素定位精度。</p><p>'+esc('；'.join(warnings))+'</p><h2>成分统计</h2><p>'+esc(STATS_NOTE)+'</p><table><tr><th>元素</th><th>均值 wt%</th><th>中位数 wt%</th><th>点间标准差</th><th>点数</th></tr>'+stat_cells+'</table><p><a href="composition.csv">成分CSV</a> · <a href="coordinates.csv">坐标CSV</a> · <a href="sources.json">来源与SHA</a></p>'
    if png:report+='<img src="'+png+'" alt="实测成分点图">'
    report+='<table><tr><th>PMC</th>'+''.join('<th>'+esc(e)+' wt% ± err</th>' for e in data['elements'])+'</tr>'+''.join(rows)+'</table></body></html>'
    from product_theme import report_style
    (out/'index.html').write_text(report_style(report),encoding='utf-8')
    paths=[obs['path']]+([data['rxl_path']] if data['rxl_path'] else [])
    if data['frames'] and data['frames'][frame]['path']:paths.append(data['frames'][frame]['path'])
    paths+=[p.with_suffix('.xml') for p in list(paths) if p.with_suffix('.xml').exists()]
    import json
    (out/'sources.json').write_text(json.dumps({'files':[{'file':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in paths],
        'axis':axis,'pixel_origin_offset':origin,'error_label':data['error_label'],'note':NOTE},ensure_ascii=False,indent=2),encoding='utf-8')
    return {'csv':'composition.csv','map':png,'report':'index.html','warnings':warnings}
