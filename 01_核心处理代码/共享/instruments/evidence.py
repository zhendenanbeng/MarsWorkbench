"""Empirical shot recurrence, never a detection probability or concentration estimate."""
from collections import OrderedDict
from copy import deepcopy
import hashlib
import threading
import numpy as np
from . import lines

_CACHE=OrderedDict();_LOCK=threading.RLock()
SCOPE='经验检出率=通过同一谱线规则的保留发数/全部保留发数；不是置信概率、独立重复实验或含量。SNR中位数包含有局部峰但未达门槛的发，无局部峰的发不参与中位数。†为截断宽度下界。'

def shot_statistics(x,shots,tol=.35):
    x=np.asarray(x,float);shots=np.asarray(shots,float)
    if shots.ndim!=2 or shots.shape[0]!=len(x) or shots.shape[1]<1 or not np.isfinite(shots).all():
        raise ValueError('逐发证据需要与波长轴匹配的有限二维发数矩阵')
    digest=hashlib.sha256(x.tobytes()+shots.tobytes()+repr((tol,lines.LIBS_LINES)).encode()).hexdigest()
    with _LOCK:
        if digest in _CACHE:
            _CACHE.move_to_end(digest);return deepcopy(_CACHE[digest])
    samples=[lines.match_lines(x,shots[:,i],lines.LIBS_LINES,tol=tol,min_lines=2,wavelength_frame='vacuum_nm') for i in range(shots.shape[1])]
    stats=[]
    for j in range(len(lines.LIBS_LINES)):
        records=[marks[j] for marks,_ in samples];elem=records[0]['elem']
        snr=[m['snr_raw'] for m in records if m.get('snr_raw') is not None]
        widths=[m['fwhm'] for m in records if m['fwhm'] is not None]
        covered=sum(m['available'] for m in records);detected=sum(m['detected'] for m in records)
        supported=sum(elem in elements for _,elements in samples)
        stats.append({'shot_detected':detected,'shot_total':len(samples),'shot_available':covered,
                      'shot_detection_fraction':detected/len(samples) if covered else None,
                      'shot_snr_median':float(np.median(snr)) if snr else None,'shot_snr_count':len(snr),
                      'shot_fwhm_median':float(np.median(widths)) if widths else None,
                      'shot_censored':sum(m.get('fwhm_censored',False) for m in records),
                      'element_shot_detected':supported,'evidence_scope':SCOPE})
    with _LOCK:
        _CACHE[digest]=stats
        while len(_CACHE)>32:_CACHE.popitem(last=False)
    return deepcopy(stats)

def enrich(x,shots,marks):
    stats=shot_statistics(x,shots)
    if len(stats)!=len(marks):raise ValueError('逐发证据与当前参考线表不一致')
    return [dict(mark,**stat) for mark,stat in zip(marks,stats)]

def element_statistics(marks,elements):
    result=[]
    for elem in elements:
        own=[m for m in marks if m['elem']==elem and 'shot_total' in m]
        if own:result.append({'element':elem,'detected':own[0]['element_shot_detected'],'total':own[0]['shot_total'],
                              'fraction':own[0]['element_shot_detected']/own[0]['shot_total']})
    return sorted(result,key=lambda r:(-r['fraction'],r['element']))

def summary(marks,elements):
    return '；'.join('%s %d/%d发'%(r['element'],r['detected'],r['total']) for r in element_statistics(marks,elements))

def line_text(m):
    if 'shot_total' not in m:return ''
    width=m.get('fwhm_observed');snr=m.get('shot_snr_median')
    return '%s %d/%d发 · SNR中位%s（%d发） · FWHM %s%s nm'%(m['label'],m['shot_detected'],m['shot_total'],
        '—' if snr is None else '%.2f'%snr,m['shot_snr_count'],'—' if width is None else '%.4g'%width,'†' if m.get('fwhm_censored') else '')+(' · %d发宽度截断†'%m['shot_censored'] if m.get('shot_censored') else '')

def append_png(path,marks,elements):
    from PIL import Image,ImageDraw
    from sherloc.render import font
    rows=[line_text(m) for m in marks if m['detected'] and 'shot_total' in m]
    if not rows:return
    with Image.open(path) as original:
        image=Image.new('RGB',(original.width,original.height+100+24*len(rows)),'white')
        image.paste(original,(0,0));y=original.height+12
    draw=ImageDraw.Draw(image)
    for text in ['逐发元素支持：'+summary(marks,elements),'经验统计，非置信概率或含量；†为截断宽度下界；SNR中位数按有局部峰的保留发计算。',*rows]:
        draw.text((28,y),text,font=font(16),fill='#334155');y+=24
    image.save(path)
