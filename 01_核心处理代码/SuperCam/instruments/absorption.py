"""Continuum-relative trough references for calibrated CP2 I/F only, not mineral identification."""
import numpy as np
from .lines import IR_REFERENCES,robust_sigma
from sherloc.peaks import half_maximum

SCOPE='反射吸收特征，快速浏览，非矿物鉴定。深度=1−R/局部连续谱；局部差分分数不是置信概率。重叠参考不构成独立物相证据；大气校正、温度及背景残余仍可造假带。2.5µm附近尤其需复核热敏感性。'

def measure(x,y,definition):
    key,label,limits,shoulders,source=definition
    row={'key':key,'label':label,'range_um':list(limits),'source':source,'depth':None,'snr':None,'center_um':None,'detected':False,'available':False}
    step=float(np.median(np.diff(x)))
    if shoulders[0]<x[0]-step/2 or shoulders[1]>x[-1]+step/2:return dict(row,reason='窗口外或缺少两侧连续谱')
    li=np.sort(np.argsort(abs(x-shoulders[0]))[:3]);ri=np.sort(np.argsort(abs(x-shoulders[1]))[:3])
    left,right=float(np.median(x[li])),float(np.median(x[ri]))
    use=np.flatnonzero((x>=x[li[0]])&(x<=x[ri[-1]]))
    if len(use)<9 or not np.isfinite(y[use]).all():return dict(row,reason='缺少足够连续有限通道')
    xx,yy=x[use],y[use];dif=np.diff(xx)
    if np.max(dif)>5*np.median(dif):return dict(row,reason='连续谱窗口跨缺失通道空档')
    a,b=float(np.median(y[li])),float(np.median(y[ri]))
    continuum=a+(b-a)*(xx-left)/(right-left)
    if np.any(continuum<=0):return dict(row,reason='连续谱非正，不能计算反射吸收深度')
    depth=1-yy/continuum
    noise=np.r_[np.diff(depth[xx<limits[0]]),np.diff(depth[xx>limits[1]])]
    sigma=max(robust_sigma(noise)/np.sqrt(2),64*np.finfo(float).eps)
    candidates=[i for i in range(1,len(xx)-1) if limits[0]<=xx[i]<=limits[1] and depth[i]>depth[i-1] and depth[i]>=depth[i+1]]
    if not candidates:return dict(row,available=True,reason='窗口内无局部吸收极小')
    i=max(candidates,key=lambda i:depth[i]);width=half_maximum(xx,depth,i)
    if yy[i]<=0:return dict(row,reason='中心反射值非正，不能解释为吸收深度')
    score=float(depth[i]/sigma);det=bool(yy[i]>0 and depth[i]>=.01 and score>=5 and width and width['fwhm']>=1.5*np.median(dif))
    return dict(row,available=True,depth=float(depth[i]),snr=score,center_um=float(xx[i]),fwhm_um=width['fwhm'] if width else None,detected=det,
                reason='达到浏览门槛（非物相确认）' if det else '未达深度/局部分数/完整峰翼浏览门槛')

def payload(x,y=None,column='',factor=1):
    x=np.asarray(x,float)
    if x.ndim!=1 or len(x)<3 or not np.isfinite(x).all() or not np.all(np.diff(x)>0):raise ValueError('吸收参考需要有限递增波长轴')
    if y is not None and np.asarray(y).shape!=x.shape:raise ValueError('吸收谱与坐标不匹配')
    rows=[]
    for definition in IR_REFERENCES:
        if y is not None:row=measure(x,np.asarray(y,float),definition)
        else:
            k,label,r,shoulders,url=definition
            row={'key':k,'label':label,'range_um':list(r),'source':url,'detected':False,'depth':None,'snr':None,'available':False,'reason':'仅参考；请开启CP2反射吸收参考' if r[0]>=x[0] else '窗口外'}
        rows.append(row)
    return {'enabled':y is not None,'column':column,'axis_factor':factor,'references':rows,'scope':SCOPE,
            'method':'3通道中位肩点局部连续谱、窗口内局部倒峰；不是论文固定中心参数的逐项复现',
            'thresholds':{'min_depth':.01,'local_score':5.,'min_width_in_channels':1.5}}

