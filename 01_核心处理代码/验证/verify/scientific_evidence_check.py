# -*- coding: utf-8 -*-
"""文献坐标制、独立谱线、重叠与失败语义的护栏；不把合成谱称作矿物盲测。"""

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout
import json
from pathlib import Path
import sys
from unittest.mock import patch

import numpy as np

ROOT=_layout.project_layout.APP
sys.path.insert(0,str(ROOT))
from instruments import lines
from sherloc import qc,minerals
from sherloc.peaks import half_maximum

TOTAL=FAIL=0
def check(ok,name):
    global TOTAL,FAIL
    TOTAL+=1
    if not ok:FAIL+=1;print(' FAIL '+name)


def wavelength():
    # SDSS公开的空气/真空H-beta和H-alpha对，独立于待测Peck-Reeder函数。
    actual=lines.air_to_vacuum_nm(np.array([486.1363,656.2801]))
    check(np.allclose(actual,[486.2721,656.4614],rtol=0,atol=.00005),'标准空气/真空波长对')
    check(actual[1]>656.45,'Halpha不能仍按空气656.28标注')
    h=next(w for e,w,label in lines.LIBS_LINES if label=='Hα 656')
    check(abs(h-656.4613)<.0001,'LIBS实际使用真空参考表')
    check(all(b[1]>a[1] for a,b in zip(lines.LIBS_LINES_AIR,lines.LIBS_LINES) if a[0] not in ('C2','CN')),'全部原子线转换且不混用分子带')
    try:lines.air_to_vacuum_nm(180);check(False,'越界波长应拒绝')
    except ValueError:check(True,'越界波长应拒绝')


def widths():
    x=np.array([0.,10.,20.]);y=np.array([0.,10.,0.])
    measured=half_maximum(x,y,1)
    check(measured['fwhm']==10.,'单像素三角尖刺的宽度是10，不能量成20')
    x=np.array([0.,1.,3.,6.,8.]);y=np.array([0.,2.,8.,2.,0.])
    check(np.isclose(half_maximum(x,y,2)['fwhm'],10/3),'不等间距半高交点正确插值')
    xx=np.linspace(-5,5,2001);yy=np.exp(-xx**2/2)
    check(abs(half_maximum(xx,yy,1000)['fwhm']-2*np.sqrt(2*np.log(2)))<1e-5,'高斯解析半高宽')
    check(half_maximum([0,1,2],[8,10,0],1) is None,'截断峰不能报告完整FWHM')
    c=half_maximum([0,1,2],[8,10,0],1,allow_censored=True)
    check(c['censored'] and c['fwhm']==1.5,'截断峰只返回显式下界')
    check(half_maximum([0,1,2],[np.nan,10,0],1) is None,'NaN不能冒充半高交点')


def emission():
    rng=np.random.default_rng(43120)
    x=np.arange(495.,525.,.02)
    def peak(center,amp=100,width=.08):return amp*np.exp(-.5*((x-center)/width)**2)
    reference=[('X',500.,'X1'),('X',520.,'X2')]
    y=1000+2*x+rng.normal(0,.3,len(x))+peak(500)
    marks,elems=lines.match_lines(x,y,reference,.35)
    check(marks[0]['detected'] and not elems,'单条强谱线仅作候选')
    marks,elems=lines.match_lines(x,y+peak(520),reference,.35)
    check(elems==['X'] and all(m['element_supported'] for m in marks),'两条独立线支持元素候选')
    check(all(m['fwhm']>0 and m['peak_x'] for m in marks),'候选须保存实测位置和完整峰宽')
    marks,elems=lines.match_lines(x,y,[('X',499.9,'X1'),('X',500.1,'X2')],.35)
    check(not elems and marks[0]['peak_index']==marks[1]['peak_index'],'同一峰不能充当两条独立线')
    envelope=1000+peak(500,width=.25)+peak(500.6,width=.25)
    em,ee=lines.match_lines(x,envelope,[('X',500.,'X1'),('X',500.6,'X2')],.35)
    check(not ee and all(m['independent_lines']==1 for m in em),
          '同一元素未分辨包络内的两个极大值也不能当两条独立线')
    marks,elems=lines.match_lines(x,y,[('X',499.9,'X1'),('Z',500.1,'Z1')],.35,min_lines=1)
    check(not elems and all(m['ambiguous_with'] for m in marks),'共享峰不能同时确认两个元素')
    check(not lines.match_lines(x,1000+2*x,reference,.35)[1],'倾斜连续背景不构成发射线')
    broad=1000+peak(500,width=4)
    marks,elems=lines.match_lines(x,broad,reference,.35,min_lines=1)
    check(not elems and not any(m['detected'] for m in marks),'宽背景包不构成窄发射线')
    gap=(x<499)|(x>502)
    marks,_=lines.match_lines(x[gap],y[gap],reference,.35)
    check(not marks[0]['available'] and not marks[0]['detected'],'不跨仪器空档取基线或外推峰')
    apx=np.arange(.5,20,.02)
    apy=2+50*np.exp(-.5*((apx-4.51)/.05)**2)+100*np.exp(-.5*((apx-4.931)/.06)**2)
    marks,elems=lines.match_lines(apx,apy,lines.APXS_LINES,.09,min_lines=1)
    check('Ti' in elems and 'V' not in elems,'Ti Kbeta/V Kalpha未解混时不得报V')
    check(any(m['ambiguous_with'] for m in marks if m['elem']=='V'),'保存V/Ti重叠原因')
    check(not any(e=='O' for e,*_ in lines.APXS_LINES),'MSL APXS不在0.525 keV直接识别氧')
    ar=next(row for row in lines.APXS_LINES if row[0]=='Ar')
    check(ar[3] is True,'大气Ar不是目标组成的独立证据')
    source_y=2+100*np.exp(-.5*((apx-14.22)/.12)**2)
    sm,se=lines.match_lines(apx,source_y,lines.APXS_LINES,.09,min_lines=1)
    check('Sr' not in se and any(m['ambiguous_with'] for m in sm if m['elem']=='Sr'),
          '源区散射峰与Sr重叠时不得当目标Sr')
    cx=np.arange(322.,330.,.01)
    cu=[r for r in lines.LIBS_LINES if r[0]=='Cu']
    cy=sum(100*np.exp(-.5*((cx-r[1])/.04)**2) for r in cu)
    cm,ce=lines.match_lines(cx,cy,cu,.35,wavelength_frame='vacuum_nm')
    check('Cu' not in ce and all(m['ambiguous_with'] for m in cm),'文献已知Cu两条混合带不能仅靠双线确认')
    false=0
    nx=np.arange(240.,900.,.05)
    for _ in range(100):
        _,elements=lines.match_lines(nx,rng.normal(0,1,len(nx)),lines.LIBS_LINES,.35,wavelength_frame='vacuum_nm')
        false+=bool(elements)
    check(false==0,'100条冻结噪声谱经完整LIBS表无多线元素误报（合成用例）')
    check('NaN' not in json.dumps(marks,allow_nan=False),'候选元数据JSON不能包含NaN')
    for xx,yy in [(x,np.r_[y[:-1],np.nan]),(x[::-1],y),(x,y[:-1])]:
        try:lines.match_lines(xx,yy,reference,.35);check(False,'坏谱不可输出无元素的成功结果')
        except ValueError:check(True,'坏谱必须报错')


def fail_closed():
    rng=np.random.default_rng(2)
    shift=np.linspace(0,4500,2148);wl=1/(1/qc.LASER_NM-shift*1e-7)
    data=rng.integers(0,20,(2,len(wl)))
    for wx,yy,kw in [(wl,data[:,:-1],{}),(wl,np.zeros((0,len(wl))),{}),
                     (np.r_[wl[:-1],np.nan],data,{}),(np.ones(len(wl))*260,data,{}),
                     (wl,data,{'noise_window':(5000,6000)}),(wl,data,{'snr_min':np.nan}),
                     (wl,np.full(data.shape,np.nan),{}),(wl,data,{'r2':np.ones((1,len(wl)))}),
                     (wl,data,{'max_points':-1})]:
        try:qc.analyse_spectra(wx,yy,**kw);check(False,'坏输入不能显示没有拉曼信号')
        except ValueError:check(True,'坏输入必须报错')
    with patch.object(minerals,'analyse_scan',side_effect=ValueError('injected science failure')):
        try:qc.analyse_spectra(wl,data);check(False,'内核异常不可吞成无矿物')
        except RuntimeError as e:check('injected science failure' in str(e),'保留真实失败原因')
    a=qc.analyse_spectra(wl,data);b=qc.analyse_spectra(wl,data.astype(float))
    check(a['sub'].dtype.kind=='f' and np.array_equal(a['sub'],b['sub']),'整型输入也必须以浮点残差分析')
    g={'n_points':2,'band':{'center':1084},'classes':[{'key':'carbonate','class':'碳酸盐'}]}
    sentence=qc.headline({'mineral_groups':[g,g]})
    check('1 类' in sentence,'同一类别的两个组不能说成两类矿物')
    import view
    class Scan:
        key='failure-case'
        def read_regions(self):return wl,data,None,None
    state={'mineral_cache':{}}
    with patch.object(qc,'analyse_spectra',side_effect=RuntimeError('injected failure')):
        result=view.mineral_summary(view.scan_minerals(state,Scan()))
    check(not result['ok'] and 'injected failure' in result['error'],'工作台显示分析失败原因，不等同未检出')
    check(not state['mineral_cache'],'失败不能缓存成无矿物或阻止重试')


def main():
    wavelength();widths();emission();fail_closed()
    print('共 %d 项断言，%d 项失败'%(TOTAL,FAIL))
    return int(FAIL>0)


if __name__=='__main__':raise SystemExit(main())
