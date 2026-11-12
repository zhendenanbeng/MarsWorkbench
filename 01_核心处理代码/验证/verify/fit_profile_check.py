# -*- coding: utf-8 -*-
"""拟合参数/RMS一致性、物理约束、边界和复核产物的行为护栏。"""

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout
import csv
import json
from pathlib import Path
import sys
import tempfile
from unittest.mock import patch

import numpy as np

ROOT=_layout.project_layout.APP;sys.path.insert(0,str(ROOT))
from sherloc import minerals as m
from sherloc.fitting import background_context,solve_nonnegative_profiles
import review_peak_fits

TOTAL=FAIL=0
def check(ok,name):
    global TOTAL,FAIL
    TOTAL+=1
    if not ok:FAIL+=1;print(' FAIL '+name)


def reconstruct(x,r):
    return m._pseudo_voigt(x,r['amp'],r['center'],r['fwhm'],r['eta'])+r['slope']*x+r['intercept']


def numeric():
    x=np.linspace(900,1100,41);g,l=m._voigt_basis(x,1000,50);p=np.column_stack([g,l])
    y=30*g-10*l+3
    fit=m.fit_voigt(x,y,1000,50,within=(970,1030))
    rms=float(np.sqrt(np.mean((y-reconstruct(x,fit))**2)))
    check(abs(rms-fit['rms'])<1e-9,'有符号候选参数必须重建其所报告RMS')
    check(fit['eta']<0 and not fit['profile_valid'],'负峰形系数不可截成合法混合比例')
    check(fit['model']=='signed_gaussian_lorentz','非法物理模型必须明确命名')
    physical=m.fit_voigt(x,y,1000,50,within=(970,1030),nonnegative=True)
    check(physical['profile_valid'] and 0<=physical['eta']<=1,'物理拟合只允许非负峰形系数')
    check(abs(np.sqrt(np.mean((y-reconstruct(x,physical))**2))-physical['rms'])<1e-9,'物理拟合RMS须对应同一条曲线')
    check(physical['rms']>1e-4,'不可用未约束的零残差冒充物理解')
    coefficient,error,model=solve_nonnegative_profiles(x,y,p)
    check(np.all(coefficient[:2]>=0),'活动面求解系数非负')
    check(abs(error-np.sum((y-model)**2))<1e-10,'求解器平方误差可重建')
    check(np.allclose(model,p@coefficient[:2]+coefficient[2]*x+coefficient[3],atol=1e-10),'背景系数须还原原始坐标')
    target=20*g+8*l-.02*x+3
    c,e,ym=solve_nonnegative_profiles(x,target,p)
    check(np.allclose(c,[20,8,-.02,3],atol=1e-9),'已知双基函数及自由背景恢复')
    check(e<1e-16,'精确合成模型残差接近零')
    yc=m._pseudo_voigt(x,40,1003,40,.3)+3
    bounded=m.fit_voigt(x,yc,1000,40,within=(990,1000))
    check(990<=bounded['center']<=1000,'细搜索不能越过显式峰心边界')
    check(bounded['at_bound'],'受约束边界解必须标识')
    same=m.fit_voigt(x,yc,1000,40,within=(990,1000),nonnegative=True)
    check(990<=same['center']<=1000,'物理拟合也遵守同一峰心边界')
    ctx=background_context(x,target);original=target.copy();target[0]+=1
    try:solve_nonnegative_profiles(x,target,p,context=ctx);check(False,'失效上下文应拒绝')
    except ValueError:check(True,'失效上下文应拒绝')
    check(np.array_equal(ctx['y'],original),'上下文不能随输入数组修改而悄悄改变')
    for mode in [False,True]:
        y=m._pseudo_voigt(x,40,1004,42,.15)+.03*(x-1000)+10
        r=m.fit_voigt(x,y,1000,40,within=(970,1030),nonnegative=mode)
        check(abs(r['center']-1004)<=1 and abs(r['fwhm']-42)<=2,'已知正峰的中心与峰宽恢复')
        check(np.allclose(y,reconstruct(x,r),atol=1e-6),'合法正峰不受约束修正破坏')
    x2=np.linspace(900,1180,57)
    y2=m._pseudo_voigt(x2,40,970,40,.2)+m._pseudo_voigt(x2,25,1090,45,.4)+3
    two=m.fit_two_voigt(x2,y2,1010,60,.1,within=(940,1120),nonnegative=True)
    check(two['ok'],'两个分离强峰应可拆分')
    if two['ok']:
        ym=sum(m._pseudo_voigt(x2,c['amp'],c['center'],c['fwhm'],c['eta']) for c in two['comps'])+two['slope']*x2+two['intercept']
        check(abs(np.sqrt(np.mean((y2-ym)**2))-two['rms'])<1e-9,'双峰参数须重建双峰RMS')
        check(two['profile_valid'] and all(c['profile_valid'] for c in two['comps']),'双峰物理状态须逐组分一致')
    ox=np.linspace(3000,4000,101);oy=10*np.exp(-.5*((ox-3400)/50)**2)
    bad={'amp':10.,'center':3400.,'fwhm':100.,'at_bound':False,'rms':.1,'profile_valid':False}
    with patch.object(m,'fit_two_voigt',return_value={'ok':False}),patch.object(m,'fit_voigt',return_value=bad):
        oh=m.hydration(ox,oy,1.)
    check(oh['ok'] and not oh['profile_valid'] and oh['requires_review'],'单水带非物理状态必须向上传递')
    pair={'ok':True,'comps':[dict(bad)],'rms':.1,'profile_valid':False}
    with patch.object(m,'fit_two_voigt',return_value=pair):oh=m.hydration(ox,oy,1.)
    check(oh['ok'] and not oh['profile_valid'] and oh['requires_review'],'双水带非物理状态必须向上传递')


def products():
    from sherloc import qc
    rng=np.random.default_rng(17);shift=np.linspace(0,4500,2148)
    wl=1/(1/qc.LASER_NM-shift*1e-7)
    r1=np.array([30+100*np.exp(-.5*((shift-1084)/25)**2)+rng.normal(0,1,len(shift)) for _ in range(6)])
    class Scan:
        sol='9999';sequence_id='synthetic';sclk='123';products=[]
        def read_regions(self):return wl,r1,None,None
    original=r1.copy()
    with tempfile.TemporaryDirectory(prefix='sherlock-fit-profile-') as folder:
        result=review_peak_fits.review(Scan(),Path(folder))
        check(bool(result['groups']),'合成强峰须生成复核组')
        check(np.array_equal(r1,original),'复核不得更改输入光谱')
        for g in result['groups']:
            if g['physical']:
                check(abs(g['physical']['rms']-g['physical']['reconstructed_rms'])<1e-9,'复核参数及重建曲线一致')
        check(all((Path(folder)/g['csv']).is_file() for g in result['groups']),'每个组都有完整曲线CSV')
        json.dumps(result,allow_nan=False)
        check(True,'复核JSON不得包含NaN/Infinity')


def main():
    numeric();products()
    print('共 %d 项断言，%d 项失败'%(TOTAL,FAIL))
    return int(FAIL>0)


if __name__=='__main__':raise SystemExit(main())
