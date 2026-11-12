# -*- coding: utf-8 -*-
"""基于估计噪声模型的峰组重复性检验；不等同矿物真值或化学确认。"""
from __future__ import annotations
import numpy as np


def peak_hits(x,spectra,threshold=3.,min_fwhm=30.,max_fwhm=200.):
    """与默认离散候选宽度口径一致；向量化初筛，保留原始跨格宽度。

    输入已按逐点判定用σ归一；返回(峰位,点号)。不拟合Voigt，不作物相判定。
    """
    x,y=np.asarray(x,float),np.asarray(spectra,float)
    selected=(y[:,1:-1]>y[:,:-2])&(y[:,1:-1]>=y[:,2:])&(y[:,1:-1]>=threshold)
    rows,columns=np.nonzero(selected);hits=[]
    for row,j0 in zip(rows,columns):
        j=int(j0)+1;lo=hi=j;half=y[row,j]/2
        while lo>0 and y[row,lo]>half:lo-=1
        while hi<len(x)-1 and y[row,hi]>half:hi+=1
        width=x[hi]-x[lo]
        if min_fwhm<=width<=max_fwhm:hits.append((float(x[j]),int(row)))
    return hits


def max_cluster_count(hits,cluster_cm=15.):
    """完整搜索范围中的最大重复点数，包含自适应挑峰/聚类的搜索效应。"""
    maximum=0;first=last=None;points=set()
    for position,point in sorted(hits):
        if first is None or position-last>cluster_cm or position-first>4*cluster_cm:
            maximum=max(maximum,len(points));first=position;points=set()
        points.add(point);last=position
    return max(maximum,len(points))


def calibrate(x,quiet_normalized,groups,simulations=199,seed=20261023,alpha=.05,model='gaussian_ar1'):
    """用Gaussian AR(1)或联合通道块重放，模拟全窗口最大组点数。

    用全窗口最大统计量作参照，避免选完最高组后当作预指定检验。
    输出平滑的模拟超越率，不把估计模型的模拟频率包装成已校准p值或物相置信度。
    """
    x,q=np.asarray(x,float),np.asarray(quiet_normalized,float)
    if (x.ndim!=1 or len(x)<10 or not np.isfinite(x).all() or not np.all(np.diff(x)>0)
            or q.ndim!=2 or not 2<=q.shape[0]<=10000 or q.shape[1]<10 or not np.isfinite(q).all()):
        raise ValueError('重复性检验需要有限递增指纹轴和至少2点×10通道的安静区')
    if not isinstance(simulations,int) or isinstance(simulations,bool) or simulations<19 or simulations>1999 or not 0<alpha<1 or model not in ('gaussian_ar1','joint_blocks'):
        raise ValueError('模拟次数需19–1999，参考尾部比例需在0–1内，模型需为gaussian_ar1或joint_blocks')
    if any(not g.get('pts') or any(not isinstance(i,(int,np.integer)) or not 0<=i<len(q) for i in g['pts']) for g in groups):
        raise ValueError('峰组必须包含有效的原始点号')
    center=np.median(q,axis=1);scale=1.4826*np.median(abs(q-center[:,None]),axis=1)
    if not np.isfinite(scale).all() or np.any(scale<=1e-12):
        return {'status':'unsupported','model':model,'simulations':simulations,'seed':seed,
                'reason':'安静区噪声尺度为零或无法估计',
                'groups':[{'group_index':i,'observed_points':len(set(g['pts'])),
                           'simulation_exceedance_fraction':None,'exceeds_simulated_reference':False,
                           'status':'噪声模型无法估计'} for i,g in enumerate(groups)]}
    z=(q-center[:,None])/scale[:,None]
    # 极端点影响诊断，但不让少数尖刺主导相关系数拟合。
    robust=np.clip(z,-3,3);robust-=robust.mean(axis=1,keepdims=True)
    denominator=np.sum(robust[:,:-1]**2)
    rho=float(np.sum(robust[:,:-1]*robust[:,1:])/denominator) if denominator else 0.
    raw_rho=rho;rho=float(np.clip(rho,-.95,.95))
    tail_fraction=float(np.mean(abs(z)>3));diagnostics=[]
    if model=='gaussian_ar1' and tail_fraction>.01:diagnostics.append('安静区厚尾/尖刺多，Gaussian AR(1)可能低估极端噪声')
    if model=='gaussian_ar1' and abs(raw_rho)>=.95:diagnostics.append('AR(1)相关系数触及稳定性限制')
    # 高阶相关与AR(1)预期严重不一致时，不自动解释模拟参考范围。
    lag2=float(np.mean(robust[:,:-2]*robust[:,2:])/max(np.mean(robust**2),1e-12))
    if model=='gaussian_ar1' and abs(lag2-rho*rho)>.15:diagnostics.append('高阶相关不符合简单AR(1)模型')
    block_length=min(q.shape[1],max(8,int(np.ceil(200/np.median(np.diff(x))))+2))
    rng=np.random.default_rng(seed);maxima=[]
    for _ in range(simulations):
        if model=='gaussian_ar1':
            noise=rng.normal(size=(len(q),len(x)))
            for col in range(1,len(x)):
                noise[:,col]=rho*noise[:,col-1]+np.sqrt(1-rho*rho)*noise[:,col]
            noise=noise*scale[:,None]+center[:,None]
        else:
            # 全部点共同抽取连续通道块，保留厚尾和安静区内的点间共同结构。
            starts=rng.integers(0,q.shape[1],size=int(np.ceil(len(x)/block_length)))
            indices=((starts[:,None]+np.arange(block_length))%q.shape[1]).ravel()[:len(x)]
            noise=q[:,indices]
        maxima.append(max_cluster_count(peak_hits(x,noise)))
    maxima=np.asarray(maxima);records=[]
    for i,group in enumerate(groups):
        count=len(set(group['pts']));p=float((1+np.count_nonzero(maxima>=count))/(simulations+1))
        supported=not diagnostics and p<=alpha
        records.append({'group_index':i,'observed_points':count,'simulation_exceedance_fraction':p,
                        'exceeds_simulated_reference':supported,
                        'status':'高于本次模拟参考范围' if supported else ('模型诊断待复核' if diagnostics else '未高于本次模拟参考范围')})
    return {'status':'diagnostics_failed' if diagnostics else 'reference_only','model':model,
            'simulations':simulations,'seed':seed,'alpha':alpha,'monte_carlo_resolution':1/(simulations+1),
            'rho':rho,'lag2':lag2,'tail_fraction':tail_fraction,'diagnostic_warnings':diagnostics,
            'block_length':block_length if model=='joint_blocks' else None,
            'detector':{'threshold_sigma':3.,'min_width_cm-1':30.,'max_width_cm-1':200.,'cluster_gap_cm-1':15.,'max_cluster_span_cm-1':60.},
            'row_centers':center.tolist(),'row_scales':scale.tolist(),
            'null_max_cluster_points':maxima.tolist(),'groups':records,
            'limits':'模拟超越率只对当前估计模型作描述，不是已校准p值或矿物置信度。安静区噪声可迁移到指纹区是一项未独立验证的假设；两者的波长响应和背景曲率可能不同。Gaussian模型未描述所有点间相关；经验块模型受有限安静区和拼接边界影响。超出噪声仍可能是仪器背景；未超出也不等于没有弱信号。'}
