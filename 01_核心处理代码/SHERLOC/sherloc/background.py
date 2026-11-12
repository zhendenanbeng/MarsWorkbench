# -*- coding: utf-8 -*-
"""供验证的对称稳健多项式背景；避免单侧裁剪造成正残差偏置。"""
import numpy as np


def symmetric_baseline(x,y,degree=5,n_iter=10,clip=1.5):
    x,y=np.asarray(x,float),np.asarray(y,float)
    if (x.ndim!=1 or x.shape!=y.shape or len(x)<degree+3 or not np.isfinite(x).all()
            or not np.isfinite(y).all() or not np.all(np.diff(x)>0)):
        raise ValueError('背景拟合需要足够、有限、递增且等长的一维数据')
    xx=(x-x.mean())/np.ptp(x);keep=np.ones(len(x),bool)
    coefficient=np.polyfit(xx,y,degree)
    for _ in range(n_iter):
        residual=y-np.polyval(coefficient,xx)
        center=float(np.median(residual))
        # 使用全部残差的稳健尺度，不能在已经截断的子集反复估计σ使其逐轮缩小。
        scale=float(1.4826*np.median(abs(residual-center)))
        if scale<=32*np.finfo(float).eps*max(1.,float(np.max(abs(y)))):break
        new=abs(residual-center)<clip*scale
        if new.sum()<degree+3 or np.array_equal(new,keep):break
        keep=new;coefficient=np.polyfit(xx[keep],y[keep],degree)
    return np.polyval(coefficient,xx)
