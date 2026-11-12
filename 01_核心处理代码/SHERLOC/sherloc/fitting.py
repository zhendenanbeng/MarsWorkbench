# -*- coding: utf-8 -*-
"""非负峰形系数与自由线性背景的有界最小二乘；仅依赖numpy。"""
from __future__ import annotations

from itertools import combinations
import numpy as np


def background_context(x, y):
    """背景用中心化坐标求解，输出仍兼容原始坐标的斜率/截距。"""
    x,y=np.array(x,dtype=float,copy=True),np.array(y,dtype=float,copy=True)
    if (x.ndim!=1 or x.shape!=y.shape or len(x)<3 or not np.isfinite(x).all()
            or not np.isfinite(y).all() or not np.all(np.diff(x)>0)):
        raise ValueError('拟合需要有限、递增且等长的一维数据')
    origin=float(np.mean(x));scale=float(np.ptp(x))
    b=np.column_stack(((x-origin)/scale,np.ones(len(x))))
    q,_=np.linalg.qr(b,mode='reduced')
    projected=y-q@(q.T@y)
    return {'x':x,'y':y,'background':b,'q':q,'projected_y':projected,
            'origin':origin,'scale':scale}


def solve_nonnegative_profiles(x, y, profiles, context=None):
    """峰形系数>=0，线性底不受符号约束。

    先投影掉自由背景，再在至多4个峰形系数的全部活动面上求解；
    不能用求解后截断系数替代有约束最小二乘。返回系数/平方误差/同一模型。
    """
    ctx=background_context(x,y) if context is None else context
    if context is not None and (not np.array_equal(np.asarray(x),ctx['x']) or not np.array_equal(np.asarray(y),ctx['y'])):
        raise ValueError('背景上下文与当前数据不一致')
    p=np.asarray(profiles,dtype=float)
    if p.ndim!=2 or p.shape[0]!=len(ctx['x']) or not 1<=p.shape[1]<=4 or not np.isfinite(p).all():
        raise ValueError('峰形设计矩阵必须有限，且有1–4列')
    q=ctx['q'];a=p-q@(q.T@p);target=ctx['projected_y']
    unconstrained=np.linalg.lstsq(a,target,rcond=None)[0]
    if np.all(unconstrained>=0):
        best=unconstrained
    else:
        n=p.shape[1];best=np.zeros(n);best_error=float(target@target)
        # 非负二次问题的最优解必在某个活动面。每个面的背景已精确消去。
        for size in range(1,n+1):
            for active in combinations(range(n),size):
                active=list(active)
                coefficient=np.linalg.lstsq(a[:,active],target,rcond=None)[0]
                if np.any(coefficient<0):continue
                residual=target-a[:,active]@coefficient
                error=float(residual@residual)
                if error<best_error:
                    best_error=error;best=np.zeros(n);best[active]=coefficient
    background=np.linalg.lstsq(ctx['background'],ctx['y']-p@best,rcond=None)[0]
    model=p@best+ctx['background']@background
    slope=float(background[0]/ctx['scale'])
    intercept=float(background[1]-slope*ctx['origin'])
    coefficient=np.r_[best,slope,intercept]
    residual=ctx['y']-model
    return coefficient,float(residual@residual),model
