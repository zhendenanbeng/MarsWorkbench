# -*- coding: utf-8 -*-
"""采样谱的半高宽：两侧半高交点线性插值，不把采样格宽当成峰宽。"""
from __future__ import annotations

import numpy as np


def half_maximum(x, y, index, allow_censored=False):
    """返回宽度及交点；默认两侧交点不完整返回 None。

    allow_censored=True 时可返回可观察宽度的下界，并标记 censored，不能
    将下界当完整 FWHM。用于筛选潜在带，完整宽度仍需拟合或更大测量窗口。

    y 必须已扣除基线；这是采样折线的宽度，不是亚采样的真实仪器峰形。
    不外推、不把 NaN/窗口边缘当成零；不提供对截断峰的完整宽度判断。
    """
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if x.ndim != 1 or x.shape != y.shape or not 0 < index < len(x)-1:
        return None
    height = y[index]
    if not np.isfinite(height) or height <= 0:
        return None
    half = height / 2
    left = right = int(index)
    while left > 0 and np.isfinite(y[left]) and y[left] > half:
        left -= 1
    while right < len(y)-1 and np.isfinite(y[right]) and y[right] > half:
        right += 1
    if (not left < index < right or not np.isfinite(y[left:right+1]).all()
            or not np.isfinite(x[left:right+1]).all() or not np.all(np.diff(x[left:right+1]) > 0)):
        return None
    censored = bool(y[left] > half or y[right] > half)
    if censored and not allow_censored:
        return None
    xl = (x[left] if y[left] > half else
          x[left] + (half-y[left])/(y[left+1]-y[left])*(x[left+1]-x[left]))
    xr = (x[right] if y[right] > half else
          x[right-1] + (half-y[right-1])/(y[right]-y[right-1])*(x[right]-x[right-1]))
    return {'fwhm':float(xr-xl), 'left':float(xl), 'right':float(xr),
            'channels':int(np.count_nonzero(y[left:right+1] > half)), 'censored':censored}
