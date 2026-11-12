# -*- coding: utf-8 -*-
"""元素谱线参考表与「检出元素」判定。

LIBS（ChemCam / SuperCam 共用，波长窗口 240–900 nm）：
  主动 LIBS 判元素看的是中性/一次电离原子的发射线。这里收录的是
  原始参考线为NIST空气波长，运行时转换为与仪器标定相同的真空波长。
  元素候选至少需要两条独立、无重叠的原子线；单线仍保留为谱线候选。

APXS（x 轴 keV）：
  X 射线荧光看 K 系线；另标注仪器/大气及Cm源Pu子体的背景参考，
  它们参与画参考线但不进「检出元素」名单。

判定口径（对每条参考线就地打分，不另做全局找峰）：
  要求连续测量段内存在局部极大、完整双侧背景与峰翼，局部差分MAD分数达标。
  重叠峰、分子带及仪器/大气线不作为目标元素的独立支持；不是定量反演。
"""

from __future__ import annotations

import numpy as np

SNR_MIN = 5.0

# Solid-state near-IR reference families, not NIST atomic lines or unique mineral labels.
# (key, label, center search interval µm, continuum shoulders µm, primary reference)
IR_REFERENCES = [
    ('water14','H₂O/OH 1.4', (1.38,1.45),(1.34,1.49),'https://doi.org/10.1029/2022JE007450'),
    ('water19','H₂O 1.9',(1.88,1.98),(1.82,2.04),'https://doi.org/10.1029/2022JE007450'),
    ('water21','结合水/单水合盐等 2.13',(2.10,2.16),(1.82,2.24),'https://doi.org/10.1029/2020JE006372'),
    ('carbonate23','CO₃相关/重叠 2.3',(2.30,2.35),(2.24,2.38),'https://doi.org/10.1029/2022JE007463'),
    ('carbonate25','CO₃相关/重叠 2.5',(2.50,2.55),(2.49,2.59),'https://doi.org/10.1029/2022JE007463'),
    ('sulfate24','硫酸盐等/非专属 2.4',(2.40,2.46),(2.36,2.51),'https://doi.org/10.1002/2014JE004707'),
    ('metaloh','Fe/Mg–OH等 2.28–2.33',(2.28,2.33),(2.24,2.38),'https://doi.org/10.1029/2022JE007450'),
    ('ferrous','Fe²⁺宽吸收约1.0（CP2窗口外）',(.85,1.15),(.75,1.25),'https://doi.org/10.1029/2022JE007463'),
]

def air_to_vacuum_nm(wavelength):
    """标准空气→真空波长；Peck & Reeder (1972)，NIST ASD 使用的色散关系。

    输入185–1700 nm。折射率自变量是真空波长，迭代反解；不修改产品的原始横轴。
    """
    air = np.asarray(wavelength, dtype=float)
    if not np.isfinite(air).all() or np.any((air < 185) | (air > 1700)):
        raise ValueError('空气波长必须有限且在185–1700 nm内')
    vac = air.copy()
    for _ in range(6):
        s2 = (1000.0 / vac) ** 2
        refractive = 1 + 1e-8 * (5791817.0 / (238.0185 - s2) + 167909.0 / (57.362 - s2))
        vac = air * refractive
    return float(vac) if vac.ndim == 0 else vac


# --- LIBS 发射线（原始参考nm，NIST标准空气波长） ----------------------------
# (元素符号, 波长, 人读标注)
LIBS_LINES_AIR = [
    ("Li", 610.35, "Li I 610"), ("Li", 670.78, "Li I 671"),
    ("C", 247.86, "C I 248"),
    ("H", 486.13, "Hβ 486"), ("H", 656.28, "Hα 656"),
    ("O", 777.19, "O I 777"), ("O", 844.64, "O I 845"),
    ("Na", 588.99, "Na I 589"), ("Na", 589.59, "Na I 590"),
    ("K", 766.49, "K I 766"), ("K", 769.90, "K I 770"),
    ("Ca", 393.37, "Ca II 393"), ("Ca", 396.85, "Ca II 397"),
    ("Ca", 422.67, "Ca I 423"), ("Ca", 612.22, "Ca I 612"), ("Ca", 616.22, "Ca I 616"),
    ("Mg", 279.55, "Mg II 280"), ("Mg", 280.27, "Mg II 280"), ("Mg", 285.21, "Mg I 285"),
    ("Al", 308.21, "Al I 308"), ("Al", 309.27, "Al I 309"),
    ("Al", 394.40, "Al I 394"), ("Al", 396.15, "Al I 396"),
    ("Si", 250.69, "Si I 251"), ("Si", 251.43, "Si I 251"), ("Si", 288.16, "Si I 288"),
    ("Si", 390.55, "Si I 391"),
    ("Fe", 238.20, "Fe II 238"), ("Fe", 259.84, "Fe II 260"),
    ("Fe", 273.96, "Fe I 274"), ("Fe", 404.58, "Fe I 405"),
    ("Ti", 334.94, "Ti II 335"), ("Ti", 336.12, "Ti II 336"), ("Ti", 338.38, "Ti II 338"),
    ("Mn", 403.08, "Mn I 403"), ("Mn", 403.45, "Mn I 403"),
    ("Ni", 341.48, "Ni I 341"), ("Ni", 352.45, "Ni I 352"),
    ("Cr", 425.43, "Cr I 425"), ("Cr", 427.48, "Cr I 427"),
    ("Zn", 330.26, "Zn I 330"), ("Zn", 334.50, "Zn I 335"),
    ("Cu", 324.75, "Cu I 325"), ("Cu", 327.40, "Cu I 327"),
    ("Sr", 407.77, "Sr II 408"), ("Sr", 421.55, "Sr II 422"),
    ("Ba", 455.40, "Ba II 455"), ("Ba", 614.17, "Ba II 614"),
    ("Rb", 780.03, "Rb I 780"),
    ("C2", 516.52, "C₂ 分子带 517"), ("CN", 388.31, "CN 分子带 388"),
]

# ChemCam/SuperCam 的波长标定都以真空线表为准。分子带旧参考位置未确认介质，
# 保留原值仅作候选，不把它们包装成NIST原子线，也不参与元素支持票数。
LIBS_LINES = [(e, air_to_vacuum_nm(w) if e not in ('C2','CN') else w, label)
              for e,w,label in LIBS_LINES_AIR]

# --- APXS 特征能量（keV） --------------------------------------------------
# instrument=True 的是仪器自身带来的线，画参考线但不算「目标元素」
APXS_LINES = [
    # MSL APXS公开工作能区约0.7–25 keV，O Kα不作为可测目标线。
    ("Na", 1.041, "Na Kα"), ("Mg", 1.254, "Mg Kα"),
    ("Al", 1.487, "Al Kα"), ("Si", 1.740, "Si Kα"), ("P", 2.014, "P Kα"),
    ("S", 2.308, "S Kα"), ("Cl", 2.621, "Cl Kα"), ("Ar", 2.957, "Ar Kα 大气", True),
    ("K", 3.312, "K Kα"), ("K", 3.589, "K Kβ"),
    ("Ca", 3.691, "Ca Kα"), ("Ca", 4.012, "Ca Kβ"),
    ("Ti", 4.510, "Ti Kα"), ("Ti", 4.931, "Ti Kβ"),
    ("V", 4.949, "V Kα"), ("Cr", 5.412, "Cr Kα"), ("Mn", 5.895, "Mn Kα"),
    ("Fe", 6.391, "Fe Kα"), ("Fe", 7.057, "Fe Kβ"),
    ("Ni", 7.471, "Ni Kα"), ("Zn", 8.615, "Zn Kα"), ("Br", 11.907, "Br Kα"),
    ("Sr", 14.140, "Sr Kα"),
    ("Zr", 15.770, "Zr Kα 仪器", True),
    # Rieder et al. 2003 §3.2/Fig9：244Cm的240Pu子体提供约14.3/18.4 keV激发光。
    # 源的散射峰不能当样品Sr等元素。这里只放近似背景参考，不作源强反演。
    ("Pu", 14.300, "Pu 子体源散射约14.3 keV", True),
    ("Pu", 18.400, "Pu 子体源散射约18.4 keV", True),
]


def robust_sigma(v):
    """MAD 口径的稳健标准差估计。"""
    v = np.asarray(v, dtype=float)
    v = v[np.isfinite(v)]
    if v.size < 8:
        return float(np.std(v)) if v.size else 0.0
    med = np.median(v)
    return float(1.4826 * np.median(np.abs(v - med))) or float(np.std(v)) or 1e-12


def _line_evidence(x, y, xline, tol):
    """连续测量段内找局部峰，双侧背景+局部差分MAD；无覆盖不猜测。

    这是浏览用的候选筛选，不是去混叠、拟合定量或已校准的检测概率。
    """
    unavailable = {'snr':None,'sigma':None,'peak_x':None,'peak_index':None,'fwhm':None,
                   'interval':None,'available':False,'reason':'无完整局部覆盖'}
    if xline < x[0] or xline > x[-1]:return unavailable
    center = int(np.argmin(abs(x-xline)))
    dx = np.diff(x)
    step = float(np.median(dx[max(0,center-12):min(len(dx),center+12)]))
    # 禁止跨光谱仪空档借另一段的边界当背景。
    breaks = np.flatnonzero(dx > 5*step)
    left = int(breaks[breaks < center][-1]+1) if np.any(breaks < center) else 0
    right = int(breaks[breaks >= center][0]+1) if np.any(breaks >= center) else len(x)
    local = np.arange(left,right)[abs(x[left:right]-xline) <= 6*tol]
    if len(local)<9:return unavailable
    xx,yy=x[local],y[local]
    before=(xx < xline-2*tol);after=(xx > xline+2*tol)
    if before.sum()<3 or after.sum()<3:return unavailable
    xa,xb=float(np.median(xx[before])),float(np.median(xx[after]))
    ya,yb=float(np.median(yy[before])),float(np.median(yy[after]))
    residual=yy-(ya+(yb-ya)*(xx-xa)/(xb-xa))
    if not np.isfinite(residual).all():
        raise ValueError('局部背景计算发生数值溢出，不能输出谱线判定')
    # 双侧噪声分别差分，不能把侧窗之间的跳跃/峰体计入噪声。
    differences=np.r_[np.diff(residual[before]),np.diff(residual[after])]
    scale=max(robust_sigma(differences)/np.sqrt(2),
              32*np.finfo(float).eps*max(1.,float(np.max(abs(yy)))))
    candidates=[j for j in range(1,len(xx)-1) if abs(xx[j]-xline)<=tol
                and residual[j]>residual[j-1] and residual[j]>=residual[j+1]]
    if not candidates:
        return dict(unavailable,available=True,sigma=float(scale),reason='窗口内无局部峰')
    j=max(candidates,key=lambda k:residual[k])
    from sherloc.peaks import half_maximum
    width=half_maximum(xx,residual,j)
    diagnostic=half_maximum(xx,residual,j,allow_censored=True)
    snr=float(residual[j]/scale)
    reason='局部峰候选'
    if width is None:reason='峰翼不完整，可能为连续背景或边界峰'
    elif width['fwhm'] > 4*tol:reason='峰过宽，未作为元素线证据'
    return {'snr':snr,'sigma':float(scale),'peak_x':float(xx[j]),'peak_index':int(local[j]),
            'fwhm':width['fwhm'] if width else None,
            'fwhm_observed':diagnostic['fwhm'] if diagnostic else None,
            'fwhm_censored':bool(diagnostic and diagnostic['censored']),
            'interval':(width['left'],width['right']) if width else None,
            'available':True,'reason':reason,
            'shape_ok':width is not None and width['fwhm']<=4*tol}


def _line_snr(x, y, xline, tol):
    evidence=_line_evidence(np.asarray(x,float),np.asarray(y,float),xline,tol)
    return evidence['snr'],evidence['sigma']


def match_lines(x, y, table, tol, snr_min=SNR_MIN, min_lines=2, wavelength_frame=None):
    """逐线打分。返回 (标注列表, 检出元素列表)。

    标注列表：[{x, label, elem, snr, detected}]，全部线都返回
    （前端只画 detected 的，但完整表能显示「参考线窗口内有没有」）。
    元素候选：去掉分子带/背景/重叠线，按独立峰数排序；LIBS默认需2条。
    """
    x,y=np.asarray(x,float),np.asarray(y,float)
    if (x.ndim!=1 or x.shape!=y.shape or len(x)<3 or not np.isfinite(x).all()
            or not np.isfinite(y).all() or not np.all(np.diff(x)>0)):
        raise ValueError('谱线判定需要有限、严格递增且与强度等长的坐标轴')
    if not np.isfinite([tol,snr_min]).all() or tol<=0 or snr_min<=0 or min_lines<1:
        raise ValueError('无效的谱线判定参数')
    marks, votes = [], {}
    for row in table:
        elem, xline, label = row[0], float(row[1]), row[2]
        instrument = bool(row[3]) if len(row) > 3 else False
        evidence=_line_evidence(x,y,xline,tol)
        snr=evidence['snr']
        det = snr is not None and snr >= snr_min and evidence.get('shape_ok',False)
        # Goetz et al. 2023, 10.1029/2021JE007101 §2：Cu324/Cu327两条带
        # 均有Ti等谱线混合。两条同时出现也不是Cu的独立确认，须联合拟合。
        known_blend = (elem=='Cu' and wavelength_frame=='vacuum_nm' and 324<xline<328)
        marks.append({"x": xline, "label": label, "elem": elem,
                      "snr": round(float(snr), 1) if snr is not None else None,
                      'snr_raw':snr,
                      "detected": bool(det), "instrument": instrument,
                      'available':evidence['available'],'peak_x':evidence['peak_x'],
                      'peak_index':evidence['peak_index'],'fwhm':evidence['fwhm'],
                      'fwhm_observed':evidence.get('fwhm_observed'),
                      'fwhm_censored':evidence.get('fwhm_censored',False),
                      'interval':evidence['interval'],'reason':evidence['reason'],
                      'ambiguous_with':(['Ti等谱线（文献Cu324/Cu327混合区，需联合拟合）'] if known_blend else []),
                      'element_supported':False,
                      'reference_frame':('legacy_unspecified' if elem in ('C2','CN') else wavelength_frame)})
    for i,mark in enumerate(marks):
        if not mark['detected']:continue
        for other in marks[i+1:]:
            if not other['detected'] or other['elem']==mark['elem']:continue
            # 同一局部极大或相互重叠的半高区间不能作为两个元素的独立证据。
            a,b=mark['interval'],other['interval']
            if mark['peak_index']==other['peak_index'] or (a and b and max(a[0],b[0])<=min(a[1],b[1])):
                mark['ambiguous_with'].append(other['label'])
                other['ambiguous_with'].append(mark['label'])
    for mark in marks:
        if mark['detected'] and not mark['instrument'] and not mark['ambiguous_with'] and mark['elem'] not in ('C2','CN'):
            votes.setdefault(mark['elem'],[]).append(mark['interval'])
    # 同一元素的两个参考位置也可能落在同一个未分辨包络上。
    # 合并半高区间相交的候选；不能把一个包络上的噪声极大值算成独立谱线。
    for elem,intervals in votes.items():
        independent=[]
        for low,high in sorted(intervals):
            if independent and low<=independent[-1][1]:
                independent[-1][1]=max(high,independent[-1][1])
            else:independent.append([low,high])
        votes[elem]=independent
    elements = [e for e,peaks in sorted(votes.items(),key=lambda kv:-len(kv[1])) if len(peaks)>=min_lines]
    for mark in marks:
        mark['independent_lines']=len(votes.get(mark['elem'],()))
        mark['element_supported']=mark['elem'] in elements and mark['detected'] and not mark['ambiguous_with'] and not mark['instrument']
        if mark['ambiguous_with']:mark['reason']='与其他参考线重叠：'+'、'.join(mark['ambiguous_with'])
        elif mark['detected'] and mark['instrument']:mark['reason']='仪器/大气背景候选，不计入目标元素'
        elif mark['detected'] and mark['elem'] not in elements:mark['reason']='单线或非独立谱线候选，未达到元素支持条件'
    return marks, elements
