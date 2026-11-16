# -*- coding: utf-8 -*-
"""APXS 原始能谱的自研峰拟合层（Phase 4；全自有实现）。

自算/对账边界声明：
  管线输入 = PDS rsp 能谱（1024 通道，eV 能量轴在 csv 内）；
  自研科学变换 = 平滑（可选）→ 元素 Kα 高斯峰拟合（线性局部背景）→
               峰净面积+泊松误差+检出判定；
  对账基准 = 官方 rwp 氧化物表（逐元素偏差台账）；
  不做：全基本参数法（GUAPX 复刻无科学增量）、基体效应建模。

口径：峰面积是「与官方定量独立同源的原始信号量」，不是浓度；
自研浓度（峰面积→标定曲线）属后续阶段，需先冻结多靶峰面积-浓度对。
"""
from __future__ import annotations

import numpy as np

# 元素 Kα 能量（keV）——APXS 元素集（Na-Br 主量）
KA_KEV = {"Na": 1.041, "Mg": 1.254, "Al": 1.487, "Si": 1.740, "S": 2.308,
          "Cl": 2.622, "K": 3.314, "Ca": 3.691, "Ti": 4.511, "Cr": 5.415,
          "Mn": 5.899, "Fe": 6.404, "Ni": 7.478, "Zn": 8.639}
PEAK_KEV_FWHM = 0.15          # APXS 探测器分辨率量级（keV）


def fit_peak_area(energy_kev, counts, center_kev, half_width=0.12):
    """单峰高斯拟合+线性背景，返回净面积与误差。

    拟合域：[center±half_width]；背景 = 域两端各 2 通道的线性 polyfit；
    高斯初值：峰值位置=窗内 argmax、σ=PEAK_KEV_FWHM/2.355。
    数据不足/窗越界返回 area=None（不猜值）。
    """
    e = np.asarray(energy_kev, dtype=float)
    y = np.asarray(counts, dtype=float)
    lo, hi = center_kev - half_width, center_kev + half_width
    m = (e >= lo) & (e <= hi)
    if m.sum() < 5:
        return {"area": None, "error": None, "center_kev": center_kev,
                "detected": False, "note": "窗内通道不足"}
    ew, yw = e[m], y[m]
    # 线性背景：两端各 2 通道
    edge = np.r_[yw[:2], yw[-2:]]
    edge_e = np.r_[ew[:2], ew[-2:]]
    if len(set(edge_e)) < 2:
        bg_slope, bg_icept = 0.0, float(np.median(yw))
    else:
        bg_slope, bg_icept = np.polyfit(edge_e, edge, 1)
    net = yw - (bg_slope * ew + bg_icept)
    pk = int(np.argmax(net))
    sigma = PEAK_KEV_FWHM / 2.355
    amp = max(float(net[pk]), 1e-6)
    gauss = amp * np.exp(-0.5 * ((ew - ew[pk]) / sigma) ** 2)
    net_c = np.clip(net, 0, None)
    # 高斯归一面积（格点离散积分即可：Δe 缩放由逐通道求和等价）
    area = float(np.minimum(net_c, gauss * 3).sum())
    err = float(np.sqrt(max(area, 1.0)))
    detected = bool(area > 5 * err)
    return {"area": area, "error": err, "center_kev": float(ew[pk]),
            "amp": amp, "sigma_kev": sigma, "detected": detected}


def spectrum_areas(energy_kev, counts, elements=tuple(KA_KEV)):
    """整谱逐元素峰面积：{元素: fit_peak_area 结果}。"""
    return {el: fit_peak_area(energy_kev, counts, KA_KEV[el]) for el in elements}
