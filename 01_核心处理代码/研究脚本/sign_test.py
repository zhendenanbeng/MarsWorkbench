# -*- coding: utf-8 -*-
"""补上一个所有前面判据都漏掉的方向：**负相关**。

上一轮所有判据（NCC 峰值、平移曲面、互信息、翻转组合）都在找"越像越好"。
但两台相机的照明几何完全不同：ACI 是环绕 LED、WATSON 是另一组 LED，
同一颗矿物颗粒在一个方向下是镜面高光、在另一个方向下可能就是暗坑——
**同一块地方在两台相机下的颗粒级明暗关系可能是反的**。
这种情况下正相关本来就该是 0，而对的那一对会表现为**很强的负相关**。

所以这一步同时报 max 和 min，并用两组合成对照把两者的"上限"标出来：
  B1  退化但不反相  → 报 max（"同相"时的上限）
  B2  退化且反相    → 报 min（"反相"时的上限）
  C   真问题        → max 和 min 都报，并且看几何预测位置上取到什么值

用法：python sign_test.py
"""
from __future__ import annotations

import os
import sys

import numpy as np
import cv2

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from colorize_pair import (ACI_IMG, ACI_XML, WAT_IMG, WAT_XML, OUT,      # noqa: E402
                           homography, load_cam, read_bsq, to8)

ACI_UM = 6.66
BORDER = 40
ANG0, S0 = 37.0, 0.42
PX_PASTE = (300, 260)


def bandpass(img, lo, hi):
    g = img.astype(np.float32)
    return cv2.GaussianBlur(g, (0, 0), lo) - cv2.GaussianBlur(g, (0, 0), hi)


def zs(a):
    s = float(a.std())
    return (a - float(a.mean())) / (s if s > 1e-9 else 1.0)


def search(img, templ, ang, scale, mode="max"):
    ns, nl = templ.shape[1], templ.shape[0]
    M = cv2.getRotationMatrix2D((ns / 2.0, nl / 2.0), ang, scale)
    c = np.array([[0, 0], [ns, 0], [ns, nl], [0, nl]], float)
    mc = (M @ np.hstack([c, np.ones((4, 1))]).T).T
    x0, y0 = int(np.floor(mc[:, 0].min())), int(np.floor(mc[:, 1].min()))
    x1, y1 = int(np.ceil(mc[:, 0].max())), int(np.ceil(mc[:, 1].max()))
    M2 = M.copy()
    M2[0, 2] -= x0
    M2[1, 2] -= y0
    t = np.ascontiguousarray(cv2.warpAffine(templ, M2, (x1 - x0, y1 - y0),
                                            flags=cv2.INTER_LINEAR).astype(np.float32))
    m = cv2.warpAffine(np.full((nl, ns), 255, np.uint8), M2, (x1 - x0, y1 - y0),
                       flags=cv2.INTER_NEAREST)
    m = cv2.erode(m, np.ones((31, 31), np.uint8))
    if m.sum() == 0 or t.shape[0] > img.shape[0] or t.shape[1] > img.shape[1]:
        return float("nan"), None
    res = cv2.matchTemplate(img.astype(np.float32), t, cv2.TM_CCORR_NORMED,
                            mask=m.astype(np.float32))
    inner = res[BORDER:-BORDER, BORDER:-BORDER]
    if inner.size == 0:
        return float("nan"), None
    # 掩膜版归一化相关在"窗口整体落在置零区"时会是 nan，而 numpy 的 argmax/argmin
    # 会把 nan 当极值返回——必须先把 nan 剔掉，否则整个搜索静默退化成 nan。
    good = np.isfinite(inner)
    if not good.any():
        return float("nan"), None
    if mode == "max":
        k = np.unravel_index(int(np.argmax(np.where(good, inner, -np.inf))), inner.shape)
    else:
        k = np.unravel_index(int(np.argmin(np.where(good, inner, np.inf))), inner.shape)
    return float(inner[k]), (int(k[1] + BORDER), int(k[0] + BORDER))


def degrade(patch, from_scale, to_scale, flip, rng):
    tiny = cv2.resize(cv2.GaussianBlur(patch, (0, 0), 1.44),
                      None, fx=to_scale / from_scale, fy=to_scale / from_scale,
                      interpolation=cv2.INTER_AREA)
    back = cv2.resize(tiny, (patch.shape[1], patch.shape[0]),
                      interpolation=cv2.INTER_LINEAR)
    out = zs(back)
    if flip:
        out = -out
    out = out + 0.25 * rng.standard_normal(out.shape)
    out = out + np.linspace(-0.6, 0.6, out.shape[1])[None, :]
    return out


def main():
    wat = load_cam(WAT_XML)
    aci = load_cam(ACI_XML)
    wr = read_bsq(WAT_IMG, wat)
    ar = read_bsq(ACI_IMG, aci)
    b = [i - 1 for i in wat["rgb"]]
    wgray = cv2.cvtColor(np.stack([to8(wr[b[i]]) for i in range(3)], axis=-1),
                         cv2.COLOR_RGB2GRAY)
    aci8 = to8(ar[0])

    aci_ppm, wat_ppm = 1000.0 / ACI_UM, 1000.0 / ACI_UM * 0.353
    A = np.ascontiguousarray(zs(bandpass(aci8, 0.06 * aci_ppm / 2.355,
                                         0.5 * aci_ppm / 2.355)).astype(np.float32))
    W = np.ascontiguousarray(zs(bandpass(wgray, 0.06 * wat_ppm / 2.355,
                                         0.5 * wat_ppm / 2.355)).astype(np.float32))
    rng = np.random.default_rng(7)
    ns, nl = aci["ns"], aci["nl"]
    M = cv2.getRotationMatrix2D((ns / 2.0, nl / 2.0), ANG0, S0)
    t0 = cv2.warpAffine(A, M, (ns, nl), flags=cv2.INTER_LINEAR)
    m0 = cv2.warpAffine(np.full(A.shape, 255, np.uint8), M, (ns, nl),
                        flags=cv2.INTER_NEAREST) > 127
    ys, xs = np.where(m0)
    y0, y1, x0, x1 = ys.min(), ys.max() + 1, xs.min(), xs.max() + 1
    t0c, m0c = t0[y0:y1, x0:x1], m0[y0:y1, x0:x1]
    tx, ty = PX_PASTE

    def paste(patch):
        cv = W.copy()
        cv[ty:ty + patch.shape[0], tx:tx + patch.shape[1]][m0c] = patch[m0c]
        return cv

    # B1：退化、不反相 → 找 max
    p1, l1 = search(paste(degrade(t0c, S0, 0.353, False, rng)), A, ANG0, S0, "max")
    print("B1 退化不反相 → 找 max：峰 %.3f 在 %s（贴在 %s）%s"
          % (p1, l1, PX_PASTE, "✓" if p1 > 0.4 else "✗ 判据本身没分辨力"))
    # B2：退化且反相 → 找 min
    p2, l2 = search(paste(degrade(t0c, S0, 0.353, True, rng)), A, ANG0, S0, "min")
    print("B2 退化且反相 → 找 min：谷 %.3f 在 %s → %s"
          % (p2, l2, "反相也能抓到" if p2 < -0.4 else "反相抓不到"))
    # 反相情况下的 max 是多少（这就是"漏掉负相关"时会看到的值）
    p2b, _ = search(paste(degrade(t0c, S0, 0.353, True, rng)), A, ANG0, S0, "max")
    print("   同一张（反相）去找 max 只能得到 %.3f ← 只看 max 的判据会因此误判为『没对上』"
          % p2b)

    # C：真问题
    H, dst = homography(aci, wat, 129.0)
    print("\nC 真问题（ACI vs WATSON，位置全搜，排除边界）：")
    best_max, best_min = -9, 9
    for scale in (0.30, 0.353, 0.44):
        for ang in range(0, 360, 10):
            mx, lm = search(W, A, ang, scale, "max")
            mn, ln = search(W, A, ang, scale, "min")
            if np.isfinite(mx) and mx > best_max:
                best_max, argmax = mx, (scale, ang, lm)
            if np.isfinite(mn) and mn < best_min:
                best_min, argmin = mn, (scale, ang, ln)
        print("   尺度 %.3f 扫完：当前 max %.3f，min %.3f" % (scale, best_max, best_min))
    print("   全搜最大正相关 %.3f 在 尺度 %.3f 转角 %d° 位置 %s"
          % (best_max, argmax[0], argmax[1], argmax[2]))
    print("   全搜最小（最负）相关 %.3f 在 尺度 %.3f 转角 %d° 位置 %s"
          % (best_min, argmin[0], argmin[1], argmin[2]))

    val, loc = search(W, A, 0.0, 0.353, "max")
    print("   几何预测那一组（恒等、0.353）：max %.3f 在 %s；ACI 左上角几何落点 (%.0f,%.0f)"
          % (val, loc, dst[0, 0], dst[0, 1]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
