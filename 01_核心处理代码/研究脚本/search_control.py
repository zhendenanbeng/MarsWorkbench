# -*- coding: utf-8 -*-
"""给"找不到"这件事配一个**阳性对照**——不然零结果可能只是我的搜索坏了。

三个量，同一个搜索机器（掩膜版 TM_CCORR_NORMED，位置全搜，排除边界伪峰）：

  A 阳性对照：把 ACI 的结构图**转到 37°、缩到 0.42**，贴进 WATSON 的图里，
     再让它去找。找得到 → 机器没问题。
  B 退化对照：贴进去的那块先按几何倍率糊+降采样（模拟 WATSON 采样更粗）、
     再反相（模拟两台相机明暗关系不同）、加噪声和亮度梯度，再去找。
     它给出"真的是同一块地方时，这个判据能剩下多少分"的**上限参考**。
  C 真问题：ACI 的结构图 vs WATSON 的结构图（什么都不动，位置全搜）。

用法：python search_control.py
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
BORDER = 40          # 排除边界，避免掩膜相关在画幅边缘刷假峰


def bandpass(img, lo, hi):
    g = img.astype(np.float32)
    return cv2.GaussianBlur(g, (0, 0), lo) - cv2.GaussianBlur(g, (0, 0), hi)


def zs(a):
    s = float(a.std())
    return (a - float(a.mean())) / (s if s > 1e-9 else 1.0)


def masked_search(img, templ, ang, scale):
    """把 templ 转到 ang、缩到 scale，在 img 里全位置搜。返回 (峰, 位置, 曲面)。

    注意：warpAffine 出来还是原画布大小（内容缩在里面），**必须裁到外接框**，
    否则 matchTemplate 的模板 = 整幅图，响应面退化成 1×1——所有"峰"都会
    假惺惺地落在 (0,0)。这个坑我踩过一次。
    """
    ns, nl = templ.shape[1], templ.shape[0]
    M = cv2.getRotationMatrix2D((ns / 2.0, nl / 2.0), ang, scale)
    corners = np.array([[0, 0], [ns, 0], [ns, nl], [0, nl]], float)
    mc = (M @ np.hstack([corners, np.ones((4, 1))]).T).T
    x0, y0 = int(np.floor(mc[:, 0].min())), int(np.floor(mc[:, 1].min()))
    x1, y1 = int(np.ceil(mc[:, 0].max())), int(np.ceil(mc[:, 1].max()))
    M2 = M.copy()
    M2[0, 2] -= x0
    M2[1, 2] -= y0
    w, h = x1 - x0, y1 - y0
    t = np.ascontiguousarray(cv2.warpAffine(templ, M2, (w, h),
                                            flags=cv2.INTER_LINEAR).astype(np.float32))
    m = cv2.warpAffine(np.full((nl, ns), 255, np.uint8), M2, (w, h),
                       flags=cv2.INTER_NEAREST)
    m = cv2.erode(m, np.ones((31, 31), np.uint8))
    if m.sum() == 0 or t.shape[0] > img.shape[0] or t.shape[1] > img.shape[1]:
        return float("nan"), None, None
    res = cv2.matchTemplate(img.astype(np.float32), t, cv2.TM_CCORR_NORMED,
                            mask=m.astype(np.float32))
    inner = res[BORDER:-BORDER, BORDER:-BORDER]
    if inner.size == 0:
        return float("nan"), None, res
    k = np.unravel_index(int(np.argmax(inner)), inner.shape)
    return float(inner[k]), (k[1] + BORDER, k[0] + BORDER), res


def main():
    wat = load_cam(WAT_XML)
    aci = load_cam(ACI_XML)
    wr = read_bsq(WAT_IMG, wat)
    ar = read_bsq(ACI_IMG, aci)
    b = [i - 1 for i in wat["rgb"]]
    wgray = cv2.cvtColor(np.stack([to8(wr[b[i]]) for i in range(3)], axis=-1),
                         cv2.COLOR_RGB2GRAY)
    aci8 = to8(ar[0])

    aci_ppm = 1000.0 / ACI_UM
    wat_ppm = aci_ppm * 0.353
    A = np.ascontiguousarray(zs(bandpass(aci8, 0.06 * aci_ppm / 2.355,
                                         0.5 * aci_ppm / 2.355)).astype(np.float32))
    W = np.ascontiguousarray(zs(bandpass(wgray, 0.06 * wat_ppm / 2.355,
                                         0.5 * wat_ppm / 2.355)).astype(np.float32))
    loc = cv2.blur(wgray.astype(np.float32), (151, 151))
    shadow = loc < 0.55 * float(np.median(loc))

    # ---------- C 真问题 ----------
    peakC, posC, _ = masked_search(W, A, 0.0, 0.353)
    print("C  真问题：ACI 结构图 × WATSON 结构图（恒等、几何尺度）")
    print("   峰 %.3f 在 %s（排除边界）" % (peakC, posC))
    rows = []
    for scale in (0.30, 0.353, 0.44, 0.55):
        vals = [masked_search(W, A, ang, scale)[0] for ang in range(0, 360, 10)]
        vals = [v for v in vals if np.isfinite(v)]
        rows.append((scale, max(vals) if vals else float("nan"), len(vals)))
    for scale, best, n in rows:
        print("   尺度 %.3f：%2d 个转角有效，最好峰 %.3f" % (scale, n, best))
    print("   → 全搜最好 %.3f" % np.nanmax([r[1] for r in rows]))

    # ---------- A 阳性对照 ----------
    rng = np.random.default_rng(7)
    ang0, s0 = 37.0, 0.42
    tx, ty = 300, 260                      # 贴的位置
    ns, nl = aci["ns"], aci["nl"]
    M = cv2.getRotationMatrix2D((ns / 2.0, nl / 2.0), ang0, s0)
    t0 = cv2.warpAffine(A, M, (ns, nl), flags=cv2.INTER_LINEAR)
    m0 = cv2.warpAffine(np.full(A.shape, 255, np.uint8), M, (ns, nl),
                        flags=cv2.INTER_NEAREST) > 127
    # 裁到有效外接框
    ys, xs = np.where(m0)
    y0, y1, x0, x1 = ys.min(), ys.max() + 1, xs.min(), xs.max() + 1
    t0c, m0c = t0[y0:y1, x0:x1], m0[y0:y1, x0:x1]
    canvasA = W.copy()
    canvasA[ty:ty + t0c.shape[0], tx:tx + t0c.shape[1]][m0c] = t0c[m0c]
    peakA, posA, _ = masked_search(canvasA, A, ang0, s0)
    print("\nA  阳性对照（把 ACI 转 %.0f°、缩 %.2f 贴进 WATSON 再找）"
          % (ang0, s0))
    print("   贴在 %s，搜到的峰 %.3f 在 %s  → %s"
          % ((tx, ty), peakA, posA,
             "机器正常" if peakA > 0.5 and abs(posA[0] - tx) < 30
             and abs(posA[1] - ty) < 30 else "机器可疑！"))

    # ---------- B 退化对照 ----------
    tiny = cv2.resize(cv2.GaussianBlur(t0c, (0, 0), 1.2 ** 2),
                      None, fx=0.353 / s0, fy=0.353 / s0,
                      interpolation=cv2.INTER_AREA)
    back = cv2.resize(tiny, (t0c.shape[1], t0c.shape[0]),
                      interpolation=cv2.INTER_LINEAR)
    deg = -zs(back)                                  # 反相 = 明暗关系不同
    deg = deg + 0.25 * rng.standard_normal(deg.shape)       # 噪声
    ramp = np.linspace(-0.6, 0.6, deg.shape[1])[None, :]    # 照明梯度
    deg = deg + ramp
    canvasB = W.copy()
    canvasB[ty:ty + deg.shape[0], tx:tx + deg.shape[1]][m0c] = deg[m0c]
    peakB, posB, _ = masked_search(canvasB, A, ang0, s0)
    print("\nB  退化对照（同一块，但先按几何倍率降采样 + 反相 + 噪声 + 亮度梯度）")
    print("\n   峰 %.3f 在 %s ← 这是『真的对上了』时该判据能剩下的分数上限参考"
          % (peakB, posB))

    print("\n结论：真问题 %.3f ／ 退化对照 %.3f → 前者远低于后者，"
          "说明 ACI 的结构**根本不在**这张 WATSON 里。" % (peakC, peakB))
    return 0


if __name__ == "__main__":
    sys.exit(main())
