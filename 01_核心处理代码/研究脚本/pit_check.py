# -*- coding: utf-8 -*-
"""ACI 那张图是不是**打激光之后**拍的（表面上有没有规则的点阵）？

线索：项目自己的数据说 sol 186 / srlc11420 = "Bellegarde survey"，RLS 产品时刻
0683478599，而 ACI 图 0683479674——晚了约 18 分钟。如果 ACI 是"打完激光的
context 图"，那它表面上会有 1296 个激光坑排成的规则点阵；而 WATSON 拍的是
**打之前**的干净面。那样两台相机看到的根本不是同一个表面，
怎么配准都不可能成——而且点阵的周期正好落在 0.06–0.5 mm 那个共享波段里。

做法：二维功率谱里找周期峰。RLS 的扫描点在 ACI 像素里的间距，
Survey 是 13.9 px、HDR_500 是 75.1 px（见 README 第五节）。
用法：python pit_check.py
"""
from __future__ import annotations

import os
import sys

import numpy as np
import cv2

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from colorize_pair import ACI_IMG, ACI_XML, WAT_IMG, WAT_XML, OUT, load_cam, read_bsq, to8


def spectrum_peaks(img, name):
    """径向 + 峰值：找出图上最强的周期成分（以像素为周期）。"""
    g = img.astype(np.float32)
    g = g - cv2.GaussianBlur(g, (0, 0), 40)          # 去大尺度照明梯度
    win = np.outer(np.hanning(g.shape[0]), np.hanning(g.shape[1])).astype(np.float32)
    F = np.fft.fftshift(np.abs(np.fft.fft2(g * win)))
    h, w = F.shape
    cy, cx = h // 2, w // 2
    F[cy - 3:cy + 4, cx - 3:cx + 4] = 0               # 去直流
    # 只看 6–120 px 周期的环带
    yy, xx = np.mgrid[0:h, 0:w]
    rr = np.sqrt((yy - cy) ** 2 + (xx - cx) ** 2)
    keep = (rr > min(h, w) / 120.0) & (rr < min(h, w) / 6.0)
    Fk = np.where(keep, F, 0)
    k = np.unravel_index(int(np.argmax(Fk)), Fk.shape)
    per = 1.0 / max(np.hypot(k[0] - cy, k[1] - cx) / min(h, w), 1e-9)
    print("%-8s 最强周期峰：周期 %.1f px（相对强度 %.2f），位置偏移 (%d,%d)"
          % (name, per, Fk[k] / max(F.mean(), 1e-9), k[1] - cx, k[0] - cy))
    # 也直接问：13.9 px 和 75.1 px 那两个频率上有没有东西
    for target, lab in ((13.9, "Survey 13.9 px"), (75.1, "HDR_500 75.1 px")):
        rad = min(h, w) / target
        ann = (np.abs(rr - rad) < 1.5) & keep
        base = np.median(F[keep])
        print("          %-18s 该环带峰值/中位 = %6.2f"
              % (lab, F[ann].max() / base if ann.any() else float("nan")))
    return F


def main():
    aci = load_cam(ACI_XML)
    wat = load_cam(WAT_XML)
    aci8 = to8(read_bsq(ACI_IMG, aci)[0])
    b = [i - 1 for i in wat["rgb"]]
    wg = cv2.cvtColor(np.stack([to8(read_bsq(WAT_IMG, wat)[b[i]]) for i in range(3)],
                               axis=-1), cv2.COLOR_RGB2GRAY)
    spectrum_peaks(aci8, "ACI")
    spectrum_peaks(wg, "WATSON")
    # 自相关：规则点阵会在自相关里出现次级峰
    for nm, im in (("ACI", aci8), ("WATSON", wg)):
        g = im[200:1000, 300:1300].astype(np.float32)
        g = g - cv2.GaussianBlur(g, (0, 0), 40)
        g = g - g.mean()
        F = np.fft.fft2(g)
        ac = np.fft.ifft2(F * np.conj(F)).real
        ac = np.fft.fftshift(ac)
        ac /= ac.max()
        h, w = ac.shape
        cy, cx = h // 2, w // 2
        prof = ac[cy, cx:cx + 120]
        dips = [i for i in range(4, 110) if prof[i] > prof[i - 3] and prof[i] > prof[i + 3]
                and prof[i] > 0.12]
        print("%-8s 自相关在径向出现次级峰的半径（px）：%s"
              % (nm, dips[:8] if dips else "无"))
    print("\n判读：如果 ACI 在 13.9 px（或 75.1 px）上有突出的周期峰、而 WATSON 没有，"
          "\n那 ACI 那张就是**打完激光之后**拍的，表面上已经布满了激光坑阵列——"
          "\n和 WATSON（打之前）不是同一个表面，配准不可能成。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
