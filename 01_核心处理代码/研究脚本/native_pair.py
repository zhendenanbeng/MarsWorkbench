# -*- coding: utf-8 -*-
"""在 WATSON 原生网格上做"同一块地方吗"的肉眼比对图。

两边都放在 **WATSON 的 21 µm/px 网格**上：WATSON 用自己的像素，
ACI 用几何单应投过来（并按 WATSON 的 PSF 补糊）。这样不存在
"把 WATSON 放大 2.83 倍"带来的假象，两张图是同一物理采样的。

出图 out/vc_pair_129.png：
  第 1 行  原始灰度（各自拉过对比度）      WATSON | ACI 投过来
  第 2 行  同一物理波段的结构图（0.06–0.5 mm） WATSON | ACI 投过来
  第 3 行  红=WATSON 绿=ACI 的结构叠加（黄色=两边都有结构）
用法：python native_pair.py [standoff_mm]
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


def stretch(g, lo=1.0, hi=99.0):
    a, b = np.percentile(g, (lo, hi))
    return np.clip((g - a) / max(b - a, 1e-9) * 255, 0, 255).astype(np.uint8)


def bandpass(img, lo, hi):
    g = img.astype(np.float32)
    return cv2.GaussianBlur(g, (0, 0), lo) - cv2.GaussianBlur(g, (0, 0), hi)


def zs(a, m=None):
    v = a[m] if m is not None else a
    s = float(v.std())
    return (a - float(v.mean())) / (s if s > 1e-9 else 1.0)


def bar(w, text):
    b = np.full((34, w, 3), 255, np.uint8)
    cv2.putText(b, text, (6, 24), cv2.FONT_HERSHEY_SIMPLEX, 0.55, (0, 0, 0), 2,
                cv2.LINE_AA)
    return b


def main():
    d = float(sys.argv[1]) if len(sys.argv) > 1 else 129.0
    wat = load_cam(WAT_XML)
    aci = load_cam(ACI_XML)
    wr = read_bsq(WAT_IMG, wat)
    ar = read_bsq(ACI_IMG, aci)
    b = [i - 1 for i in wat["rgb"]]
    wgray = cv2.cvtColor(np.stack([to8(wr[b[i]]) for i in range(3)], axis=-1),
                         cv2.COLOR_RGB2GRAY)
    aci8 = to8(ar[0])

    Hm, dst = homography(aci, wat, d)
    aci_w = cv2.warpPerspective(aci8.astype(np.float32), Hm,
                                (wat["ns"], wat["nl"]), flags=cv2.INTER_LINEAR)
    cov = aci_w > 0
    # ACI 的 PSF 到 WATSON 网格只剩 ~1.2 px，WATSON 自己是 ~3 px → 补糊对齐
    aci_w = cv2.GaussianBlur(aci_w, (0, 0), 2.75)

    # 物理尺度：WATSON 一像素多少 µm（用几何比 0.353 反推）
    wat_um = ACI_UM / 0.353
    wat_ppm = 1000.0 / wat_um
    sl = 0.06 * wat_ppm / 2.355
    sh = 0.50 * wat_ppm / 2.355
    print("WATSON %.2f µm/px；带通 σ = %.1f–%.1f px" % (wat_um, sl, sh))

    w_band = zs(bandpass(wgray, sl, sh))
    a_band = zs(bandpass(aci_w, sl, sh), cov)
    a_band[~cov] = 0

    pad = 18
    x0 = max(0, int(dst[:, 0].min()) - pad)
    x1 = min(wat["ns"], int(dst[:, 0].max()) + pad)
    y0 = max(0, int(dst[:, 1].min()) - pad)
    y1 = min(wat["nl"], int(dst[:, 1].max()) + pad)
    print("比对区：WATSON [%d..%d, %d..%d] = %d×%d px，即岩面上 %.1f×%.1f mm"
          % (x0, x1, y0, y1, x1 - x0, y1 - y0,
             (x1 - x0) * wat_um / 1000, (y1 - y0) * wat_um / 1000))

    wraw = cv2.cvtColor(stretch(wgray[y0:y1, x0:x1]), cv2.COLOR_GRAY2RGB)
    araw = aci_w[y0:y1, x0:x1]
    araw8 = cv2.cvtColor(stretch(araw), cv2.COLOR_GRAY2RGB)
    araw8[~(cov[y0:y1, x0:x1])] = 0

    def g3(a):
        return cv2.cvtColor(np.clip(128 + 60 * a, 0, 255).astype(np.uint8),
                            cv2.COLOR_GRAY2RGB)

    wst = g3(w_band[y0:y1, x0:x1])
    ast = g3(a_band[y0:y1, x0:x1])
    ast[~(cov[y0:y1, x0:x1])] = 0

    rg = np.full((y1 - y0, x1 - x0, 3), 128, np.uint8)
    rg[..., 0] = np.clip(128 + 70 * w_band[y0:y1, x0:x1], 0, 255)
    rg[..., 1] = np.clip(128 + 70 * a_band[y0:y1, x0:x1], 0, 255)
    rg[~(cov[y0:y1, x0:x1])] = 0

    W = x1 - x0
    rows = [
        np.hstack([bar(W, "1  WATSON raw (inverted gray)"), bar(W, "1  ACI raw, mapped to WATSON grid")]),
        np.hstack([wraw, araw8]),
        np.hstack([bar(W, "2  WATSON structure 0.06-0.5 mm"), bar(W, "2  ACI structure, same band")]),
        np.hstack([wst, ast]),
        np.hstack([bar(W, "3  red = WATSON structure, green = ACI structure (yellow = agree)"),
                   bar(W, "    standoff %.0f mm" % d)]),
        np.hstack([rg, rg * 0]),
    ]
    out = np.vstack(rows)
    cv2.imwrite(os.path.join(OUT, "vc_pair_%d.png" % int(d)),
                cv2.cvtColor(out, cv2.COLOR_RGB2BGR))
    big = cv2.resize(out, None, fx=1.7, fy=1.7, interpolation=cv2.INTER_NEAREST)
    cv2.imwrite(os.path.join(OUT, "vc_pair_%d_big.png" % int(d)),
                cv2.cvtColor(big, cv2.COLOR_RGB2BGR))
    print("写到 out/vc_pair_%d.png (%d×%d) 与 _big.png" % (int(d), out.shape[1],
                                                          out.shape[0]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
