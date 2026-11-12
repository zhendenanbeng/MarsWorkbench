# -*- coding: utf-8 -*-
"""上色 v2：把"阴影/高光处的色相不可信"这件事做进流程里。

v1（colorize_pair.py）直接把 WATSON 的 a/b 搬到 ACI 网格上。现在知道：
  * WATSON 是**阳光**照明，画面 37% 是投影阴影；ACI 是**自己的环光**、几乎无阴影；
  * 按几何算，WATSON 的阴影恰好压在 ACI 视场的右半边；
  * 阴影和高光像素的色相不是"岩石的颜色"，是照明条件。
所以 v2 做三件事：
  1. 用 RLS 锚定的 standoff（108 mm，而不是几何交点 129 mm）；
  2. 把 WATSON 的**阴影区与镜面高光**掩掉，用周围的 a/b 内插补回来；
  3. 把 a/b 平滑到 ~0.4 mm 尺度——逐像素色相本来就没有可验证的依据，
     平滑后至少不敢宣称比这更细的颜色分布。
亮度仍然完全用 ACI 自己的，一个像素都不动。

用法：python colorize_v2.py [standoff_mm]
输出：out/v2_colorized.png, out/v2_compare.png, out/v2_mask.png
"""
from __future__ import annotations

import os
import sys

import numpy as np
import cv2

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from colorize_pair import (ACI_IMG, ACI_XML, OUT, WARP_FLAGS, WAT_IMG, WAT_XML,
                           homography, load_cam, read_bsq, sample_watson, to8)


def main():
    d = float(sys.argv[1]) if len(sys.argv) > 1 else 108.0
    wat = load_cam(WAT_XML)
    aci = load_cam(ACI_XML)
    wr = read_bsq(WAT_IMG, wat)
    ar = read_bsq(ACI_IMG, aci)
    b = [i - 1 for i in wat["rgb"]]
    wcol = np.stack([to8(wr[b[i]]) for i in range(3)], axis=-1)     # RGB
    wgray = cv2.cvtColor(wcol, cv2.COLOR_RGB2GRAY)
    aci8 = to8(ar[0])

    Hm, dst = homography(aci, wat, d)
    warped = sample_watson(wcol, Hm, (aci["ns"], aci["nl"]))
    wgray_w = sample_watson(wgray, Hm, (aci["ns"], aci["nl"]))
    valid = sample_watson(np.full((wat["nl"], wat["ns"]), 255, np.uint8), Hm,
                          (aci["ns"], aci["nl"])) > 127

    # --- 掩膜：在 WATSON 自己的网格上判阴影/高光，再投到 ACI 网格 ---
    loc = cv2.blur(wgray.astype(np.float32), (151, 151))
    shadow = loc < 0.55 * float(np.median(loc))                 # 投影阴影
    hi = wgray >= np.percentile(wgray[wgray > 0], 99.5)          # 镜面高光
    bad = (shadow | hi).astype(np.uint8) * 255
    bad_w = sample_watson(bad, Hm, (aci["ns"], aci["nl"]),
                          ) if False else cv2.warpPerspective(
        bad, Hm, (aci["ns"], aci["nl"]),
        flags=cv2.INTER_NEAREST | cv2.WARP_INVERSE_MAP)
    badm = bad_w > 127
    print("WATSON 自身：阴影 %.0f%%、高光 %.1f%%；投到 ACI 网格后被判不可信的占 %.0f%%"
          % (100 * shadow.mean(), 100 * hi.mean(), 100 * badm[valid].mean()))

    # --- a/b：先把不可信处用内插补上，再平滑到 0.4 mm ---
    lab = cv2.cvtColor(warped, cv2.COLOR_RGB2LAB).astype(np.float32)
    ab = lab[..., 1:3].copy()
    um_aci = 6.66
    for k in range(2):
        ch = ab[..., k]
        fill = ch.copy()
        fill[badm] = 0
        w = (~badm).astype(np.float32)
        num = cv2.GaussianBlur(fill * w, (0, 0), 8)
        den = cv2.GaussianBlur(w, (0, 0), 8)
        ch = np.where(badm, num / np.maximum(den, 1e-6), ch)
        ch[~valid] = ch[valid].mean() if valid.any() else 128
        ab[..., k] = ch
    s = 0.4 * (1000.0 / um_aci) / 2.355          # 0.4 mm → px
    ab = cv2.GaussianBlur(ab, (0, 0), s)
    print("a/b 已平滑到 ~0.4 mm（σ=%.1f px）" % s)

    out = cv2.cvtColor(aci8, cv2.COLOR_GRAY2RGB)
    lab_a = cv2.cvtColor(out, cv2.COLOR_RGB2LAB).astype(np.float32)
    lab_a[..., 1] = ab[..., 0]
    lab_a[..., 2] = ab[..., 1]
    col = cv2.cvtColor(np.clip(lab_a, 0, 255).astype(np.uint8), cv2.COLOR_LAB2RGB)
    col[~valid] = 0
    cv2.imwrite(os.path.join(OUT, "v2_colorized.png"), cv2.cvtColor(col, cv2.COLOR_RGB2BGR))
    m = np.zeros_like(aci8)
    m[badm & valid] = 255
    m[valid & ~badm] = 90
    cv2.imwrite(os.path.join(OUT, "v2_mask.png"), m)

    def panel(img, txt):
        bar = np.full((34, img.shape[1], 3), 255, np.uint8)
        cv2.putText(bar, txt, (8, 24), cv2.FONT_HERSHEY_SIMPLEX, 0.6, (0, 0, 0), 2,
                    cv2.LINE_AA)
        return np.vstack([bar, img])

    h = 720
    def rs(im):
        f = h / im.shape[0]
        return cv2.resize(im, (int(im.shape[1] * f), h))
    sides = [panel(rs(cv2.cvtColor(aci8, cv2.COLOR_GRAY2RGB)), "ACI gray (L)"),
             panel(rs(cv2.cvtColor(warped, cv2.COLOR_RGB2BGR)),
                   "WATSON colour warped (%d mm) — 阴影不可用" % d),
             panel(rs(cv2.cvtColor(col, cv2.COLOR_RGB2BGR)),
                   "v2: ab 去阴影+平滑到 0.4 mm, L=ACI")]
    C = np.hstack([sides[0], sides[2]])
    cv2.imwrite(os.path.join(OUT, "v2_compare.png"), cv2.cvtColor(C, cv2.COLOR_RGB2BGR))
    cv2.imwrite(os.path.join(OUT, "v2_compare3.png"),
                cv2.cvtColor(np.hstack(sides), cv2.COLOR_RGB2BGR))
    print("写到 out/v2_colorized.png, out/v2_mask.png, out/v2_compare.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
