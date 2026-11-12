# -*- coding: utf-8 -*-
"""把光谱画成 PNG。用 PIL 手绘，不依赖 matplotlib。"""

from __future__ import annotations

import os

import numpy as np
from PIL import Image, ImageDraw, ImageFont

FONT_CANDIDATES = (r"C:\Windows\Fonts\msyh.ttc", r"C:\Windows\Fonts\simhei.ttf",
                   r"C:\Windows\Fonts\arial.ttf", "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf")
_FONTS = {}


def font(size):
    if size not in _FONTS:
        f = None
        for path in FONT_CANDIDATES:
            try:
                f = ImageFont.truetype(path, size)
                break
            except Exception:                       # noqa: BLE001
                continue
        _FONTS[size] = f or ImageFont.load_default()
    return _FONTS[size]


def _nice_ticks(lo, hi, n=6):
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return [lo], lo, hi
    raw = (hi - lo) / n
    mag = 10.0 ** np.floor(np.log10(raw))
    for m in (1, 2, 2.5, 5, 10):
        step = m * mag
        if step >= raw:
            break
    start = np.ceil(lo / step) * step
    ticks = np.arange(start, hi + step * 0.5, step)
    return [float(t) for t in ticks], float(lo), float(hi)


def render_png(path, series, x, xlabel, title, subtitle="", notes=None,
               size=(1440, 860), ylabel="强度 (计数)"):
    """series: [(标签, 数值数组, 颜色), ...]；x 与各数值数组等长。"""
    W, H = size
    img = Image.new("RGB", (W, H), "white")
    d = ImageDraw.Draw(img)
    L, R, T, B = 108, 40, 108, 136
    pw, ph = W - L - R, H - T - B

    d.text((36, 26), title, fill=(20, 22, 26), font=font(28))
    if subtitle:
        d.text((36, 64), subtitle, fill=(110, 116, 124), font=font(16))

    x = np.asarray(x, dtype=float)
    vals = [v for _, v, _ in series if v is not None and len(v)]
    if not vals or len(x) == 0:
        d.text((L, T + 40), "没有可绘制的数据", fill=(150, 60, 60), font=font(20))
        img.save(path)
        return path

    ylo = float(min(np.nanmin(v) for v in vals))
    yhi = float(max(np.nanmax(v) for v in vals))
    if yhi == ylo:
        yhi, ylo = ylo + 1, ylo - 1
    pad = (yhi - ylo) * 0.08
    ylo, yhi = ylo - pad, yhi + pad

    xlo, xhi = float(np.nanmin(x)), float(np.nanmax(x))
    yticks, _, _ = _nice_ticks(ylo, yhi, 7)
    xticks, _, _ = _nice_ticks(xlo, xhi, 8)

    fx = lambda v: L + (v - xlo) / (xhi - xlo) * pw
    fy = lambda v: T + (yhi - v) / (yhi - ylo) * ph

    for t in yticks:
        y = fy(t)
        d.line([L, y, L + pw, y], fill=(234, 237, 241), width=1)
        d.text((L - 12, y - 9), "%.4g" % t, fill=(120, 126, 134), font=font(15), anchor="ra")
    for t in xticks:
        xx = fx(t)
        d.line([xx, T, xx, T + ph], fill=(243, 245, 248), width=1)
        d.text((xx, T + ph + 12), "%.0f" % t, fill=(120, 126, 134), font=font(15), anchor="ma")

    d.rectangle([L, T, L + pw, T + ph], outline=(150, 156, 164), width=2)
    if ylo < 0 < yhi:
        d.line([L, fy(0), L + pw, fy(0)], fill=(190, 196, 204), width=2)

    for label, y, color in series:
        if y is None or not len(y):
            continue
        pts = [(fx(xx), fy(yy)) for xx, yy in zip(x, y) if np.isfinite(yy)]
        if len(pts) > 1:
            d.line(pts, fill=color, width=2, joint="curve")

    d.text((L + pw / 2, T + ph + 42), xlabel, fill=(30, 34, 40), font=font(17), anchor="ma")
    # 纵向轴标题独立排版，不能把整行文字居中在 x=30（会越过画布左边）。
    label_font = font(17)
    bbox = d.textbbox((0, 0), ylabel, font=label_font)
    label = Image.new("RGBA", (bbox[2] - bbox[0] + 4, bbox[3] - bbox[1] + 4))
    ImageDraw.Draw(label).text((2 - bbox[0], 2 - bbox[1]), ylabel,
                               fill=(30, 34, 40, 255), font=label_font)
    label = label.rotate(90, expand=True)
    img.paste(label, (12, round(T + ph / 2 - label.height / 2)), label)

    lx = L + 12
    for label, y, color in series:
        if y is None or not len(y):
            continue
        d.line([lx, T + 18, lx + 26, T + 18], fill=color, width=3)
        d.text((lx + 32, T + 10), label, fill=(40, 46, 54), font=font(16))
        lx += 40 + d.textlength(label, font=font(16)) + 26

    if notes:
        yy = H - 30
        for line in reversed(notes[-2:]):
            d.text((36, yy), line, fill=(120, 126, 134), font=font(14))
            yy -= 22

    if isinstance(path, (str, bytes, os.PathLike)):
        img.save(path)
    else:
        img.save(path, format="PNG")
    return path
