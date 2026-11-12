#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""矿物分布图：把逐点矿物类别画到 RLS 真实坐标上。

两种模式：
  * 有底图——给 --aci-dir 指一个放着 ACI 灰度图的目录，就把点叠在图像上。
    RLS 里的 samp/line 本来就是 ACI 图像的像素坐标，不需要任何配准。
  * 无底图——不给 --aci-dir，就直接按坐标画在白底上，带比例尺和图例。
    没有 ACI 图像也能出图，这一点是有意设计的：图像下载还没打通，
    但空间分布（谁挨着谁、哪一类聚成一团）本身不需要底图就能看。

用法：
    python make_mineral_map.py <数据文件夹> [-o 输出目录] [--aci-dir 图片目录]
                               [--sol 257] [--seq srlc11360]

不给 --aci-dir 时，会自动去找 `<文件夹>\_aci`（工作台的「下载缺失的 ACI」
默认下到这里）以及项目里的 `images\aci`；都没有就按无底图模式出图。

依赖：pillow、numpy（和质控台一样，不需要额外装东西）
"""

from __future__ import annotations

import source_paths  # source roots; no-op for bundled imports

import argparse
import math
import os
import re
import sys

import numpy as np

from PIL import Image, ImageDraw, ImageFont

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from sherloc import minerals, qc                   # noqa: E402
from sherloc.organize import scan_folder           # noqa: E402
from sherloc.vicar import load_aci                 # noqa: E402

# ACI 是 10.1 µm/像素（Corpolongo et al., 2023 第 2.1 节）。
# 2026-09-20 独立定标验证：PIXL rxl 臂坐标(米) vs MCC 像素线性拟合
# （R²=0.9998）解出 MCC≈50.7 µm/px，除以 MCC↔ACI 配准尺度 4.87-4.97
# 得 ACI≈10.3 µm/px，与 10.1 在 3% 内吻合。曾被 RLS 反推质疑为 6.66
# （06_算法实验与依据/README.md 线索2），该推导错在把 HDR 扫描步长当成 500 µm
# （实测 ~758 µm）。10.1 维持不变。
ACI_UM_PER_PX = 10.1
BG_LIGHTEN = 1.25          # 灰度底图提亮一点，彩色圆点才跳得出来
MAP_TARGET_PX = 1200       # 裁剪后把画面放大到长边这么多像素


def crop_box(xs, ys, size, margin=0.12, min_pad_px=30):
    """按扫描点的包围盒算裁剪框（含边距），并夹到图像范围以内。

    为什么要裁：一张 ACI 是 16.6 × 12.1 mm，而一次扫描往往只有 5×5 mm（Survey）
    甚至 0.9×0.9 mm（1 mm 的 Detail）。不裁的话扫描区在图上只占中间一小块，
    为了让点不重叠只能把圆点缩得很小——效果还是不好。
    """
    W, H = size
    x0, x1 = min(xs), max(xs)
    y0, y1 = min(ys), max(ys)
    pad = max(margin * max(x1 - x0, y1 - y0), min_pad_px)
    cx0 = int(max(0, math.floor(x0 - pad)))
    cy0 = int(max(0, math.floor(y0 - pad)))
    cx1 = int(min(W, math.ceil(x1 + pad)))
    cy1 = int(min(H, math.ceil(y1 + pad)))
    # 只有一个点（或几个点重合）时给一个最小窗口，不然裁出来是 1 像素
    if cx1 - cx0 < 60:
        c = (cx0 + cx1) // 2
        cx0, cx1 = max(0, c - 30), min(W, c + 30)
    if cy1 - cy0 < 60:
        c = (cy0 + cy1) // 2
        cy0, cy1 = max(0, c - 30), min(H, c + 30)
    return cx0, cy0, cx1, cy1


def median_spacing(pts):
    """一批点的最近邻距离中位数（单位：像素）。用来决定圆点画多大。

    为什么要这个：圆的半径原来是写死的 15 px，只适合火山口底那种
    「100 个点铺满 7×7 mm」的稀疏扫描（点间距约 75 px）。
    后来的数据密度高得多——1296 点的 Survey 点间距只有约 14 px，
    1×1 mm 的 Detail 只有约 9.5 px——半径 15 的圆点会重叠 2–3 倍，
    整张图糊成一坨（2026-09-17 用户反馈）。所以改成按疏密自适应。
    """
    a = np.asarray(pts, dtype=float)
    if a.shape[0] < 2:
        return None
    d = np.sqrt((a[:, None, 0] - a[None, :, 0]) ** 2
                + (a[:, None, 1] - a[None, :, 1]) ** 2)
    np.fill_diagonal(d, np.inf)          # 自己不算自己的邻居
    return float(np.median(d.min(axis=1)))


def get_font(size):
    for name in (r"C:\Windows\Fonts\msyh.ttc", r"C:\Windows\Fonts\simhei.ttf",
                 "msyh.ttc", "simhei.ttf", "arial.ttf"):
        try:
            return ImageFont.truetype(name, size)
        except Exception:                          # noqa: BLE001
            continue
    return ImageFont.load_default()


def hex_rgb(h):
    h = h.lstrip("#")
    return tuple(int(h[i:i + 2], 16) for i in (0, 2, 4))


def index_aci(aci_dirs):
    """给 ACI 目录建两个索引，用来把 RLS 里的图像名对上实际图片文件。

    坑：RLS 里写的是 `..._0000LMJ03.IMG`，而手上下到的图可能是
    `..._0000LMJ01.png`——尾部的两位数字是产品**版本号**，会不一样。
    所以主索引要先把版本号剥掉；再留一个按 SCLK 的兜底索引。

    aci_dirs 可以是一个目录，也可以是多个（图片可能分散在几处，比如
    images\\aci 和 images\\sol_00181）。后面的目录不会覆盖前面的。
    """
    if isinstance(aci_dirs, str):
        aci_dirs = [aci_dirs]
    by_stem, by_sclk = {}, {}
    for aci_dir in aci_dirs or []:
        if not aci_dir or not os.path.isdir(aci_dir):
            continue
        for root, _dirs, files in os.walk(aci_dir):
            for f in files:
                # .img 是 PDS 上 ACI 的原始格式（VICAR），得自己解，见 sherloc/vicar.py
                if not f.lower().endswith((".png", ".tif", ".tiff", ".jpg", ".img")):
                    continue
                stem = os.path.splitext(f)[0]
                full = os.path.join(root, f)
                by_stem.setdefault(re.sub(r"\d+$", "", stem).lower(), full)
                parts = stem.split("_")
                if len(parts) > 2:
                    by_sclk.setdefault(parts[2], []).append(full)
    return by_stem, by_sclk


def guess_aci_dirs(folder):
    r"""没给 --aci-dir 时，猜图片可能放在哪几个目录。返回列表（可能为空）。

    按这个顺序找：
      1. `<数据文件夹>\_aci`        —— 工作台「下载缺失的 ACI」的默认位置
      2. `<数据文件夹>\images\aci`  —— 本项目放 12 张官方图的地方
      3. 往上一级 / 上两级的 images\aci（数据放在 raw\sol_xxxx 这种子目录时用得上）
      4. `<数据文件夹>\aci`
      5. 各个 sol 子目录下的 `_aci`（2026-09-17 加）——工作台的当前文件夹是
         单个 sol 时，图就下在 `processed\sol_xxxx\_aci` 里；跑顶层文件夹出图时
         原来找不到它们，那张图就退回白底版了。

    第 2 条是必须的：工作台默认指向的就是 02_练习数据与验证基准 这个顶层目录，
    而图片在它下面的 images\aci 里。曾经漏了这条，于是顶层文件夹出图时
    底图永远找不到，只能出白底版。
    """
    parent = os.path.dirname(folder)
    cands = [os.path.join(folder, "_aci"),
             os.path.join(folder, "images", "aci"),
             os.path.join(parent, "images", "aci"),
             os.path.join(os.path.dirname(parent), "images", "aci"),
             os.path.join(folder, "aci"),
             os.path.join(folder, "images")]        # 图可能按 sol 分子目录放
    # 第 5 条：把各级子目录里的 _aci 也收进来。只往下走 3 层，
    # 不在几十 GB 的数据树上乱翻。
    root_depth = folder.rstrip("\\/").count(os.sep)
    for dirpath, dirs, _files in os.walk(folder):
        if dirpath.rstrip("\\/").count(os.sep) - root_depth >= 3:
            dirs[:] = []                              # 不再往下
            continue
        for d in dirs:
            if d.lower() in ("_aci", "aci"):
                cands.append(os.path.join(dirpath, d))
    out = []
    for c in cands:
        if not os.path.isdir(c) or c in out:
            continue
        has = False
        for _root, _dirs, files in os.walk(c):
            if any(f.lower().endswith((".img", ".png", ".tif", ".jpg")) for f in files):
                has = True
                break
        if has:
            out.append(c)
    return out


def find_aci(aci_index, image_name):
    """按 RLS 的图像名找对应图片。先按剥离版本号的主干，再按 SCLK 兜底。"""
    by_stem, by_sclk = aci_index
    stem = os.path.splitext(image_name)[0]
    hit = by_stem.get(re.sub(r"\d+$", "", stem).lower())
    if hit:
        return hit
    parts = stem.split("_")
    if len(parts) > 2:
        cands = by_sclk.get(parts[2], [])
        if len(cands) == 1:                 # 不唯一就别猜
            return cands[0]
    return None


def point_classes(det):
    """点号 -> (主类别 key, 字母串)。主类别取第一个（内部已按 SNR 排序）。"""
    out = {}
    for d in det:
        i = int(d["point"])
        if d["classes"]:
            key = d["classes"][0]["key"]
            letters = "+".join(minerals.SHORT.get(b["key"], "?") for b in d["classes"])
        elif d["n_bands"]:
            key, letters = None, "·"
        else:
            key, letters = None, ""
        out[i] = (key, letters, d)
    return out


def layout(pts, canvas_w, canvas_h, pad_ratio=0.10):
    """无底图模式：把 (samp, line) 线性映射到画布内部，保持长宽比。"""
    xs = [p["samp"] for p in pts]
    ys = [p["line"] for p in pts]
    x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
    dx, dy = max(x1 - x0, 1.0), max(y1 - y0, 1.0)
    inner_w = canvas_w * (1 - 2 * pad_ratio)
    inner_h = canvas_h - 150          # 底下留给图例和比例尺
    scale = min(inner_w / dx, inner_h / dy)
    ox = (canvas_w - dx * scale) / 2
    oy = 90 + (inner_h - dy * scale) / 2
    return lambda p: (ox + (p["samp"] - x0) * scale, oy + (p["line"] - y0) * scale), scale


def draw_scale_bar(draw, x, y, px_per_mm, font, color=(40, 40, 40)):
    """比例尺：优先画 1 mm，太短就退到 5 mm。"""
    for mm in (1.0, 2.0, 5.0):
        length = mm * px_per_mm
        if 60 <= length <= 260:
            draw.line([(x, y), (x + length, y)], fill=color, width=3)
            for xx in (x, x + length):
                draw.line([(xx, y - 7), (xx, y + 7)], fill=color, width=3)
            draw.text((x, y + 12), "%g mm" % mm, fill=color, font=font)
            return True
    return False


def draw_legend(draw, counts, x, y, font, font_small, line_h=22):
    """图例：出现过的类别 + 点数。"""
    draw.text((x, y), "矿物类别", fill=(30, 30, 30), font=font)
    y += line_h + 4
    order = [k for k in minerals.LABELS if k in counts]
    order += [k for k in counts if k not in order]
    for key in order:
        n = counts[key]
        if not n:
            continue
        color = hex_rgb(minerals.CLASS_COLORS.get(key, "#6b7683"))
        draw.ellipse([x, y + 4, x + 14, y + 18], fill=color,
                     outline=(255, 255, 255))
        label = minerals.LABELS.get(key, key)
        draw.text((x + 22, y), "%s  %d 点" % (label, n), fill=(40, 40, 40),
                  font=font_small)
        y += line_h
    return y


def colorized_base(aci_path, verbose=True):
    """给一张 ACI 灰度图找/生成上色底图。

    返回 (PIL RGB 图, 说明文字)；失败返回 (None, 原因文字)——原因要
    一路带到界面日志，"无配对或判据未过"这种笼统话帮不了用户
    （sol 865 的反馈实测踩过）。

    优先用缓存（<aci 目录>/colorized/<同名>.png）；没有且装了 opencv 就
    现场配：watson.find_pair 找共视的彩色图 → 下载 → colorize.colorize。
    任何一步不成（无配对、NO-GO、下载失败）都安静地回退灰度——上色是
    锦上添花，不是必需品。
    """
    try:
        from sherloc import colorize as _colorize
        from sherloc import watson as _watson
    except ImportError:
        return None, "sherloc 模块导入失败"
    if not _colorize.HAS_CV2:
        return None, "未安装 opencv"

    import numpy as np

    stem = os.path.splitext(os.path.basename(aci_path))[0]
    aci_dir = os.path.dirname(aci_path)
    cache = os.path.join(aci_dir, "colorized", stem + ".png")
    if os.path.exists(cache):
        try:
            return Image.open(cache).convert("RGB"), "上色底图（缓存）"
        except Exception:                            # noqa: BLE001
            pass

    logs = []
    pair = _watson.find_pair(stem, folder=aci_dir,
                             log=logs.append)
    if not pair or not pair.get("url"):
        # find_pair 的日志里有具体原因（无清单 / 无同定位候选 /
        # 有候选但都不共视），取最后一条作结
        why = logs[-1].strip() if logs else "找不到共视的 WATSON 彩色图"
        if verbose:
            print("    上色：%s，用灰度底图" % why)
        return None, why

    wat_path = os.path.join(aci_dir, "colorized", pair["stem"] + ".IMG")
    if not os.path.exists(wat_path):
        try:
            from sherloc.aci import download as _dl
            os.makedirs(os.path.dirname(wat_path), exist_ok=True)
            _dl(pair["url"], wat_path)
            if verbose:
                print("    上色：已下载 WATSON 彩图 %s" % pair["stem"][-24:])
        except Exception as exc:                     # noqa: BLE001
            if verbose:
                print("    上色：WATSON 下载失败（%s），用灰度底图" % exc)
            return None, "WATSON 下载失败（%s）" % exc

    try:
        aci_cam = _cam_from_xml_text(
            _watson._aci_cam_text(stem, aci_dir))
        wat_xml_text = _watson.get_cam_xml(pair["stem"], aci_dir)
        wat_cam = _cam_from_xml_text(wat_xml_text)
        if not aci_cam or not wat_cam:
            raise ValueError("相机模型不可用")
        from sherloc.vicar import read_vicar
        wat_rgb, wlab = _read_rgb(wat_path, wat_cam)
        aci_gray = np.asarray(load_aci(aci_path).convert("L"))
        res = _colorize.colorize(aci_gray, wat_rgb, aci_cam, wat_cam,
                                 pair.get("standoff") or 120.0,
                                 log=logs.append)
    except Exception as exc:                         # noqa: BLE001
        if verbose:
            print("    上色：失败（%s: %s），用灰度底图"
                  % (type(exc).__name__, exc))
        return None, "配准失败（%s: %s）" % (type(exc).__name__, exc)
    if not res.get("ok"):
        why = res.get("reason") or "配准判据未过"
        if verbose:
            print("    上色：%s，用灰度底图" % why)
        return None, why

    os.makedirs(os.path.dirname(cache), exist_ok=True)
    Image.fromarray(res["color"]).save(cache)
    if verbose:
        print("    上色：成功（%d 锚点、残差 %.0fpx、受光 %.0f%%），已缓存"
              % (res["n_anchors"], res["resid"], 100 * res["lit"]))
    return (Image.fromarray(res["color"]),
            "上色底图（L=ACI 纹理，颜色=WATSON；%d 锚点）" % res["n_anchors"])


def _cam_from_xml_text(text):
    """xml 文本 -> colorize 用的相机 dict；解析不了返回 None。"""
    import re as _re
    from sherloc import watson as _watson
    if not text:
        return None
    ns = nl = None
    m = _re.search(r"<img:samples>\s*(\d+)", text)
    n2 = _re.search(r"<img:lines>\s*(\d+)", text)
    if m and n2:
        ns, nl = int(m.group(1)), int(n2.group(1))
    return _watson.load_cam_xml(text, ns, nl)


def _read_rgb(path, cam):
    """读 WATSON 的 3 波段 IMG 成 RGB uint8（大端 16 位 BSQ）。"""
    import numpy as _np
    from sherloc.vicar import _parse_label
    raw = open(path, "rb").read()
    m = raw.find(b"END\r\n")
    if m < 0:
        m = raw.find(b"END\n")
    lab = _parse_label(raw[:m].decode("ascii", "replace"))
    lines, samps, bands = lab["LINES"], lab["LINE_SAMPLES"], lab["BANDS"]
    if bands != 3:
        raise ValueError("WATSON 产品不是 3 波段（%d）" % bands)
    rec = int(lab["RECORD_BYTES"])
    start = (int(lab["^IMAGE"]) - 1) * rec
    dt = _np.dtype(">i2") if "MSB" in lab.get("SAMPLE_TYPE", "") \
        else _np.dtype("<i2")
    arr = _np.frombuffer(raw, dtype=dt, offset=start,
                         count=lines * samps * bands).reshape(
        bands, lines, samps).astype(_np.float64)
    # 逐通道拉伸（照 spike 的 to8）：三通道合一起拉会把颜色调灰
    chan = []
    for i in range(bands):
        lo, hi = _np.percentile(arr[i], (0.5, 99.5))
        chan.append(_np.clip((arr[i] - lo) / max(hi - lo, 1e-9) * 255.0, 0, 255))
    rgb = _np.stack(chan, -1).astype(_np.uint8)
    # 通道顺序照标签的 RGB 通道声明
    order = [cam["rgb"][i] - 1 for i in range(3)] if cam.get("rgb") else [0, 1, 2]
    return rgb[..., order], lab


def render(ob, det, pts, aci_index, out_path, verbose=True):
    """把一次观测里「同一张 ACI 图像」的点画成一张图。"""
    by_image = {}
    for p in pts:
        by_image.setdefault(p["image"], []).append(p)
    pcls = point_classes(det)

    made = []
    all_counts = {}
    for image_name, group in sorted(by_image.items()):
        stem = os.path.splitext(image_name)[0]
        aci_path = find_aci(aci_index, image_name)
        label = stem.split("_")[2] if len(stem.split("_")) > 2 else stem[:24]

        if aci_path:
            # 上色底图优先（缓存过就直接用），不成回退灰度——降级是常态路径
            crop_note = ""        # 两条路径都要有初值，下面统一追加
            col_img, col_note = colorized_base(aci_path, verbose=verbose)
            if col_img is not None:
                base = col_img.point(lambda v: min(255, int(v * BG_LIGHTEN)))
            else:
                base = load_aci(aci_path).convert("RGB")
                base = base.point(lambda v: min(255, int(v * BG_LIGHTEN)))
            # 裁到这次扫描的范围（2026-09-17 改）。
            #
            # 原来整张 ACI 都画出来：那是 16.6 × 12.1 mm，而一次扫描往往只有
            # 5×5 mm（Survey）、甚至 0.9×0.9 mm（1 mm 间距的 Detail）。扫描区
            # 在图上只占中间一小块，为了不重叠只能把圆点缩到 7 px，看着还是不行。
            # 正确做法是裁底图、让画面去适配扫描区，圆点该多大就多大。
            xs = [p["samp"] for p in group if p.get("samp") is not None]
            ys = [p["line"] for p in group if p.get("line") is not None]
            if len(xs) >= 2:
                x0, y0, x1, y1 = crop_box(xs, ys, base.size)
                base = base.crop((x0, y0, x1, y1))
                cw, ch = base.size
                sc = max(1.0, min(10.0, MAP_TARGET_PX / float(max(cw, ch))))
                sc = min(12.0, max(sc, 560.0 / ch))     # 太扁就再放大，留出图例的位置
                if sc != 1.0:
                    base = base.resize((max(1, int(cw * sc)), max(1, int(ch * sc))),
                                       Image.LANCZOS)
                crop_note = crop_note + "（已裁到扫描范围 %.1f × %.1f mm，放大 %.1f 倍）" % (
                    (x1 - x0) * ACI_UM_PER_PX / 1000.0,
                    (y1 - y0) * ACI_UM_PER_PX / 1000.0, sc)
            else:
                x0 = y0 = 0
                sc = 1.0
                crop_note = "（只有不到 2 个点，未裁剪）"
            img = base
            W, H = img.size
            to_px = lambda p: ((p["samp"] - x0) * sc,                # noqa: E731
                               (p["line"] - y0) * sc)
            px_per_mm = 1000.0 / ACI_UM_PER_PX * sc
            mode = "底图：%s%s" % (os.path.basename(aci_path), crop_note)
        else:
            W, H = 1100, 900
            img = Image.new("RGB", (W, H), (250, 250, 248))
            to_px, scale = layout(group, W, H)
            px_per_mm = scale * 1000.0 / ACI_UM_PER_PX
            mode = "无底图（坐标按 RLS 相对位置绘制）"

        draw = ImageDraw.Draw(img)
        font_big = get_font(26)
        font = get_font(19)
        font_small = get_font(16)
        font_tiny = get_font(13)

        # 圆点半径跟着点的疏密走（2026-09-17 修）。
        # 原来的 15/16 px 是写死的，只适合火山口底那种「100 点铺满 7×7 mm」
        # 的稀疏扫描（点间距约 75 px）。后来的数据密度高得多：1296 点的 Survey
        # 点间距只有约 14 px、1×1 mm 的 Detail 约 9.5 px——半径 15 的圆点
        # 会重叠 2–3 倍，整张图糊成一坨。现在取最近邻间距中位数的 40% 作半径，
        # 圆点刚好挨着不叠；稀疏时封顶 15/16 px，老图观感不变。
        visible = [to_px(p) for p in group if int(p["number"]) in pcls]
        visible = [(x, y) for x, y in visible if 0 <= x < W and 0 <= y < H]
        sp = median_spacing(visible)
        cap = 20.0 if aci_path else 16.0
        # 半径要减掉描边宽度：圆点外面还有一圈白描边，不减的话相邻同色圆点
        # 会连成一片。裁底图放大之后点间距变大了，圆点也能跟着大起来。
        radius = max(2.5, min(cap, 0.40 * sp - 2.0)) if sp else cap
        # 点小的时候描边也得细，不然描边比圆还粗
        ow = 2 if radius >= 9.0 else 1
        # 点太密时字母根本读不出来（字母串比圆点宽得多），只画色点、靠图例认类别
        show_letters = radius >= 9.0
        # 字母大小跟着圆点走，圆点大了字也大一点
        font_letters = get_font(int(max(13, min(18, round(radius * 0.95)))))
        counts = {}
        for p in group:
            i = int(p["number"])
            if i not in pcls:
                continue
            key, letters, d = pcls[i]
            x, y = to_px(p)
            if not (0 <= x < W and 0 <= y < H):
                continue
            if key:
                counts[key] = counts.get(key, 0) + 1
                color = hex_rgb(minerals.CLASS_COLORS.get(key, "#6b7683"))
                faint = d["snr"] is not None and d["snr"] < 3.0
                draw.ellipse([x - radius, y - radius, x + radius, y + radius],
                             fill=color, outline=(255, 255, 255), width=ow)
                if faint:
                    # 未达 3σ 的加一圈虚线感的浅环
                    ring = max(2.0, radius * 0.25)
                    draw.ellipse([x - radius - ring, y - radius - ring,
                                  x + radius + ring, y + radius + ring],
                                 outline=(255, 255, 255), width=1)
                if letters and show_letters:
                    tw = draw.textlength(letters, font=font_letters)
                    # 字母用黑色（2026-09-17 改）：底色是无底图模式的近白 (250,250,248)
                    # 或提亮过的 ACI 灰度图，而字母串比这个半径 15-16px 的圆宽得多，
                    # 大部分笔画落在圆**外面**的浅色底上——白字在那里等于没写。
                    # 白色描边是为了压在深色圆上时也读得清。
                    draw.text((x - tw / 2, y - radius * 0.45), letters,
                              fill=(17, 17, 17), font=font_letters,
                              stroke_width=2, stroke_fill=(255, 255, 255))
            else:
                draw.ellipse([x - 6, y - 6, x + 6, y + 6],
                             outline=(150, 150, 150), width=2)

        # 标题
        img_hdr = Image.new("RGB", (W, 78), (255, 255, 255))
        img.paste(img_hdr, (0, 0))
        draw.text((16, 12), "sol %s · %s · %s" % (ob.sol, ob.sequence_id, label),
                  fill=(20, 20, 20), font=font_big)
        draw.text((16, 46), "%s ｜ %d 点 ｜ %s" % (mode, len(group), ob.sclk),
                  fill=(110, 110, 110), font=font_small)

        # 图例（放在底部一条白带上，避免压到图）
        band_h = 150
        band = Image.new("RGB", (W, band_h), (255, 255, 255))
        img.paste(band, (0, H - band_h))
        draw = ImageDraw.Draw(img)
        n_lab = sum(v for v in counts.values())
        draw.text((16, H - band_h + 10),
                  "判出类别 %d / %d 点　空心小圈＝有波段但未归类或没波段"
                  % (n_lab, len(group)),
                  fill=(90, 90, 90), font=font_small)
        draw_legend(draw, counts, 16, H - band_h + 40, font, font_small)
        draw_scale_bar(draw, W - 260, H - band_h + 40, px_per_mm, font_small)

        img.save(out_path if len(by_image) == 1
                 else out_path.replace(".png", "_%s.png" % label))
        made.append(out_path if len(by_image) == 1
                    else out_path.replace(".png", "_%s.png" % label))
        if verbose:
            for k, v in counts.items():
                all_counts[k] = all_counts.get(k, 0) + v
            detail = "、".join("%s %d" % (minerals.LABELS.get(k, k), v)
                               for k, v in sorted(counts.items(),
                                                  key=lambda kv: -kv[1]))
            print("    %s：%d 点，判出类别 %d（%s）"
                  % (label, len(group), n_lab, detail or "无"))
    return made, all_counts


def run(folder, out_dir=None, aci_dir=None, sol=None, seq=None, verbose=True):
    """给一个数据文件夹出一套矿物分布图。返回一个结果字典。

    界面上的「生成矿物分布图」按钮调的就是这个函数，所以它不打印、
    也不退出——把信息装在返回值里交给调用方。
    """
    folder = os.path.abspath(folder)
    out_dir = out_dir or os.path.join(folder, "_maps")
    os.makedirs(out_dir, exist_ok=True)
    aci_dirs = [aci_dir] if aci_dir else guess_aci_dirs(folder)
    aci_index = index_aci(aci_dirs)
    if verbose:
        if aci_dirs:
            print("ACI 目录 %s，索引到 %d 张图"
                  % ("；".join(aci_dirs), len(aci_index[0])))
        else:
            print("没找到 ACI 图像目录，按无底图模式出图"
                  "（可以在工作台里点「下载数据…→ 下载缺失的 ACI」补齐）")

    _products, observations = scan_folder(folder)
    n, made, skipped, errors = 0, [], 0, []
    for ob in observations:
        if ob.spectra is None:
            continue
        if sol and ob.sol != sol:
            continue
        if seq and ob.sequence_id != seq:
            continue
        wl, r1, r2, r3 = ob.read_regions()
        if wl is None or r1 is None:
            continue
        res = qc.analyse_spectra(wl, r1, r2, r3)
        det = res.get("point_detail") or []
        pts = ob.read_shot_positions()
        if not det or not pts:
            skipped += 1
            if verbose:
                print("  跳过 sol %s %s（没有逐点结果或没有 RLS 坐标）"
                      % (ob.sol, ob.sequence_id))
            continue
        # 文件名必须带上 SCLK：一个 sol 上同一个序列号可能有好几次扫描
        # （sol 162 的 srlc10500 就有 4 次），只用 sol+序列号会互相覆盖。
        base = "map_sol%s_%s_%s.png" % (
            ob.sol, ob.sequence_id,
            "_".join(s.replace(".", "_") for s in ob.all_sclks))
        if verbose:
            print("  sol %s %s：%d 点" % (ob.sol, ob.sequence_id, len(pts)))
        # 坏一个观测只废一张图，不中止整批——和 run_qc.analyse_folder 同口径
        try:
            got, _counts = render(ob, det, pts, aci_index,
                                  os.path.join(out_dir, base), verbose=verbose)
        except Exception as exc:                        # noqa: BLE001
            errors.append("sol %s %s：%s: %s"
                          % (ob.sol, ob.sequence_id,
                             type(exc).__name__, exc))
            if verbose:
                print("  出图失败 sol %s %s：%s: %s"
                      % (ob.sol, ob.sequence_id,
                         type(exc).__name__, exc))
            continue
        made.extend(got)
        n += 1
    if verbose:
        print("完成：%d 次观测 -> %s" % (n, out_dir))
        if errors:
            print("其中 %d 次观测出图失败（其余不受影响）" % len(errors))
    return {"count": n, "files": len(made), "out": out_dir,
            "aci": "；".join(aci_dirs) if aci_dirs else None,
            "maps": made, "skipped": skipped, "errors": errors}


def main(argv=None):
    ap = argparse.ArgumentParser(description="画矿物分布图")
    ap.add_argument("folder", help="含 PDS4 label 的数据文件夹")
    ap.add_argument("-o", "--out", default=None, help="输出目录")
    ap.add_argument("--aci-dir", default=None, help="ACI 灰度图所在目录（可选）")
    ap.add_argument("--sol", default=None, help="只处理这个 sol")
    ap.add_argument("--seq", default=None, help="只处理这个序列号")
    args = ap.parse_args(argv)
    res = run(args.folder, out_dir=args.out, aci_dir=args.aci_dir,
              sol=args.sol, seq=args.seq)
    return 0 if res["count"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
