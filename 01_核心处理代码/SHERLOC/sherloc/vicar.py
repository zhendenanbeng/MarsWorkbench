# -*- coding: utf-8 -*-
"""读 Mars 2020 imgops 的 ACI 图像（VICAR/ODL 格式）。

PDS 上 `mars2020_imgops` 数据集里的 ACI 图像是 `.IMG`，不是 PNG——
文件开头是一段 ASCII 的 ODL 标签，后面跟二进制栅格。标签里写着
LINES / LINE_SAMPLES / SAMPLE_BITS / BANDS，据此就能把栅格取出来。

PIL 不认这个格式，所以这里自己解。只用 numpy + PIL。
"""

from __future__ import annotations

import numpy as np
from PIL import Image


def _parse_label(text):
    """把 ODL 标签里 `KEY = VALUE` 的行抽成字典。忽略注释和多行 GROUP。"""
    out = {}
    for line in text.splitlines():
        s = line.strip()
        if not s or s.startswith("/*") or "=" not in s:
            continue
        key, val = s.split("=", 1)
        key = key.strip()
        if key.startswith("^"):          # 指针，单独保留
            out[key] = val.strip().strip('"')
            continue
        val = val.strip()
        if val.startswith('"'):
            val = val.strip('"')
        else:
            val = val.split("<")[0].split("(")[0].strip()
        try:
            out[key] = int(val)
        except ValueError:
            try:
                out[key] = float(val)
            except ValueError:
                out[key] = val
    return out


def read_vicar(path):
    """读一个 VICAR/ODL 图像，返回 (numpy 二维数组, 标签字典)。"""
    raw = open(path, "rb").read()

    # ODL 标签以 END 结尾；找到它，同时记住文本长度
    m = raw.find(b"END\r\n")
    cut = m + 5 if m >= 0 else raw.find(b"END\n") + 4
    if cut <= 0:
        raise ValueError("找不到 ODL 标签的 END：%s" % path)
    label = _parse_label(raw[:cut].decode("ascii", "replace"))

    lines = int(label["LINES"])
    samples = int(label["LINE_SAMPLES"])
    bits = int(label.get("SAMPLE_BITS", 8))
    bands = int(label.get("BANDS", 1))
    if bands != 1:
        raise ValueError("只支持单波段，这个文件是 %d 波段" % bands)

    # 字节序照标签的 SAMPLE_TYPE 来：MSB_INTEGER = 大端。ACI 的 16 位
    # 产品（FDR/EDR）实际是大端——用小端读会把高低字节颠倒（中位数从
    # 408 变 31745，图像完全坏掉，2026-09-18 修）。8 位没有字节序。
    st = label.get("SAMPLE_TYPE", "")
    if bits == 8:
        dtype = np.uint8
    elif "MSB" in st:
        dtype = np.dtype(">u" + str(bits // 8))
    elif "LSB" in st:
        dtype = np.dtype("<u" + str(bits // 8))
    else:
        dtype = {16: np.dtype(">u2"), 32: np.dtype(">u4")}.get(bits)
    if dtype is None or np.dtype(dtype).itemsize * 8 != bits:
        raise ValueError("不支持的位深 %d" % bits)
    need = lines * samples * np.dtype(dtype).itemsize
    if len(raw) < need:
        raise ValueError("文件比标签声明的还短")

    # 栅格在文件末尾（前面是 ODL 标签 + 原始 VICAR2 头）。
    # 这是 Mars 2020 imgops 的实际排布；用长度反推最稳，不依赖指针的正负号约定。
    buf = np.frombuffer(raw[-need:], dtype=dtype)
    img = buf.reshape((lines, samples))
    return img, label


def load_aci(path):
    """把 ACI 图像读成 PIL 灰度图。PNG/TIF/JPG 直接交给 PIL，.IMG 走 VICAR。"""
    if path.lower().endswith((".img", ".vic", ".vicar")):
        arr, label = read_vicar(path)
        if arr.dtype != np.uint8:
            # 16 位按最大值归一化到 8 位，只为显示
            hi = float(arr.max()) or 1.0
            arr = (arr.astype(np.float64) / hi * 255.0).astype(np.uint8)
        return Image.fromarray(arr, mode="L")
    return Image.open(path).convert("L")
