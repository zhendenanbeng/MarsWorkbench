#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
按 label 的字节偏移读取数据表。

传统写法要数表头行数、猜 skiprows，换个产品就崩。
label 里的 offset 字段直接给出每张表第一条数据记录的起始字节，
seek 过去读 records 行即可。光谱数、通道数也全部从 label 取。
"""

from __future__ import annotations

import io
import os
import re
from typing import Optional

import numpy as np

from .label import Label, decode_product_id, read_label, label_path_for


class Product:
    """一个数据文件 + 它的 label。"""

    def __init__(self, data_path: str, label: Optional[Label] = None):
        self.data_path = data_path
        self.label = label or read_label(label_path_for(data_path) or "")
        self.name = os.path.basename(data_path)
        self._pid_cache = decode_product_id(self.name)

    @property
    def pid(self):
        return self._pid_cache

    # -- 基本信息：一律优先取 label --------------------------------------
    @property
    def sol(self) -> str:
        return self.label.sol or self.pid.sol

    @property
    def sequence_id(self) -> str:
        return (self.label.sequence_id or self.pid.sequence_id).lower()

    @property
    def sclk(self) -> str:
        return self.label.sclk_start or self.pid.sclk

    @property
    def product_type(self) -> str:
        return self.pid.product_type

    @property
    def n_spectra(self) -> Optional[int]:
        return self.label.n_spectra

    @property
    def has_label(self) -> bool:
        return bool(self.label.path and not self.label.errors and self.label.product_id)

    @property
    def size_mb(self) -> float:
        return os.path.getsize(self.data_path) / 1048576.0

    # -- 读表 ------------------------------------------------------------
    def tables(self):
        return [t for t in self.label.tables if t.offset is not None]

    def find_table(self, pattern: str):
        for t in self.tables():
            if re.search(pattern, t.name, re.I):
                return t
        return None

    def read_table(self, table, max_rows: Optional[int] = None) -> np.ndarray:
        if table is None or table.offset is None:
            raise ValueError("该产品没有可用的表定义")
        n = table.records if max_rows is None else min(table.records or 0, max_rows)
        if not n:
            raise ValueError("label 未给出 records")
        with open(self.data_path, "rb") as f:
            f.seek(table.offset)
            raw = b"".join(f.readline() for _ in range(n))
        return np.loadtxt(io.BytesIO(raw), delimiter=",", ndmin=2)

    def wavelength_nm(self):
        t = self.find_table(r"WAVELENGTH")
        return None if t is None else self.read_table(t).ravel()

    def read_text_table(self, table, max_rows: Optional[int] = None):
        """读混合类型表（如 RLS 的坐标表），返回字符串二维列表。"""
        if table is None or table.offset is None:
            raise ValueError("该产品没有可用的表定义")
        n = table.records if max_rows is None else min(table.records or 0, max_rows)
        with open(self.data_path, "rb") as f:
            f.seek(table.offset)
            rows = []
            for _ in range(n):
                line = f.readline()
                if not line:
                    break
                rows.append(line.decode("utf-8", "replace").strip().split(","))
        return rows

    def shot_positions(self):
        """RLS：返回 [{number, image, samp, line}, ...]；没有则返回 None。"""
        t = self.find_table(r"SHOT POSITION|LASER SHOT")
        if t is None:
            return None
        out = []
        for r in self.read_text_table(t):
            if len(r) < 5:
                continue
            try:
                out.append({"number": int(r[0]), "spec": r[1], "image": r[2],
                            "samp": float(r[3]), "line": float(r[4])})
            except ValueError:
                continue
        return out or None

    def region(self, idx: int, max_rows: Optional[int] = None) -> np.ndarray:
        t = self.find_table(r"REGION[\s_]*%d\b" % idx)
        if t is None:
            spec = [x for x in self.tables() if "SPECTRA" in x.name.upper()
                    or "INTENSITY" in x.name.upper()]
            if len(spec) >= idx:
                t = spec[idx - 1]
        if t is None:
            raise ValueError("找不到 Region %d 的表" % idx)
        return self.read_table(t, max_rows=max_rows)

    def regions(self, max_rows: Optional[int] = None):
        out = []
        for i in (1, 2, 3):
            try:
                out.append(self.region(i, max_rows=max_rows))
            except Exception:                        # noqa: BLE001
                out.append(None)
        return out


SKIP_SUFFIXES = (".xml", ".txt", ".lbl", ".pdf", ".html", ".htm", ".jpg", ".png", ".md", ".py")
SKIP_DIRS = ("_qc", ".git", "__pycache__")


def collect_products(folder: str, recursive: bool = True, skip_suffixes=SKIP_SUFFIXES):
    """扫描文件夹（默认递归），返回所有数据文件；label 随 Product 一起加载。"""
    paths = []
    if recursive:
        for root, dirs, files in os.walk(folder):
            dirs[:] = [d for d in dirs if d not in SKIP_DIRS]
            for name in files:
                paths.append(os.path.join(root, name))
    else:
        paths = [os.path.join(folder, n) for n in os.listdir(folder)
                 if os.path.isfile(os.path.join(folder, n))]
    out = []
    for p in sorted(paths):
        name = os.path.basename(p)
        if os.path.splitext(name)[1].lower() in skip_suffixes:
            continue
        if name.startswith("~$") or name.startswith("."):
            continue
        out.append(Product(p))
    return out
