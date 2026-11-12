#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PDS4 label 读取 —— SHERLOC 数据产品的权威信息来源
====================================================

为什么以 label 为准：

* ``<records>``   直接给出每张表的光谱条数，不用去数 CSV 行数
* ``<offset>``    直接给出每张表的字节偏移，不用猜 skiprows
* ``<name>``      表名里写着做过什么处理（例如 "LASER-NORMALIZED SPECTRA REGION 1"）
* ``<mars2020:sol_number>`` / ``spacecraft_clock_start`` / 时间 / 版本号都是官方值

例外：**光谱产品的处理标记（W/C/G/N）在 label 里没有**，只在文件名里。
PDS4 把文件名作为 ``<alternate_id>``（VICAR PRODUCT_ID）记录在案，所以解析
文件名是有官方依据的，但仅限于这一处，其余一律走 label。
"""

from __future__ import annotations

import os
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from typing import Optional


def _strip_ns(tag: str) -> str:
    """去掉 XML 命名空间前缀，方便按名字取值。"""
    return tag.rsplit("}", 1)[-1] if "}" in tag else tag


def _text(root, name, default=None):
    for el in root.iter():
        if _strip_ns(el.tag) == name and el.text and el.text.strip():
            return el.text.strip()
    return default


def _all(root, name):
    return [el for el in root.iter() if _strip_ns(el.tag) == name]


# --------------------------------------------------------------------------
# 文件名解码（附录 A）—— 处理标记的唯一来源
# --------------------------------------------------------------------------

PRODUCT_TYPES = {
    "ERA": "科学光谱·主动帧（原始）",
    "ERB": "科学光谱·暗背景帧（原始）",
    "ECA": "定标光谱·主动帧（原始）",
    "ECB": "定标光谱·暗背景帧（原始）",
    "ESH": "仪器健康状况（科学模式）",
    "ECH": "仪器健康状况（定标模式）",
    "ESP": "扫描镜位置（毫米）",
    "EXH": "扫描镜位置（像素）",
    "EPA": "激光光电二极管阵列",
    "ERP": "车载 process_data 光谱",
    "RRS": "科学光谱·已扣暗背景",
    "RRA": "科学光谱·处理后主动帧",
    "RRB": "科学光谱·处理后暗背景帧",
    "RCS": "定标光谱·已扣暗背景",
    "RCA": "定标光谱·处理后主动帧",
    "RCB": "定标光谱·处理后暗背景帧",
    "RCC": "内部 AlGaN 定标拟合结果",
    "RLS": "激光点在 ACI 图像中的坐标",
    "RLI": "激光强度图",
    "RLP": "扫描镜位置表（方位/俯仰/目标空间）",
    "RMO": "激光点位置与光谱强度汇总",
    "RAD": "辐射校正图像",
    "EDR": "常规图像",
    "ECM": "companded 图像",
    "ECZ": "车载 Z-stack 最佳焦图",
    "EZS": "地面 Z-stack 最佳焦图",
    "EDM": "对焦合并深度图",
}

RAMAN_TYPES = {"ERA", "ERB", "RRS", "RRA", "RRB", "ERP", "ECA", "ECB", "RCS", "RCA", "RCB"}


@dataclass
class ProductId:
    """文件名 55 位的逐位解读。"""

    raw: str
    instrument: str = ""
    instrument_name: str = ""
    configuration: str = ""
    sol: str = ""
    sclk: str = ""
    product_type: str = ""
    product_name: str = ""
    site: str = ""
    drive: str = ""
    sequence_id: str = ""
    processing_raw: str = ""
    flags: dict = field(default_factory=dict)
    provider: str = ""
    version: str = ""
    extension: str = ""
    is_image: bool = False

    @property
    def flag_text(self) -> str:
        f = self.flags
        if self.is_image:
            return "图像产品，处理标记不适用"
        out = []
        out.append("波长校正" if f.get("wavelength") else "未做波长校正")
        out.append("去宇宙射线" if f.get("cosmic") == "C" else
                   ("跳过去宇宙射线" if f.get("cosmic") == "Z" else "未做去宇宙射线"))
        out.append("增益校正" if f.get("gain") == "G" else
                   ("车载 process_data" if f.get("gain") == "P" else "未做增益校正"))
        out.append("激光归一化" if f.get("normalize") == "N" else
                   ("跳过归一化" if f.get("normalize") == "Z" else "未做激光归一化"))
        return " · ".join(out)


def decode_product_id(name: str) -> ProductId:
    """按 PDS4 附录 A 解码文件名（去掉扩展名后的 54 位主体）。"""
    stem = name.strip()
    ext = ""
    if "." in stem:
        stem, ext = stem.rsplit(".", 1)
    pid = ProductId(raw=name, extension=ext.upper())
    s = stem.upper()
    if len(s) < 26:
        return pid

    pid.instrument = s[0:2]
    pid.instrument_name = {"SS": "SHERLOC 光谱", "SC": "SHERLOC ACI",
                           "SI": "SHERLOC WATSON", "SE": "SHERLOC 工程"}.get(pid.instrument, pid.instrument)
    pid.is_image = pid.instrument in ("SC", "SI")
    pid.configuration = s[2]
    pid.sol = s[4:8]
    pid.sclk = "%s.%s" % (s[9:19], s[20:23])
    pid.product_type = s[23:26].rstrip("_")
    pid.product_name = PRODUCT_TYPES.get(pid.product_type, pid.product_type or "未知")
    pid.site = s[28:31]
    pid.drive = s[31:35]
    pid.sequence_id = s[35:44]
    pid.processing_raw = s[44:51]
    pid.provider = s[51] if len(s) > 51 else ""
    pid.version = s[52:54] if len(s) > 53 else ""

    if not pid.is_image:
        p = (pid.processing_raw + "_______")[:7]
        pid.flags = {
            "wavelength": p[0] == "W",
            "expid": p[1],
            "aci_count": p[2:4],
            "cosmic": p[4],
            "gain": p[5],
            "normalize": p[6],
        }
    return pid


# --------------------------------------------------------------------------
# label 读取
# --------------------------------------------------------------------------


@dataclass
class Table:
    name: str
    offset: Optional[int]
    records: Optional[int]
    channels: Optional[int]
    raw_name: str = ""


@dataclass
class Label:
    path: str
    logical_identifier: str = ""
    product_id: str = ""          # <alternate_id>，即官方 VICAR PRODUCT_ID
    title: str = ""
    data_file: str = ""
    start_time: str = ""
    stop_time: str = ""
    sol: str = ""
    sclk_start: str = ""
    sclk_stop: str = ""
    release: str = ""
    sequence_id: str = ""
    version: str = ""
    target: str = ""
    instrument: str = ""
    tables: list = field(default_factory=list)
    errors: list = field(default_factory=list)

    @property
    def pid(self) -> ProductId:
        return decode_product_id(self.product_id or self.data_file or os.path.basename(self.path))

    @property
    def n_spectra(self):
        """光谱条数（取所有光谱类表里最大的 records）。"""
        vals = [t.records for t in self.tables
                if t.records and re.search(r"SPECTRA|INTENSITY|PHOTODIODE|POSITION", t.name, re.I)]
        return max(vals) if vals else None

    @property
    def processing_summary(self) -> str:
        """从表名读出处理状态 —— 这是 label 里唯一能反映处理的地方。"""
        names = " ".join(t.name.upper() for t in self.tables)
        hits = []
        if "LASER-NORMALIZED" in names or "LASER_NORMALIZED" in names:
            hits.append("激光归一化")
        if "DARK-SUBTRACTED" in names or "DARK_SUBTRACTED" in names:
            hits.append("已扣暗背景")
        if "ACTIVE" in names:
            hits.append("主动帧")
        if "BACKGROUND" in names or "DARK" in names:
            hits.append("暗背景帧")
        return " · ".join(hits) if hits else "（表名未标明处理）"


def label_path_for(data_path: str) -> Optional[str]:
    """给定数据文件，找到同名的 .xml label。"""
    stem = os.path.splitext(data_path)[0]
    for cand in (stem + ".xml", stem + ".XML", stem + "_label.xml"):
        if os.path.exists(cand):
            return cand
    return None


def read_label(path: str) -> Label:
    """解析一个 PDS4 label。"""
    lab = Label(path=path)
    try:
        tree = ET.parse(path)
    except Exception as exc:                     # noqa: BLE001
        lab.errors.append("XML 解析失败: %s" % exc)
        return lab
    root = tree.getroot()

    lab.logical_identifier = _text(root, "logical_identifier", "") or ""
    lab.product_id = _text(root, "alternate_id", "") or ""
    lab.title = _text(root, "title", "") or ""
    lab.data_file = _text(root, "file_name", "") or ""
    lab.start_time = _text(root, "start_date_time", "") or ""
    lab.stop_time = _text(root, "stop_date_time", "") or ""
    lab.sol = _text(root, "sol_number", "") or ""
    lab.sclk_start = _text(root, "spacecraft_clock_start", "") or ""
    lab.sclk_stop = _text(root, "spacecraft_clock_stop", "") or ""
    lab.release = _text(root, "release_number", "") or ""
    lab.sequence_id = _text(root, "sequence_id", "") or ""
    lab.version = _text(root, "version_id", "") or ""

    for el in _all(root, "Target_Identification"):
        n = _text(el, "name")
        if n:
            lab.target = n
            break
    for el in _all(root, "Observing_System_Component"):
        n = _text(el, "name")
        if n and "sherloc" in n.lower():
            lab.instrument = n
            break

    # 表：name / offset / records / repetitions
    for td in _all(root, "Table_Delimited"):
        name = _text(td, "name", "") or ""
        off = _text(td, "offset")
        rec = _text(td, "records")
        reps = None
        for el in td.iter():
            if _strip_ns(el.tag) == "repetitions" and el.text:
                try:
                    reps = int(el.text)
                except ValueError:
                    reps = None
                break
        lab.tables.append(
            Table(
                name=name.replace("_", " ").strip(),
                offset=int(off) if off and off.isdigit() else None,
                records=int(rec) if rec and rec.isdigit() else None,
                channels=reps,
                raw_name=name,
            )
        )
    return lab
