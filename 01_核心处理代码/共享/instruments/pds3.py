# -*- coding: utf-8 -*-
"""PDS3 detached label（.lbl）解析 + 配套表格读取。

ChemCam / APXS 的 RDR 都是这个格式：一个文本 .lbl 描述结构（键 = 值），
数据在同名 .csv/.tab 里。这两个仪器的 label 都把列定义写在自己身上
（没有 ^STRUCTURE 指到别的文件），所以这里只做一个足够诚实的子集：

- ``KEY = 值``，值可以是带引号字符串、数字、(a, b, c) 元组
- ``/* 注释 */`` 跨行
- 行尾续行（下一行缩进开头且没有 = 号）
"""

from __future__ import annotations

import re


class PDS3Label(dict):
    """dict 子类，多一个 get_str/get_int 的宽容取值。"""

    def get_int(self, key, default=None):
        v = self.get(key)
        if v is None:
            return default
        try:
            return int(str(v).strip().strip('"'))
        except ValueError:
            return default


_ASSIGN = re.compile(r"^([A-Z0-9_:]+)\s*=\s*(.*)$", re.I)


def parse_lbl(path):
    """解析 .lbl，返回 PDS3Label。解析不了结构也没关系——只要能取出
    PLANET_DAY_NUMBER / SEQUENCE_ID 这类标识键，扫描整理就够用了。"""
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            text = f.read()
    except OSError:
        return PDS3Label()

    # 去掉跨行注释
    text = re.sub(r"/\*.*?\*/", " ", text, flags=re.S)
    out = PDS3Label()
    for line in text.splitlines():
        line = line.rstrip()
        m = _ASSIGN.match(line.strip())
        if not m:
            continue
        key, val = m.group(1).upper(), m.group(2).strip()
        if val.startswith('"'):
            out[key] = val.strip('"').strip()
        elif val.startswith("("):
            # 元组：本行截断了也无所谓，正则抓行内所有数字/字符串
            items = re.findall(r'"[^"]*"|[-+0-9.eE]+', val)
            out[key] = [i.strip('"') for i in items]
        else:
            out[key] = val
    return out


def read_table(path, delimiter=",", skiprows=0):
    """读 .csv/.tab 数字表为 numpy 数组。'#' 开头的注释行自动跳过。"""
    import numpy as np
    return np.loadtxt(path, delimiter=delimiter, skiprows=skiprows,
                      comments="#", ndmin=2)


# --------------------------------------------------------------------------
# 跨仪器对照用的元数据：目标名 + 停车位（site/drive）
# --------------------------------------------------------------------------

import io
import os


def target_from_chemcam_csv(path):
    """ChemCam 光谱 csv 的 '#' 注释头里读 `# TARGET = <名>`。

    侦察确认：真目标名只在这里（label 的 TARGET_NAME 一律 MARS）。
    没有 TARGET 行或文件打不开时返回 None——旧数据不炸。
    """
    try:
        with io.open(path, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                if not line.startswith("#"):
                    break                     # 注释头只在整个文件开头
                m = re.match(r"#\s*TARGET\s*=\s*(.+?)\s*$", line, re.I)
                if m:
                    # TARGET 行出现在 csv 表头段，同行后面可能跟着空值列
                    # （Newland , , , , …）：目标名只取第一个逗号前的部分
                    return m.group(1).split(",")[0].strip()
    except OSError:
        pass
    return None


def target_from_apxs_header(header_line):
    """APXS rsp csv 首行 `Target: <名>, Contains DAT files: ...` 里的目标名。"""
    if not header_line:
        return None
    m = re.match(r'\s*"?\s*Target\s*:\s*([^,]+?)\s*,?\s*"?$', header_line.split(",")[0] + ",")
    if m:
        return m.group(1)
    m = re.match(r'\s*"?\s*Target\s*:\s*(.+?)(?:\s*,\s*Contains.*)?\s*"?$', header_line.strip().strip('"'))
    return m.group(1).strip() if m else None


def site_drive_from_lbl(label):
    """PDS3 label 的 ROVER_MOTION_COUNTER 前两分量 → (site, drive) int 元组。

    两仪器（ChemCam/APXS）的 label 都带这个 10 元组；SITE+DRIVE 相同 =
    同一停车位置，是跨仪器同址判据。解析失败返回 None。
    """
    if label is None:
        return None
    rmc = label.get("ROVER_MOTION_COUNTER")
    if not isinstance(rmc, (list, tuple)) or len(rmc) < 2:
        return None
    try:
        return int(float(str(rmc[0]))), int(float(str(rmc[1])))
    except ValueError:
        return None


def normalize_target(name):
    """目标名规范化：小写、下划线/空格统一、去常见前后缀（_rp 之类）。

    用于「同名匹配」的宽松比较；返回 None 表示无法规范化。
    """
    if not name:
        return None
    s = str(name).strip().strip('"').strip("'").lower()
    s = re.sub(r"[\s\-]+", "_", s)
    s = re.sub(r"_(rp|residue|brushed|rb|na|nosed|raster)$", "", s)
    s = re.sub(r"^(the_)?", "", s)
    return s or None
