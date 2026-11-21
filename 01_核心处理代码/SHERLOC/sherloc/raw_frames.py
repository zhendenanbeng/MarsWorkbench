# -*- coding: utf-8 -*-
"""SHERLOC 原始帧（EDR）读取与自行扣暗（ERA−ERB）。

官方 RRS = ERA（激光主动帧）− ERB（暗背景帧），再经激光能量归一化、
波长校正、宇宙射线/增益合并的产物。本模块只做其中**自洽可验证**的
一步：逐点位 ERA−ERB 暗帧相减；波长轴借用同序列（或同 sol 任意）
RRS 的官方波长校正表——对原始通道是近似，峰位可能有微小系统差。
激光能量归一化不在本层：公式含仪器内标细节，自行近似会产出
"看起来官方"的错误数字——边界如实声明，输出文件/文件夹均带标记。

文件结构（实测 ss__*ERA/ERB*J06.CSV）：
  ERA 主动帧段 = `ACTIVE_SPECTRA_REGION_1/2/3:`，ERB 暗帧段 =
  `DARK_SPECTRA_REGION_1/2/3:`——各一张表：列头行（R1_Channel_0…2147）
  + 逐激光点数据行（每行 2148 列）；另含 WAVELENGTH (NM):（仅 RRS 有）。
  ERB 暗帧与对应 ERA 同构、行数相同（暗采集镜像扫描几何），
  按行号一一对应相减。

输出：RRS 同布局 CSV（WAVELENGTH 行 + DARK_SUBTRACTED_SPECTRA:
_REGION_N 段）+ 适配 offset/records 的配对 XML——现有 SHERLOC 台
可直接打开输出文件夹（organize/reader 按名字匹配 REGION 表）。
"""
from __future__ import annotations

import io
import re
from pathlib import Path

import numpy as np

SECTION_RE = re.compile(r"^([A-Z][A-Z0-9_ ()]*):$")   # 含 WAVELENGTH (NM): 的括号空格
NAME_RE = re.compile(
    r"^[A-Za-z]+__(\d+)_(\d+)_(\d{3})([A-Za-z]{2,4})__(\w+?)_(\d+)_{3,}([A-Z]\d+)\.")


def parse_sections(path):
    """按 `大写标记:` 切分 SHERLOC CSV，返回 {标记: [数据行字符串]}。

    每段的第一行是列头（R1_Channel_…），这里一并保留，由调用方跳过。
    """
    sections = {}
    cur = None
    for raw in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        ln = raw.rstrip("\r")
        m = SECTION_RE.match(ln.strip())
        if m:
            cur = m.group(1)
            sections.setdefault(cur, [])
        elif cur is not None and ln.strip():
            sections[cur].append(ln)
    return sections


def read_region_frames(path):
    """读原始帧的三个 Region 矩阵：{1: (n,2148) 数组, 2: …, 3: …}。

    ERA 主动帧段名 = ACTIVE_SPECTRA_REGION_N；ERB 暗帧段名 =
    DARK_SPECTRA_REGION_N——两者都接受。缺段不报错（该段不在字典里）。
    """
    sec = parse_sections(path)
    out = {}
    for n in (1, 2, 3):
        for prefix in ("ACTIVE_SPECTRA_REGION_%d" % n, "DARK_SPECTRA_REGION_%d" % n):
            rows = sec.get(prefix, [])
            data = [ln for ln in rows if not ln.upper().startswith("R%d_CHANNEL_" % n)]
            if not data:
                continue
            try:
                out[n] = np.loadtxt(io.StringIO("\n".join(data)), delimiter=",", ndmin=2)
            except (ValueError, OSError):
                continue
            break
    return out


def parse_name(path):
    """文件名 → (类型, 序列, expID, sol)；不匹配返回 None。"""
    m = NAME_RE.match(Path(path).name)
    if not m:
        return None
    _mission, sol, _sclk, ptype, seq, exp, _jn = m.groups()
    return ptype.upper(), seq, exp, sol


def find_raw_products(folder):
    """目录里 ERA/ERB/RRS 三类产品 → {类型: [(seq, exp, 路径), …]}。"""
    found = {"ERA": [], "ERB": [], "RRS": []}
    unparsed = 0
    for f in sorted(Path(folder).iterdir()):
        if f.suffix.lower() != ".csv":
            continue
        if re.search(r"_\d{3}rrs__", f.name, re.I):
            # RRS 文件名存在 w108cgn/b108zpz/w1 等多种尾变体，只认类型位
            found["RRS"].append((None, None, None, f))
            continue
        info = parse_name(f)
        if info and info[0] in found:
            found[info[0]].append((info[1], info[2], info[3], f))
        else:
            unparsed += 1
    found["unparsed"] = unparsed
    return found


def pair_dark(era_list, erb_list):
    """按 (序列, expID) 配对主动/暗帧。

    返回 (配对列表 [(seq, exp, era_path, erb_path)], 未配对 [(seq, exp, 原因)])。
    """
    darks = {(seq, exp): (sol, p) for seq, exp, sol, p in erb_list}
    pairs, missing = [], []
    for seq, exp, sol, p in era_list:
        d = darks.get((seq, exp))
        if d is None:
            missing.append((seq, exp, "无同序列同 expID 的 ERB 暗帧"))
        else:
            pairs.append((seq, exp, sol, p, d[1]))
    return pairs, missing


def read_wavelength(rrs_path):
    """从 RRS 读官方波长轴（nm）：WAVELENGTH 段标记后的第二行是数值行。

    段键名按原文捕获（含 "(NM)"），按前缀 WAVELENGTH 查找。
    """
    sec = parse_sections(rrs_path)
    key = next((k for k in sec if k.upper().startswith("WAVELENGTH")), None)
    rows = sec.get(key, []) if key else []
    if len(rows) < 2:
        return None
    try:
        wl = np.array([float(x) for x in rows[1].split(",") if x.strip()],
                      dtype=float)
        return wl if wl.size else None
    except ValueError:
        return None


def dark_subtract(active, dark):
    """逐点位相减；行数不一致时按较少数量的行截断并如实返回信息。"""
    a = np.asarray(active, dtype=float)
    d = np.asarray(dark, dtype=float)
    if a.shape[1] != d.shape[1]:
        raise ValueError("ERA 与 ERB 通道数不一致（%d vs %d）"
                         % (a.shape[1], d.shape[1]))
    n = min(a.shape[0], d.shape[0])
    return a[:n] - d[:n]


def write_scientific_csv(out_path, wavelength, regions):
    """写 RRS 同布局的科学谱 CSV（三段 DARK_SUBTRACTED_SPECTRA）。

    返回 (offsets, records)：{表名: offset} 与 {表名: records} 供 XML
    使用（波长表 1 行、Region 表 n_points 行）；wavelength 为 None 时
    省略波长段（此时台里无法出图，仅存档）。
    """
    offsets = {}
    records = {}
    n_ch = next(iter(regions.values())).shape[1]
    with open(out_path, "wb") as f:
        pos = 0

        def w(line):
            nonlocal pos
            f.write((line + "\r\n").encode("utf-8"))
            pos += len((line + "\r\n").encode("utf-8"))

        if wavelength is not None:
            w("WAVELENGTH (NM):")
            w(",".join("Channel_%d" % i for i in range(n_ch)))
            offsets["WAVELENGTH NM"] = pos
            records["WAVELENGTH NM"] = 1
            w(",".join("%.6f" % v for v in wavelength))
        for n in sorted(regions):
            mat = regions[n]
            w("DARK_SUBTRACTED_SPECTRA:_REGION_%d" % n)
            w(",".join("R%d_Channel_%d" % (n, i) for i in range(mat.shape[1])))
            offsets["DARK SUBTRACTED SPECTRA REGION %d" % n] = pos
            records["DARK SUBTRACTED SPECTRA REGION %d" % n] = mat.shape[0]
            for row in mat:
                w(",".join("%.3f" % v for v in row))
    return offsets, records


LABEL_TEMPLATE = """<?xml version="1.0" encoding="UTF-8"?>
<!-- 此文件为 MarsWorkbench 自扣暗派生产品（ERA−ERB 逐点位相减）的配对标签，
     非官方 PDS 产品；波长轴借用同 sol 官方 RRS 的波长校正表。 -->
<Product_Observational xmlns="http://pds.nasa.gov/pds4/pds/v1">
  <Identification_Area>
    <logical_identifier>urn:nasa:pds:mars2020_sherloc:marsworkbench_derived</logical_identifier>
    <alternate_id>__FILE_STEM__</alternate_id>
    <version_id>1.0</version_id>
    <title>MarsWorkbench self dark-subtracted spectra (ERA-ERB)</title>
  </Identification_Area>
  <Observing_System>
    <Observing_System_Component>
      <name>Sherloc Deep Uv Resonance Raman And Flourescence Spectrometer</name>
      <type>Instrument</type>
    </Observing_System_Component>
  </Observing_System>
  <Target_Identification>
    <name>Mars surface</name>
  </Target_Identification>
  <File_Area_Observational>
    <File>
      <file_name>__FILE_NAME__</file_name>
    </File>
__TABLES__
  </File_Area_Observational>
</Product_Observational>
"""

TABLE_TEMPLATE = """    <Table_Delimited>
      <name>__NAME__</name>
      <offset unit="byte">__OFFSET__</offset>
      <records unit="count">__RECORDS__</records>
      <Record_Delimited>
        <fields>1</fields>
        <groups>1</groups>
        <Group_Field_Delimited>
          <repetitions>__REPS__</repetitions>
        </Group_Field_Delimited>
      </Record_Delimited>
    </Table_Delimited>"""


def write_label(out_xml, csv_name, offsets, records_map, n_ch):
    """为自扣暗 CSV 生成配对 XML（offset/records 为实测字节布局）。"""
    blocks = []
    names = {}
    if "WAVELENGTH NM" in offsets:
        names["WAVELENGTH NM"] = "WAVELENGTH NM"
    for k in sorted(offsets):
        if k.startswith("DARK"):
            names[k] = k
    for key, name in names.items():
        if key not in offsets:
            continue
        blocks.append(TABLE_TEMPLATE
                      .replace("__NAME__", name)
                      .replace("__OFFSET__", str(offsets[key]))
                      .replace("__RECORDS__", str(records_map.get(key, 1)))
                      .replace("__REPS__", str(n_ch)))
    xml = (LABEL_TEMPLATE
           .replace("__FILE_NAME__", csv_name)
           .replace("__FILE_STEM__", Path(csv_name).stem)
           .replace("__TABLES__", "\n".join(blocks)))
    Path(out_xml).write_text(xml, encoding="utf-8")


def process_folder(folder, out_dir=None, region_list=(1, 2, 3)):
    """整个目录的自扣暗流程：配对 → 相减 → 写 CSV+XML → 返回报告。

    报告字段：pairs（成功明细）、skipped（未配对及原因）、out_dir、
    wavelength_source、boundary（口径声明）。不动原文件。
    """
    folder = Path(folder)
    out_dir = Path(out_dir) if out_dir else folder / "自扣暗_输出"
    found = find_raw_products(folder)
    pairs, missing = pair_dark(found["ERA"], found["ERB"])

    wl = None
    wl_src = None
    for _seq, _exp, _sol, p in found["RRS"]:
        wl = read_wavelength(p)
        if wl is not None:
            wl_src = p.name
            break

    out_dir.mkdir(parents=True, exist_ok=True)
    pairs_out, skipped = [], list(missing)
    for seq, exp, sol, era_p, erb_p in pairs:
        try:
            act = read_region_frames(era_p)
            dk = read_region_frames(erb_p)
            regions = {}
            for n in region_list:
                if n in act and n in dk:
                    regions[n] = dark_subtract(act[n], dk[n])
            if not regions:
                skipped.append((seq, exp, "ERA/ERB 均无可解析的 Region 帧"))
                continue
            n_points = next(iter(regions.values())).shape[0]
            # 输出名沿用 ERA 名做类型/版本替换（保留 3 位后缀，pid 可解析）
            out_name = era_p.name.replace("ERA__", "rrs__")
            out_name = re.sub(r"___J\d+\.", "___J96.", out_name)
            csv_path = out_dir / out_name
            offsets, records_map = write_scientific_csv(csv_path, wl, regions)
            write_label(str(csv_path).replace(".CSV", ".xml").replace(".csv", ".xml"),
                        csv_path.name, offsets, records_map,
                        len(wl) if wl is not None else next(iter(regions.values())).shape[1])
            pairs_out.append({"seq": seq, "exp": exp, "points": n_points,
                              "regions": sorted(regions), "csv": csv_path.name})
        except (ValueError, OSError) as exc:
            skipped.append((seq, exp, str(exc)[:120]))

    (out_dir / "README.txt").write_text(
        "MarsWorkbench 自扣暗派生产品（ERA−ERB 逐点位相减）\n"
        "来源文件夹: %s\n"
        "边界: 未做激光能量归一化/宇宙射线/增益合并——绝对强度与官方 RRS 不可比；\n"
        "波长轴借用官方 RRS 波长校正表（对原始通道是近似）。判读看峰位与相对形态。\n"
        % folder.name, encoding="utf-8")
    return {"pairs": pairs_out, "skipped": skipped, "out_dir": str(out_dir),
            "wavelength_source": wl_src,
            "boundary": "仅暗帧相减；无激光能量归一化"}


