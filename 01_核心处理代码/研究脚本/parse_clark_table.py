# -*- coding: utf-8 -*-
"""把 Clark et al. (2026) 的 Table 1 解析成结构化峰位表。

论文 = Clark, J.V. et al. (2026), "Deep-Ultraviolet (DUV) Raman and Fluorescence
of Perchlorates and Chlorates: Implications for the Detection of Oxychlorines by
the SHERLOC Instrument in Jezero Crater, Mars", Earth and Space Science,
DOI 10.1029/2025EA004512。表 1 = ACRONM 实测的各氧氯酸盐拉曼峰参数
（峰位 ω、归一化强度、观测峰宽 Γ_obs、仪器展宽 Γ_IRB、**换算到 SHERLOC 的峰宽
Γ_SRLC**）。

为什么要脚本解析而不是手抄：这张表有上百行、峰位精确到小数，
手抄迟早出错；而且以后换版本可以重跑。

用法：python parse_clark_table.py
输出：oxychlorine_table1.json（样品 → [(峰位, 归一化强度), ...]）
"""

from __future__ import annotations

import json
import os
import re
import sys

from pypdf import PdfReader

HERE = os.path.dirname(os.path.abspath(__file__))
PDF_GLOBS = [
    os.path.join(os.path.expanduser("~"), "Desktop"),
]


def find_pdf():
    for d in PDF_GLOBS:
        for name in os.listdir(d):
            if name.lower().endswith(".pdf") and "Perchlorate" in name:
                return os.path.join(d, name)
    raise SystemExit("找不到那篇 PDF（文件名里应含 Perchlorate）")


def load_text(pdf):
    r = PdfReader(pdf)
    # Table 1 跨第 7–9 页（0 基 6–8）
    return "\n".join((r.pages[i].extract_text() or "") for i in range(6, 10))


HEAD = re.compile(r"^([A-Z][a-z]?(?: [IVX]+)? (?:perchlorate|chlorate));\s*(.*)$")
# 行里是「模式 峰位(误差) 强度 ...」，模式可能是一两个词（ClO4 νs、H2O ν、Overtone）
ROW = re.compile(r"^\s*((?:\S+\s+){0,2}?)(\d{3,4}\.\d)\s*\(\s*[\d.]+\s*\)\s+([\d.]+)")


def parse(text):
    text = text.replace("\u2010", "-").replace("\u2019", "'")
    data = {}
    cur = None
    for line in text.split("\n"):
        m = HEAD.match(line.strip())
        if m:
            cur = m.group(1)
            data[cur] = {"vendor": m.group(2).strip(), "bands": []}
            continue
        if cur is None:
            continue
        m2 = ROW.match(line)
        if m2:
            data[cur]["bands"].append(
                {"mode": " ".join(m2.group(1).split()),
                 "pos": float(m2.group(2)), "intensity": float(m2.group(3))})
    return data


def mode_key(mode):
    """把论文里的模式标签翻成 library.py 的键名。

    用标签而不是用位置猜：论文里 1563 是"VW"、1621 是"H2O β"，
    两者位置很近，按位置分不开。
    """
    m = mode.replace(" ", "")
    if "δas" in m:
        return "ClO_delta_as"
    if "δs" in m:
        return "ClO_delta_s"
    if "νas" in m:
        return "ClO_nu_as"
    if "νs" in m:
        return "ClO_nu_s"
    if "Overtone" in m:
        return "overtone"
    if m.startswith("H2O") and ("β" in m):
        return "H2O_bend"
    if m.startswith("H2O"):
        return "H2O_stretch"
    return "VW_unassigned"


def emit_python(data):
    """输出可以直接粘进 library.py 的条目。"""
    lines = ['    "oxychlorine": {',
             '        "label": "高氯酸盐与氯酸盐",',
             '        # 论文原话：对称伸缩 925-943（Na perchlorate 例外，950）、',
             '        # 反对称伸缩 971-1190、对称 Cl-O 变形 600-650、',
             '        # 反对称 Cl-O 变形 400-500、H2O 弯曲 1600-1630。',
             '        "modes": ["ClO_nu_s", "ClO_nu_as", "ClO_delta_s", "ClO_delta_as",',
             '                  "overtone", "H2O_bend", "H2O_stretch", "VW_unassigned"],',
             '        "samples": [']
    for name, info in data.items():
        groups = {}
        for b in info["bands"]:
            groups.setdefault(mode_key(b["mode"]), []).append(b["pos"])
        parts = []
        for k in ("ClO_nu_s", "ClO_nu_as", "ClO_delta_s", "ClO_delta_as",
                  "overtone", "H2O_bend", "H2O_stretch", "VW_unassigned"):
            if k not in groups:
                continue
            v = sorted(groups[k])
            parts.append('"%s": %s' % (k, "%.1f" % v[0] if len(v) == 1
                                       else "(" + ", ".join("%.1f" % x for x in v) + ")"))
        vendor = info["vendor"].split(";")[0].strip() or "?"
        lines.append('            ("%s", "%s",\n             {%s}),'
                     % (name, vendor, ", ".join(parts)))
    lines += ["        ],", "    },"]
    return "\n".join(lines)


def main():
    pdf = find_pdf()
    data = parse(load_text(pdf))
    out_path = os.path.join(HERE, "oxychlorine_table1.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=1)
    if "--emit-python" in sys.argv:
        print(emit_python(data))
        return 0
    print("来源：%s" % os.path.basename(pdf))
    print("解析到 %d 个样品\n" % len(data))
    print("%-22s %4s %-34s %s" % ("样品", "峰数", "指纹区 800–1800（峰位(强度)）", "最强峰"))
    for k, v in data.items():
        bands = v["bands"]
        if not bands:
            print("%-22s %4d %s" % (k, 0, "（没解析到峰位——表格跨页或版式特殊）"))
            continue
        strong = max(bands, key=lambda b: b["intensity"])
        fp = ", ".join("%.1f(%.3f)" % (b["pos"], b["intensity"])
                       for b in bands if 800 <= b["pos"] <= 1800)
        print("%-22s %4d %-34s %.1f" % (k, len(bands), fp[:34], strong["pos"]))
    print("\n已写 %s" % out_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
