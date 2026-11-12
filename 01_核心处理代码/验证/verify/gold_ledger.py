# -*- coding: utf-8 -*-
"""金标准差异台账：把《SHERLOC 对照报告》的逐靶点差异表机器化。

产出 verify/gold_ledger.json——第五期所有「差异收账」动作的基线快照：
  {sol, seq, target, unit, 类别, 论文点数, 工具点数, 论文峰位, 工具峰位,
   峰位差, 性质}
性质 ∈ {一致, 点数差, 峰位差, 缺口, 模糊, 假阳性}。

用法：
    python verify/gold_ledger.py            # 生成并校验（与 regression.GOLD 对账）
    python verify/gold_ledger.py --check    # 只校验现有台账（护栏用）

数据来源是对照报告的人类判读——台账里的「工具点数」是 2026-09 生成报告时的
快照，不是当下判定结果；当下漂移检测由 R5 的 ledger_guard 负责（允许区间
沿用 regression 的区间，不引用这里的快照数字）。
"""

from __future__ import annotations

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout

import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = str(_layout.project_layout.APP)
sys.path.insert(0, ROOT)

REPORT = os.path.join(os.path.dirname(ROOT), "04_使用文档", "SHERLOC_对照报告.md")
LEDGER = os.path.join(HERE, "gold_ledger.json")

# 对照报告第三节各靶点的小节标题行 → (sol, seq, target, unit)
TARGETS = [
    ("293", "srlc11360", "Quartier HDR_500_1", "Séítah"),
    ("257", "srlc11421", "Dourbes HDR_500_1", "Séítah"),
    ("207", None, "Garde HDR_500_2", "Séítah"),          # 对照报告只有比例行
    ("162", None, "Guillaumes HDR_250_1", "Máaz"),
    ("186", "srlc10500", "Bellegarde HDR_250_2", "Máaz"),
    ("349", "srlc11420", "Montpezat Survey", "Máaz"),
    ("370", "srlc11360", "Rimplas Survey", "Máaz"),
]

ROW_RX = re.compile(
    r"^\|\s*([^|]+?)\s*\|\s*(\d+)\s*\|\s*\*{0,2}(\d+)[^|]*?\|\s*(\S+)[^|]*\|"
    r"\s*(\S+)[^|]*\|\s*(.+?)\s*\|\s*$")


def parse_report():
    """Read the frozen numeric transcription; historical human reports are retired."""
    source = _layout.project_layout.VERIFY/'fixtures/历史数值基准/gold_report_rows.json'
    return json.loads(source.read_text(encoding='utf-8'))['rows']


def classify(rows):
    """给每行打性质标签。"""
    out = []
    for r in rows:
        kind = "点数差"
        if "假阳性" in r["note"] or "假阳" in r["note"]:
            kind = "假阳性"
        elif "缺口" in r["note"]:
            kind = "缺口"
        elif "模糊" in r["note"]:
            kind = "模糊"
        elif r["paper_points"] == r["tool_points"]:
            kind = "一致"
        r["kind"] = kind
        out.append(r)
    return out


def check_against_regression(rows):
    """与 regression.GOLD 对账：台账里 sol293/257 的论文点数与峰位必须一致。"""
    from verify.regression import GOLD          # noqa: PLC0415
    problems = []
    for r in rows:
        g = GOLD.get(r["sol"], {}).get(r["category"])
        if g is None:
            continue
        if g["n"][0] != r["paper_points"]:
            problems.append("sol%s %s 论文点数 台账%d vs GOLD%d"
                            % (r["sol"], r["category"], r["paper_points"], g["n"][0]))
        pos = str(int(g["pos"][0])) if float(g["pos"][0]).is_integer() else str(g["pos"][0])
        if pos not in (r["paper_pos"], r["tool_pos"]):
            problems.append("sol%s %s 峰位对不上: 台账%s/%s vs GOLD %s"
                            % (r["sol"], r["category"], r["paper_pos"], r["tool_pos"], pos))
    return problems


def main(argv):
    rows = classify(parse_report())
    if not rows:
        print("FAIL 解析不到差异行——对照报告结构变了？")
        return 1
    problems = check_against_regression(rows)
    if "--check" in argv:
        # 护栏模式：台账文件必须存在、条数不缩、对账通过
        if not os.path.exists(LEDGER):
            print("FAIL gold_ledger.json 不存在，先跑生成模式")
            return 1
        saved = json.load(open(LEDGER, encoding="utf-8"))
        saved_rows = len(saved.get("rows", []))          # 数行数，不是 dict 键数
        checks = fails = 0
        for name, ok in [
            ("台账条数不缩（%d ≥ %d）" % (len(rows), saved_rows),
             len(rows) >= saved_rows),
            ("台账与 regression.GOLD 对账", not problems),
        ]:
            checks += 1
            print("  %s  %s" % ("OK  " if ok else "FAIL", name))
            fails += 0 if ok else 1
        print("共 %d 项断言，%d 项失败" % (checks, fails))
        return 1 if fails else 0

    with open(LEDGER, "w", encoding="utf-8") as f:
        json.dump({"source": os.path.relpath(REPORT, ROOT),
                   "generated": "phase5-R1",
                   "rows": rows}, f, ensure_ascii=False, indent=1)
    print("台账写入 %d 行 → %s" % (len(rows), LEDGER))
    if problems:
        print("与 regression.GOLD 对账发现问题：")
        for p in problems:
            print("  -", p)
        return 1
    by_kind = {}
    for r in rows:
        by_kind[r["kind"]] = by_kind.get(r["kind"], 0) + 1
    print("性质分布：", by_kind)
    print("对账通过（sol293/257 的论文点数与峰位与 regression.GOLD 一致）")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
