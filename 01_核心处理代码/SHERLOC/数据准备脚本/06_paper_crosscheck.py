#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SHERLOC 数据准备脚本 06：软件分类输出 vs 已发表结论的逐靶对比检查
====================================================================
目的：把工作台的自主分类链（qc.analyse_spectra + minerals 规则候选）放到
行业标准（同行评审论文的靶点矿物结论）下检验——这是"用软件处理的数据
支撑科学结论"之前必须先过的对账关。

输入：
  * published_targets.json（已发表靶点矿物声明台账）
  * 本地已下载的 SHERLOC sol 数据目录（--pair 参数指定 靶名=sol 目录）

判定口径（诚实边界，不宣称准确率）：
  * 复现 reproduced：论文声明检出某类，软件在同一 sol 的观测中有该类
    （或其模糊键）候选组 ≥1 点；
  * 未复现 not_reproduced：论文检出而软件无候选——记录为差异，可能是
    数据子集/处理口径/声明粒度（靶点 vs sol 级观测）差异，不判对错；
  * 排除 excluded：论文类目在 SHERLOC 深紫外拉曼类别体系中无对应
    （如铁氧化物拉曼惰性、halite、天然表面观测缺失）；
  * 额外候选 extra：软件有候选而论文该靶未声明——记录，不判错。

用法：
    python 06_paper_crosscheck.py --pair Garde=03_下载数据/SHERLOC/练习数据/raw/sol_00207 \
                                   --pair Alfalfa=03_下载数据/SHERLOC/练习数据/raw/sol_00370
"""

from __future__ import annotations

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0, str(_BootPath(__file__).resolve().parents[3] / "02_图形界面代码"))
import source_paths as _layout            # noqa: E402

import argparse
import datetime
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

import project_layout                      # noqa: E402
from sherloc import organize, qc           # noqa: E402

# 论文类目 → 软件类别键（minerals.py RANGES/LABELS 的英文键）。
# 模糊键（"sulfate_or_pyroxene" 等）同时计入两个论文类目，但标 weak。
PAPER_TO_KEYS = {
    "olivine": {"strong": {"olivine"}, "weak": set()},
    "pyroxene": {"strong": {"pyroxene"}, "weak": {"sulfate_or_pyroxene"}},
    "sulfate": {"strong": {"sulfate"}, "weak": {"sulfate_or_pyroxene"}},
    "feldspar": {"strong": {"silicate"}, "weak": set()},
    "phyllosilicates": {"strong": {"silicate"}, "weak": set()},
    "chlorite": {"strong": {"silicate"}, "weak": set()},
    "carbonate": {"strong": {"carbonate"}, "weak": set()},
    "phosphate": {"strong": {"phosphate"}, "weak": {"perchlorate_or_phosphate"}},
    "perchlorate": {"strong": {"perchlorate_chlorate", "perchlorate"},
                    "weak": {"perchlorate_or_phosphate"}},
}
EXCLUDED_PAPER_CLASSES = {"fe-oxides", "feti-oxides", "halite", "organics",
                          "quartz", "ca-sulfate", "mg-ca-sulfate", "hydrated-phases",
                          "fe-oxide", "amorphous-silicate"}


def normalise(cls):
    return cls.strip().lower().replace("_", "-")


def analyse_sol(sol_dir):
    """跑一个 sol 目录的全部光谱扫描，返回 {扫描key: {class_keys: 点数并集}}。"""
    _, scans = organize.scan_folder(str(sol_dir))
    scans = [s for s in scans if s.parts and s.spectra
             and s.spectra.product_type in ("RRS", "RRA")]
    out = {}
    for sc in scans:
        try:
            wl, r1, r2, r3 = sc.read_regions()
            if wl is None or r1 is None:
                continue
            result = qc.analyse_spectra(wl, r1, r2, r3)
        except Exception as exc:                         # noqa: BLE001
            out[sc.key] = {"error": str(exc)[:80]}
            continue
        class_points = {}
        for p in result.get("point_detail", []):
            for c in p.get("classes", []):
                key = c.get("key") or c.get("class")
                if key:
                    class_points[key] = class_points.get(key, 0) + 1
        groups = {}
        for g in result.get("mineral_groups", []):
            key = g.get("class")
            if key:
                groups[key] = groups.get(key, 0) + 1
        out[sc.key] = {"class_points": class_points, "groups": groups,
                       "n_points": int(sc.n_points or 0)}
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description="软件分类 vs 已发表结论逐靶对比")
    ap.add_argument("--pair", action="append", required=True,
                    help="靶名=sol数据目录（可重复），靶名须与 published_targets.json 的 target 一致")
    ap.add_argument("--ledger", default=None, help="台账路径（默认 fixtures 内 published_targets.json）")
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)

    ledger_path = Path(args.ledger) if args.ledger else (
        _layout.project_layout.VERIFY / "fixtures/历史数值基准/published_targets.json")
    ledger = json.loads(ledger_path.read_text(encoding="utf-8"))
    rows_by_target = {}
    for r in ledger["rows"]:
        rows_by_target.setdefault(r["target"], []).append(r)

    pairs = []
    for item in args.pair:
        target, _, sol_dir = item.partition("=")
        if target not in rows_by_target:
            print("台账里没有靶 %r；核对 target 名。" % target)
            return 1
        if not Path(sol_dir).is_dir():
            print("sol 目录不存在：%s" % sol_dir)
            return 1
        pairs.append((target, Path(sol_dir)))

    results = {}
    for target, sol_dir in pairs:
        print("== 处理靶 %s（%s）==" % (target, sol_dir))
        scans = analyse_sol(sol_dir)
        paper_classes = set()
        for r in rows_by_target[target]:
            paper_classes |= {normalise(c) for c in r.get("published_classes", [])}
        # 汇总该 sol 全部扫描的类别点数
        union_points = {}
        scan_errors = []
        for key, res in scans.items():
            if "error" in res:
                scan_errors.append({"scan": key, "error": res["error"]})
                continue
            for k, n in res["class_points"].items():
                union_points[k] = max(union_points.get(k, 0), n)
        present = set(union_points)

        verdict = {}
        for pc in sorted(paper_classes):
            entry = {"paper_class": pc}
            if pc in EXCLUDED_PAPER_CLASSES:
                entry["status"] = "excluded"
                entry["note"] = "SHERLOC 深紫外拉曼类别体系无对应（拉曼惰性/无单矿物键）"
            else:
                mapping = PAPER_TO_KEYS.get(pc)
                if mapping is None:
                    entry["status"] = "unmapped"
                    entry["note"] = "映射表未覆盖（请人工核对后补充 PAPER_TO_KEYS）"
                else:
                    strong = sorted(mapping["strong"] & present)
                    weak = sorted(mapping["weak"] & present)
                    if strong:
                        entry["status"] = "reproduced"
                        entry["software_keys"] = strong
                        entry["max_points"] = max(union_points[k] for k in strong)
                    elif weak:
                        entry["status"] = "weak"
                        entry["software_keys"] = weak
                        entry["note"] = "仅模糊键候选（缺次级峰区分），计部分复现"
                        entry["max_points"] = max(union_points[k] for k in weak)
                    else:
                        entry["status"] = "not_reproduced"
                        entry["note"] = "本 sol 观测中无该类候选（记录差异，不判对错）"
            verdict[pc] = entry
            print("  %-18s %s %s" % (pc, entry["status"],
                                     entry.get("software_keys", "")))
        extra = sorted(present - {k for m in PAPER_TO_KEYS.values()
                                  for k in (m["strong"] | m["weak"])})
        results[target] = {
            "sol_dir": str(sol_dir), "scans": {k: v for k, v in scans.items()},
            "paper_classes": sorted(paper_classes), "verdict": verdict,
            "extra_software_candidates": extra,
            "scan_errors": scan_errors,
        }
        if extra:
            print("  额外候选（论文该靶未声明）：%s" % ", ".join(extra))

    n_repro = sum(1 for t in results.values() for e in t["verdict"].values()
                  if e["status"] == "reproduced")
    n_all = sum(1 for t in results.values() for e in t["verdict"].values()
                if e["status"] in ("reproduced", "weak", "not_reproduced"))
    print("可判定类目 %d 项，复现（含弱）%d 项" % (n_all, n_repro))

    stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = Path(args.out) if args.out else (
        Path(__file__).resolve().parent / "out" / "paper_crosscheck" / stamp)
    out.mkdir(parents=True, exist_ok=True)
    hashes = project_layout.source_hashes(include_verification=False)
    data_files = sorted({f for _, p in pairs for f in Path(p).rglob("*.csv")
                         if "rrs" in f.name.lower()})
    doc = {
        "scope": "软件分类链 vs 已发表靶点结论的对账检查（冻结测量）。判定口径是"
                 "类别候选复现与否，不是准确率：样本是 sol 级观测、声明是靶点级，"
                 "粒度差异如实记录。不构成矿物鉴定正确率证据。",
        "ledger": str(ledger_path),
        "ledger_sha256": hashlib.sha256(ledger_path.read_bytes()).hexdigest(),
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "pairs": {t: {"sol_dir": str(p)} for t, p in pairs},
        "data_sha256": {f.name: hashlib.sha256(f.read_bytes()).hexdigest()
                        for f in data_files},
        "code_sha256": hashes,
        "mapping": {k: {"strong": sorted(v["strong"]), "weak": sorted(v["weak"])}
                    for k, v in PAPER_TO_KEYS.items()},
        "results": results,
        "summary": {"n_judgeable": n_all, "n_reproduced": n_repro},
    }
    (out / "numbers.json").write_text(json.dumps(doc, ensure_ascii=False, indent=1),
                                      encoding="utf-8")
    lines = ["# 软件分类输出 vs 已发表结论（逐靶对比检查）", "",
             "- 台账：%s" % ledger_path.name,
             "- 判定口径：类别候选复现性（非准确率）", ""]
    for t, r in results.items():
        lines.append("## %s（%s）" % (t, r["sol_dir"]))
        lines.append("")
        lines.append("| 论文类目 | 判定 | 软件键 | 说明 |")
        lines.append("|---|---|---|---|")
        for pc, e in r["verdict"].items():
            lines.append("| %s | %s | %s | %s |" % (
                pc, e["status"], ", ".join(e.get("software_keys", [])),
                e.get("note", "检出点数 %s" % e.get("max_points", "-"))))
        if r["extra_software_candidates"]:
            lines.append("")
            lines.append("额外候选：%s" % ", ".join(r["extra_software_candidates"]))
        lines.append("")
    (out / "report.md").write_text("\n".join(lines), encoding="utf-8")
    print("报告与数字写入 %s" % out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
