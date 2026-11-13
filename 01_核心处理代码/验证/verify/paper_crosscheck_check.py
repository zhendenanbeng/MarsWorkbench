# -*- coding: utf-8 -*-
"""SHERLOC 逐靶对账摘要护栏（发行自包含）：复现数只增不缩 + 判定结构完整。"""
import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0, str(_BootPath(__file__).resolve().parents[3] / "02_图形界面代码"))
import source_paths as _layout

import json
import sys

SUMMARY = (_layout.project_layout.VERIFY
           / "fixtures/历史数值基准/paper_crosscheck_summary.json")
STATES = {"reproduced", "weak", "not_reproduced", "excluded", "unmapped"}
MIN_JUDGEABLE = 14

TOTAL = FAIL = 0


def check(ok, name):
    global TOTAL, FAIL
    TOTAL += 1
    FAIL += not ok
    print(("OK " if ok else "FAIL ") + name, flush=True)


def main():
    if not SUMMARY.is_file():
        print("FAIL 摘要文件缺失")
        print("共 %d 项断言，%d 项失败" % (TOTAL, FAIL))
        return 1
    doc = json.loads(SUMMARY.read_text(encoding="utf-8"))
    n_judge = doc.get("n_judgeable", 0)
    n_repro = doc.get("n_reproduced", 0)
    check(n_judge >= MIN_JUDGEABLE, "可判定类目不缩（现 %d，下限 %d）" % (n_judge, MIN_JUDGEABLE))
    check(n_repro >= 7, "复现数不缩（现 %d，下限 7）" % n_repro)
    check(n_repro <= n_judge, "复现数 ≤ 可判定数")
    bad = []
    for tgt, r in doc.get("per_target", {}).items():
        for pc, e in r.get("verdict", {}).items():
            if e.get("status") not in STATES:
                bad.append("%s/%s" % (tgt, pc))
    check(not bad, "判定状态只用约定枚举（%d 类目）" % (n_judge,))
    check(bool(doc.get("source_numbers_sha256")), "来源 numbers.json 指纹在场")
    check(bool(doc.get("pairs")), "对账靶对清单在场")
    print("共 %d 项断言，%d 项失败" % (TOTAL, FAIL))
    return int(FAIL > 0)


if __name__ == "__main__":
    raise SystemExit(main())
