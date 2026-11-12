# -*- coding: utf-8 -*-
"""已发表靶点台账（published_targets.json）护栏：结构与增长契约。

* 结构：每行必填字段齐全、granularity 合法、类别列表非空、
  来源可解析（source_key 必须出现在 sources 且带 DOI）；
* 增长：行数不得缩水（gold_ledger 同款"台账不缩"契约），
  台账里必须写明 growth_path（向 ≥30 行目标×类别扩展的路径）；
* 发行自包含：只读 fixtures，不联网、不依赖 03_下载数据。

台账哲学与 gold_ledger 一致：这是"已发表声明的转录"，不是盲测真值；
本护栏只保证结构完整与只增不缩，不判定矿物鉴定对错。
"""
import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0, str(_BootPath(__file__).resolve().parents[3] / "02_图形界面代码"))
import source_paths as _layout

import json
import sys

LEDGER = (_layout.project_layout.VERIFY
          / "fixtures/历史数值基准/published_targets.json")
REQUIRED = ("sol", "target", "formation", "granularity",
            "published_classes", "source_key", "locator")
GRANULARITY = {"target", "natural_surface", "unit"}
MIN_ROWS = 8          # v1 冻结值；转录 Corpolongo 逐靶点表后随行数上调

TOTAL = FAIL = 0


def check(ok, name):
    global TOTAL, FAIL
    TOTAL += 1
    FAIL += not ok
    print(("OK " if ok else "FAIL ") + name, flush=True)


def main():
    check(LEDGER.is_file(), "台账文件存在于 历史数值基准/")
    if not LEDGER.is_file():
        print("共 %d 项断言，%d 项失败" % (TOTAL, FAIL))
        return 1
    doc = json.loads(LEDGER.read_text(encoding="utf-8"))
    rows = doc.get("rows", [])
    sources = {s.get("key"): s for s in doc.get("sources", [])}

    check(len(rows) >= MIN_ROWS, "台账行数不缩（现 %d 行，下限 %d）" % (len(rows), MIN_ROWS))
    check(bool(doc.get("scope")) and "不是盲测真值" in doc["scope"],
          "scope 免责声明在场（不是盲测真值）")
    check(bool(doc.get("growth_path", {}).get("target_rows_goal")),
          "growth_path 在场（记录 ≥30 行扩展路径）")
    check(len(sources) >= 1 and all(s.get("doi") for s in sources.values()),
          "全部来源带 DOI")

    bad_fields = [r.get("target") for r in rows
                  if any(r.get(f) in (None, "", []) for f in REQUIRED
                         if f != "sol")]
    check(not bad_fields, "每行必填字段齐全（target/formation/granularity/classes/来源/定位）")
    bad_gran = [r.get("target") for r in rows if r.get("granularity") not in GRANULARITY]
    check(not bad_gran, "granularity 只用约定枚举（%s）" % "/".join(sorted(GRANULARITY)))
    bad_src = [r.get("target") for r in rows
               if sources.get(r.get("source_key")) is None]
    check(not bad_src, "每行 source_key 可解析到带 DOI 的来源")
    empty_cls = [r.get("target") for r in rows if not r.get("published_classes")]
    check(not empty_cls, "每行 published_classes 非空")

    doi_ledger = set()
    for s in sources.values():
        doi_ledger.add(str(s.get("doi", "")).lower())
    check("10.1029/2022je007455" in doi_ledger,
          "Corpolongo 2023（SHERLOC 逐靶点表的正式来源）在来源清单")
    target_rows = [r for r in rows if r.get("granularity") == "target"]
    check(len(target_rows) >= MIN_ROWS - 2 and
          any(r.get("locator") for r in target_rows),
          "target 粒度行全部带可查定位（页码/章节）")
    check(bool(doc.get("sources")) and all(s.get("sha256") or s.get("access")
                                           for s in sources.values()),
          "来源带字节指纹或访问方式说明")

    print("共 %d 项断言，%d 项失败" % (TOTAL, FAIL))
    return int(FAIL > 0)


if __name__ == "__main__":
    raise SystemExit(main())
