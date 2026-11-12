# -*- coding: utf-8 -*-
"""区间数值基准与标定台账护栏（四仪器自主定量管线的对账契约）。

* gold_intervals.json：每锚点每氧化物必须有 (reference, lower, upper) 且
  lower < reference < upper；sources 带 DOI 或 access；行数只增不缩。
* published_calibrations.json：sources 带 DOI；行必填字段齐全；
  granularity 枚举；scope 免责在场；growth_path 在场。
* 发行自包含：只读 fixtures，不联网、不依赖 03_下载数据。
"""
import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0, str(_BootPath(__file__).resolve().parents[3] / "02_图形界面代码"))
import source_paths as _layout

import json
import sys

BASE = _layout.project_layout.VERIFY / "fixtures/历史数值基准"
GRANULARITY = {"target_inventory", "metric", "method"}
MIN_CAL_ROWS = 14
MIN_ANCHOR_OXIDES = 8

TOTAL = FAIL = 0


def check(ok, name):
    global TOTAL, FAIL
    TOTAL += 1
    FAIL += not ok
    print(("OK " if ok else "FAIL ") + name, flush=True)


def main():
    gi_path = BASE / "gold_intervals.json"
    cal_path = BASE / "published_calibrations.json"
    check(gi_path.is_file() and cal_path.is_file(), "两个基准 JSON 在场")
    if not (gi_path.is_file() and cal_path.is_file()):
        print("共 %d 项断言，%d 项失败" % (TOTAL, FAIL))
        return 1
    gi = json.loads(gi_path.read_text(encoding="utf-8"))
    cal = json.loads(cal_path.read_text(encoding="utf-8"))

    # --- gold_intervals.json ---
    anchors = gi.get("anchors", {})
    check(set(anchors) >= {"chemcam", "apxs", "supercam"},
          "三仪器锚点组在场（chemcam/apxs/supercam）")
    n_ox = 0
    bad_ranges = []
    for inst, group in anchors.items():
        for anchor_name, anchor in group.items():
            for oxide, spec in anchor.get("oxides", {}).items():
                n_ox += 1
                if not (spec["lower"] < spec["reference"] < spec["upper"]):
                    bad_ranges.append("%s/%s/%s" % (inst, anchor_name, oxide))
    check(n_ox >= MIN_ANCHOR_OXIDES * 3, "区间锚点氧化物总数 ≥24（现 %d）" % n_ox)
    check(not bad_ranges, "全部区间满足 lower < reference < upper")
    check(bool(gi.get("scope")) and "区间护栏" in gi["scope"],
          "gold_intervals scope 免责在场")
    check(any(s.get("doi") or s.get("access") for s in gi.get("sources", [])),
          "gold_intervals 来源可追溯（DOI 或访问方式）")

    # --- published_calibrations.json ---
    rows = cal.get("rows", [])
    sources = {s.get("key"): s for s in cal.get("sources", [])}
    check(len(rows) >= MIN_CAL_ROWS, "标定台账行数不缩（现 %d，下限 %d）" % (len(rows), MIN_CAL_ROWS))
    check(bool(cal.get("scope")) and "不是盲测真值" in cal["scope"],
          "标定台账 scope 免责在场")
    check(bool(cal.get("growth_path", {}).get("rows_goal")),
          "growth_path 在场（成分数值表转录路径）")
    check(len(sources) >= 3 and all(s.get("doi") for s in sources.values()),
          "全部来源带 DOI")
    bad = [r.get("item") for r in rows
           if not (r.get("instrument") and r.get("item") and r.get("granularity")
                   and r.get("detail") and r.get("source_key") and r.get("locator"))]
    check(not bad, "每行必填字段齐全（instrument/item/granularity/detail/source/locator）")
    check(not [r.get("item") for r in rows if r.get("granularity") not in GRANULARITY],
          "granularity 只用约定枚举")
    check(not [r.get("item") for r in rows if sources.get(r.get("source_key")) is None],
          "每行 source_key 可解析到带 DOI 的来源")
    check(any(r.get("item") == "PLS_loo_rmse_ground_5m" for r in rows),
          "Wiens 2012 Table 8c 留一 RMSE 基线在册（Phase 1 CV 对账锚点）")
    check(any(r.get("item") == "inflight_cal_target_rmsep" for r in rows),
          "Ytsma 2023 在轨靶精度条目在册（待转录，growth_path 记录）")

    print("共 %d 项断言，%d 项失败" % (TOTAL, FAIL))
    return int(FAIL > 0)


if __name__ == "__main__":
    raise SystemExit(main())
