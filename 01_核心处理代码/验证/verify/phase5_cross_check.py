# -*- coding: utf-8 -*-
"""T3 跨仪器对照护栏：匹配判据、并排负载与不支持口径（探针缓存不在则跳过）。"""

from __future__ import annotations

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = str(_layout.project_layout.APP)
sys.path.insert(0, ROOT)

from instruments import REGISTRY, cross           # noqa: E402

checks = fails = 0


def chk(name, ok, detail=""):
    global checks, fails
    checks += 1
    print("  %s  %s %s" % ("OK  " if ok else "FAIL", name, detail))
    fails += 0 if ok else 1


def main():
    cprobe = os.path.normpath(os.path.join(str(_layout.project_layout.legacy_cache_root()),"phase3_probe","sol01000"))
    aprobe = os.path.normpath(os.path.join(str(_layout.project_layout.legacy_cache_root()),"phase5_probe","apxs_sol999"))
    if not (os.path.isdir(cprobe) and os.path.isdir(aprobe)):
        print("  （探针缓存不在本机，T3 对照护栏跳过——跑过一次 phase5 探针下载后生效）")
        print("\n共 0 项断言，0 项失败")
        return 0

    c = REGISTRY["chemcam"]
    a = REGISTRY["apxs"]
    cobs = [o for o in c.scan_folder(cprobe) if not getattr(o, "passive", False)][0]
    aobs = a.scan_folder(aprobe)[0]

    # 1) 匹配判据：site/drive 主判据命中，目标名不同也匹配
    ms = cross.find_matches(cobs, [aobs])
    chk("site48/drive1570 命中同址判据", len(ms) == 1 and ms[0]["basis"] == "site_drive",
        "（%s）" % (ms[0]["detail"] if ms else "无匹配"))

    # 2) 判据排序：伪造一个同名观测，site/drive 仍排最前
    import types
    fake = types.SimpleNamespace(sol="999", seq="x", site=99, drive=99,
                                  target=cobs.target, points=[])
    ms2 = cross.find_matches(cobs, [fake, aobs])
    chk("同名不敌同址（site/drive 排序最前）",
        ms2[0]["basis"] == "site_drive" and any(m["basis"] == "target" for m in ms2))

    # 3) 负载：元素并排两来源，映射正确
    pl = cross.payload(c, cobs, a, aobs)
    chk("Newland×Wallace_DRT 负载生成（9 点）",
        pl["chemcam"]["target"] == "Newland" and pl["chemcam"]["points_loaded"] == 9)
    fe = next(r for r in pl["rows"] if r["element"] == "Fe")
    chk("Fe 行：LIBS 证据 + APXS FeO 定量并排",
        fe["libs"] and fe["apxs"]["oxide"] == "FeO" and fe["apxs"]["wt"] > 0)
    si = next(r for r in pl["rows"] if r["element"] == "Si")
    chk("高硅靶一致（SiO2>50 wt% 且 LIBS 9/9 点）",
        si["apxs"]["wt"] > 50 and si["libs"]["points"] == "9/9")
    chk("口径说明随负载下发", "不合成数值" in pl["note"])

    print("\n共 %d 项断言，%d 项失败" % (checks, fails))
    return 1 if fails else 0


if __name__ == "__main__":
    raise SystemExit(main())
