# -*- coding: utf-8 -*-
"""护栏：**processed 数据不许被二次处理**（2026-09-16）。

起因：用户发现"用 processed 数据时，工具还是给它做了一遍预处理"，问是不是重复了。

查清的结论分两半：

  A. **真的重复了，已修**。判定"PDS 有没有做过激光归一"原来只看文件名第 51 位。
     实测 sol 293 的 RRS 叫 `...srlc11360w108cgzj03`——第 51 位是 `z`，
     按文件名读就是"没归一"，工具于是又归一了一次，还顺手做了二次去宇宙射线
     （改了 522 个点）。可它的 label 里三张表都明明白白写着
     `LASER-NORMALIZED SPECTRA REGION 1/2/3`，PDS 早就归一过了。
     **改法：以 label 表名为准**（`LASER-NORMALIZED` ⇒ 已归一 + 已去宇宙射线）。

  B. **不是重复，必须保留**。
     * 死区裁剪：PDS 不做。处理后产品里首尾那 100 个不感光通道照样躺在里面。
     * 基线扣除：PDS 从不扣基线，荧光背景原样留着；不扣就看不见几波数宽的拉曼带。

这个脚本就是钉住这两条结论的。它会：
  1. 遍历全库，把每份光谱产品的 label 表名和实际跑的步骤打出来（--list）；
  2. 断言 label 写 LASER-NORMALIZED 的产品，在"开关全开"（连 despike 都勾上）
     的情况下也必须 laser=skip / despike=skip，而 deadpixel 和 baseline 照做；
  3. 断言勾上「强制二次清理」之后 despike 才真的执行（这条开关得有用）；
  4. 断言没归一的产品（表名 PROCESS DATA / 原始帧）反过来必须被归一——
     别把"不要重复"修成"该做也不做"。

用法：python verify\\duplicate_steps.py [--list]
"""

from __future__ import annotations

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = str(_layout.project_layout.APP)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from sherloc import preprocess                          # noqa: E402
from sherloc.organize import scan_folder                # noqa: E402

DATA = os.path.join(os.path.dirname(ROOT), "03_下载数据/SHERLOC/练习数据")

# 界面默认勾上的四个（despike 默认不勾）；这里故意连 despike 也勾上，
# 因为要验的正是"就算你勾了，已处理产品也不该被再动一遍"。
ALL_ON = {k: True for k in preprocess.STEP_KEYS}

OK = FAIL = 0


def check(cond, msg):
    global OK, FAIL
    if cond:
        OK += 1
        print("  OK  " + msg)
    else:
        FAIL += 1
        print("  !!  " + msg)


def label_tables(scan):
    if not scan.parts:
        return ""
    return " | ".join((t.name or "") for t in scan.parts[0].tables())


def by_step(log, key):
    for e in log:
        if e["step"] == key:
            return e
    return None


def dump(which):
    """打印全库每份光谱产品的"label 怎么说 / 工具怎么做"。"""
    for area in ("processed", "raw"):
        base = os.path.join(DATA, area)
        if not os.path.isdir(base):
            continue
        print("\n########## %s ##########" % area)
        for sol in sorted(os.listdir(base)):
            folder = os.path.join(base, sol)
            if not os.path.isdir(folder):
                continue
            try:
                _, scans = scan_folder(folder)
            except Exception as exc:                    # noqa: BLE001
                print("  %-12s 读不了：%s" % (sol, exc))
                continue
            for sc in scans:
                sp = sc.spectra
                if sp is None:
                    continue
                names = label_tables(sc)
                tag = ("LASER-NORMALIZED" if "LASER-NORMALIZED" in names.upper()
                       else "PROCESS DATA" if "PROCESS DATA" in names.upper()
                       else "其他")
                try:
                    res = preprocess.run_pipeline(sc, ALL_ON, points=[0])
                    lg = {e["step"]: e["kind"] for e in res["log"]}
                except Exception as exc:                # noqa: BLE001
                    lg = {"err": str(exc)[:40]}
                print("  %-12s %-5s %-28s label=%-16s laser=%-5s despike=%-5s "
                      "dead=%-5s base=%s"
                      % (sol, sp.product_type, sp.name[-28:], tag,
                         lg.get("laser", "-"), lg.get("despike", "-"),
                         lg.get("deadpixel", "-"), lg.get("baseline", "-")))


def main(argv):
    if "--list" in argv:
        dump(argv)
        return 0

    for area in ("processed", "raw"):
        base = os.path.join(DATA, area)
        if not os.path.isdir(base):
            continue
        print("\n=== %s ===" % area)
        for sol in sorted(os.listdir(base)):
            folder = os.path.join(base, sol)
            if not os.path.isdir(folder):
                continue
            _, scans = scan_folder(folder)
            for sc in scans:
                sp = sc.spectra
                if sp is None:
                    continue
                names = label_tables(sc)
                rdr = "LASER-NORMALIZED" in names.upper()
                tag = "%s %s %s" % (sol, sp.product_type, sp.name[-24:])

                res = preprocess.run_pipeline(sc, ALL_ON, points=[0])
                laser = by_step(res["log"], "laser")
                desp = by_step(res["log"], "despike")
                dead = by_step(res["log"], "deadpixel")
                base_ = by_step(res["log"], "baseline")

                if rdr:
                    check(laser and laser["kind"] == "skip",
                          "%s：已归一 ⇒ 激光归一步骤必须 skip" % tag)
                    check(desp and desp["kind"] == "skip",
                          "%s：已去宇宙射线 ⇒ despike 必须 skip" % tag)
                    # 勾上「强制二次清理」就该真的动手（开关得有用）
                    forced = dict(ALL_ON, despike_force=True)
                    res2 = preprocess.run_pipeline(sc, forced, points=[0])
                    d2 = by_step(res2["log"], "despike")
                    check(d2 and d2["kind"] == "done",
                          "%s：勾「强制二次清理」后 despike 必须真的执行" % tag)
                else:
                    check(laser and laser["kind"] in ("done", "warn"),
                          "%s：没归一的产品必须做激光归一（别修过头）" % tag)

                check(dead and dead["kind"] == "done",
                      "%s：死区裁剪必须照做（PDS 不做这一步）" % tag)
                # 基线要有波长表才做得成；原始帧（ECA/ERA）没波长解，
                # 这一步本来就跑不了，不是被"跳过重复"跳掉的。
                if sp.product_type in ("RRS", "RCS", "RRA", "RCA"):
                    check(base_ and base_["kind"] == "done",
                          "%s：基线扣除必须照做（PDS 从不扣基线）" % tag)

    print("\n" + "=" * 62)
    print("共 %d 项断言，%d 项失败" % (OK + FAIL, FAIL))
    return 1 if FAIL else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
