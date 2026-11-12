# -*- coding: utf-8 -*-
"""集成路径护栏：盖住 2026-10-08 实测抓到、但旧护栏没覆盖的那批 bug。

旧的 regression / duplicate_steps 两条护栏盯的是**判类数学**，全部通过也
拦不住下面这些——它们坏在报告渲染和新集成的路径上：

* `make_mineral_map.render` 的裁剪代码块被复制粘贴了两份：上色路径
  UnboundLocalError 必崩、灰度路径裁两次（出过 7977×560 的失控图）。
  09-18 生效的上色缓存存在后没人再点过按钮，所以一直没被发现。
* `pixl.arbitrate` 返回了 minerals 里不存在的类别键（perchlorate_chloride），
  点阵透明、图例英文。这类"键漂移"用一个断言就能整类堵死。
* 导出 PNG 的文件名只用 sol+序列号，同号多次扫描互相覆盖
  （实测导 999 张磁盘只剩 631 张）。
* 前端导出请求不带 steps/mode，界面勾选形同虚设。

用法（一般不用单跑，跑总检即可）：

    python verify\\integration_guards.py
"""

from __future__ import annotations

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout

import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = str(_layout.project_layout.APP)
sys.path.insert(0, ROOT)

DATA = os.environ.get("SHERLOC_DATA",
                      os.path.join(os.path.dirname(ROOT), "03_下载数据/SHERLOC/练习数据"))

_n = [0, 0]      # [断言数, 失败数]


def check(cond, msg):
    _n[0] += 1
    if not cond:
        _n[1] += 1
        print("  FAIL %s" % msg)
        return False
    return True


def main():
    from sherloc import minerals, pixl

    # ---------- 1. 类别键合法性：所有能写进报告的键必须真实存在 ----------
    # （perchlorate_chloride 那次就是 pixl 和 minerals 的键漂移，
    #   点阵出透明格、图例显示英文原键）
    for k in pixl._CHEM_RULES:
        check(k in minerals.LABELS,
              "_CHEM_RULES 的键 %r 不在 minerals.LABELS 里" % k)
    for k in minerals.LABELS:
        check(k in minerals.CLASS_COLORS, "类别 %r 没有配色" % k)
        check(k in minerals.SHORT, "类别 %r 没有缩写" % k)

    # arbitrate 的四个分支逐个触发，返回的新键必须都真实存在
    med = {"SO3": 2.0, "Cl": 1.0, "P2O5": 2.0}
    cases = [
        # (class_keys, elems, 预期返回非空)
        (["sulfate_or_pyroxene"], {"SO3": (1.0, 0.1)}, "SO3 否决 → 辉石"),
        (["sulfate_or_pyroxene"], {"SO3": (4.0, 0.1)}, "SO3 确认 → 硫酸盐"),
        (["perchlorate_or_phosphate"], {"Cl": (3.0, 0.1), "P2O5": (2.0, 0.1)},
         "Cl 富集 → 高氯酸盐"),
        (["perchlorate_or_phosphate"], {"Cl": (1.0, 0.1), "P2O5": (4.0, 0.1)},
         "P2O5 富集 → 磷酸盐"),
    ]
    for keys, elems, what in cases:
        old, new, note = pixl.arbitrate(keys, elems, med)
        ok = check(new is not None, "arbitrate 应触发分支：%s" % what)
        if ok:
            for k in new:
                check(k in minerals.LABELS,
                      "arbitrate 返回的键 %r 不在 minerals.LABELS 里（%s）"
                      % (k, what))

    # ---------- 2. 导出文件名：同 sol+序列号的多次扫描不许撞名 ----------
    if os.path.isdir(DATA):
        import view
        scans = view.scan_listing(DATA)[1]
        if scans:
            names = [view._export_name(sc, 0) for sc in scans]
            check(len(names) == len(set(names)),
                  "导出文件名有撞名：同 sol+序列号的扫描会互相覆盖")
            pairs = [(sc.sol, sc.sequence_id) for sc in scans]
            # 数据里确实存在同号多次扫描（sol 141 ×3）才真正考验命名
            if len(set(pairs)) < len(pairs):
                dup = [sc for sc in scans if pairs.count(
                    (sc.sol, sc.sequence_id)) > 1]
                a, b = dup[0], dup[1]
                check(view._export_name(a, 0) != view._export_name(b, 0),
                      "同 sol+序列号的两次扫描导出名相同（%s）"
                      % view._export_name(a, 0))
        else:
            print("  （数据文件夹没有扫描，导出命名检查跳过）")

    # ---------- 3. 地图冒烟：上色缓存存在时 render 不许崩 ----------
    # （P0：裁剪块重复导致 UnboundLocalError；只要 colorized/ 里有一张
    #   缓存就会触发。跑 sol 293 一条观测，几十秒。）
    if os.path.isdir(DATA) and os.path.isdir(
            os.path.join(DATA, "_aci", "colorized")):
        import make_mineral_map
        tmp = tempfile.mkdtemp(prefix="sherloc_mapguard_")
        try:
            res = make_mineral_map.run(DATA, out_dir=tmp, sol="293",
                                       seq="srlc11360", verbose=False)
            check(not res.get("errors"),
                  "地图生成报错：%s" % (res.get("errors") or "")[:200])
            check(res.get("count", 0) >= 1, "sol 293 应至少出 1 次观测的图")
            from PIL import Image
            for p in (res.get("maps") or []):
                w, h = Image.open(p).size
                check(max(w, h) <= 2500,
                      "地图尺寸失控：%s 是 %dx%d（裁剪/缩放坏了）"
                      % (os.path.basename(p), w, h))
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
    else:
        print("  （没有上色缓存或数据文件夹，地图冒烟跳过）")

    # ---------- 4. 前端契约：导出请求必须带 steps/mode，且可取消 ----------
    src = open(os.path.join(ROOT, "view.py"), encoding="utf-8").read()
    check("&mode='+$('mode').value+'&steps='" in src,
          "导出请求没带界面的 mode/steps（所见非所得回归了）")
    check("/api/export/cancel" in src, "导出没有取消入口")
    check("var CLASS_COLORS=__CLASS_COLORS__;" in src,
          "前端类别配色必须是注入的单一来源，不许再抄一份硬编码")

    print()
    print("共 %d 项断言，%d 项失败" % (_n[0], _n[1]))
    return 1 if _n[1] else 0


if __name__ == "__main__":
    raise SystemExit(main())
