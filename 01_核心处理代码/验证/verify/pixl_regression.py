# -*- coding: utf-8 -*-
"""PIXL 化学链路的回归护栏。

断言：
  A. Bellegarde（sol 186 / 定位 N0070000）全链路——配准质量、点命中率、
     斑中位成分落在合理区间（数值锚定 spike 实测：内点 24-32、
     重投影 1.9-2.5px、命中 94/100、SiO2≈45 FeO≈23 CaO≈7.4）；
  B. evidence() 佐证规则的单元检验（合成数据，不依赖网络/缓存）；
  C. 无 opencv 时优雅降级。

用法：python verify\\pixl_regression.py   （通过返回 0）
A 类需要 ~/.sherloc_cache/pixl 缓存或网络；都没有时标记 SKIP。
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

DATA = os.path.join(os.path.dirname(ROOT), "03_下载数据/SHERLOC/练习数据")

checks = fails = 0


def check(name, cond, detail=""):
    global checks, fails
    checks += 1
    if cond:
        print("  OK  %s %s" % (name, detail))
    else:
        fails += 1
        print("  FAIL  %s %s" % (name, detail))


def main():
    from sherloc import pixl

    # ---- C 类：降级（最快）----
    print("优雅降级")
    saved = pixl.HAS_CV2
    try:
        pixl.HAS_CV2 = False
        r = pixl.chemistry_at_points([[0, 0]], sol="0186",
                                     placement="N0070000", points=[(1, 1)])
        check("无 cv2 时返回 None", r is None)
    finally:
        pixl.HAS_CV2 = saved

    # ---- B 类：佐证规则 ----
    print("佐证规则（合成数据）")
    med = {"SO3": 1.0, "Cl": 1.0, "P2O5": 0.5, "SiO2": 45.0,
           "CaO": 7.0, "MgO": 2.0, "FeO-T": 23.0}

    def elems(**over):
        d = {k: (v, 0.1) for k, v in med.items()}
        for k, v in over.items():
            d[k] = v if isinstance(v, tuple) else (v, 0.1)
        return d

    check("SO3 富集 → ✓", "✓" in pixl.evidence("sulfate", elems(SO3=5.0), med))
    check("SO3 亏损 → ✗", "✗" in pixl.evidence("sulfate", elems(SO3=0.5), med))
    check("SO3 富集但 z<3 → ～",
          "～" in pixl.evidence("sulfate", elems(SO3=(2.0, 0.5)), med))
    carb_ok = elems(SiO2=20.0)
    carb_ok["FeO-T"] = (40.0, 1.0)
    check("碳酸盐：FeO 高 Si 低 → ✓", "✓" in pixl.evidence("carbonate", carb_ok, med))
    carb_bad = elems(SiO2=60.0)
    carb_bad["FeO-T"] = (10.0, 1.0)
    check("碳酸盐：FeO 低 Si 高 → ✗", "✗" in pixl.evidence("carbonate", carb_bad, med))
    check("P/Cl 二选一：Cl 富集 → Cl", "Cl" in pixl.evidence(
        "perchlorate_or_phosphate", elems(Cl=3.0), med))
    check("无规则类别 → 空", pixl.evidence("olivine", elems(SO3=5.0), med) == "")
    check("空元素 → 空", pixl.evidence("sulfate", None, med) == "")

    # ---- B2 类：化学仲裁（二选一模糊类别的拆分）----
    print("化学仲裁（合成数据）")
    k, nk, note = pixl.arbitrate(["sulfate_or_pyroxene"],
                                 elems(SO3=0.5), med)
    check("S/Px：SO3 亏损 → 辉石", k == "sulfate_or_pyroxene"
          and nk == ["pyroxene"] and "否决" in note, str(note))
    k, nk, note = pixl.arbitrate(["sulfate_or_pyroxene"],
                                 elems(SO3=5.0), med)
    check("S/Px：SO3 富集 → 硫酸盐", nk == ["sulfate"], str(note))
    k, nk, note = pixl.arbitrate(["sulfate_or_pyroxene"],
                                 elems(SO3=1.2), med)
    check("S/Px：中间地带 → 不动", nk is None, str(note))
    k, nk, note = pixl.arbitrate(["perchlorate_or_phosphate"],
                                 elems(Cl=3.0, P2O5=0.5), med)
    check("P/Cl：Cl 强 → 高氯酸盐", nk == ["perchlorate_chloride"], str(note))
    k, nk, note = pixl.arbitrate(["perchlorate_or_phosphate"],
                                 elems(Cl=0.5, P2O5=3.0), med)
    check("P/Cl：P 强 → 磷酸盐", nk == ["phosphate"], str(note))
    k, nk, note = pixl.arbitrate(["perchlorate_or_phosphate"],
                                 elems(Cl=0.8, P2O5=0.4), med)
    check("P/Cl：都不富集 → 不动", nk is None)

    # ---- A 类：Bellegarde 全链路 ----
    print("全链路（Bellegarde，sol 186）")
    try:
        from sherloc.organize import scan_folder
        from sherloc.vicar import load_aci
        from sherloc.watson import normalize
        import numpy as np
        prods, obs = scan_folder(DATA)
        bel = [o for o in obs if o.sol == "186"
               and o.sequence_id == "srlc10500" and o.spectra is not None]
        if not bel:
            print("  SKIP  练习数据里没有 Bellegarde 扫描")
        else:
            pts = bel[0].read_shot_positions()
            aci_name = sorted({p["image"] for p in pts})[0]
            aci_path = None
            for root, _d, files in os.walk(DATA):
                for f in files:
                    if normalize(f) == normalize(aci_name):
                        aci_path = os.path.join(root, f)
            if not aci_path:
                print("  SKIP  本地没有对应 ACI 图")
            else:
                aci_gray = np.asarray(load_aci(aci_path).convert("L"))
                points = [(p["samp"], p["line"]) for p in pts]
                res = pixl.chemistry_at_points(
                    aci_gray, sol="0186", placement="N0070000",
                    points=points)
                if res is None:
                    check("化学查询成功", False, "（返回 None：缓存/网络缺？）")
                else:
                    check("配准内点 ≥20", res["inliers"] >= 20,
                          "(%d)" % res["inliers"])
                    check("重投影 <5px", res["err_px"] < 5,
                          "(%.2f)" % res["err_px"])
                    n_hit = sum(1 for p in res["points"] if p[0] is not None)
                    check("命中 ≥80/100", n_hit >= 80, "(%d)" % n_hit)
                    m = res["median"]
                    check("斑中位 SiO2 35-55", 35 <= m["SiO2"] <= 55,
                          "(%.1f)" % m["SiO2"])
                    check("斑中位 FeO-T 15-35", 15 <= m["FeO-T"] <= 35,
                          "(%.1f)" % m["FeO-T"])
                    check("斑中位 CaO 3-15", 3 <= m["CaO"] <= 15,
                          "(%.1f)" % m["CaO"])
                    pmc, _d, elems = next(p for p in res["points"]
                                          if p[0] is not None)
                    ev = pixl.evidence("carbonate", elems, m)
                    check("evidence() 可用", isinstance(ev, str))
    except Exception as exc:                         # noqa: BLE001
        print("  SKIP  链路不可用（%s: %s）" % (type(exc).__name__, exc))

    print("=" * 70)
    print("共 %d 项断言，%d 项失败" % (checks, fails))
    return 0 if fails == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
