# -*- coding: utf-8 -*-
"""R2 实验终版：非负峰形复核能否否决 sol349 Survey 的假阳性——结论是不能。

背景：峰形复核（2026-10-23）发现 sol349 的 14/14 组候选需要负峰形系数，
本实验检验「负系数/非负幅度」能否作为否决判据接入默认判定。

两项检验（2026-10-24 实测）：

  1. 组平均谱非负复核（fitting.solve_nonnegative_profiles，G+L 非负 + 线性底，
     峰高 = G+L 系数和）：sol349 的矿物组并不塌缩——碳酸盐组 @1083 峰高/σ=5.46、
     硫酸盐或辉石组 @991 =3.51，均超过 3σ 门槛。用「非负幅度 <3σ」否决清不掉
     sol349（第一版实验只计 G 分量得出 15/16 否决是口径 bug，已修正）。
  2. 冻结证据对比（峰形复核交付的 review.json）：sol293 的**真检出**
     （碳酸盐 48 点、硫酸盐 22 点）同样有 6/8 组带「候选需要负峰形系数」。
     负系数在真/假检出之间不具判别力——强荧光背景上的真弱拉曼带同样需要
     负分量吸收背景形状。

结论（闸门不通过，按计划不动默认判定）：
  KNOWN_FALSE_POSITIVE 上限维持 60；sol349 的 18 点次假阳性留档不改。
  默认判定已带的 profile_valid / requires_review 标记继续作为提示，
  不升级为否决。

用法：python verify/phase5_r2_experiment.py   （约 1 分钟，读冻结复核 + 跑 sol349）
"""

from __future__ import annotations

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = str(_layout.project_layout.APP)
sys.path.insert(0, ROOT)

REVIEW_293 = str(_layout.project_layout.VERIFY/'fixtures/历史数值基准/fit_review/交付示例_Sol293/review.json')
REVIEW_349 = str(_layout.project_layout.VERIFY/'fixtures/历史数值基准/fit_review/交付示例_Sol349/review.json')
DATA = os.path.join(os.path.dirname(ROOT), "03_下载数据/SHERLOC/练习数据")

checks = fails = 0


def chk(name, ok, detail=""):
    global checks, fails
    checks += 1
    print("  %s  %s %s" % ("OK  " if ok else "FAIL", name, detail))
    fails += 0 if ok else 1


def main():
    # ---- 检验 1：冻结证据——负系数不判别真假 ----
    for tag, path, expect_named in (("sol293 真检出", REVIEW_293, ("碳酸盐", "硫酸盐")),
                                    ("sol349 假阳性", REVIEW_349, None)):
        d = json.load(open(path, encoding="utf-8"))
        groups = d["scans"][0]["groups"]
        neg = sum(1 for g in groups if "候选需要负峰形系数" in (g.get("review_reasons") or []))
        chk("%s 负系数组存在（%d/%d）——非负标记本身不能当否决闸" % (tag, neg, len(groups)),
            neg > 0)
        if expect_named:
            named = [g for g in groups
                     if any(c in (g.get("default_classes") or []) for c in expect_named)]
            named_neg = sum(1 for g in named
                            if "候选需要负峰形系数" in (g.get("review_reasons") or []))
            chk("%s 的真检出组也带负系数（%d/%d）→ 判别失败成立" % (tag, named_neg, len(named)),
                named_neg >= 1)

    # ---- 检验 2：sol349 矿物组在非负幅度复核下存活 ----
    from sherloc import fitting as F
    from sherloc import minerals as M, qc
    from sherloc.organize import scan_folder
    _, scans = scan_folder(DATA)
    hit = [s for s in scans if s.sol == "349" and s.sequence_id == "srlc11420"
           and s.spectra is not None][0]
    wl, r1, _r2, _r3 = hit.read_regions()
    shift = qc.nm_to_raman(wl)
    order = np.argsort(shift)
    x, y = shift[order], np.asarray(r1, float)[:, order]
    res = M.analyse_scan(x, y)
    pc = res["point_classes"]
    survived = []
    for g in res.get("groups", []):
        pts = [p for p in g["pts"] if p < len(pc)]
        cls = {c for p in pts for c in pc[p]
               if c and c not in ("未归类", "未定硅酸盐")}
        if not cls:
            continue
        band = g["band"]
        center, fwhm = float(band["center"]), float(band.get("fwhm", 60.0))
        mean = y[pts].mean(axis=0)
        resid = mean - M.baseline(x, mean)
        win = (x >= center - fwhm - 40.0) & (x <= center + fwhm + 40.0)
        design = np.column_stack((
            np.exp(-0.5 * ((x[win] - center) / (fwhm / 2.355)) ** 2),
            1.0 / (1.0 + ((x[win] - center) / (fwhm / 2.0)) ** 2)))
        coef, _e, model = F.solve_nonnegative_profiles(x[win], resid[win], design)
        r = resid[win] - model
        sig = 1.4826 * np.median(np.abs(r - np.median(r)))
        amp = float(coef[0] + coef[1])
        if amp / max(sig, 1e-12) >= 3.0:
            survived.append((sorted(cls), round(center), len(pts)))
    chk("sol349 矿物组有 %d 组在非负幅度复核下存活（≥3σ）→ 幅度否决也清不掉"
        % len(survived), len(survived) >= 1, str(survived[:4]))

    print("\n共 %d 项断言，%d 项失败" % (checks, fails))
    print("结论：两项判据都不能分离 sol349 假阳性与真弱检出，默认判定保持不变。")
    return 1 if fails else 0


if __name__ == "__main__":
    raise SystemExit(main())
