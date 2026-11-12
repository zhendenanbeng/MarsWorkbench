# -*- coding: utf-8 -*-
"""R3：KNOWN_GAP 四条缺口的逐条复核断言（结论=四条均不可自动补，各有依据）。

诊断（2026-10-24，口径与 regression.class_summary 对齐）：

  sol293 未定硅酸盐（论文 3 点 @~1050）
    → 组级判类口径差异：生产管线（qc.analyse_spectra）逐点与组级均无
      未定硅酸盐——~1050 的弱硅酸盐带被并入碳酸盐主导组或落为未归类
      （9 点）。强行拆组把未定硅酸盐从碳酸盐组里剥出来，会威胁已通过
      的碳酸盐金标准区间（48 ∈ [26,62]），不做。
  sol293 高氯酸盐或磷酸盐（论文 1 点 @952）
    → 类别命名差异：工具检出同族 5 点判为「高氯酸盐/氯酸盐 @939」
      （更细的族名）。不强行改类名——那是把分类往论文词汇表上凑。
  sol257 辉石（论文 9 点 @993）
    → 工具的模糊类「硫酸盐或辉石（缺次级峰）」16 点 @987 有意保留：
      次级峰低于检测限时无法自动拆分，论文靠目视判读。已有 PIXL
      化学仲裁通道（pixl.arbitrate）可对模糊类给出化学倾向。
  sol257 硫酸盐（论文 1 点）→ 同上，在模糊组内。

本护栏断言的是**诊断事实的稳定性**：若将来分类器改动导致这些检出结构
变化（模糊类消失、氯酸盐点丢失等），此护栏会红，提示重新审视 KNOWN_GAP。

用法：python verify/phase5_r3_gaps.py
"""

from __future__ import annotations

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout

import os
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = str(_layout.project_layout.APP)
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

from sherloc import qc                                   # noqa: E402
from sherloc.organize import scan_folder                 # noqa: E402

DATA = os.path.join(os.path.dirname(ROOT), "03_下载数据/SHERLOC/练习数据")

checks = fails = 0


def chk(name, ok, detail=""):
    global checks, fails
    checks += 1
    print("  %s  %s %s" % ("OK  " if ok else "FAIL", name, detail))
    fails += 0 if ok else 1


def main():
    _p, scans = scan_folder(DATA)

    def get(sol, seq):
        return [s for s in scans
                if s.sol == sol and s.sequence_id == seq and s.spectra is not None]

    # ---- sol293：两条缺口的诊断事实 ----
    s293 = get("293", "srlc11360")[0]
    wl, r1, r2, r3 = s293.read_regions()
    res293 = qc.analyse_spectra(wl, r1, r2, r3)
    point_cnt = Counter()
    for d in res293.get("point_detail") or []:
        for b in d["classes"]:
            point_cnt[b["class"]] += 1
    chk("sol293 生产管线确无未定硅酸盐（缺口如实），论文 3 点落在未归类/碳酸盐组",
        point_cnt.get("未定硅酸盐", 0) == 0 and point_cnt.get("未归类", 0) >= 3,
        "（未归类 %d 点；拆组将威胁碳酸盐金标准 48∈[26,62]）"
        % point_cnt.get("未归类", 0))
    chk("sol293 高氯酸盐族已检出（命名为高氯酸盐/氯酸盐，族名差异非缺口）",
        point_cnt.get("高氯酸盐/氯酸盐", 0) >= 1,
        "（%d 点次）" % point_cnt.get("高氯酸盐/氯酸盐", 0))

    # ---- sol257：两条缺口都落在有意设计的模糊类 ----
    s257 = get("257", "srlc11360")[0]
    wl, r1, r2, r3 = s257.read_regions()
    res257 = qc.analyse_spectra(wl, r1, r2, r3)
    fuzzy = 0
    fuzzy_pos = None
    for g in res257.get("mineral_groups") or []:
        if g.get("classes") and "硫酸盐或辉石" in g["classes"][0]["class"]:
            fuzzy += g["n_points"]
            fuzzy_pos = g["band"]["center"]
    chk("sol257 辉石/硫酸盐缺口由模糊类承载（缺次级峰不可自动拆分）",
        fuzzy >= 10 and fuzzy_pos is not None and 970 <= fuzzy_pos <= 1000,
        "（模糊类 %d 点 @%.0f，论文辉石 9 @993 + 硫酸盐 1）"
        % (fuzzy, fuzzy_pos or -1))
    # PIXL 仲裁通道存在且可对模糊类表态（不依赖数据在线，仅检函数契约）
    from sherloc import pixl
    chk("模糊类的 PIXL 化学仲裁通道存在（人工可选的第二证据）",
        hasattr(pixl, "arbitrate"))

    print("\n共 %d 项断言，%d 项失败" % (checks, fails))
    print("结论：四条 KNOWN_GAP 均为口径/命名/判读差异，不可自动补，依据见")
    print("04_使用文档/第五期验证_2026-10-24/R3_缺口复核结论.md")
    return 1 if fails else 0


if __name__ == "__main__":
    raise SystemExit(main())
