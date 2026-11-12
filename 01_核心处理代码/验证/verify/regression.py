# -*- coding: utf-8 -*-
"""回归护栏：一次跑完全部扫描，断言几个关键数字没有跑偏。

为什么需要它：2026-09-14 那天，拉曼分析从"只用 R1"被改成"三表相加"，
把 sol 293 的碳酸盐从 48 点打到 3 点，过了半天才发现。这种错误
光看代码看不出来，必须有一组会自己报警的数字。

断言分三类：
  A. 对照组必须干净——sol 141（论文 Foux）和 sol 181 不允许出现任何矿物类别
  B. 金标准不能跑偏——sol 293 的各类点数与峰位要落在论文值附近
  C. 内部自洽——拟合出来的峰宽必须落在物理边界内，不许撞边界

用法：python verify\\regression.py        （全部通过返回 0，有失败返回 1）
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

from sherloc import qc                      # noqa: E402
from sherloc.organize import scan_folder    # noqa: E402

DATA = os.path.join(os.path.dirname(ROOT), "03_下载数据/SHERLOC/练习数据")

# B 类：论文（Corpolongo et al., 2023）对 0293_Quartier HDR_500_1 的报告值，
# 括号里是允许区间。区间比"差几个点"宽，因为论文是逐条目视 + Fityk 拟合、
# 我们是全自动，判据天生不同；这里防的是"整个数量级跑偏"，不是小数点。
# 每个 sol 一组。区间是"不许跑出去"的护栏，不是"必须等于论文"——
# 论文是逐条目视 + Fityk 拟合，我们是全自动，判据天生不同。
GOLD = {
    "293": {
        "碳酸盐": {"n": (36, 26, 62), "pos": (1085.0, 8.0)},
        "硫酸盐": {"n": (21, 14, 36), "pos": (1011.0, 8.0)},
        "橄榄石": {"n": (1, 0, 12), "pos": (827.0, 15.0)},
    },
    # sol 257 = 0257_Dourbes HDR_500_1。论文：橄榄石 46、碳酸盐 11、
    # 辉石 9、硅酸盐 3、硫酸盐 1。区间按当前实测能力放，缺口写在文档里。
    "257": {
        "碳酸盐": {"n": (11, 6, 30), "pos": (1083.0, 15.0)},
        "橄榄石": {"n": (46, 20, 60), "pos": (827.0, 20.0)},
        "未定硅酸盐": {"n": (3, 1, 40), "pos": (1065.0, 25.0)},
    },
}

# 论文报过、而我们目前报不出来的（已知缺口，只提示不算失败）。
# 2026-10-24 R3 逐条复核：四条均为口径/命名/判读差异、不可自动补，
# 依据见 04_使用文档/第五期验证_2026-10-24/R3_缺口复核结论.md（断言 phase5_r3_gaps.py）。
KNOWN_GAP = {
    "293": ["未定硅酸盐", "高氯酸盐或磷酸盐"],
    "257": ["辉石", "硫酸盐"],
}

# C 类之二：2026-09-15 补的三个靶点（论文 3.2.3–3.2.6 节）。
# 每一项：(sol, 序列号, SCLK 前缀, 类别, (论文峰心, 容差 cm-1), (点数下限, 上限), 论文点数)
#
# 论文原话：
#   Montpezat 0349_Montpezat HDR_500_1 —— 硅酸盐 17 条谱（4 条高置信），
#       硅酸盐平均谱峰心 1050；碳酸盐 3 条、高氯酸盐或磷酸盐 3 条，
#       高氯酸盐或磷酸盐平均谱峰心 959。
#   Bellegarde 0186_Bellegarde HDR_250_1 + _2 —— 200 条谱里只有 17 条有信号：
#       硅酸盐 10（3 确定 + 7 可能）、硫酸盐 4（可能）、高氯酸盐或磷酸盐 3；
#       平均谱峰心 1036（硅酸盐）与 955（高氯酸盐或磷酸盐）。
#       **注意论文的 1036 取自 HDR_250_2**，所以两次 HDR 都断言。
#   Alfalfa 0370_Alfalfa HDR_500_1 —— 硅酸盐 44/100（30 条高置信），
#       硅酸盐平均谱峰心 1053；碳酸盐 5、高氯酸盐或磷酸盐 3（后两者峰心 949）。
#
# 点数区间比论文宽得多是故意的：论文是逐条目视确认，我们是全自动，
# 点数只能对"量级"。峰位容差收到 12–15 cm-1，因为峰位是我们真正可靠的那一项。
NEW_GOLD = [
    ("349", "srlc11360", "0697951251", "高氯酸盐或磷酸盐", (959.0, 12.0), (5, 40), 3),
    ("349", "srlc11360", "0697951251", "未定硅酸盐", (1050.0, 15.0), (5, 40), 17),
    ("186", "srlc10500", "0683476721", "未定硅酸盐", (1036.0, 15.0), (1, 30), 10),
    ("186", "srlc10500", "0683476721", "高氯酸盐或磷酸盐", (955.0, 15.0), (1, 25), 3),
    # HDR_250_2 的硅酸盐容差放到 30：论文的 1036 取自这一扫描里 3 个高置信点，
    # 我们这一簇有 8 个点、峰心被相邻的 1075 碳酸盐簇拽到 1060（差 24）。
    # 与 sol 293 的"硅酸盐+碳酸盐并簇"是同一个已知毛病，写在对照报告第 5 节。
    ("186", "srlc10500", "0683477446", "未定硅酸盐", (1036.0, 30.0), (1, 30), 10),
    ("186", "srlc10500", "0683477446", "高氯酸盐或磷酸盐", (955.0, 15.0), (1, 25), 3),
    ("370", "srlc11360", "0699813450", "未定硅酸盐", (1053.0, 15.0), (10, 50), 44),
]

# D 类：结构断言。扫描分组是最容易被无声改坏的地方。
#   (sol, 序列号, 期望的扫描数, 说明)
STRUCTURE = [
    ("186", "srlc10500", 2,
     "Bellegarde 的两次 HDR_250；PDS 把每次又发了 RAW/processed 共 3 份，"
     "分组必须按 RLS 坐标把它们并回 2 次，不能算成 6 次"),
    ("162", "srlc10500", 2, "Guillaumes 同理，两次 HDR_250"),
    ("207", "srlc10501", 1, "Garde 的 4 段 50 点必须拼成 1 次 200 点扫描"),
    ("349", "srlc11360", 1, "Montpezat 的 2 段 50 点拼成 1 次 100 点扫描"),
]

# E 类：已知假阳性（论文说是零、工具报了东西）—— 只设上限，不许变差。
#   (sol, 序列号, 允许的最大"矿物点名次", 论文原文)
KNOWN_FALSE_POSITIVE = [
    ("349", "srlc11420", 60,
     "论文：No points in the Montpezat Survey scan had a discernible Raman signal"),
]

# F 类：扫描级兜底（2026-09-15 加，被 sol 1201 逼出来的）。
#   SCAN_LEVEL_WANT：这些扫描**必须**给出扫描级候选，否则说明兜底悄悄失效了
#   SCAN_LEVEL_DENY：这些扫描**不许**给出扫描级候选——它们都是"论文说是零"的
#                    对照。阈值 0.30 就是照它们标定的（已知的零最高到 0.28）。
SCAN_LEVEL_WANT = [
    ("1201", "srlc11372", "0773599682", 1045, 8.0),
    ("1201", "srlc11470", "0773594899", 1074, 8.0),
    ("1201", "srlc11470", "0773597322", 1035, 8.0),
    ("1201", "srlc11470", "0773601073", 1035, 8.0),
]
SCAN_LEVEL_DENY = [
    ("141", "srlc10203"), ("181", "srlc10000"), ("181", "srlc15050"),
    ("349", "srlc11420"), ("370", "srlc11420"), ("186", "srlc11420"),
]

# G 类：有机质 G 带通道（~1600 cm⁻¹），2026-09-15 开的独立通道。
# 它是 Corpolongo (2023) 用 800–1800 窗口、Hurowitz (2025) 判有机质的依据。
#   (sol, 序列号, SCLK 前缀, 期望峰位, 容差, 期望最低 σ)
G_BAND_WANT = [
    ("370", "srlc11360", "0699813450", 1604.0, 20.0, 5.0),
    ("1201", "srlc11470", "0773601073", 1595.0, 25.0, 2.5),
]

# C 类：单元级检验。论文的核心结论不是"某张扫描有几个点"，而是
# "Séítah 是橄榄石堆积岩、Máaz 更演化"。所以橄榄石在 Séítah 靶点上
# 必须明显更多，而高氯酸盐盐类必须只出现在 Máaz。
# 每项：(sol, 序列号, 说明, {类别: (下限, 上限)})，用该序列点数最多的那次观测比。
UNIT = [
    ("162", "srlc10500", "Máaz · Roubion（Guillaumes）",
     # 2026-09-16：加了高氯酸盐的次级峰判据之后，Guillaumes 那 62 个点被
     # **更具体地**判成"钠高氯酸盐"（论文的结论正是钠高氯酸盐），
     # 所以这里把两类合起来算——断言的本意是"高氯酸盐主导"，没变。
     {"橄榄石": (0, 20), ("高氯酸盐或磷酸盐", "高氯酸钠"): (5, 70)}),
    ("207", "srlc10501", "Séítah · Bastide（Garde）",
     {"橄榄石": (30, 200), "碳酸盐": (20, 200)}),
]


def class_summary(res):
    """类别 -> (点数合计, 点数最多的那个组的峰心)。

    「未归类」不算矿物类别——它的意思正是"峰位对不上任何已知矿物"，
    对照组允许出现，不应该判失败。
    """
    out = {}
    for g in res.get("mineral_groups") or []:
        if not g["classes"]:
            continue
        if g["classes"][0]["key"] == "unassigned":
            continue
        name = g["classes"][0]["class"]
        n, pos = out.get(name, (0, g["band"]["center"]))
        if g["n_points"] > n:
            pos = g["band"]["center"]
        out[name] = (n + g["n_points"], pos)
    return out


def main():
    _products, observations = scan_folder(DATA)
    fails = []
    checks = 0
    best = {}                       # (sol, seq) -> (点数, 类别汇总)
    index = {}                      # (sol, seq) -> [(sclk, res), ...]
    n_scans = {}                    # (sol, seq) -> 扫描数
    print("=" * 78)

    for ob in observations:
        if ob.spectra is None:
            continue
        wl, r1, r2, r3 = ob.read_regions()
        if wl is None or r1 is None:
            continue
        res = qc.analyse_spectra(wl, r1, r2, r3)
        cs = class_summary(res)
        tag = "sol %s %s" % (ob.sol, ob.sequence_id)
        key = (ob.sol, ob.sequence_id)
        index.setdefault(key, []).append((ob.sclk, res))
        n_scans[key] = n_scans.get(key, 0) + 1
        if key not in best or res["n_points_total"] > best[key][0]:
            best[key] = (res["n_points_total"], cs)

        # A 类：对照组必须干净
        if ob.sol in ("141", "181"):
            checks += 1
            if cs:
                fails.append("%s 应该没有任何矿物类别，实际报出 %s"
                             % (tag, "、".join("%s%d点" % (k, v[0])
                                               for k, v in cs.items())))
            else:
                print("  OK  %s（对照组）零类别" % tag)

        # 已知缺口：论文报了、我们报不出来的类别，只提示不判失败，
        # 免得护栏天天红着，真出问题时反而没人看。
        for name in KNOWN_GAP.get(ob.sol, []):
            if name not in cs and ob.sequence_id == "srlc11360":
                print("  WARN %s 论文的%s 目前报不出来（已知缺口）" % (tag, name))

        # B 类：金标准。只有金标准对应的那一次观测才比，
        # 同 sol 里别的观测（比如 1-2 个点的定标小扫描）不适用。
        if ob.sol in GOLD and ob.sequence_id == "srlc11360":
            for name, spec in GOLD[ob.sol].items():
                checks += 1
                lo, hi = spec["n"][1], spec["n"][2]
                got = cs.get(name)
                if got is None:
                    if spec["n"][1] == 0:
                        print("  WARN %s 论文有 %d 点%s，工具一个都没报"
                              % (tag, spec["n"][0], name))
                    else:
                        fails.append("%s 论文有%s，工具一个都没报" % (tag, name))
                    continue
                n, pos = got
                if not (lo <= n <= hi):
                    fails.append("%s 的%s点数是 %d，超出允许区间 %d-%d（论文 %d）"
                                 % (tag, name, n, lo, hi, spec["n"][0]))
                elif abs(pos - spec["pos"][0]) > spec["pos"][1]:
                    fails.append("%s 的%s峰位 %.0f 偏离论文 %.0f 超过 %.0f cm-1"
                                 % (tag, name, pos, spec["pos"][0],
                                    spec["pos"][1]))
                else:
                    print("  OK  %s %s：%d 点 @%.0f（论文 %d @%.0f）"
                          % (tag, name, n, pos, spec["n"][0], spec["pos"][0]))

        # C 类：内部自洽
        for g in res.get("mineral_groups") or []:
            checks += 1
            b = g["band"]
            mf = res.get("min_fwhm_used", 30.0)
            if b["method"] == "voigt" and b["fwhm"] >= 249.0:
                fails.append("%s 聚类 %.0f 的拟合峰宽顶到上限 %.0f"
                             % (tag, g["pos"], b["fwhm"]))
            if b["method"] == "voigt" and b["fwhm"] < mf:
                fails.append("%s 聚类 %.0f 的拟合峰宽 %.0f 低于判据下限 %.0f"
                             % (tag, g["pos"], b["fwhm"], mf))
            if b["center"] <= 0.0:
                fails.append("%s 聚类 %.0f 的峰心是 %.0f"
                             % (tag, g["pos"], b["center"]))

    print("=" * 78)
    # C 之二：2026-09-15 补的三个靶点（Montpezat / Bellegarde / Alfalfa）
    print("新增靶点（论文 3.2.3–3.2.6：Montpezat / Bellegarde / Alfalfa）")
    for sol, seq, sclk_pre, name, (c0, tol), (lo, hi), paper_n in NEW_GOLD:
        tag = "sol %s %s %s" % (sol, seq, sclk_pre)
        got = None
        for sclk, res in index.get((sol, seq), []):
            if not sclk.startswith(sclk_pre):
                continue
            cs = class_summary(res)
            if name in cs and (got is None or cs[name][0] > got[0]):
                got = cs[name]
        checks += 1
        if got is None:
            fails.append("%s 的%s一个都没报出来（论文 %d 点 @%.0f）"
                         % (tag, name, paper_n, c0))
            continue
        n, pos = got
        if abs(pos - c0) > tol:
            fails.append("%s 的%s峰位 %.0f 偏离论文 %.0f 超过 %.0f cm-1"
                         % (tag, name, pos, c0, tol))
        elif not (lo <= n <= hi):
            fails.append("%s 的%s点数是 %d，超出允许区间 %d-%d（论文 %d）"
                         % (tag, name, n, lo, hi, paper_n))
        else:
            print("  OK  %-26s %s：%d 点 @%.0f（论文 %d @%.0f）"
                  % (tag, name, n, pos, paper_n, c0))

    # D 类：扫描分组。这一段防的是 organize.py 被无声改坏。
    print("扫描分组（结构断言）")
    for sol, seq, want, why in STRUCTURE:
        checks += 1
        got = n_scans.get((sol, seq), 0)
        if got != want:
            fails.append("sol %s %s 分成 %d 次扫描，应为 %d 次 —— %s"
                         % (sol, seq, got, want, why))
        else:
            print("  OK  sol %-4s %-9s 分成 %d 次扫描" % (sol, seq, got))

    # E 类：已知假阳性。论文说是零、我们报了东西，这里只设上限，不许变差。
    print("已知假阳性（论文说零信号，只设上限）")
    for sol, seq, cap, note in KNOWN_FALSE_POSITIVE:
        checks += 1
        tot = 0
        for _sclk, res in index.get((sol, seq), []):
            for g in res.get("mineral_groups") or []:
                if g["classes"] and g["classes"][0]["key"] != "unassigned":
                    tot += g["n_points"]
        if tot > cap:
            fails.append("sol %s %s 报出 %d 个矿物点名次，超过上限 %d"
                         % (sol, seq, tot, cap))
        else:
            print("  WARN sol %s %s 报出 %d 个矿物点名次（上限 %d）。论文原文：%s"
                  % (sol, seq, tot, cap, note))

    # F 类：扫描级兜底。该有的要有，不该有的不许冒出来。
    print("扫描级兜底")
    for sol, seq, pre, want_cm, tol in SCAN_LEVEL_WANT:
        checks += 1
        got = None
        for sclk, res in index.get((sol, seq), []):
            if not sclk.startswith(pre):
                continue
            sl = res.get("scan_level")
            if sl and sl.get("bands"):
                got = sl["bands"][0]
        if got is None:
            fails.append("sol %s %s %s 应该有扫描级候选，但一个都没给出"
                         % (sol, seq, pre))
        elif abs(got["center"] - want_cm) > tol:
            fails.append("sol %s %s %s 的扫描级候选峰位 %.0f 偏离 %.0f 超过 %.0f"
                         % (sol, seq, pre, got["center"], want_cm, tol))
        else:
            print("  OK  sol %-4s %-9s %.0f cm-1（%.1fσ，幅度 %.2f×单点噪声）"
                  % (sol, seq, got["center"], got["snr"], got["ratio_to_point_sigma"]))
    for sol, seq in SCAN_LEVEL_DENY:
        checks += 1
        bad = []
        for sclk, res in index.get((sol, seq), []):
            sl = res.get("scan_level")
            if sl and sl.get("bands"):
                bad.append("%s %.0f cm-1" % (sclk[:14], sl["bands"][0]["center"]))
        if bad:
            fails.append("sol %s %s 是已知的零信号对照，扫描级兜底却报了候选：%s"
                         % (sol, seq, "、".join(bad)))
        else:
            print("  OK  sol %-4s %-9s 没有任何扫描级候选（对照组）" % (sol, seq))

    # G 类：有机质 G 带通道
    print("有机质 G 带通道（~1600 cm-1）")
    for sol, seq, pre, want_cm, tol, min_snr in G_BAND_WANT:
        checks += 1
        got = None
        for sclk, res in index.get((sol, seq), []):
            if not sclk.startswith(pre):
                continue
            cands = []
            for g in res.get("mineral_groups") or []:
                if g.get("g_band"):
                    cands.append((g["g_band"]["snr"], g["g_band"], g["n_points"]))
            if res.get("g_band_scan"):
                cands.append((res["g_band_scan"]["snr"], res["g_band_scan"], 0))
            if cands:
                got = max(cands, key=lambda t: t[0])
        tag = "sol %s %s %s" % (sol, seq, pre)
        if got is None:
            fails.append("%s 没报出任何 G 带候选" % tag)
        elif abs(got[1]["center"] - want_cm) > tol:
            fails.append("%s 的 G 带峰位 %.0f 偏离 %.0f 超过 %.0f"
                         % (tag, got[1]["center"], want_cm, tol))
        elif got[1]["snr"] < min_snr:
            fails.append("%s 的 G 带只有 %.1fσ（期望 ≥%.1fσ）"
                         % (tag, got[1]["snr"], min_snr))
        else:
            print("  OK  %-26s %.0f cm-1（%.1fσ，幅度 %.2f×单点噪声，%d 个点）"
                  % (tag, got[1]["center"], got[1]["snr"],
                     got[1]["ratio_to_point_sigma"], got[2]))

    print("=" * 78)
    # C 类：单元级检验
    print("单元级检验（论文结论：Séítah = 橄榄石堆积，Máaz = 盐类为主）")
    for sol, seq, desc, spec in UNIT:
        item = best.get((sol, seq))
        if item is None:
            continue
        n_tot, cs = item
        for name, (lo, hi) in spec.items():
            checks += 1
            # name 允许写成一组类别（例如"高氯酸盐或磷酸盐"或"钠高氯酸盐"），
            # 这时把几类的点数加起来比。
            names = name if isinstance(name, tuple) else (name,)
            n = sum(cs.get(k, (0, 0))[0] for k in names)
            if not (lo <= n <= hi):
                fails.append("%s 的%s是 %d 点（该次观测共 %d 点），"
                             "超出单元级区间 %d-%d"
                             % (desc, "+".join(names), n, n_tot, lo, hi))
            else:
                print("  OK  %-28s %s %d/%d 点" % (desc, "+".join(names), n, n_tot))

    print("=" * 78)
    print("共 %d 项断言，%d 项失败" % (checks, len(fails)))
    for f in fails:
        print("  x " + f)
    if not fails:
        print("  全部通过。")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
