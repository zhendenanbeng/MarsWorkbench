# -*- coding: utf-8 -*-
"""参数扫描：σ 窗口、FWHM 下限、稳健 σ 开关，各组合下点数怎么变。

动因：手册（A-5）和论文（2.3 节）给的口径不一样——
  手册：σ 取 2000-2100，FWHM 下限 25，异常光谱改用扫描中位 σ
  论文：σ 取 2500-3000，FWHM 下限 30，同样建议用中位数
本脚本跑遍组合，看"点位数"这个最不稳的量到底被哪一个参数支配。

数值全部来自工具自己的 qc.analyse_spectra，不是另写一套实现，
所以看到的就是工具真跑出来的结果。

用法：
    python verify\\param_sweep.py
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

COMBOS = [
    # (标签, σ窗口, FWHM下限, 稳健σ, SNR门槛)
    ("论文口径  s2500-3000 FWHM30 SNR3", (2500, 3000), 30.0, True, 3.0),
    ("手册口径  s2000-2100 FWHM25 SNR3", (2000, 2100), 25.0, True, 3.0),
    ("论文口径  s2500-3000 FWHM30 SNR0", (2500, 3000), 30.0, True, 0.0),
    ("手册口径  s2000-2100 FWHM25 SNR0", (2000, 2100), 25.0, True, 0.0),
    ("手册口径  s2000-2100 FWHM30 SNR0", (2000, 2100), 30.0, True, 0.0),
    ("论文口径  s2500-3000 FWHM30 SNR1", (2500, 3000), 30.0, True, 1.0),
]

# 论文里报出的逐类点数（Corpolongo et al., 2023），用来对照
PAPER = {
    "293": {"碳酸盐": 36, "硫酸盐": 21, "未定硅酸盐": 3, "橄榄石": 1,
            "高氯酸盐或磷酸盐": 1},
    "141": {"碳酸盐": 0, "硫酸盐": 0, "未定硅酸盐": 0, "橄榄石": 0,
            "高氯酸盐或磷酸盐": 0},
}


def class_counts(res):
    """按类别统计点数。一个点归多类就各算一次，和论文的口径一致。"""
    out = {}
    for g in res.get("mineral_groups") or []:
        if not g["classes"]:
            continue
        name = g["classes"][0]["class"]
        out[name] = out.get(name, 0) + g["n_points"]
    return out


def class_peak(res):
    """每个类别的代表峰位：取点数最多的那个组。"""
    out = {}
    for g in res.get("mineral_groups") or []:
        if not g["classes"]:
            continue
        name = g["classes"][0]["class"]
        if name not in out or g["n_points"] > out[name][1]:
            out[name] = (g["band"]["center"], g["n_points"])
    return out


def main():
    _products, observations = scan_folder(DATA)
    todo = [ob for ob in observations
            if ob.sol in ("293", "141") and ob.spectra is not None]
    if not todo:
        print("没找到 sol 293 / 141 的观测")
        return 1

    for ob in todo:
        wl, r1, r2, r3 = ob.read_regions()
        if wl is None or r1 is None:
            print("sol %s %s：读不到波长表或 Region 1，跳过"
                  % (ob.sol, ob.sequence_id))
            continue
        print("=" * 92)
        head = "sol %s  %s  %s" % (ob.sol, ob.sequence_id, ob.sclk)
        paper = PAPER.get(ob.sol)
        if paper:
            head += "    论文：" + "  ".join("%s %d" % (k, v)
                                            for k, v in paper.items())
        print(head)
        print("%-34s %6s %6s %6s %6s %6s %8s" %
              ("配置", "有波段", "碳酸盐", "硫酸盐", "硅酸盐", "橄榄石", "高氯酸盐"))
        for label, nw, mf, robust, snr in COMBOS:
            res = qc.analyse_spectra(wl, r1, r2, r3, noise_window=nw,
                                     min_fwhm=mf, robust=robust, snr_min=snr)
            cc = class_counts(res)
            print("%-34s %6d %6d %6d %6d %6d %8d" % (
                label, res.get("n_points_with_band") or 0,
                cc.get("碳酸盐", 0), cc.get("硫酸盐", 0),
                cc.get("未定硅酸盐", 0), cc.get("橄榄石", 0),
                cc.get("高氯酸盐或磷酸盐", 0)))

        print("  峰位与 σ 替换情况（稳健σ开）：")
        for label, nw, mf, robust, snr in (COMBOS[0], COMBOS[1], COMBOS[2]):
            res = qc.analyse_spectra(wl, r1, r2, r3, noise_window=nw,
                                     min_fwhm=mf, robust=robust, snr_min=snr)
            pk = class_peak(res)
            txt = "、".join("%s %.0f(%d点)" % (k, v[0], v[1])
                           for k, v in sorted(pk.items(),
                                              key=lambda x: -x[1][1])[:5])
            print("    %-30s %s   σ替 %d"
                  % (label[:30], txt, len(res.get("sigma_replaced") or [])))

        # 第二个实验：信号与噪声到底该取「三表相加」还是「只用 Region 1」。
        # 拉曼峰只在 R1 里有，而噪声是三个表都在贡献——相加会把噪声按平方和抬高，
        # 信号却不变，所以相加谱的 SNR 天然比 R1 单表低。
        print("  三表相加 vs 只用 Region 1（论文口径 σ2500-3000 FWHM30 SNR3）：")
        for tag, args in (("三表相加", (r1, r2, r3)), ("只用 R1", (r1, None, None))):
            res = qc.analyse_spectra(wl, *args)
            cc = class_counts(res)
            print("    %-8s 有波段 %3d  碳酸盐 %3d  硫酸盐 %3d  硅酸盐 %3d  "
                  "橄榄石 %3d  高氯酸盐 %3d"
                  % (tag, res.get("n_points_with_band") or 0,
                     cc.get("碳酸盐", 0), cc.get("硫酸盐", 0),
                     cc.get("未定硅酸盐", 0), cc.get("橄榄石", 0),
                     cc.get("高氯酸盐或磷酸盐", 0)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
