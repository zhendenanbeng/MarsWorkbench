#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""PIXL 数据准备脚本 09：火星尘厚度与成分估计（VanBommel 2026 方法移植）。
=====================================================

来源：WUSTL Library DOI 10.7936/6rxs-108366 的 PIXL_CT_dust_thickness.m
（PNAS 投稿 "Fine and Ultrafine Martian Dust Present Additional Risks for
Future Crewed Missions to Mars"）。物理模型：

    测量谱 = LF(T) × 尘成分信号 + (1 − LF(T)) × 底材信号

LF(T) 为官方预计算的采样分数表（质量厚度 T 时信号中来自尘层的比例，
逐元素不同）。两个 PIXL 定标靶提供约束：PTFE 靶（8 元素全来自尘）
约束尘成分；BHVO-2 靶用 P、S 两元素（岩石几乎不含、尘富含）约束厚度。
蒙特卡洛逐步增减扰动搜索 8 元素尘成分，每轮取 χ² 最优 T，中位数报告。

边界（如实声明）：
  * 支持数据（仪器灵敏度、PTFE/BHVO 飞行计数、初始成分）按 .m 原文
    硬编码转录；
  * 结果与 .m 在统计上同类但随机路径不同（Python/ MATLAB 随机数不同）——
    中位厚度与成分应稳定（先验与数据共同约束），不逐位一致；
  * 官方论文的正式数值在 PNAS 文章中，本脚本不与其数值做断言。

用法：
    python 09_dust_thickness.py                    # 默认使用随包下载的官方 CSV
    python 09_dust_thickness.py --numiter 200 --seed 42
    python 09_dust_thickness.py --csv <自定义采样分数表>
"""

from __future__ import annotations

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0, str(_BootPath(__file__).resolve().parents[3] / "02_图形界面代码"))
import source_paths as _layout            # noqa: E402

import argparse
import sys
from pathlib import Path

from instruments import dust_model as dm   # noqa: E402

DR = _layout.project_layout.data_root()
DEFAULT_CSV = DR / "PIXL" / "下载" / "pixl_dust_model" / "PIXL_Sampling_Fractions.csv"


def main(argv=None):
    ap = argparse.ArgumentParser(description="PIXL 火星尘厚度与成分估计（VanBommel 2026 移植）")
    ap.add_argument("--csv", default=str(DEFAULT_CSV),
                    help="采样分数表 CSV（默认随包下载的官方表）")
    ap.add_argument("--numiter", type=int, default=100,
                    help="蒙特卡洛轮数（默认 100；论文量级 1e5 步，取值越大越稳）")
    ap.add_argument("--seed", type=int, default=20261121,
                    help="随机种子（默认固定，保证结果可复现）")
    args = ap.parse_args(argv)

    csv_p = Path(args.csv)
    if not csv_p.is_file():
        print("采样分数表不存在：%s\n先运行下载（WUSTL Library DOI 10.7936/6rxs-108366，"
              "文件 PIXL_Sampling_Fractions.csv）" % csv_p)
        return 1

    T, LF = dm.load_fractions(csv_p)
    print("采样分数表: %s" % csv_p.name)
    print("  质量厚度 %g – %g g/cm²（%d 档）× 8 元素" % (T.min(), T.max(), LF.shape[1]))

    inp_e = dm.load_support()
    inp = dm.build_input(inp_e, LF)

    import time
    t0 = time.time()
    summary, comp_log = dm.mcmc(inp, T, numiter=args.numiter, seed=args.seed)
    print("MCMC: %d 轮 × %d 步，用时 %.1f s（seed=%d，确定性可复现）"
          % (args.numiter, 270, time.time() - t0, args.seed))

    res = dm.summarize(summary, comp_log)
    print("\n--- 中位尘厚度（与 .m 输出段同口径）---")
    print("PTFE:  T = %.2e g/cm² | %.3f µm @ %.1f g/cm³（χ² %.1f）"
          % (res["ptfe_T_gcm2"], res["ptfe_T_um"], 2.0, res["ptfe_chi2"]))
    print("BHVO2: T = %.2e g/cm² | %.3f µm @ %.1f g/cm³（χ² %.1f）"
          % (res["bhvo_T_gcm2"], res["bhvo_T_um"], 2.0, res["bhvo_chi2"]))

    print("\n--- 尘成分中位（wt%，95% CI）---")
    for e, v in res["dust_composition"].items():
        print("  %-3s %6.2f [%6.2f, %6.2f]"
              % (e, v["median_wt"], v["ci95"][0], v["ci95"][1]))

    print("\n边界：本移植逐式对照 .m 实现；与 MATLAB 的随机路径不同，中位结果"
          "应稳定但步级轨迹不逐位一致。官方论文的正式数值以 PNAS 文章为准。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
