#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SHERLOC 数据准备脚本 07：原始帧自行扣暗（ERA−ERB）→ 科学谱输出。
=====================================================

用途：从 PDS 原始件（data_raw 的 EDR）出发做自己的数据处理。
官方 RRS = ERA（主动帧）− ERB（暗帧）+ 激光能量归一化 + 波长校正 +
宇宙射线/增益合并；本脚本只复现其中**自洽可验证**的第一步：

  * 逐激光点 ERA − ERB 相减（暗采集镜像扫描几何，按行号对齐）；
  * 波长轴借用同 sol 官方 RRS 的波长校正表（对原始通道是近似）；
  * 输出 RRS 同布局 CSV + 配对 XML（J96 版本位）到「自扣暗_输出」
    子文件夹——SHERLOC 台可直接打开该文件夹，全部既有功能可用。

边界（如实声明）：
  * 未做激光能量归一化/宇宙射线/增益合并——绝对强度与官方 RRS
    不可比，判读看峰位与相对形态；
  * Region 2/3（荧光段）同样相减输出，供有机荧光证据自查；
  * 无同序列同 expID ERB 的 ERA 段如实跳过并列原因。

用法：
    python 07_raw_dark_subtract.py --dir <下载文件夹>
    python 07_raw_dark_subtract.py --dir <文件夹> --out <输出目录> --region 1,2,3
"""

from __future__ import annotations

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0, str(_BootPath(__file__).resolve().parents[3] / "02_图形界面代码"))
import source_paths as _layout            # noqa: E402

import argparse
import json
import sys
from pathlib import Path

from sherloc import raw_frames            # noqa: E402

DR = _layout.project_layout.data_root()


def main(argv=None):
    ap = argparse.ArgumentParser(description="SHERLOC 原始帧自行扣暗（ERA−ERB）")
    ap.add_argument("--dir", required=True, help="含原始件的文件夹（ERA/ERB，可含 RRS/RLS）")
    ap.add_argument("--out", default=None, help="输出目录（默认 <dir>/自扣暗_输出）")
    ap.add_argument("--region", default="1,2,3", help="处理的 Region（默认 1,2,3）")
    args = ap.parse_args(argv)

    folder = Path(args.dir)
    if not folder.is_dir():
        print("文件夹不存在：%s" % folder)
        return 1
    regions = tuple(int(x) for x in re_split(args.region))
    out_dir = Path(args.out) if args.out else folder / "自扣暗_输出"

    report = raw_frames.process_folder(folder, out_dir=out_dir, region_list=regions)

    print("自扣暗输出目录：%s" % report["out_dir"])
    print("波长轴来源：%s" % (report["wavelength_source"] or "（无 RRS——输出仅通道空间，台里无法出图）"))
    for p in report["pairs"]:
        print("  ✓ %s exp%s：%d 点，Region %s → %s"
              % (p["seq"], p["exp"], p["points"],
                 "/".join(str(r) for r in p["regions"]), Path(p["csv"]).name))
    for seq, exp, why in report["skipped"]:
        print("  ✗ %s exp%s 跳过：%s" % (seq, exp, why))
    print("\n口径边界：%s——判读看峰位与相对形态，绝对强度与官方 RRS 不可比。"
          % report["boundary"])
    if not report["pairs"]:
        print("没有成功配对——检查文件夹里是否有 ERA 与对应 ERB。")
        return 1
    print("\n下一步：SHERLOC 台「打开本地数据」选输出文件夹（%s）即可正常浏览/判读。"
          % report["out_dir"])
    return 0


def re_split(s):
    return [x.strip() for x in s.split(",") if x.strip()]


if __name__ == "__main__":
    sys.exit(main())
