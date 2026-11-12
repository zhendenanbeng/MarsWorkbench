# -*- coding: utf-8 -*-
"""一次性脚本：查清"用 processed 数据还做预处理"这件事的真实影响面。

两件事：
 1. 全库对账 —— 文件名的 normalize/cosmic 标记 vs label 表名，谁跟谁不一致；
 2. 量化 —— 那些被"文件名规则"误判的产品，实际被改动了多少。
"""
import os
import sys
from pathlib import Path

import numpy as np

# 旧开发树位置：重做实验时改成你的本地路径
DEV_ROOT = Path.home() / "Desktop" / "Sherlock"
sys.path.insert(0, str(DEV_ROOT / "01_软件工作台"))
DATA = str(DEV_ROOT / "02_练习数据与验证基准")

from sherloc import preprocess                     # noqa: E402
from sherloc.organize import scan_folder           # noqa: E402


def items():
    for area in ("processed", "raw"):
        base = os.path.join(DATA, area)
        if not os.path.isdir(base):
            continue
        for sol in sorted(os.listdir(base)):
            folder = os.path.join(base, sol)
            if not os.path.isdir(folder):
                continue
            _, scans = scan_folder(folder)
            for sc in scans:
                if sc.spectra is not None:
                    yield area, sol, sc


print("=" * 100)
print("对账：文件名标记  vs  label 表名")
print("=" * 100)
bad = []
for area, sol, sc in items():
    sp = sc.spectra
    f = dict(sp.pid.flags)
    names = " | ".join((t.name or "") for t in sc.parts[0].tables()).upper()
    rdr = "LASER-NORMALIZED" in names
    pdr = "PROCESS DATA" in names
    tag = "LASER-NORMALIZED" if rdr else ("PROCESS DATA" if pdr else "其他")
    old_skip = (f.get("normalize") == "N")           # 旧规则：只看文件名
    new_skip = rdr or old_skip                       # 新规则：label 优先
    if rdr and not old_skip:                         # 明明归一过，却被再归一
        bad.append((area, sol, sc))
    print("%-10s %-12s %-5s %-30s norm=%s cosmic=%s | label=%-16s 旧=%-8s 新=%-8s %s"
          % (area, sol, sp.product_type, sp.name[-30:],
             f.get("normalize"), f.get("cosmic"), tag,
             "跳过" if old_skip else "会重做", "跳过" if new_skip else "会重做",
             "  <== 文件名读错" if (rdr and not old_skip) else ""))

print()
print("=" * 100)
print("被旧规则误判、会做二次归一的产品：%d 个" % len(bad))
print("=" * 100)

gui = {"deadpixel": True, "dark": True, "laser": True,
       "despike": False, "baseline": False}

for area, sol, sc in bad:
    sp = sc.spectra
    stage0 = preprocess.run_pipeline(sc, dict(gui, laser=False), points=None,
                                     region=1)["r1"]
    new = preprocess.run_pipeline(sc, gui, points=None, region=1)["r1"]
    e = preprocess.laser_energy(sc)
    n = stage0.shape[0]
    e = np.asarray(e, dtype=float)[:n]
    gain = (e[0] if e[0] else e.mean()) / np.where(e == 0, np.nan, e)
    old = stage0 * gain[:, None]                      # 旧规则：再归一一次
    rel = np.abs(old - new) / np.maximum(np.abs(new), 1e-12)
    print("%s %s" % (sol, sp.name))
    print("  二次归一倍数 %.4f ~ %.4f（中位 %.4f）"
          % (gain.min(), gain.max(), np.median(gain)))
    print("  改到 %d / %d 个元素（%.2f%%），相对差异中位 %.3f%%、最大 %.2f%%"
          % ((np.abs(old - new) > 1e-9).sum(), old.size,
             100.0 * (np.abs(old - new) > 1e-9).sum() / old.size,
             100 * np.median(rel), 100 * rel.max()))
    print("  每个点的谱型被整体乘了不同的数 ⇒ 峰的相对高度会变")

    withsp = preprocess.run_pipeline(sc, dict(gui, despike=True, despike_force=True),
                                     points=None, region=1)["r1"]
    print("  二次去宇宙射线（勾上才发生）：改了 %d 个元素"
          % (np.abs(withsp - new) > 1e-9).sum())
    print()
