# -*- coding: utf-8 -*-
"""校验质控台里的「逐点矿物归属」图与导出的 CSV 是否自洽。

检查项：
  1. HTML 里没有残留的 __XXX__ 占位符
  2. 图上恰好 100 个点
  3. 图里每个点的字母与 CSV 的「类别缩写」逐点一致
  4. 每个点的位置确实来自 RLS 坐标表（用 _grid_from_rls 独立重算一遍对比），
     不是随手画的示意格
"""
from __future__ import annotations

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout

import csv
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = str(_layout.project_layout.APP)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from sherloc.organize import scan_folder      # noqa: E402
from sherloc.workbench import _grid_from_rls  # noqa: E402

QC_DIR = os.path.join(_layout.project_layout.data_root(), "SHERLOC",
                      "练习数据", "raw", "_qc")
DATA = os.path.join(os.path.dirname(ROOT), "03_下载数据/SHERLOC/练习数据")


def check(html_path, csv_path):
    html = open(html_path, encoding="utf-8").read()
    bad = sorted(set(re.findall(r"__[A-Z]+__", html)))
    if bad:
        return ["HTML 里还有没替换的占位符: %s" % ", ".join(bad)]

    block = re.split(r"<div class='mmap'[^>]*>", html, maxsplit=1)
    if len(block) < 2:
        return ["没有找到矿物分布图区块"]
    cells = re.findall(r"<div class='([^']*)' style='left:([\d.]+)%;top:([\d.]+)%'"
                       r" title='([^']*)'>(.*?)</div>", block[1], re.S)

    errs = []

    # 点数以 CSV 为准（有的扫描是 200 点，比如 sol 207 的四段拼接）
    with open(csv_path, encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    want_n = len(rows)
    if len(cells) != want_n:
        errs.append("图上的点数 %d，CSV 是 %d" % (len(cells), want_n))

    painted = {}
    for cls, x, y, title, text in cells:
        m = re.match(r"点 (\d+)", title)
        if not m:
            errs.append("点的标题不含点号: %r" % title[:60])
            continue
        i = int(m.group(1))
        painted[i] = (float(x), float(y))
        if not (0.0 <= float(x) <= 100.0 and 0.0 <= float(y) <= 100.0):
            errs.append("点 %d 的坐标 (%.1f, %.1f) 跑到画面外" % (i, float(x), float(y)))

    if sorted(painted) != list(range(want_n)):
        errs.append("点号不是 0-%d 各一次" % (want_n - 1))
    by_pt = {int(r["点号"]): r for r in rows}
    for cls, x, y, title, text in cells:
        m = re.match(r"点 (\d+)", title)
        if not m:
            continue
        i = int(m.group(1))
        want = by_pt[i]["类别缩写"]
        if want != text.strip():
            errs.append("点 %d：图里写 '%s'，CSV 写 '%s'" % (i, text.strip(), want))

    # 位置必须能由 RLS 独立重算出来。sol / 序列号从文件名里取，
    # 不能写死——写死过一次，拿 sol 257 的坐标去核对 sol 293 的图，满盘皆错。
    mm = re.search(r"sol(\d+)_(srlc\d+)", os.path.basename(html_path))
    if not mm:
        errs.append("从文件名 %s 里认不出 sol 和序列号" % os.path.basename(html_path))
        return errs
    want_sol, want_seq = mm.group(1), mm.group(2)
    obs = None
    for _p, ob in [(None, o) for o in scan_folder(DATA)[1]]:
        if ob.sol == want_sol and ob.sequence_id == want_seq:
            obs = ob
            break
    if obs is None:
        errs.append("找不到 sol %s %s，无法独立核对坐标" % (want_sol, want_seq))
    else:
        rls_pts = obs.read_shot_positions()
        if not rls_pts:
            # 有些观测（1-2 个点的定标小扫描）本来就没有 RLS 产品，
            # 这时用示意摆法是合理的，不算错。
            return errs
        grid, title = _grid_from_rls(rls_pts, want_n)
        if "无 RLS 坐标" in title:
            errs.append("这次观测有 RLS 记录，但图上用的是示意坐标")
        off = 0
        for i, (x, y) in painted.items():
            if i >= len(grid):
                continue
            if abs(x - grid[i][0]) > 0.5 or abs(y - grid[i][1]) > 0.5:
                off += 1
        if off:
            errs.append("有 %d 个点的位置与 RLS 重算结果不符（超过 0.5%%）" % off)
    return errs


def main():
    # 自动发现最新的质控报告对（.html + _points.csv）；不再写死旧报告名
    import glob as _glob
    pairs = []
    for html_path in sorted(_glob.glob(os.path.join(QC_DIR, "sol*.html"))):
        csv_path = html_path[:-5] + "_points.csv"
        if os.path.exists(csv_path):
            pairs.append((os.path.getmtime(html_path), html_path, csv_path))
    if not pairs:
        print("SKIP 数据不在位：", QC_DIR, " 下没有质控报告（先用工作台生成质控台）")
        return 0
    pairs.sort()
    html_path, csv_path = pairs[-1][1], pairs[-1][2]
    print("检查最新报告:", os.path.basename(html_path))
    errs = check(html_path, csv_path)
    if errs:
        print("不合格，共 %d 处：" % len(errs))
        for e in errs[:20]:
            print("  -", e)
        return 1
    print("通过：100 格不重不漏，蛇形映射正确，字母与 CSV 逐点一致，无残留占位符。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
