# -*- coding: utf-8 -*-
"""APXS 峰拟合层护栏（发行自包含，合成能谱）。"""
import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0, str(_BootPath(__file__).resolve().parents[3] / "02_图形界面代码"))
import source_paths as _layout

import sys

import numpy as np

ROOT = _layout.project_layout.APP
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from instruments import apxs_fit                  # noqa: E402

TOTAL = FAIL = 0


def check(ok, name):
    global TOTAL, FAIL
    TOTAL += 1
    FAIL += not ok
    print(("OK " if ok else "FAIL ") + name, flush=True)


def main():
    rng = np.random.default_rng(23)
    e = np.linspace(0.5, 12.0, 1024)
    counts = 100 + 0.5 * e
    for el, amp in (("Fe", 8000.0), ("Ca", 2500.0), ("Si", 12000.0)):
        c = apxs_fit.KA_KEV[el]
        sigma = apxs_fit.PEAK_KEV_FWHM / 2.355
        counts = counts + amp * np.exp(-0.5 * ((e - c) / sigma) ** 2)
    counts = rng.poisson(counts).astype(float)

    areas = apxs_fit.spectrum_areas(e, counts, elements=("Fe", "Ca", "Si", "Ni"))
    for el, expected in (("Fe", 8000.0), ("Ca", 2500.0), ("Si", 12000.0)):
        r = areas[el]
        check(r.get("amp") is not None and 0.6 * expected <= r["amp"] <= 1.6 * expected,
              "%s 峰高量级恢复（%.0f vs 合成 %.0f）" % (el, r["amp"], expected))
        check(r["detected"], "%s 检出" % el)
    r = areas["Ni"]
    check(not r["detected"], "无 Ni 合成峰时 Ni 不检出")

    # 窗越界：Na(1.041)-0.12 在 0.5 起始窗内可拟合；Br 无 Kα 表则跳过
    out = apxs_fit.fit_peak_area(e, counts, 0.0)
    check(out["area"] is None, "窗越界返回 area=None（不猜值）")

    # 数据不足
    out2 = apxs_fit.fit_peak_area(np.array([6.40, 6.42]), np.array([1.0, 2.0]), 6.404)
    check(out2["area"] is None, "窗内通道不足返回 area=None")

    print("共 %d 项断言，%d 项失败" % (TOTAL, FAIL))
    return int(FAIL > 0)


if __name__ == "__main__":
    raise SystemExit(main())
