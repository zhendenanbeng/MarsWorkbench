# -*- coding: utf-8 -*-
"""仪器本底扣除（instbkg）护栏：合成谱 + 隔离数据目录，离线可跑。

契约（对应 preprocess._load_stowed_background / run_pipeline 2b 步）：
* 纯本底合成扫描：M 落在 [0,2] 且接近真值，800 cm-1 主带减到接近 0；
* 含峰合成扫描：831 cm-1 真实峰（橄榄石口径）减后仍立得住、峰位不跑偏；
* 开关默认关；缺文件、坏 JSON、激光波长不符时一律 warn 跳过且不改变数据；
* 减除后基线步骤照常执行；重复运行结果逐位一致；管线不改输入数据。
"""
import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0, str(_BootPath(__file__).resolve().parents[3] / "02_图形界面代码"))
import source_paths as _layout

import json
import os
import re
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

ROOT = _layout.project_layout.APP
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from sherloc import preprocess                     # noqa: E402

TOTAL = FAIL = 0


def check(ok, name):
    global TOTAL, FAIL
    TOTAL += 1
    FAIL += not ok
    print(("OK " if ok else "FAIL ") + name, flush=True)


LASER_NM = 248.6
SHIFT = np.linspace(214.0, 4124.0, 2048)
M_TRUE = 0.21


def _gauss(x, center, amp, fwhm):
    return amp * np.exp(-4.0 * np.log(2.0) * ((x - center) / fwhm) ** 2)


def _true_background(x):
    """熔石英本底形状：四带幅度对齐 Phase 1 在收纳臂谱上的实测量级。"""
    return (10.0 + _gauss(x, 490.0, 140.0, 36.0) + _gauss(x, 800.0, 250.0, 44.0)
            + _gauss(x, 1055.0, 105.0, 50.0) + _gauss(x, 1190.0, 70.0, 56.0))


def make_background_file(folder, shift, spec, laser_nm=248.6, raw=None):
    d = Path(folder) / "data" / "SHERLOC" / "缓存"
    d.mkdir(parents=True, exist_ok=True)
    target = d / "stowed_background.json"
    target.write_text(raw if raw is not None else json.dumps({
        "scope": "护栏合成数据", "laser_nm": laser_nm,
        "shift_cm1": [float(v) for v in shift],
        "background": [float(v) for v in spec],
        "scans_used": [{"scan": "synthetic"}],
    }), encoding="utf-8")
    return target


def make_scan(rng, scale, with_peak):
    """合成扫描：真本底 × scale（+831 cm-1 峰）+ 缓坡 + 噪声。"""
    wl = 1.0e7 / (1.0e7 / LASER_NM - SHIFT)          # 升序波长轴
    bg = _true_background(SHIFT) * scale
    peak = _gauss(SHIFT, 831.0, 60.0, 24.0) if with_peak else 0.0
    slope = 0.002 * SHIFT
    rows = [bg + peak + slope + rng.normal(0.0, 2.0, SHIFT.size) for _ in range(6)]
    r1 = np.vstack(rows)

    class FakeProduct:
        product_type = "RRS"
        name = "ss__9999_synthetic__0090000srlc99999w108cgzj03.csv"
        n_spectra = r1.shape[0]
        pid = SimpleNamespace(sclk="123", flags={})

        def __init__(self):
            self._wl, self._r1 = wl, r1

        def wavelength_nm(self):
            return self._wl

        def regions(self):
            return (self._r1, None, None)

        def tables(self):
            return []

    class Scan:
        sol = "9999"
        sequence_id = "synthetic"
        sclk = "123"

        def __init__(self):
            self.parts = [FakeProduct()]
            self.products = []

    return Scan(), wl, r1


def w8_mask(shift):
    return (shift >= 775.0) & (shift <= 825.0)


def instbkg_log(log):
    for entry in log:
        if entry.get("step") == "instbkg":
            return entry
    return None


STEPS_ON = {"deadpixel": True, "dark": True, "laser": True,
            "instbkg": True, "baseline": True}


def main():
    rng = np.random.default_rng(41)
    scan_a, wl, raw_a = make_scan(rng, scale=M_TRUE, with_peak=False)
    scan_b, _, raw_b = make_scan(rng, scale=M_TRUE, with_peak=True)
    original_a = raw_a.copy()
    on = dict(STEPS_ON)
    off = {k: v for k, v in on.items() if k != "instbkg"}

    with tempfile.TemporaryDirectory(prefix="sherlock-instbkg-") as folder:
        shift_axis = 1.0e7 / LASER_NM - 1.0e7 / np.asarray(wl, dtype=float)
        make_background_file(folder, shift_axis, _true_background(shift_axis))
        env = {"MARSWORKBENCH_DATA_HOME": folder}
        with patch.dict(os.environ, env):
            res_a = preprocess.run_pipeline(scan_a, dict(on))
            res_a2 = preprocess.run_pipeline(scan_a, dict(on))
            res_b = preprocess.run_pipeline(scan_b, dict(on))
            res_off = preprocess.run_pipeline(scan_a, dict(off))
            entry = instbkg_log(res_a["log"])

        check(entry is not None and entry["kind"] == "done", "本底在位时 instbkg 记 done 日志")
        m = re.search(r"M=([0-9.]+)", entry["text"] or "") if entry else None
        scale = float(m.group(1)) if m else -1.0
        check(0.0 <= scale <= 2.0, "缩放 M 落在 [0,2]")
        check(abs(scale - M_TRUE) < 0.03, "纯本底扫描的 M 恢复真值 %.2f（±0.03，实得 %.3f）"
              % (M_TRUE, scale))

        mask = w8_mask(np.asarray(res_a["shift"], dtype=float))
        after = float(np.median(np.asarray(res_a["r1"], dtype=float)[:, mask]))
        before = float(np.median(raw_a[:, preprocess.DEAD_LO:preprocess.DEAD_HI][:, mask]))
        check(abs(after) < 15.0, "800 主带减后中位接近 0（实得 %.1f）" % after)
        check(before > 30.0, "减前主带确有本底（实得 %.1f）" % before)

        peak_on = np.asarray(res_b["r1"], dtype=float).mean(axis=0)
        peak_off = np.asarray(preprocess.run_pipeline(scan_b, dict(off))["r1"],
                              dtype=float).mean(axis=0)
        sh = np.asarray(res_b["shift"], dtype=float)
        ol = (sh >= 760.0) & (sh <= 900.0)
        pos_on = float(sh[ol][int(np.argmax(peak_on[ol]))])
        pos_off = float(sh[ol][int(np.argmax(peak_off[ol]))])
        check(abs(pos_on - 831.0) <= 5.0 and abs(pos_off - 831.0) <= 5.0,
              "831 cm-1 真实峰位在减除前后都不跑偏（%.1f/%.1f）" % (pos_on, pos_off))
        check(float(peak_on[ol].max()) > 30.0,
              "减除后真实峰仍立得住（峰高 %.1f）" % float(peak_on[ol].max()))

        entry_off = instbkg_log(res_off["log"])
        check(entry_off is not None and entry_off["kind"] == "skip",
              "开关关闭时 instbkg 记 skip 日志")
        check(np.isfinite(np.asarray(res_a["r1"], dtype=float)).all(),
              "减除不得引入 NaN/Inf")
        check(all(e["step"] == "baseline" and e["kind"] == "done"
                  for e in res_a["log"] if e["step"] == "baseline"),
              "减除后基线步骤照常执行")
        check(np.array_equal(np.asarray(res_a["r1"], dtype=float),
                             np.asarray(res_a2["r1"], dtype=float)),
              "重复运行逐位一致")
        check(np.array_equal(raw_a, original_a), "管线不得修改输入扫描数据")

        # --- 缺文件 / 坏 JSON / 激光波长不符 / 缺省开关：一律安全侧 ---
        with tempfile.TemporaryDirectory(prefix="sherlock-instbkg-missing-") as empty:
            with patch.dict(os.environ, {"MARSWORKBENCH_DATA_HOME": empty}):
                res_missing = preprocess.run_pipeline(scan_a, dict(on))
        e2 = instbkg_log(res_missing["log"])
        check(e2 is not None and e2["kind"] == "warn", "缺本底文件时记 warn")
        check(np.allclose(np.asarray(res_missing["r1"], dtype=float),
                          np.asarray(res_off["r1"], dtype=float)),
              "缺本底文件时结果与关闭开关完全一致")

        with tempfile.TemporaryDirectory(prefix="sherlock-instbkg-bad-") as bad:
            make_background_file(bad, SHIFT, _true_background(SHIFT), raw="{not json")
            with patch.dict(os.environ, {"MARSWORKBENCH_DATA_HOME": bad}):
                res_bad = preprocess.run_pipeline(scan_a, dict(on))
        e3 = instbkg_log(res_bad["log"])
        check(e3 is not None and e3["kind"] == "warn", "坏 JSON 时记 warn 不崩溃")
        check(np.allclose(np.asarray(res_bad["r1"], dtype=float),
                          np.asarray(res_off["r1"], dtype=float)),
              "坏 JSON 时结果与关闭开关完全一致")

        with tempfile.TemporaryDirectory(prefix="sherlock-instbkg-laser-") as wl_mm:
            make_background_file(wl_mm, SHIFT, _true_background(SHIFT), laser_nm=300.0)
            with patch.dict(os.environ, {"MARSWORKBENCH_DATA_HOME": wl_mm}):
                res_laser = preprocess.run_pipeline(scan_a, dict(on))
        e4 = instbkg_log(res_laser["log"])
        check(e4 is not None and e4["kind"] == "warn" and "248.6" in (e4["text"] or ""),
              "激光波长不符时 warn 并给出双方波长")

        with patch.dict(os.environ, {"MARSWORKBENCH_DATA_HOME": folder}):
            res_default = preprocess.run_pipeline(scan_a, {})
        e5 = instbkg_log(res_default["log"])
        check(e5 is not None and e5["kind"] == "skip", "steps 缺省时 instbkg 默认跳过")

    # 部分覆盖契约：本底谱只覆盖部分波段时，覆盖外通道必须逐位不变
    with tempfile.TemporaryDirectory(prefix="sherlock-instbkg-partial-") as folder2:
        cov = SHIFT <= 2600.0          # 覆盖 ~61%（裁边后仍 >50%，不触发整步跳过）
        make_background_file(folder2, SHIFT[cov], (_true_background(SHIFT) * M_TRUE)[cov])
        env2 = {"MARSWORKBENCH_DATA_HOME": folder2}
        on_nb = {"deadpixel": True, "dark": True, "laser": True,
                 "instbkg": True, "baseline": False}   # 关基线：隔离本底步效应
        off_nb = {k: v for k, v in on_nb.items() if k != "instbkg"}
        with patch.dict(os.environ, env2):
            res_on = preprocess.run_pipeline(scan_a, dict(on_nb))
            res_off = preprocess.run_pipeline(scan_a, dict(off_nb))
        # 掩码用管线自身返回的 shift 轴（deadpixel 已裁边，通道数与 r1 一致）
        shift_out = np.asarray(res_on["shift"], dtype=float)
        r_on_m = np.asarray(res_on["r1"], dtype=float)
        r_off_m = np.asarray(res_off["r1"], dtype=float)
        out_mask = shift_out > 2800.0                 # 覆盖外（>2600 cm-1 留裕量）
        in_mask = (shift_out >= 400.0) & (shift_out <= 1200.0)   # 覆盖内含 800 主带
        same_out = bool(np.array_equal(r_on_m[:, out_mask], r_off_m[:, out_mask]))
        check(same_out, "部分覆盖：覆盖外通道逐位不变（instbkg 不越界）")
        diff_in = not np.array_equal(r_on_m[:, in_mask], r_off_m[:, in_mask])
        check(diff_in, "部分覆盖：覆盖内通道确有减除（对照有效）")

    print("共 %d 项断言，%d 项失败" % (TOTAL, FAIL))
    return int(FAIL > 0)


def r1_raw_of(scan):
    return scan.parts[0]._r1


if __name__ == "__main__":
    raise SystemExit(main())
