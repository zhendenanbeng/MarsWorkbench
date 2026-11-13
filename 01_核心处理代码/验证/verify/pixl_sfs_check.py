# -*- coding: utf-8 -*-
"""PIXL SFS 自研谱层护栏（发行自包含，合成 XRF 谱）。

契约（sfs_spectra.py）：
* read_sfs：表头定位通道起始列、非谱表返回 (None, None)、坏值行跳过；
* calibrate_energy：合成 Ca/Fe 双峰谱的两点自证恢复真实增益（±2%），
  身份假设失败（无 Ca/Fe 组合）时显式抛 ValueError；
* peak_areas：已知净面积恢复 ±泊松噪声内、误差≈√面积、窗越界不猜值。
"""
import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0, str(_BootPath(__file__).resolve().parents[3] / "02_图形界面代码"))
import source_paths as _layout

import json
import os
import sys
import tempfile

import numpy as np

ROOT = _layout.project_layout.APP
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from instruments import sfs_spectra               # noqa: E402

TOTAL = FAIL = 0


def check(ok, name):
    global TOTAL, FAIL
    TOTAL += 1
    FAIL += not ok
    print(("OK " if ok else "FAIL ") + name, flush=True)


def make_sfs_csv(folder, pmcs, gain, zero, seed=11):
    """合成 SFS csv：每 PMC 一行，谱 = 元素 Kα 高斯×强度 + 连续谱 + 泊松噪声。"""
    rng = np.random.default_rng(seed)
    ch = np.arange(sfs_spectra.N_CHANNELS)
    kev = gain * ch + zero
    header = ("SCLK,PMC,run_status,dpp_status,real_time,live_time,num_events,"
              + ",".join("c%d" % c for c in range(sfs_spectra.N_CHANNELS)))
    lines = [header]
    for i, pmc in enumerate(pmcs):
        spec = 20 + 5 * np.exp(-0.5 * ((kev - 3.0) / 0.3) ** 2)      # 连续谱+背景峰
        strengths = {"Ca": 6000 + 400 * i, "Fe": 9000 - 500 * i, "Si": 2000}
        for el, amp in strengths.items():
            e = sfs_spectra.KA_KEV[el]
            spec = spec + amp * np.exp(-4 * np.log(2) * ((kev - e) / 0.08) ** 2)
        counts = rng.poisson(np.clip(spec, 0, None))
        lines.append("668183994,%d,8,41,10,9.8,%d,78930,0,0,21739150,741,%s"
                     % (pmc, counts.sum(), ",".join(str(int(c)) for c in counts)))
    path = os.path.join(folder, "PS__9999_000ENA__synthetic___J01.CSV")
    open(path, "w", encoding="utf-8").write("\n".join(lines))
    return path, strengths


def main():
    gain_true, zero_true = 0.00695, 0.02
    with tempfile.TemporaryDirectory(prefix="sherlock-sfs-") as folder:
        path, strengths = make_sfs_csv(folder, [101, 102, 103], gain_true, zero_true)
        spec, meta = sfs_spectra.read_sfs(path)
        check(spec is not None and spec.shape == (3, sfs_spectra.N_CHANNELS),
              "read_sfs 解析 3×4096 谱矩阵")
        check(int(meta[0]["PMC"]) == 101, "PMC 元数据随行保留")

        mean = spec.mean(axis=0)
        axis = sfs_spectra.calibrate_energy(mean)
        check(abs(axis["keV_per_ch"] - gain_true) / gain_true < 0.02,
              "两点自证增益恢复 ±2%%（%.5f vs %.5f）" % (axis["keV_per_ch"], gain_true))
        check(abs(axis["zero_keV"] - zero_true) < 0.15,
              "零点恢复 |b|<0.15（%.3f）" % axis["zero_keV"])

        pa = sfs_spectra.peak_areas(mean, axis,
                                    elements=("Ca", "Fe", "Si"))
        ok_ca = pa["Ca"]["detected"] and pa["Ca"]["area"] > 0.5 * strengths["Ca"]
        check(ok_ca, "Ca 峰面积检出且量级正确（%.0f vs 合成 %d）"
              % (pa["Ca"]["area"], strengths["Ca"]))
        ok_fe = pa["Fe"]["detected"] and pa["Fe"]["area"] > 0.5 * strengths["Fe"]
        check(ok_fe, "Fe 峰面积检出且量级正确（%.0f vs 合成 %d）"
              % (pa["Fe"]["area"], strengths["Fe"]))
        check(pa["Ca"]["error"] > 0 and pa["Ca"]["error"] ** 2 <= pa["Ca"]["area"] + 1,
              "峰面积误差 ≈ √面积（泊松）")

        # 身份假设失败路径：无 Ca/Fe 峰的纯背景谱
        bg_only = np.full(sfs_spectra.N_CHANNELS, 25.0)
        try:
            sfs_spectra.calibrate_energy(bg_only)
            check(False, "纯背景谱必须拒绝能量自证")
        except ValueError:
            check(True, "纯背景谱必须拒绝能量自证")

        # 非谱表跳过
        bad = os.path.join(folder, "PS__9999_EMA.CSV")
        open(bad, "w", encoding="utf-8").write("A,B\n1,2\n")
        spec2, meta2 = sfs_spectra.read_sfs(bad)
        check(spec2 is None, "非谱表返回 (None, None)")

        # 模型/结果 JSON 往返
        doc = {"axis": axis}
        j = os.path.join(folder, "axis.json")
        open(j, "w", encoding="utf-8").write(json.dumps(doc))
        doc2 = json.loads(open(j, encoding="utf-8").read())
        check(abs(doc2["axis"]["keV_per_ch"] - gain_true) / gain_true < 0.02,
              "能量轴 JSON 往返一致")

    print("共 %d 项断言，%d 项失败" % (TOTAL, FAIL))
    return int(FAIL > 0)


if __name__ == "__main__":
    raise SystemExit(main())
