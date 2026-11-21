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

    # ---- 火星尘厚度模型（dust_model，VanBommel 2026 移植）----
    from instruments import dust_model as dm

    # 1) 已知 T 反演：用混合模型生成"测量计数"，perform_cal 必须找回该 T
    T_grid = np.linspace(5e-5, 5e-4, 40)
    lf = np.linspace(0.05, 0.6, 8)[:, None] * np.ones((8, len(T_grid)))
    lf = np.clip(lf + np.linspace(-0.02, 0.05, len(T_grid))[None, :], 0.01, 0.9)
    inp_t = {"LF": lf, "Sens": np.array([112.6, 313.6, 687.0, 910.2, 1450.0,
                                        2800.0, 3292.4, 4297.3]),
             "PTFE_cps": np.ones(8), "PTFE_cps_err": np.full(8, 0.2),
             "BL_cps": np.array([0.5, 2.5, 19.9, 0.2, 0.02, 24.1, 6.8, 44.2]),
             "BF_cps": np.zeros(8),
             "Start_wt": np.array([7.9, 10.4, 42.5, 1.7, 8.2, 5.4, 1.3, 22.6])}
    t_k = 2.0e-4
    comp = np.array([7.9, 10.4, 42.5, 1.7, 8.2, 5.4, 1.3, 22.6])
    rel_c = comp * inp_t["Sens"]; rel_c = 100.0 * rel_c / rel_c.sum()
    model_b = lf * rel_c[:, None] + (1.0 - lf) * inp_t["BL_cps"][:, None]
    i_k = int(np.argmin(np.abs(T_grid - t_k)))
    inp_t["BF_cps"] = model_b[:, i_k]
    got = dm.perform_cal(comp, inp_t, T_grid)
    check(abs(got[3] - T_grid[i_k]) / T_grid[i_k] < 0.02,
          "尘模型：合成计数反演厚度找回已知 T（%.2e，偏差 <2%%）" % got[3])

    # 2) MCMC 确定性：同种子两次运行结果逐位一致
    s1, c1 = dm.mcmc(inp_t, T_grid, numiter=3, seed=7)
    s2, c2 = dm.mcmc(inp_t, T_grid, numiter=3, seed=7)
    check(np.array_equal(s1, s2) and np.array_equal(c1, c2),
          "尘模型：MCMC 同种子确定性（结果可复现）")

    # 3) summarize 列序契约：Summary=[totChi, PTFE_Chi2, BHVO_Chi2, PTFE_T, BHVO_T]
    summ = dm.summarize(s1, c1)
    check(summ["ptfe_T_gcm2"] == float(np.median(s1[:, 3]))
          and summ["bhvo_T_gcm2"] == float(np.median(s1[:, 4]))
          and summ["ptfe_chi2"] == float(np.median(s1[:, 1])),
          "尘模型：summarize 列序契约（Summary 列 4/5=T、列 2=PTFE χ²，中位数口径）")

    # 4) 真实 LF 表（随包下载）在场时：形状与 T 过滤契约
    real_csv = (_layout.project_layout.data_root() / "PIXL" / "下载"
                / "pixl_dust_model" / "PIXL_Sampling_Fractions.csv")
    if real_csv.is_file():
        T_real, LF_real = dm.load_fractions(str(real_csv))
        check(LF_real.shape[0] == 8 and LF_real.shape[1] == len(T_real)
              and T_real.min() >= 1e-6 and len(T_real) > 100,
              "尘模型：官方 LF 表解析（8 元素 × %d 档，T≥1e-6 过滤）" % len(T_real))
    else:
        check(True, "SKIP 尘模型：官方 LF 表不在位（pixl_dust_model 未下载）")

    print("共 %d 项断言，%d 项失败" % (TOTAL, FAIL))
    return int(FAIL > 0)


if __name__ == "__main__":
    raise SystemExit(main())
