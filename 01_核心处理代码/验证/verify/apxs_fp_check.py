# -*- coding: utf-8 -*-
"""APXS 基本参数法-lite 护栏：物理表、标定冻结、确定性复算与引擎单元自检。

对账口径（03_fp_calibration.py 台账，fixtures/apxs/fp_calibration_ledger.json
冻结副本离线可复算）：
  * 69 自然靶（45 apa + 24 apb，sol 46-4216）标定 u_i（逐元素有效灵敏度，
    α 激发+荧光产额+探测器效率+XRF 折入），70/30 切分；
  * 留出验证：Sherman-lite（Kβ 串扰剥离 + 基体吸收修正 + 化学计量 O）
    预测氧化物 wt% vs 官方 rwp，主量五氧化物中位相对误差 1.9-9.5%，
    全部 ≤ 无吸收修正基线（冻结事实，逐项断言）；
  * 负结果留档：逐元素灵敏度时代漂移模型 b/σ≤1.7 不显著，弃用
    （实验见 03 号脚本试跑记录）；Na2O/K2O/Zn 留出误差 FP 略劣于基线，
    按适用边界如实入册，不作成功声明。
"""
import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0, str(_BootPath(__file__).resolve().parents[3] / "02_图形界面代码"))
import source_paths as _layout            # noqa: E402

import json
import sys
from pathlib import Path

import numpy as np

import project_layout                      # noqa: E402
from instruments import apxs_fp            # noqa: E402

HERE = Path(__file__).resolve().parent
LEDGER = HERE / "fixtures/apxs/fp_calibration_ledger.json"
INST = Path(apxs_fp.__file__).resolve().parent
TOTAL = FAIL = 0


def check(ok, name):
    global TOTAL, FAIL
    TOTAL += 1
    FAIL += not ok
    print(("OK  " if ok else "FAIL ") + name, flush=True)


def skip(name):
    print("SKIP " + name, flush=True)


def matrix_from_official(official):
    els, o_wt = apxs_fp.elements_from_oxides(official)
    comp = {k: v / 100.0 for k, v in els.items()}
    comp["O"] = o_wt / 100.0
    return comp, els


def refit_sensitivities(dataset, fit_idx):
    """独立复算 u_i（与 03 号脚本同口径，实现独立成对账）。"""
    sens = {}
    for el in apxs_fp.FP_ELEMENTS:
        e_kev = apxs_fp._tables()["line_energies_kev"][el]
        uvals = []
        for i in fit_idx:
            rec = dataset[i]
            ox = rec["official"].get(apxs_fp.OXIDE_OF[el][0])
            a = apxs_fp.correct_kbeta(rec["areas"]).get(el)
            if not ox or not a or a <= 0:
                continue
            comp, els = matrix_from_official(rec["official"])
            if els.get(el, 0) <= 0:
                continue
            uvals.append(a * apxs_fp.matrix_mass_atten(comp, e_kev) / els[el])
        if len(uvals) >= 8:
            sens[el] = float(np.median(uvals))
    return sens


def main():
    # ========== A. 冻结资产在册 ==========
    t = apxs_fp._tables()
    grid = t["energies_kev"]
    check(len(grid) == 14 and t["matrix_elements"] == list(apxs_fp._tables()["matrix_elements"])
          and len(t["cs_total"]) == 17,
          "物理表：14 能量网格 + 17 基体元素 CS_Total 在册")
    kb = t.get("kbeta") or {}
    ok_kb = all(el in kb and 0.11 <= kb[el]["ratio"] <= 0.15
                for el in ("K", "Ca", "Cr", "Mn", "Fe"))
    check(ok_kb, "Kβ 分支比在册且物理合理（K/Ca/Cr/Mn/Fe 0.11-0.15）")

    cal_p = INST / "apxs_fp_calibration.json"
    if not cal_p.is_file():
        skip("运行期标定不在册（apxs_fp_calibration.json）——其余 A/C 组跳过")
        print("共 %d 项断言，%d 项失败" % (TOTAL, FAIL))
        return int(FAIL > 0)
    cal = json.loads(cal_p.read_text(encoding="utf-8"))
    sens_f = cal.get("sensitivities") or {}
    check(len(sens_f) == 14 and all(v > 0 for v in sens_f.values()),
          "运行期标定：14 元素灵敏度在册且为正")
    check("sol" in cal.get("calib_era", "") and cal.get("n_targets") == 69,
          "运行期标定：时代口径与靶数（69 靶）在册")
    check(len(cal.get("holdout") or {}) == 13,
          "留出验证 13 氧化物误差台账在册")

    # ========== B. 引擎单元自检（合成，恒跑） ==========
    # 合成成分（不含 K/Cr/Mn——避开 Kβ 串扰链对合成前向的反演干扰）
    # → 前向面积（u,μ 同一模型）→ solve 应不动点恢复
    comp0 = {"Si": 24.0, "Fe": 12.0, "Ca": 6.0, "Mg": 8.0, "Al": 8.5,
             "Na": 2.0, "S": 3.0, "Cl": 0.8, "Ti": 0.7, "Zn": 0.1}
    fsum = sum(c / apxs_fp.OXIDE_OF[el][1] for el, c in comp0.items())
    ox_true = {el: 100.0 * c / apxs_fp.OXIDE_OF[el][1] / fsum
               for el, c in comp0.items()}
    el_wt_true = {el: 100.0 * c / fsum for el, c in comp0.items()}
    o_true = 100.0 - sum(el_wt_true.values())
    comp_m = dict(el_wt_true)
    comp_m["O"] = o_true
    comp_m = {k: v / 100.0 for k, v in comp_m.items()}
    areas0 = {}
    for el, c in comp0.items():
        e_kev = t["line_energies_kev"][el]
        areas0[el] = (100.0 * ox_true[el] * apxs_fp.OXIDE_OF[el][1]
                      * sens_f[el] / apxs_fp.matrix_mass_atten(comp_m, e_kev))
    res = apxs_fp.solve(areas0, sens_f)
    ok_rec = res.get("available") and all(
        abs(res["oxides"].get(apxs_fp.OXIDE_OF[el][0], -99) - ox_true[el]) <= 0.5
        for el in comp0)
    check(ok_rec, "solve 合成自洽恢复：10 元素氧化物 ≤0.5 wt%% 不动点复现")
    check(res.get("converged") and res.get("iterations", 999) <= 30,
          "Sherman-lite 迭代收敛（合成用例 %s 轮）" % res.get("iterations"))
    tot = sum(res["oxides"].values())
    # 输出按 3 位小数舍入：契约按 ±0.01 断言（14 行舍入残差上界 ~0.014）
    check(abs(tot - 100.0) <= 0.02 and res["o_wt"] > 0,
          "归一化契约：Σ氧化物=100±0.02（实测 %.6f，3 位小数舍入输出）且 O>0" % tot)
    res2 = apxs_fp.solve({"Si": areas0["Si"]}, sens_f)
    check(not res2.get("available") and "拒绝" in res2.get("reason", ""),
          "安全降级：Si 缺席拒绝求解（不硬算）")

    # Kβ 剥离方向与幅度：Cr 主导 → Mn 净面积 -= Kβ/Cα×Cr；
    # Mn 保持正值时其 Kβ 继续污染 Fe 窗（链式）；Mn 修负后 Kβ 贡献≈0 跳过
    a_raw = {"Cr": 10.0, "Mn": 5.0, "Si": 20.0, "Fe": 8.0, "Ca": 5.0,
             "Mg": 4.0, "Al": 6.0, "Na": 1.0, "Ti": 0.4,
             "S": 2.0, "Cl": 0.5, "Zn": 0.1, "Ni": 0.1}
    a_cor = apxs_fp.correct_kbeta(a_raw)
    check(a_cor["Mn"] < a_raw["Mn"] and abs(a_raw["Mn"] - a_cor["Mn"]
          - kb["Cr"]["ratio"] * a_raw["Cr"]) < 1e-9,
          "Kβ 剥离方向与幅度：Mn 净面积 -= Kβ/Cα×Cr（%.3f 剥离 %.3f）"
          % (kb["Cr"]["ratio"], a_raw["Mn"] - a_cor["Mn"]))
    d_fe = a_raw["Fe"] - a_cor["Fe"]
    mn_cor = a_raw["Mn"] - kb["Cr"]["ratio"] * a_raw["Cr"]
    check(abs(d_fe - kb["Mn"]["ratio"] * mn_cor) < 1e-9 and d_fe > 0,
          "Kβ 剥离链：Fe 窗 -= Kβ/Cα×(Mn 修正后面积)（剥离 %.4f）" % d_fe)
    a_neg = apxs_fp.correct_kbeta({"Cr": 10.0, "Mn": 0.5, "Si": 20.0,
                                   "Fe": 8.0, "Ca": 5.0, "Mg": 4.0,
                                   "Al": 6.0, "Na": 1.0, "Ti": 0.4,
                                   "S": 2.0, "Cl": 0.5, "Zn": 0.1, "Ni": 0.1})
    check(a_neg["Fe"] == 8.0,
          "Kβ 修负语义：Mn 修正后 ≤0 → 其 Kβ 贡献≈0，Fe 窗不再剥离")
    # 基体物理方向：Fe 富集基体在 1.74 keV 的 μ 应显著大于长英质
    mu_femrich = apxs_fp.matrix_mass_atten(
        {"Si": 0.20, "Fe": 0.55, "O": 0.25}, 1.74)
    mu_felsic = apxs_fp.matrix_mass_atten(
        {"Si": 0.47, "Al": 0.15, "K": 0.05, "Na": 0.03, "O": 0.30}, 1.74)
    check(mu_femrich > mu_felsic * 1.3,
          "基体吸收物理方向：Fe 富集 μ(1.74keV)=%.0f > 1.3×长英质 %.0f"
          % (mu_femrich, mu_felsic))

    # ========== C. 冻结台账独立复算 ==========
    if not LEDGER.is_file():
        skip("冻结台账不在位（fixtures/apxs/fp_calibration_ledger.json）")
    else:
        led = json.loads(LEDGER.read_text(encoding="utf-8"))
        ds = led["dataset"]
        ds.sort(key=lambda r: (r["sol"], r["prefix"]))
        n = len(ds)
        fit_idx = [i for i in range(n) if i % 3 != 2]
        val_idx = [i for i in range(n) if i % 3 == 2]
        sens_r = refit_sensitivities(ds, fit_idx)
        ok_s = set(sens_r) == set(sens_f) and all(
            abs(sens_r[el] / sens_f[el] - 1) <= 1e-9 for el in sens_f)
        check(ok_s, "确定性复算：台账数据集独立重拟 u_i == 冻结灵敏度（14 元素 ±1e-9）")

        errs = {}
        n_conv, n_ok = 0, 0
        for i in val_idx:
            rec = ds[i]
            p = apxs_fp.solve(rec["areas"], sens_f)
            if not p.get("available"):
                continue
            n_conv += p.get("converged", False)
            n_ok += 1
            for ox, pv in p["oxides"].items():
                tv = rec["official"].get(ox)
                if tv and tv > 0.2:
                    errs.setdefault(ox, []).append(abs(pv - tv) / tv)
        h_f = cal.get("holdout") or {}
        ok_h = True
        for ox, rows in errs.items():
            med = float(np.median(rows))
            frozen = (h_f.get(ox) or {}).get("median_rel_err")
            if frozen is None or abs(med - frozen) > 1e-9:
                ok_h = False
        check(ok_h and n_ok == len(val_idx),
              "留出复算：%d/%d 靶求解，逐氧化物中位误差 == 冻结台账（±1e-9）"
              % (n_ok, len(val_idx)))
        check(n_conv == n_ok, "Sherman 迭代：留出靶全部收敛（%d/%d）"
              % (n_conv, n_ok))
        majors = ("SiO2", "Al2O3", "FeO", "MgO", "CaO")
        ok_major = all((h_f.get(ox) or {}).get("median_rel_err", 1) <= 0.10
                       for ox in majors)
        check(ok_major, "科学界：主量五氧化物留出中位误差 ≤10%%（冻结：%s）"
              % "、".join("%.1f%%" % (100 * h_f[ox]["median_rel_err"])
                          for ox in majors))
        base_f = led.get("holdout_baseline_no_absorption") or {}
        ok_gain = all((h_f.get(ox) or {}).get("median_rel_err", 9)
                      <= (base_f.get(ox) or {}).get("median_rel_err", 0) + 1e-12
                      for ox in majors)
        check(ok_gain, "增量契约：主量五氧化物 FP 全部 ≤ 无吸收修正基线（吸收修正净增量）")
        rej = led.get("rejected") or []
        check(len(rej) >= 1 and any("live time" in r.get("reason", "")
                                    for r in rej),
              "拒收台账诚实留档：live-time 异常靶显式拒收（%d 条）" % len(rej))
        rl = led.get("holdout_ratio_layer") or {}
        check(bool(rl), "比值层对照在册（%d 元素，比值口径；FP 为氧化物绝对口径）"
              % len(rl))

    print("共 %d 项断言，%d 项失败" % (TOTAL, FAIL))
    return int(FAIL > 0)


if __name__ == "__main__":
    sys.exit(main())
