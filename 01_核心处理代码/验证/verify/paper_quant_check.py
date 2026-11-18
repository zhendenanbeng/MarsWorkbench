# -*- coding: utf-8 -*-
"""文献定量对账护栏：已发表定量数字 vs 盘上数据/管线的独立复算。

对账锚点（全部来自可开放获取的已发表定量结果）：
  1. PIXL Bellegarde（sol 187 PIXL 扫描）——Schmidt et al. 2025 Science Advances
     （doi:10.1126/sciadv.adr2613；NTRS 20250000314 接收稿含 SOM）Table 1 低盐校正
     火成岩主体成分。复算：盘上 RQC（2809 PMC，与论文 Bulk #PMC 一致）→ 论文判据
     SO3<4.85/Cl<3.0 选低盐（复得 2404 vs 论文 2411）→ 论文盐校正算法（减 CaSO4/
     NaCl 后归一）→ 与 Table 1 ±σ 对比。
     冻结结果：8/10 氧化物在 ±σ 内；FeO-T（+1.21σ）与 Na2O（−0.92σ）系统性超差——
     论文按双探测器择低加谱再定量（非逐 PMC 均值）且 RQC 版本演进，属诚实记录的
     版本差，不是数据损坏（见 data_audit 的逐 PMC 面积精确复现）。
  2. APXS Rocknest3_rp（sol 102 基岩）——Schmidt et al. 2013 LPSC（NTRS 20130009734）
     Table 1：SO3 4.0 wt%；LPSC 2014 #1504：Rocknest3 类 K2O 1.7-1.9 wt%。
     附目标身份防混淆：Morris 2015 LPSC（NTRS 20150001942）Table 1 "RN bulk" 是
     Rocknest 风成土壤（SiO2 41.18/MgO 8.33/K2O 0.47），不是 rocknest3_rp 岩石。

数据不在位时逐项 SKIP（发行机不强制）；需联网的文献值不在此处重下载，
数字以台账/本脚注为准，来源 sha256 见 published_calibrations.json。
"""
import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0, str(_BootPath(__file__).resolve().parents[3] / "02_图形界面代码"))
import source_paths as _layout

import csv
import json
import sys
from pathlib import Path

import numpy as np

import project_layout

BASE = _layout.project_layout.VERIFY / "fixtures/历史数值基准"
TOTAL = FAIL = 0


def check(ok, name):
    global TOTAL, FAIL
    TOTAL += 1
    FAIL += not ok
    print(("OK  " if ok else "FAIL ") + name, flush=True)


def skip(name):
    print("SKIP " + name, flush=True)


# Schmidt 2025 Table 1（Bellegarde 低盐校正成分 wt%±σ）
LIT_BELLEGARDE = {
    "SiO2": (48.5, 2.3), "Al2O3": (7.9, 0.4), "FeO-T": (22.9, 1.1),
    "MgO": (2.0, 0.6), "CaO": (7.2, 0.4), "Na2O": (3.9, 0.3),
    "K2O": (1.2, 0.4), "TiO2": (2.8, 0.7), "MnO": (0.4, 0.2),
    "P2O5": (3.1, 0.6),
}
# 版本差超差项的冻结偏差与允许漂移带（诚实记录：|dev| ≤ 2.0 绝对 wt%）
LIT_EXCEPTION_BAND = 2.0
N_PMC_BULK = 2809
N_PMC_LOW_SALT = 2411
LOW_SALT_COUNT_TOL = 20


def rqc_records(rqc_path):
    return list(csv.DictReader(open(rqc_path, encoding="utf-8-sig", newline="")))


def mean_over(recs, cols):
    out = {}
    for c in cols:
        arr = []
        for r in recs:
            v = r.get(c + "_wt%")
            if v:
                try:
                    arr.append(float(v))
                except ValueError:
                    pass
        if arr:
            out[c] = float(np.mean(arr))
    return out


def salt_correct(mean_raw):
    """Schmidt 2025 Table 1 口径：减 SO3 为 CaSO4、Cl 为 NaCl，去列后归一到 100%。"""
    so3, cl = mean_raw.get("SO3", 0.0), mean_raw.get("Cl", 0.0)
    out = dict(mean_raw)
    out["CaO"] = out["CaO"] - so3 * 56.08 / 80.06
    out["Na2O"] = out["Na2O"] - cl * 61.98 / 70.90
    out.pop("SO3", None)
    out.pop("Cl", None)
    s = sum(out.values())
    return {k: x / s * 100.0 for k, x in out.items()}


def main():
    data_root = Path(project_layout.data_root())
    rqc_p = next(data_root.glob("PIXL/随包示例/*rqc*.csv"), None) if data_root.is_dir() else None

    # ========== 1. PIXL Bellegarde（Schmidt 2025）==========
    if not (rqc_p and rqc_p.is_file()):
        skip("PIXL Bellegarde 对账（RQC 样例不在位）")
    else:
        rows = rqc_records(rqc_p)
        check(len(rows) == N_PMC_BULK,
              "RQC PMC 总数 %d 与论文 Bulk #PMC %d 一致（Bellegarde 扫描集）"
              % (len(rows), N_PMC_BULK))
        low = []
        for r in rows:
            try:
                so3 = float(r["SO3_wt%"])
                cl = float(r["Cl_wt%"])
            except (ValueError, KeyError):
                continue
            if so3 < 4.85 and cl < 3.0:
                low.append(r)
        check(abs(len(low) - N_PMC_LOW_SALT) <= LOW_SALT_COUNT_TOL,
              "低盐 PMC 复得 %d vs 论文 %d（差 ≤%d）"
              % (len(low), N_PMC_LOW_SALT, LOW_SALT_COUNT_TOL))
        cols = list(LIT_BELLEGARDE) + ["SO3", "Cl"]
        corr = salt_correct(mean_over(low, cols))
        n_in = n_out_band_ok = 0
        for ox, (ref, sig) in LIT_BELLEGARDE.items():
            dev = corr[ox] - ref
            if abs(dev) <= sig:
                n_in += 1
            elif abs(dev) <= LIT_EXCEPTION_BAND:
                n_out_band_ok += 1
                print("   %-6s 校正复算 %.2f vs 文献 %.2f（%+.2f，σ=%.1f，版本差记录）"
                      % (ox, corr[ox], ref, dev, sig))
            else:
                print("   %-6s 校正复算 %.2f vs 文献 %.2f（%+.2f，超版本差带！）"
                      % (ox, corr[ox], ref, dev))
        check(n_in >= 7, "盐校正主体成分 ≥7/10 氧化物落在文献 ±σ（现 %d）" % n_in)
        check(n_in + n_out_band_ok == len(LIT_BELLEGARDE),
              "全部 10 氧化物在冻结容差内（±σ 或版本差带 ±%.1f）" % LIT_EXCEPTION_BAND)
        # Mg# 交叉核验（论文 13.6）
        mg = corr["MgO"] / 40.304
        fe = corr["FeO-T"] / 71.844
        mgp = mg / (mg + fe) * 100.0
        check(abs(mgp - 13.6) <= 2.0,
              "Mg# 复算 %.1f vs 论文 13.6（±2.0）" % mgp)

    # ========== 1b. PIXL PIQUANT × RQC 对账（官方开源定量层，缓存）==========
    pq_p = data_root / "PIXL/缓存/piquant_vs_rqc.json" if data_root.is_dir() else None
    if not (pq_p and pq_p.is_file()):
        skip("PIXL PIQUANT×RQC 对账（piquant_vs_rqc.json 不在位——02 号脚本未跑或环境未部署）")
    else:
        pq = json.loads(pq_p.read_text(encoding="utf-8"))
        st = pq.get("stats", {})
        check(pq.get("detector") == "ENB" and pq.get("n_failed") == 0
              and pq.get("n_quantified", 0) == len(pq.get("rows", [])) >= 250,
              "PIQUANT 对账规模：ENB 逐 PMC %d 点零失败（步长 %s）"
              % (pq.get("n_quantified", 0), pq.get("stride")))
        majors = ("FeO-T", "SiO2", "CaO", "Al2O3")
        ok_m = all(o in st and st[o]["corr"] >= 0.97
                   and abs(st[o]["median_rel_dev_pct"]) <= 5.0 for o in majors)
        check(ok_m, "四大主量（FeO-T/SiO2/CaO/Al2O3）：相关 ≥0.97 且中位偏差 ≤5%%"
                    "（%s）" % "; ".join("%s %+.1f%%/%.3f" % (o, st[o]["median_rel_dev_pct"],
                                                            st[o]["corr"])
                                        for o in majors if o in st))
        ok_all = all(s["corr"] >= 0.75 for s in st.values()) and len(st) == 10
        check(ok_all, "全部 10 氧化物相关 ≥0.75（最弱 %s）"
              % min(st.items(), key=lambda kv: kv[1]["corr"])[0] if st else "空")
        check(pq.get("scan_dir") == "sfs_sol0186" and "rqc" in pq,
              "对账扫描与基准留痕（%s × %s）" % (pq.get("scan_dir"), pq.get("rqc")))
        # 溯源一致性：JSON 里的 exe 指纹 = piquant_path.json 注册指纹
        ppj = data_root / "PIXL/缓存/piquant_path.json"
        if ppj.is_file():
            reg = json.loads(ppj.read_text(encoding="utf-8")).get("sha256", {})
            check(reg.get("Piquant.exe") == pq.get("piquant_bin_sha256"),
                  "PIQUANT 可执行指纹与注册表一致（防静默换版）")
        else:
            skip("PIQUANT 指纹核对（piquant_path.json 不在位）")

    # ========== 2. APXS Rocknest3_rp（Schmidt 2013 / LPSC 2014）==========
    fix = _layout.project_layout.VERIFY / "fixtures/apxs/sol00102"
    rsp = fix / "apa_406545731rsp01020050178_______p1.csv"
    rwp = fix / "apa_406545731rwp01020050178_______p1.csv"
    if not (rsp.is_file() and rwp.is_file()):
        skip("APXS Rocknest3_rp 文献身份对账（fixture 不在位）")
    else:
        header = rsp.open(encoding="utf-8", errors="replace").readline()
        check("rocknest3_rp" in header.lower(),
              "rsp 首行目标名含 rocknest3_rp（基岩身份）")
        ox = {}
        for line in rwp.read_text(encoding="utf-8").splitlines()[1:]:
            parts = [p.strip() for p in line.split(",")]
            if len(parts) >= 2 and parts[0]:
                try:
                    ox[parts[0]] = float(parts[1])
                except ValueError:
                    continue
        so3 = ox.get("SO3")
        k2o = ox.get("K2O")
        check(so3 is not None and 3.5 <= so3 <= 4.5,
              "SO3 %.2f 在 Schmidt 2013 Table 1 带 [3.5, 4.5]" % (so3 or -1))
        check(k2o is not None and 1.7 <= k2o <= 1.9,
              "K2O %.2f 在 LPSC 2014 Rocknest3 类带 [1.7, 1.9]" % (k2o or -1))
        # 身份防混淆：rocknest3_rp 是高碱基岩，不是 K2O≈0.5 的 Rocknest 土壤
        check(k2o is not None and k2o > 1.2,
              "K2O %.2f 与 Morris 2015 RN 土壤 bulk 0.47 明确区分（>1.2）" % (k2o or -1))
        check(abs(ox.get("MgO", -1) - 8.33) > 2.0,
              "MgO %.2f 与土壤 bulk 8.33 明确区分（|差|>2.0）" % ox.get("MgO", -1))

    # ========== 3. PIXL 比值矩阵（01_sfs_pipeline 生成，真轴 v4）==========
    rec = json.loads((BASE / "pixl_ratio_reconciliation.json").read_text(encoding="utf-8"))
    mat = rec.get("ratio_matrix") or {}
    pairs = {p["pair"]: p["loglog_corr"] for p in rec.get("ratio_matrix", {}).get("top8", [])}
    check(mat.get("n_pairs") == 66,
          "比值矩阵 66 对在册（现 %s，12 元素两两）" % mat.get("n_pairs"))
    TOP8 = {"Si/Fe": 0.9941, "Si/S": 0.9888, "Ca/S": 0.9860, "Ca/Fe": 0.9843,
            "Al/S": 0.9834, "Si/Ca": 0.9831, "Al/Fe": 0.9812, "Ti/S": 0.9796}
    bad_top = [k for k, v in TOP8.items()
               if pairs.get(k) is None or pairs[k] < v - 0.01]
    check(not bad_top, "Top8 比值相关不低于冻结值-0.01（%s）"
          % ("；".join("%s %.4f<%.4f" % (k, pairs.get(k, -1), TOP8[k] - 0.01)
                       for k in bad_top) or "全部达标"))
    check((mat.get("min_corr") or 0) >= 0.45,
          "全矩阵最差相关 ≥0.45（现 %s，弱尾=Na 低能边/远 Z 对）" % mat.get("min_corr"))
    # 扫描内半定量标定（wt=k×area，过原点）：强元素 R² 冻结下限
    inscan = rec.get("in_scan_calibration") or {}
    INSCAN_MIN = {"Ca": 0.975, "S": 0.975, "Ti": 0.97, "K": 0.95,
                  "Fe": 0.92, "Si": 0.89}
    bad_cal = [(el, (inscan.get(el) or {}).get("r2", 0))
               for el, lo in INSCAN_MIN.items()
               if (inscan.get(el) or {}).get("r2", 0) < lo]
    check(not bad_cal and len(inscan) == 12,
          "扫描内标定 R^2 达标（Ca %.4f/S %.4f/Ti %.4f/Fe %.4f/Si %.4f；%d 元素，Na 不可用如实记录）"
          % tuple([(inscan.get(el) or {}).get("r2", 0)
                   for el in ("Ca", "S", "Ti", "Fe", "Si")] + [len(inscan)]))
    check(not bad_cal, "扫描内标定 R^2 低于冻结下限: %s"
          % ("；".join("%s %.4f<%.2f" % (el, v, INSCAN_MIN[el]) for el, v in bad_cal)
             or "无"))
    check(all((inscan.get(el) or {}).get("k", 0) > 0 for el in INSCAN_MIN),
          "扫描内标定斜率全部为正（面积-浓度方向一致）")

    # ========== 3b. APXS Si 内标比值层（02_ratio_calibration 生成，数据在位时）==========
    rc_p = _layout.project_layout.data_root() / "APXS" / "缓存" / "apxs_ratio_calib.json"
    if not rc_p.is_file():
        skip("APXS 比值层断言（apxs_ratio_calib.json 不在位）")
    else:
        rc = json.loads(rc_p.read_text(encoding="utf-8"))
        els = rc.get("elements") or {}
        RATIO_MIN = {"Ca": 0.97, "Fe": 0.94, "K": 0.92, "Mg": 0.88}
        bad_rc = [(el, (els.get(el) or {}).get("corr_all", 0))
                  for el, lo in RATIO_MIN.items()
                  if (els.get(el) or {}).get("corr_all", 0) < lo]
        ca = els.get("Ca") or {}
        check(not bad_rc,
              "APXS Si 内标比值相关达标（Ca %.4f/Fe %.4f/K %.4f/Mg %.4f）"
              % tuple((els.get(el) or {}).get("corr_all", 0)
                      for el in ("Ca", "Fe", "K", "Mg")))
        check(not bad_rc, "APXS 比值相关低于冻结下限: %s"
              % ("；".join("%s %.4f" % (el, v) for el, v in bad_rc) or "无"))
        check((ca.get("val_rel_err_median") or 1) <= 0.05
              and 0.5 <= (ca.get("power_b") or 0) <= 1.5,
              "APXS Ca 留出中位相对误差 <=5%%（%.1f%%）且幂律指数 b∈[0.5,1.5]（%.3f）"
              % ((ca.get("val_rel_err_median") or 1) * 100, ca.get("power_b") or 0))
        check(abs((rc.get("n_targets") or 0) - 25) == 0
              and (rc.get("n_fit") or 0) + (rc.get("n_val") or 0) == 25,
              "APXS 比值层 25 靶齐（拟合 %s + 留出 %s）"
              % (rc.get("n_fit"), rc.get("n_val")))

    # ========== 4. APXS 峰位自证（能量轴一致性）==========
    from instruments import apxs as apxs_backend
    rsp2 = _layout.project_layout.VERIFY / "fixtures/apxs/sol00102/apa_406545731rsp01020050178_______p1.csv"
    if not rsp2.is_file():
        skip("APXS 峰位自证（fixture 不在位）")
    else:
        from instruments import apxs_fit
        data = apxs_backend.pds3.read_table(str(rsp2), skiprows=1)
        energy = data[:, 0] / 1000.0
        counts = data[:, 1] if data.shape[1] == 2 else np.median(data[:, 1:-1], axis=1)
        re = apxs_fit.spectrum_areas(energy[1:-1], counts[1:-1])
        dev = []
        for el in apxs_fit.KA_KEV:
            r = re.get(el)
            if r and r.get("center_kev") is not None:
                dev.append((el, abs(r["center_kev"] - apxs_fit.KA_KEV[el])))
        # 峰位=argmax（道分辨率 ~24 eV）；Mg 受相邻 Na Kα 拖拽 ~2 道属物理合理。
        # 覆盖全部有拟合结果的元素（主量 7 + S/Cl/Ti/Cr/Mn/Ni/Zn，实测最大 0.053）
        elmax = max(dev, key=lambda x: x[1]) if dev else ("-", -1)
        check(len(dev) >= 14 and elmax[1] <= 0.06,
              "APXS %d 元素拟合峰位与理论 Kα 偏差 ≤0.06 keV（最大 %s %.3f，道分辨 ~24eV）"
              % (len(dev), elmax[0], elmax[1]))

    # ========== 5. 台账锚点在册 ==========
    gi = json.loads((BASE / "gold_intervals.json").read_text(encoding="utf-8"))
    check("bellegarde_schmidt2025_corrected" in gi.get("anchors", {}).get("pixl", {}),
          "gold_intervals.pixl.bellegarde_schmidt2025_corrected 锚点在册")
    check("rocknest3_schmidt2013_identity" in gi.get("anchors", {}).get("apxs", {}),
          "gold_intervals.apxs.rocknest3_schmidt2013_identity 锚点在册")
    cal = json.loads((BASE / "published_calibrations.json").read_text(encoding="utf-8"))
    keys = {s.get("key") for s in cal.get("sources", [])}
    check({"schmidt2025", "schmidt2014"} <= keys,
          "published_calibrations 来源含 schmidt2025/schmidt2014（带 DOI+sha256）")

    print("共 %d 项断言，%d 项失败" % (TOTAL, FAIL))
    return int(FAIL > 0)


if __name__ == "__main__":
    sys.exit(main())
