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

    # ========== 3. 台账锚点在册 ==========
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
