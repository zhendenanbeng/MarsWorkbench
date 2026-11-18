# -*- coding: utf-8 -*-
"""PIQUANT 桥接层护栏（官方开源定量层）。

契约（piquant_bridge.py）：
* status/位置：可执行文件+FM 配置+FM 标定三者齐备才 available；缺件显式
  给 reason，绝不猜路径；
* stdout wt% 解析：氧化物行（含数字名 SiO2/K2O/TiO2/Cr2O3）与元素行分别
  归类，Element sum 与化学计量 O 单独抓取；
* sfs_to_msa：NCOLUMNS=1 每行只写计数值（read_EMSA_PIXL 逐行取前
  numColumns 个值——带通道列会把通道号当计数读入，契约锁定此坑）；
* 官方测试谱全链：PIQUANT 仓库自带 BHVO-2G 面包板谱 × FM 配置/标定，
  经本桥 quant 应复现冻结基线 wt%（确定性构建，±容差）。

环境依赖：PIQUANT 可执行文件 + FM 配置/标定 + 官方测试谱（数据根
PIXL/缓存/piquant*，见 piquant_path.json）。缺件→整体 SKIP 不计失败。
"""
import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0, str(_BootPath(__file__).resolve().parents[3] / "02_图形界面代码"))
import source_paths as _layout

import os
import sys

ROOT = _layout.project_layout.APP
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from instruments import piquant_bridge          # noqa: E402

TOTAL = FAIL = 0
SKIPPED = False


def check(ok, name):
    global TOTAL, FAIL
    TOTAL += 1
    FAIL += not ok
    print(("OK " if ok else "FAIL ") + name, flush=True)


def skip_all(reason):
    global SKIPPED
    SKIPPED = True
    print("SKIP piquant_check —— %s" % reason, flush=True)


# 官方 BHVO-2G 测试谱（FM 配置+FM 标定配对）的冻结基线——本机 MSVC/NMake
# 构建实测值（确定性拟合，两次运行逐位一致）。断言 ±10% 相对容差：
# 换编译器/平台时拟合应高度接近；超出说明构建或桥接层被改坏。
BHVO_BASELINE = {"SiO2": 17.90, "K2O": 0.44, "P2O5": 0.08, "CaO": 11.15,
                 "TiO2": 3.55, "Cr2O3": 0.10, "MnO": 0.28, "FeO-T": 19.29,
                 "Sr": 0.1119}

BHVO_ELEMENTS = "Si_K K_K P_K Ca_K Ti_K Cr_K Mn_K Fe_K Sr_K Ar_I"


def main():
    binp = piquant_bridge.piquant_bin()
    if not binp:
        skip_all("PIQUANT 可执行文件未配置（piquant_path.json / MARSWORKBENCH_PIQUANT）")
        return 0
    td = piquant_bridge.testdata_dir()
    if not td:
        skip_all("PIQUANT 官方测试谱目录缺失（piquant_path.json 的 testdata_dir）")
        return 0
    st = piquant_bridge.status()
    if not st["available"]:
        skip_all(st["reason"])
        return 0

    # 1) status 契约
    check(st["available"] and st.get("config") and st.get("calibration"),
          "status：可执行+FM 配置+FM 标定齐备")

    # 2) stdout 解析契约（官方输出形态片段，含数字名氧化物与元素行）
    sample = """Fit results after 7 iterations, reduced chi sq = 17817.91         live time 9.37 sec.
    SiO2     38.10 %      Si_K   int 7869626.0  coeff   0.9587   re_c      0.0%   abs_err(1s)   4.5655%
     K2O      1.03 %       K_K   int 862610.4  coeff   0.5987   re_c      0.1%   abs_err(1s)   0.4452%
    TiO2      2.02 %      Ti_K   int 7648916.0  coeff   0.9692   re_c      0.0%   abs_err(1s)   0.6188%
   Cr2O3      0.09 %      Cr_K   int 332404.0  coeff   1.2006   re_c      0.3%   abs_err(1s)   0.1051%
   FeO-T     19.10 %      Fe_K   int 56932104.0  coeff   0.9961   re_c      0.0%   abs_err(1s)   1.5010%
      Sr    0.0986 %      Sr_K   int  30661.8  coeff   1.0884   re_c      1.3%   abs_err(1s)   0.0510%
Matrix elements
    O   33.7300 %     (33.73 % from stoichiometry)

    Element sum 84.24 %
"""
    q = piquant_bridge.parse_quant_stdout(sample)
    check(q["oxides"].get("SiO2") == 38.10 and q["oxides"].get("TiO2") == 2.02
          and q["oxides"].get("Cr2O3") == 0.09,
          "解析：数字名氧化物行（SiO2/TiO2/Cr2O3）全归位")
    check(q["oxides"].get("FeO-T") == 19.10,
          "解析：FeO-T 连字名归位")
    check(q["elements"].get("Sr") == 0.0986,
          "解析：元素行（Sr）与氧化物分开")
    check(q.get("element_sum_pct") == 84.24 and q.get("o_stoich_pct") == 33.73,
          "解析：Element sum 与化学计量 O")
    check("Element" not in q["oxides"] and "Matrix" not in q["oxides"]
          and len(q["oxides"]) == 5,
          "解析：非定量行（Matrix/Element sum 行首）不误入")

    # 3) MSA 导出契约：NCOLUMNS=1 每行单值
    import tempfile
    import json as _json
    axis = {"keV_per_ch": 0.007863768115942029, "zero_keV": -0.01283478260869586}
    spec = [100.0] * 4096
    with tempfile.TemporaryDirectory(prefix="mw-piquant-") as folder:
        msa = piquant_bridge.sfs_to_msa(spec, 9.37, axis,
                                        os.path.join(folder, "t.msa"), title="契约测试")
        lines = open(msa, encoding="ascii").read().splitlines()
        data = [ln for ln in lines if ln and not ln.startswith("#")]
        check(len(data) == 4096 and data[0] == "100.0",
              "sfs_to_msa：4096 行单列计数值（无通道列）")
        check("#NPOINTS : 4096" in "\n".join(lines[:16]),
              "sfs_to_msa：NPOINTS/头关键字在位")
        hdr = "\n".join(lines[:20])
        check("XPERCHAN : 7.86377" in hdr and "OFFSET : -12.8348" in hdr,
              "sfs_to_msa：能量轴写为 eV/ch（keV×1000）")
        try:
            piquant_bridge.sfs_to_msa([1.0] * 100, 1.0, axis,
                                      os.path.join(folder, "bad.msa"))
            check(False, "非 4096 通道谱必须拒绝导出")
        except ValueError:
            check(True, "非 4096 通道谱必须拒绝导出")

        # 4) 官方测试谱全链：BHVO-2G × FM 配置/标定 → 冻结基线
        bhvo = os.path.join(td, "Calibration_box_BHVO-2G_28kV_230uA_03_28_2019_bulk_sum.msa")
        if not os.path.isfile(bhvo):
            skip_all("官方 BHVO-2G 测试谱缺失（testdata 目录不完整）")
            return 0
        work = os.path.join(folder, "work")
        r, why = piquant_bridge.run_quant(bhvo, work, BHVO_ELEMENTS)
        check(r is not None, "官方 BHVO 谱全链 quant 成功（%s）"
              % ("通过" if r else why))
        if r is None:
            print("共 %d 项断言，%d 项失败" % (TOTAL, FAIL))
            return 1
        bad = []
        for k, ref in BHVO_BASELINE.items():
            v = r["oxides"].get(k, r["elements"].get(k))
            if v is None or abs(v - ref) > 0.10 * max(abs(ref), 0.01):
                bad.append("%s: %s vs 基线 %s" % (k, v, ref))
        check(not bad, "wt%% 复现冻结基线 ±10%%（%s）" % ("；".join(bad) if bad else "全部吻合"))
        check(r.get("plot_csv") and os.path.isfile(r["plot_csv"]),
              "quant 落盘 plot csv（官方曲线输出）")

    print("共 %d 项断言，%d 项失败" % (TOTAL, FAIL))
    return int(FAIL > 0)


if __name__ == "__main__":
    raise SystemExit(main())
