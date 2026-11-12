# -*- coding: utf-8 -*-
"""新仪器后端自检：用 verify/fixtures 里的 PDS 样例走一遍
扫描整理 → 读点 → 预处理 → 元素线判定。

用法：  python verify/instruments_check.py
退出码非零 = 有失败项。改动 instruments/ 后先跑这个。
"""

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = str(_layout.project_layout.APP)
sys.path.insert(0, ROOT)

from instruments import REGISTRY                                   # noqa: E402

FIX = os.path.join(HERE, "fixtures")


def check(name, fn):
    try:
        fn()
        print("[PASS] %s" % name)
        return True
    except Exception as exc:                                       # noqa: BLE001
        import traceback
        print("[FAIL] %s\n       %s" % (name, traceback.format_exc()[-600:]))
        return False


def check_chemcam():
    b = REGISTRY["chemcam"]
    obs = b.scan_folder(os.path.join(FIX, "chemcam"))
    assert obs, "没有扫描到 ChemCam 观测"
    print("       观测: %s" % [(o.sol, o.seq, o.n_points) for o in obs])
    o = next(o for o in obs if not getattr(o,'passive',False))
    d = b.load_point(o.points[0], None, {"dust": True, "despike": True})
    assert len(d["x"]) > 6000, "通道数不对: %d" % len(d["x"])
    assert d["n_shots"] == 25, "除尘后应剩 25 发，实际 %d" % d["n_shots"]
    assert 240 < d["x"][0] < 250, "波长轴起点不对: %.1f" % d["x"][0]
    assert 890 < d["x"][-1] < 910, "波长轴终点不对: %.1f" % d["x"][-1]
    marks, elems = b.annotate(d["x"], d["y"])
    assert any(m["detected"] for m in marks), "没检出任何线（火星岩应有 Fe/Ti 等）"
    print("       波长 %.1f–%.1f nm, 检出元素: %s" % (d["x"][0], d["x"][-1], elems))
    d1 = b.load_point(o.points[0], 0, {})
    assert d1["n_shots"] == 30 and len(d1["x"]) == len(d1["y"])


def check_supercam():
    b = REGISTRY["supercam"]
    obs = b.scan_folder(os.path.join(FIX, "supercam"))
    assert obs, "没有扫描到 SuperCam 观测"
    print("       观测: %s" % [(o.sol, o.seq, o.product, o.n_points) for o in obs])
    o = next(x for x in obs if x.product == "CL1")   # 被动 CP2/CP3 也可能排前面
    d = b.load_point(o.points[0], None, {"dust": True, "despike": True})
    assert len(d["x"]) > 6000, "通道数不对: %d" % len(d["x"])
    assert 240 < d["x"][0] < 260, "波长轴起点不对: %.1f" % d["x"][0]
    assert 700 < d["x"][-1] < 900, "波长轴终点不对: %.1f" % d["x"][-1]
    marks, elems = b.annotate(d["x"], d["y"])
    print("       波长 %.1f–%.1f nm, 检出元素: %s" % (d["x"][0], d["x"][-1], elems))
    d1 = b.load_point(o.points[0], 3, {})
    assert d1["n_shots"] >= 25


def check_apxs():
    b = REGISTRY["apxs"]
    obs = b.scan_folder(os.path.join(FIX, "apxs"))
    assert obs, "没有扫描到 APXS 观测"
    print("       观测: %s" % [(o.sol, o.seq, o.n_points) for o in obs])
    o = obs[0]
    d = b.load_point(o.points[0], None, {})
    assert 900 < len(d["x"]) < 1100, "通道数不对: %d" % len(d["x"])
    assert 0.4 < d["x"][0] < 0.7, "能量轴起点不对: %.3f" % d["x"][0]
    marks, elems = b.annotate(d["x"], d["y"])
    assert "Fe" in elems or "Ca" in elems, "火星岩应检出 Fe/Ca: %s" % elems
    print("       能量 %.2f–%.2f keV, 检出元素: %s" % (d["x"][0], d["x"][-1], elems))
    ox = b.oxides(o)
    assert ox and len(ox["rows"]) >= 15, "成分表行数不够: %s" % (ox and len(ox["rows"]))
    fe = [r for r in ox["rows"] if r["oxide"].strip() == "FeO"]
    assert fe and 5 < fe[0]["wt"] < 30, "FeO wt% 不合理: %s" % fe


def check_supercam_moc():
    """官方成分预测表：解析、键匹配、面板负载、缺表/缺点的明确口径。"""
    from instruments import supercam as sc
    fixture = os.path.join(FIX, "supercam", "_moc", "supercam_moc.csv")
    rows, meta = sc.load_moc(fixture)
    assert len(rows) == 3, "fixture 应有 3 行成分记录，实际 %d" % len(rows)
    assert meta.get("RMSEP", {}).get("SiO2") == 6.1, "RMSEP 元数据解析不对: %s" % meta
    # 本地主名（含两位版本号）→ 表键（去掉版本号）
    stem = "scam_0998_0755522849_085_cl1_scam02998_aegis_0998a__________06p01"
    assert sc.moc_key(stem) in rows, "版本号剥离匹配失败: %s" % sc.moc_key(stem)
    b = REGISTRY["supercam"]
    # 用 fixture 路径直接构造数据（不依赖软件目录缓存状态）
    rows2, meta2 = b.moc_table(fixture)
    rec = rows2[sc.moc_key(stem)]
    assert abs(rec["oxides"]["SiO2"][0] - 10.27) < 1e-6, "SiO2 值不对"
    assert abs(rec["oxides"]["MgO"][0] - 30.82) < 1e-6, "MgO 值不对"
    assert rec["distance"] == 3035, "距离列解析不对: %s" % rec["distance"]
    # 缺表 / 缺点两种口径
    missing = b.moc_payload("scam_9999_000000000_000_cl1_scam99999_none___01p01")
    assert missing == {"have_table": False} or missing.get("matched") is False, \
        "缓存缺表时的口径不对: %s" % missing


def check_report():
    """质控报告生成器：三个仪器都能出 HTML 汇总 + 逐点 CSV，内容对得上。"""
    import tempfile
    from instruments import report as rpt
    for bid in ("chemcam", "supercam", "apxs"):
        b = REGISTRY[bid]
        obs = b.scan_folder(os.path.join(FIX, bid))
        assert obs, "%s 没有观测" % bid
        with tempfile.TemporaryDirectory() as td:
            res = rpt.build_report(b, obs, td)
            page = open(res["index"], encoding="utf-8").read()
            assert "质控报告" in page and "观测汇总" in page, "%s 报告缺章节" % bid
            assert ("Sol %s" % obs[0].sol) in page, "%s 报告缺观测行" % bid
            assert "论文对照" in page, "%s 报告缺口径说明" % bid
            csv_path = os.path.join(td, res["csv"][0])
            head = open(csv_path, encoding="utf-8-sig").readline()
            assert "检出元素" in head and "SCLK" in head, "%s CSV 表头不对" % bid
            if bid == "supercam":
                assert "官方预测 SiO2 wt%" in head, "SuperCam CSV 缺成分列"
                assert "supercam_libs_moc" in page, "SuperCam 报告缺官方表口径说明"
            if bid == "apxs":
                assert "氧化物成分" in page, "APXS 报告缺成分表"
    # S2：MOC 趋势图（合成数据单测，不依赖多点 fixture）
    from instruments.report import moc_trend_svg
    pairs = [(i, {"SiO2": (40.0 + 5 * i, 1.0), "FeOT": (12.0 - i, 0.5)})
             for i in range(1, 5)]
    tsvg = moc_trend_svg(pairs)
    assert tsvg and tsvg.count("polyline") >= 2 and "查表展示" in tsvg, "趋势图缺内容"
    assert moc_trend_svg(pairs[:1]) == "" and moc_trend_svg([(1, {}), (2, {})]) == ""
    print("       三仪器报告均生成，SuperCam 含官方成分列与趋势图")


def check_supercam_passive():
    """SuperCam 被动/红外光谱（cp2/cp3）：能看图、µm/nm 轴正确、无发射线判定。"""
    b = REGISTRY["supercam"]
    obs = b.scan_folder(os.path.join(FIX, "supercam"))
    kinds = {o.product: o for o in obs}
    assert "CP2" in kinds and "CP3" in kinds, "应扫到 CP2/CP3 观测: %s" % list(kinds)
    d2 = b.load_point(kinds["CP2"].points[0], None, {})
    assert 1.2 < min(d2["x"]) < 1.4 and 2.5 < max(d2["x"]) < 2.7, \
        "CP2 应是 1.31–2.59µm 的被动红外谱"
    assert d2["x"][0] < d2["x"][-1], "波长轴应已排序为升序"
    assert d2["passive"] and d2["no_annotate"]
    assert d2["xlabel"] == "波长 (µm)" and len(d2["x"]) == 256
    assert any("被动" in a for a, _t in d2["log"])
    d3 = b.load_point(kinds["CP3"].points[0], None, {})
    assert d3["passive"] and 370 < min(d3["x"]) < 390, "CP3 应是 379nm 起的被动谱"
    assert len(d3["x"]) > 5000
    assert "CL1" in kinds, "CL1 观测仍在: %s" % list(kinds)
    d1 = b.load_point(kinds["CL1"].points[0], None, {})
    assert not d1.get("passive") and d1["n_shots"] >= 25
    print("       CP2 %.2f–%.2f µm · CP3 %.0f–%.0f nm · CL1 正常"
          % (d2["x"][0], d2["x"][-1], d3["x"][0], d3["x"][-1]))


def check_cross_metadata():
    """T1 跨仪器元数据：TARGET 名与 site/drive 解析进观测；缺字段降级不炸。"""
    from instruments import pds3 as P

    # ChemCam：csv 注释头 TARGET + label RMC（phase3_probe 的 sol01000 是 Newland）
    b = REGISTRY["chemcam"]
    probe = os.path.normpath(os.path.join(str(_layout.project_layout.legacy_cache_root()),"phase3_probe","sol01000"))
    have_probe = os.path.isdir(probe)
    obs = b.scan_folder(probe if have_probe else os.path.join(FIX, "chemcam"))
    assert obs, "ChemCam 无观测"
    hit = [o for o in obs if getattr(o, "target", "")]
    if have_probe:
        # sol01000 探针数据在缓存里才校验真值（不在缓存的机器走降级路径）
        assert hit and hit[0].target == "Newland",             "sol01000 目标名应为 Newland: %s" % [getattr(o, "target", "") for o in obs]
        assert getattr(hit[0], "site", None) == 48 and hit[0].drive == 1570,             "sol01000 site/drive 应为 48/1570"
    # fixtures 真值：sol00013 = Coronation（好奇号首个激光靶）site/drive 3/4；
    # PSV（tab 无 TARGET 行）为空串——缺字段降级不炸
    fx = b.scan_folder(os.path.join(FIX, "chemcam"))
    fx_by = {(str(o.sol), o.seq): o for o in fx}
    cor = fx_by[("13", "ccam02013")]
    assert cor.target == "Coronation", "sol00013 目标应为 Coronation: %r" % cor.target
    assert (cor.site, cor.drive) == (3, 4), "sol00013 site/drive 应为 3/4"
    assert fx_by[("10", "ccam04010")].target == "", "PSV 无 TARGET 行应为空串"

    # APXS：rsp 首行 Target + label RMC（fixtures sol00102 = rocknest3_rp, 5/178）
    a = REGISTRY["apxs"]
    aobs = a.scan_folder(os.path.join(FIX, "apxs"))
    assert aobs and getattr(aobs[0], "target", "") == "rocknest3_rp",         "APXS 目标名解析失败: %r" % getattr(aobs[0], "target", None)
    assert (aobs[0].site, aobs[0].drive) == (5, 178), "APXS site/drive 应为 5/178"

    # 规范化与纯函数
    assert P.target_from_apxs_header('"Target: rocknest3_rp, Contains DAT files: X"') == "rocknest3_rp"
    assert P.normalize_target("RockNest3_RP") == P.normalize_target("rocknest3 rp")
    assert P.target_from_apxs_header("") is None
    assert P.site_drive_from_lbl(None) is None
    print("       ChemCam/APXS 目标名与 site/drive 解析通过（含降级路径）")


def check_coLocated_pair():
    """T2 同址对：Newland(ChemCam) × Wallace_DRT(APXS) 同 site/drive 48/1570。

    探针数据在 .cache/phase5_probe 与 phase3_probe 时才校验（缓存不在则跳过，
    不让护栏依赖必须联网下载过的机器）。
    """
    from instruments import pds3 as P
    cprobe = os.path.normpath(os.path.join(str(_layout.project_layout.legacy_cache_root()),"phase3_probe","sol01000"))
    aprobe = os.path.normpath(os.path.join(str(_layout.project_layout.legacy_cache_root()),"phase5_probe","apxs_sol999"))
    if not (os.path.isdir(cprobe) and os.path.isdir(aprobe)):
        print("       （探针缓存不在本机，同址对校验跳过）")
        return
    c = REGISTRY["chemcam"]
    lib = [o for o in c.scan_folder(cprobe) if not getattr(o, "passive", False)]
    assert lib, "sol01000 应有 LIBS 观测"
    lib = lib[0]
    assert lib.target == "Newland" and lib.n_points == 9
    a = REGISTRY["apxs"]
    aobs = a.scan_folder(aprobe)
    assert aobs, "apb 前缀应可扫描"
    o = aobs[0]
    assert o.target == "Wallace_DRT", "APXS 目标名: %r" % o.target
    assert (lib.site, lib.drive) == (o.site, o.drive) == (48, 1570),         "同址判据失败: %s vs %s" % ((lib.site, lib.drive), (o.site, o.drive))
    ox = a.oxides(o)
    assert len(ox["rows"]) == 16 and 95 < ox["total"] < 102
    # 目标名不同但同址 → 匹配逻辑必须以 site/drive 为主、目标名为辅
    assert P.normalize_target(lib.target) != P.normalize_target(o.target)
    print("       Newland(LIBS 9点) × Wallace_DRT(APXS 16氧化物) 同址 48/1570 判据通过")


def main():
    ok = True
    ok &= check("ChemCam 扫描/读点/预处理/线判定", check_chemcam)
    ok &= check("SuperCam FITS 读取/线判定", check_supercam)
    ok &= check("APXS 谱/峰判定/成分表", check_apxs)
    ok &= check("SuperCam 官方成分表解析/匹配", check_supercam_moc)
    ok &= check("三仪器质控报告生成", check_report)
    ok &= check("SuperCam 被动/红外光谱 cp2/cp3", check_supercam_passive)
    ok &= check("跨仪器元数据（目标名/停车位）", check_cross_metadata)
    ok &= check("同址配对闸门（Newland×Wallace_DRT）", check_coLocated_pair)
    n = 6
    print("\n共 %d 项断言，%d 项失败" % (n, 0 if ok else 1))
    print("%s" % ("全部通过" if ok else "存在失败项"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
