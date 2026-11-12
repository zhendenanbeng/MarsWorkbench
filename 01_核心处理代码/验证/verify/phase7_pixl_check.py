# -*- coding: utf-8 -*-
"""第七期护栏：SuperCam×PIXL 跨仪器对照（离线断言 + 探针缓存加测）。

防三类回归：
  1. RMC 解析（site/drive）——包括「无帧索引返回 None」和健全性钳制；
  2. 氧化物↔元素映射必须覆盖 PIXL 的 16 种氧化物（漏一个就是并排行丢失）；
  3. 匹配判据与负载口径（site/drive 主判据、不合成数值的 note 随负载下发、
     缺 label 的降级文案）。
"""

from __future__ import annotations

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout

import os
import sys
import types

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = str(_layout.project_layout.APP)
sys.path.insert(0, ROOT)

from instruments import cross, supercam               # noqa: E402

checks = fails = 0


def chk(name, ok, detail=""):
    global checks, fails
    checks += 1
    print("  %s  %s %s" % ("OK  " if ok else "FAIL", name, detail))
    fails += 0 if ok else 1


RMC_XML = """<?xml version="1.0"?>
<Product_Observational xmlns:geom="http://pds.nasa.gov/pds4/geom/v1">
 <geom:Coordinate_Space_Indexed>
  <geom:coordinate_space_frame_type>ROVER_NAV_FRAME</geom:coordinate_space_frame_type>
  <geom:Coordinate_Space_Index><geom:index_id>SITE</geom:index_id>
   <geom:index_value_number>7</geom:index_value_number></geom:Coordinate_Space_Index>
  <geom:Coordinate_Space_Index><geom:index_id>DRIVE</geom:index_id>
   <geom:index_value_number>0</geom:index_value_number></geom:Coordinate_Space_Index>
 </geom:Coordinate_Space_Indexed>
</Product_Observational>
"""


def main():
    import tempfile

    # -- RMC 解析 -----------------------------------------------------------
    with tempfile.TemporaryDirectory() as td:
        good = os.path.join(td, "good.xml")
        open(good, "w", encoding="utf-8").write(RMC_XML)
        chk("site_drive_from_xml 正常解析 (7,0)",
            supercam.site_drive_from_xml(good) == (7, 0))
        open(good, "w", encoding="utf-8").write(
            RMC_XML.replace("17040", "17040").replace(
                "<geom:index_value_number>7<",
                "<geom:index_value_number>17040<"))
        chk("健全性钳制：site=17040 的帧索引判 None",
            supercam.site_drive_from_xml(good) is None)
        chk("label 缺文件返回 None",
            supercam.site_drive_from_xml(os.path.join(td, "nope.xml")) is None)
        open(good, "w", encoding="utf-8").write("<Product/>")
        chk("label 无帧索引返回 None", supercam.site_drive_from_xml(good) is None)

    # -- 氧化物映射完备性 ---------------------------------------------------
    from sherloc import pixl as ipixl
    missing = [e for e in ipixl.ELEMENTS if e not in cross.OXIDE_TO_ELEMENT]
    chk("氧化物↔元素映射覆盖 PIXL 全部 %d 种氧化物" % len(ipixl.ELEMENTS),
        not missing, "缺 %s" % missing if missing else "")
    chk("PIXL 全铁口径 FeO-T 映射 Fe", cross.OXIDE_TO_ELEMENT["FeO-T"] == "Fe")

    # -- 行构造（纯函数） ---------------------------------------------------
    rows = cross._match_rows({"Si": 9, "S": 0}, {"Si": [1.0, 1.0]}, 9,
                             {"Si": {"oxide": "SiO2", "wt": 45.2, "err": 2.3},
                              "S": {"oxide": "SO3", "wt": 1.3, "err": 0.4}})
    si = next(r for r in rows if r["element"] == "Si")
    so = next(r for r in rows if r["element"] == "S")
    chk("并排行：LIBS 证据 + 定量侧同排、未检出侧为 None",
        si["libs"]["points"] == "9/9" and si["apxs"]["wt"] == 45.2
        and so["libs"] is None and so["apxs"]["oxide"] == "SO3")

    # -- 降级：缺 label / 无 site -------------------------------------------
    bare = types.SimpleNamespace(sol="187", seq="s", target="t",
                                 site=None, drive=None, points=[])
    r = cross.pixl_cross(None, bare)
    chk("缺 site/drive 给明确降级说明（不猜）",
        "note" in r and "site/drive" in r["note"] and "payload" not in r)

    # -- 探针缓存加测：Bellegarde 全链路（sol187 × PIXL sol186） -----------
    probe = os.path.normpath(os.path.join(str(_layout.project_layout.legacy_cache_root()),"phase7_probe","sol_00187_libs"))
    if os.path.isdir(probe):
        from instruments import REGISTRY
        b = REGISTRY["supercam"]
        obs = [o for o in b.scan_folder(probe)
               if not getattr(o, "passive", False)][0]
        chk("探针观测读出 site7/drive0（Bellegarde）",
            (obs.site, obs.drive) == (7, 0),
            "site%s drive%s" % (obs.site, obs.drive))
        try:
            res = cross.pixl_cross(b, obs)
            pl = res.get("payload")
            chk("Bellegarde 匹配 PIXL sol186（site/drive 主判据）",
                pl is not None and pl["pixl"]["sol"] == "0186"
                and res["matches"][0]["basis"] == "site_drive")
            if pl:
                si = next(r for r in pl["rows"] if r["element"] == "Si")
                chk("Si 行：LIBS 9/9 检出 × PIXL SiO2>40wt% 并排",
                    si["libs"] and si["libs"]["points"] == "9/9"
                    and si["apxs"]["wt"] > 40)
                chk("口径说明随负载下发（不合成数值）",
                    "不合成数值" in pl["note"])
        except Exception as exc:                          # noqa: BLE001
            chk("Bellegarde 全链路", False, "异常：%s" % str(exc)[:60])
        # -- 报告页（第八期 B 支线） ---------------------------------------
        try:
            from instruments import report as ireport
            import tempfile, shutil
            td = tempfile.mkdtemp()
            try:
                cx = ireport.build_pixl_cross_page(b, [obs], td)
                ok_page = bool(cx) and cx["pairs"] == 1
                text = open(cx["path"], encoding="utf-8").read() if cx else ""
                chk("SuperCam×PIXL 报告页生成（Bellegarde 1 对）", ok_page)
                chk("报告页口径齐（扫描中位/不合成数值/SiO2 行）",
                    "扫描中位" in text and "不合成数值" in text and "SiO2" in text)
            finally:
                shutil.rmtree(td, ignore_errors=True)
        except Exception as exc:                          # noqa: BLE001
            chk("SuperCam×PIXL 报告页生成", False, "异常：%s" % str(exc)[:60])
    else:
        print("  （探针缓存不在本机，Bellegarde 全链路 4 项未测——跑过 phase7 探针下载后生效）")

    print("\n共 %d 项断言，%d 项失败" % (checks, fails))
    return 1 if fails else 0


if __name__ == "__main__":
    raise SystemExit(main())
