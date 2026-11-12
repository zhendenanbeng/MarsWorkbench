# -*- coding: utf-8 -*-
"""上色链路的回归护栏（独立于拉曼判定的 regression.py）。

断言三件事：
  A. Alfalfa 配对的发现路径——从 ACI 图名出发，能走完
     清单 → 候选 → 几何共视，且找到的 standoff 在物理窗口内；
  B. 配准判据——对 Alfalfa 的已知 GO 对，锚点数/残差/旋转角
     必须落在验证时的区间里（跑了就是几分钟）；
  C. 优雅降级——没装 opencv 时 colorize 必须报"未安装"，而不是抛异常。

用法：python verify\\colorize_regression.py   （通过返回 0）
网络不可用时 A 类自动跳过（标记 SKIP），B/C 仍执行。
"""

from __future__ import annotations

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = str(_layout.project_layout.APP)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

DATA = os.path.join(os.path.dirname(ROOT), "03_下载数据/SHERLOC/练习数据")
ALFALFA_ACI = "SC3_0370_0699816293_742FDR_N0110108SRLC10600_0000LMJ04.IMG"

checks = fails = 0


def check(name, cond, detail=""):
    global checks, fails
    checks += 1
    if cond:
        print("  OK  %s %s" % (name, detail))
    else:
        fails += 1
        print("  FAIL  %s %s" % (name, detail))


def main():
    from sherloc import watson
    from sherloc import colorize as C

    # ---- C 类先跑（最快，不依赖网络）----
    print("优雅降级")
    saved = C.HAS_CV2
    try:
        C.HAS_CV2 = False
        import numpy as np
        r = C.colorize(np.zeros((100, 100), np.uint8),
                       np.zeros((100, 100, 3), np.uint8),
                       {"C": [0, 0, 0], "A": [0, 0, 1], "H": [1, 0, 0],
                        "V": [0, 1, 0], "ns": 100, "nl": 100},
                       {"C": [0, 0, 1], "A": [0, 0, 1], "H": [1, 0, 0],
                        "V": [0, 1, 0], "ns": 100, "nl": 100}, 100.0)
        check("无 cv2 时", r is not None and not r.get("ok")
              and "opencv" in r.get("reason", ""))
    finally:
        C.HAS_CV2 = saved

    # ---- B 类：已知 GO 对的判据区间（用 spike 缓存的图像）----
    print("配准判据（Alfalfa 已知 GO 对）")
    sp_img = os.path.join(os.path.dirname(ROOT), "06_算法实验与依据", "out", "img")
    sp_xml = os.path.join(os.path.dirname(ROOT), "06_算法实验与依据", "out", "xml")
    wat_xml_p = os.path.join(
        sp_xml, "SIF_0367_0699544277_375FDR_N0110108SRLC08029_0000LMJ04.xml")
    aci_xml_p = os.path.join(
        sp_xml, "SC3_0370_0699816293_742FDR_N0110108SRLC10600_0000LMJ04.xml")
    wat_img_p = os.path.join(
        sp_img, "SIF_0367_0699544277_375FDR_N0110108SRLC08029_0000LMJ04.IMG")
    aci_img_p = os.path.join(
        sp_img, "SC3_0370_0699816293_742FDR_N0110108SRLC10600_0000LMJ04.IMG")
    if not all(os.path.exists(p) for p in
               (wat_xml_p, aci_xml_p, wat_img_p, aci_img_p)):
        print("  SKIP  spike 素材不全（B 类跳过）")
    else:
        import numpy as np
        from make_mineral_map import _cam_from_xml_text, _read_rgb
        from sherloc.vicar import load_aci

        wat_cam = _cam_from_xml_text(open(wat_xml_p, encoding="utf-8",
                                          errors="replace").read())
        aci_cam = _cam_from_xml_text(open(aci_xml_p, encoding="utf-8",
                                          errors="replace").read())
        wat_rgb, _lab = _read_rgb(wat_img_p, wat_cam)
        aci_gray = np.asarray(load_aci(aci_img_p).convert("L"))
        res = C.colorize(aci_gray, wat_rgb, aci_cam, wat_cam, 140.0)
        check("colorize GO", res.get("ok") is True,
              "" if res.get("ok") else res.get("reason"))
        if res.get("ok"):
            check("锚点数 ≥4", res["n_anchors"] >= 4,
                  "(%d)" % res["n_anchors"])
            check("残差 <30px", res["resid"] < 30, "(%.0f)" % res["resid"])
            check("旋转 ≤1°", abs(res["rotation"]) <= 1.0,
                  "(%.2f°)" % res["rotation"])
            check("尺度比 0.2–0.5", 0.2 <= res["scale"] <= 0.5,
                  "(%.2f)" % res["scale"])
            check("受光 ≥50%", res["lit"] >= 0.5, "(%.0f%%)" % (100 * res["lit"]))
            check("峰噪 ≥0.3", res["peak_noise"] >= 0.3,
                  "(%.2f)" % res["peak_noise"])
            check("上色图尺寸同 ACI",
                  res["color"].shape[:2] == aci_gray.shape)

    # ---- A 类：发现路径（需要网络）----
    print("配对发现路径（Alfalfa）")
    inv = watson.load_inventory(DATA)
    if not inv:
        print("  SKIP  拿不到清单（网络不可用？）")
    else:
        check("清单非空", len(inv) > 1000, "(%d stems)" % len(inv))
        cands = watson.candidates(ALFALFA_ACI, inv)
        check("同定位候选 >0", len(cands) > 0, "(%d)" % len(cands))
        pair = watson.find_pair(ALFALFA_ACI, folder=DATA)
        if pair is None:
            check("找到共视配对", False, "（find_pair 返回 None）")
        else:
            check("找到共视配对", True, pair["stem"][-24:])
            check("standoff 在物理窗口",
                  60 <= pair.get("standoff", 0) <= 400,
                  "(%.0f mm)" % pair.get("standoff", 0))

    print("=" * 70)
    print("共 %d 项断言，%d 项失败" % (checks, fails))
    return 0 if fails == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
