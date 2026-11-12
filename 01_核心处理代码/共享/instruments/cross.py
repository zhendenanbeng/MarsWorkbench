# -*- coding: utf-8 -*-
"""跨仪器对照（好奇号 ChemCam LIBS × APXS）。

匹配判据（T2 闸门实证：同址的两仪器目标名不同——Newland vs Wallace_DRT
都在 site48/drive1570——所以 site/drive 为主、目标名为辅、sol 相邻兜底）：

  1. site+drive 相同 = 同一停车位置（最强判据）；
  2. 规范化目标名一致（normalize_target：大小写/下划线/前后缀归一）；
  3. sol 相邻（|Δsol| ≤ 3）且都有目标名（最弱，仅列表提示）。

并排口径（诚实边界，写进界面）：
  - LIBS 侧是**检出证据**（逐发重复率，定性口径），APXS 侧是 **wt% 定量**；
  - 元素↔氧化物映射表只用于「同一元素两来源并排显示」，**不合成新数值**，
    不做 LIBS 定量；
  - SuperCam（毅力号）没有 APXS，它的对照即谱图下方的官方成分预测面板（MOC），
    本模块明确返回不支持而不是硬凑。
"""

from __future__ import annotations

import os
from collections import defaultdict

import numpy as np

from . import pds3

# APXS 氧化物 → 元素（并排显示映射；FeOT 为全铁口径）
OXIDE_TO_ELEMENT = {
    "FeO": "Fe", "FeOT": "Fe", "Fe2O3": "Fe", "FeO-T": "Fe",
    "SiO2": "Si", "TiO2": "Ti", "Al2O3": "Al", "MgO": "Mg",
    "CaO": "Ca", "Na2O": "Na", "K2O": "K", "SO3": "S", "P2O5": "P",
    "Cl": "Cl", "Cr2O3": "Cr", "MnO": "Mn", "Ni": "Ni", "NiO": "Ni",
    "Zn": "Zn", "ZnO": "Zn", "Br": "Br",
}

SOL_ADJACENT = 3


def find_matches(chem_obs, apxs_observations):
    """给一个 ChemCam 观测找 APXS 候选。返回按判据强度排序的列表：

    [{"index": i, "basis": "site_drive"|"target"|"sol_adjacent",
      "detail": "site 48 · drive 1570", "obs": …}, …]
    """
    out = []
    c_site, c_drive = getattr(chem_obs, "site", None), getattr(chem_obs, "drive", None)
    c_tgt = pds3.normalize_target(getattr(chem_obs, "target", ""))
    for i, a in enumerate(apxs_observations):
        a_site, a_drive = getattr(a, "site", None), getattr(a, "drive", None)
        a_tgt = pds3.normalize_target(getattr(a, "target", ""))
        if (c_site is not None and a_site is not None
                and (c_site, c_drive) == (a_site, a_drive)):
            out.append({"index": i, "basis": "site_drive",
                        "detail": "site %s · drive %s（同一停车位置）" % (a_site, a_drive),
                        "obs": a})
        elif c_tgt and a_tgt and c_tgt == a_tgt:
            out.append({"index": i, "basis": "target",
                        "detail": "目标名一致：%s" % getattr(a, "target", ""), "obs": a})
        elif (c_tgt and a_tgt and abs(int(chem_obs.sol) - int(a.sol)) <= SOL_ADJACENT):
            out.append({"index": i, "basis": "sol_adjacent",
                        "detail": "sol 相邻（%s vs %s），目标名不同" % (chem_obs.sol, a.sol),
                        "obs": a})
    rank = {"site_drive": 0, "target": 1, "sol_adjacent": 2}
    out.sort(key=lambda m: rank[m["basis"]])
    return out


def payload(chemcam_backend, chem_obs, apxs_backend, apxs_obs, max_points=9):
    """构建对照负载数据。

    LIBS 侧：观测内逐点（封顶 max_points）取默认预处理的平均谱与逐发证据，
    元素聚合为「出现点数 / 加载点数 · 平均逐发检出率」。
    APXS 侧：官方 RWP 氧化物 wt%（映射到元素并排）。
    """
    from . import evidence as ievidence

    elem_points = defaultdict(int)
    elem_frac = defaultdict(list)
    loaded = 0
    default_steps = {k: on for k, _l, on in chemcam_backend.steps}
    for pt in chem_obs.points[:max_points]:
        try:
            d = chemcam_backend.load_point(pt, None, dict(default_steps))
            marks, elements = chemcam_backend.annotate_data(d)
        except Exception:                                  # noqa: BLE001
            continue
        loaded += 1
        stats = ievidence.element_statistics(marks, elements)
        for s in stats:
            elem_points[s["element"]] += 1 if s["detected"] else 0
            if s["fraction"] is not None:
                elem_frac[s["element"]].append(s["fraction"])

    ox = apxs_backend.oxides(apxs_obs) or {}
    apxs_by_elem = {}
    for r in ox.get("rows", []):
        el = OXIDE_TO_ELEMENT.get(r["oxide"].strip())
        if el is not None:
            apxs_by_elem[el] = {"oxide": r["oxide"], "wt": r["wt"], "err": r["err"]}

    elements = sorted(set(elem_points) | set(apxs_by_elem))
    rows = []
    for el in elements:
        a = apxs_by_elem.get(el)
        n_pts = elem_points.get(el, 0)
        fr = elem_frac.get(el, [])
        rows.append({
            "element": el,
            "libs": ({"points": "%d/%d" % (n_pts, loaded),
                      "mean_fraction": round(float(np.mean(fr)), 2)} if n_pts else None),
            "apxs": a,                                     # None = APXS 未测/未映射
        })
    return {
        "chemcam": {"sol": chem_obs.sol, "seq": chem_obs.seq,
                    "target": getattr(chem_obs, "target", ""),
                    "site": c_sd(chem_obs)[0], "drive": c_sd(chem_obs)[1],
                    "points_loaded": loaded},
        "apxs": {"sol": apxs_obs.sol, "seq": apxs_obs.seq,
                 "target": getattr(apxs_obs, "target", ""),
                 "site": c_sd(apxs_obs)[0], "drive": c_sd(apxs_obs)[1],
                 "file": ox.get("file"), "total": ox.get("total")},
        "rows": rows,
        "note": "LIBS 侧为逐发检出证据（定性+重复率），APXS 侧为官方 wt% 定量；"
                "元素↔氧化物映射仅用于并排显示，不合成数值、不做 LIBS 定量。",
    }


def c_sd(obs):
    return getattr(obs, "site", None), getattr(obs, "drive", None)


def apxs_source_folder(settings):
    """APXS 观测从哪读：记忆的 APXS 文件夹（没有则 None）。"""
    folder = (settings.get("folders") or {}).get("apxs")
    return folder if folder and os.path.isdir(folder) else None


# ---------------------------------------------------------------------------
# SuperCam（毅力号 LIBS）× PIXL（毅力号 XRF）
# ---------------------------------------------------------------------------

def _libs_side(backend, obs, max_points=9):
    """LIBS 侧逐点证据聚合（ChemCam/SuperCam 通用；失败点跳过）。"""
    from . import evidence as ievidence

    elem_points = defaultdict(int)
    elem_frac = defaultdict(list)
    loaded = 0
    default_steps = {k: on for k, _l, on in backend.steps}
    for pt in obs.points[:max_points]:
        try:
            d = backend.load_point(pt, None, dict(default_steps))
            marks, elements = backend.annotate_data(d)
        except Exception:                                  # noqa: BLE001
            continue
        loaded += 1
        stats = ievidence.element_statistics(marks, elements)
        for s in stats:
            elem_points[s["element"]] += 1 if s["detected"] else 0
            if s["fraction"] is not None:
                elem_frac[s["element"]].append(s["fraction"])
    return elem_points, elem_frac, loaded


def _match_rows(elem_points, elem_frac, loaded, other_by_elem):
    """并排行构造（APXS/PIXL 共用；other_by_elem: {元素: {oxide, wt, err}}）。"""
    elements = sorted(set(elem_points) | set(other_by_elem))
    rows = []
    for el in elements:
        a = other_by_elem.get(el)
        n_pts = elem_points.get(el, 0)
        fr = elem_frac.get(el, [])
        rows.append({
            "element": el,
            "libs": ({"points": "%d/%d" % (n_pts, loaded),
                      "mean_fraction": round(float(np.mean(fr)), 2)} if n_pts else None),
            "apxs": a,                                     # None = 另一侧未测/未映射
        })
    return rows


def pixl_candidates_near(sol, window=None):
    """PIXL rxl 扫描候选（处理 sol 窗口内），带 site/drive。

    返回 [{"info", "site", "drive"}]；网络/解析失败抛异常（调用方给明确错误）。
    清单与产品都走 pixl.py 的全局缓存（~/.sherloc_cache/pixl）。
    """
    import re as _re
    from sherloc import pixl as ipixl

    if window is None:
        window = ipixl.SOL_WINDOW
    out = []
    for lid in ipixl.load_inventory(ipixl.COLL_POS):
        info = ipixl.parse_lid(lid)
        if not info or info["type"] != "rxl":
            continue
        try:
            s = int(info["sol"])
        except ValueError:
            continue
        if not (int(sol) - window <= s <= int(sol) + window):
            continue
        csv_path = ipixl.get_product(info)
        xml = csv_path[:-4] + ".xml"
        text = open(xml, encoding="utf-8", errors="replace").read()
        idx = dict(_re.findall(
            r"index_id>([A-Z]+)</\w+:index_id>\s*<\w+:index_value_number>(\d+)",
            text))
        if "SITE" not in idx or "DRIVE" not in idx:
            continue
        site, drive = int(idx["SITE"]), int(idx["DRIVE"])
        if not (0 <= site < 1000 and 0 <= drive < 100000):   # 帧索引健全性钳制
            continue
        out.append({"info": info, "site": site, "drive": drive})
    out.sort(key=lambda c: abs(int(c["info"]["sol"]) - int(sol)))
    return out


def pixl_oxides(rxl_info):
    """一次 PIXL 扫描的逐点成分 → 各氧化物「扫描中位 wt%/err」+ 点数。

    rqc 是处理 sol 产品（可与采集 sol 差几天）；同一 placement 的多个 rqc
    取行数最多的（其余多为中间版本）。返回 None = 该扫描没有可用 rqc。
    """
    from sherloc import pixl as ipixl

    scans = ipixl.find_scans(rxl_info["pos"], rxl_info["sol"],
                             inventories={ipixl.COLL_POS: [rxl_info["lid"]],
                                          ipixl.COLL_CHEM: ipixl.load_inventory(ipixl.COLL_CHEM)})
    best, best_n = None, -1
    for cinfo in scans[ipixl.COLL_CHEM]:
        try:
            chem = ipixl.read_rqc(ipixl.get_product(cinfo))
        except Exception:                                  # noqa: BLE001
            continue
        if len(chem) > best_n:
            best, best_n = chem, len(chem)
    if not best:
        return None
    rows = {}
    for oxide in ipixl.ELEMENTS:
        vals = [best[p][oxide][0] for p in best if oxide in best[p]
                and np.isfinite(best[p][oxide][0])]
        errs = [best[p][oxide][1] for p in best if oxide in best[p]
                and np.isfinite(best[p][oxide][1])]
        if vals:
            rows[oxide] = {"wt": round(float(np.median(vals)), 2),
                           "err": round(float(np.median(errs)), 2) if errs else 0.0}
    return {"rows": rows, "n_points": len(best), "sol": rxl_info["sol"],
            "pos": rxl_info["pos"]}


def pixl_cross(supercam_backend, scam_obs, max_points=9):
    """SuperCam×PIXL 对照入口（/api/cross?app=supercam 用）。

    匹配判据：site/drive 相同（双端都在 label 的 RMC 里；PIXL label 无目标名，
    「目标名一致」这一档在 PIXL 侧不可用）。返回结构与 ChemCam×APXS 同构，
    side 键用 "pixl"。
    """
    s_site, s_drive = c_sd(scam_obs)
    if s_site is None:
        return {"supported": True,
                "note": "本地 SuperCam 数据缺 xml label，读不出 site/drive，"
                        "无法做同址对照（重新下载数据会带 label）。"}
    cands = pixl_candidates_near(scam_obs.sol)
    matches = []
    for c in cands:
        if (c["site"], c["drive"]) == (s_site, s_drive):
            matches.append(c)
    if not matches:
        near = cands[0] if cands else None
        return {"supported": True, "pixl_count": len(cands),
                "note": "PIXL 清单里 sol %s±15 没有同址扫描（site%s/drive%s）；"
                        "%s" % (scam_obs.sol, s_site, s_drive,
                                ("最近的扫描在 sol %s（site%s/drive%s）"
                                 % (near["info"]["sol"], near["site"], near["drive"])
                                 if near else "该窗口内没有 PIXL 扫描"))}
    c = matches[0]
    ox = pixl_oxides(c["info"])
    if not ox:
        return {"supported": True, "error": "同址扫描（sol %s）没有可用的 rqc 成分表"
                % c["info"]["sol"]}
    other_by_elem = {}
    for oxide, r in ox["rows"].items():
        el = OXIDE_TO_ELEMENT.get(oxide)
        if el is not None:
            other_by_elem[el] = {"oxide": oxide, "wt": r["wt"], "err": r["err"]}
    elem_points, elem_frac, loaded = _libs_side(supercam_backend, scam_obs, max_points)
    rows = _match_rows(elem_points, elem_frac, loaded, other_by_elem)
    return {
        "supported": True, "pixl_count": len(cands),
        "matches": [{"basis": "site_drive",
                     "detail": "site %s · drive %s（同一停车位置）" % (s_site, s_drive)}],
        "payload": {
            "supercam": {"sol": scam_obs.sol, "seq": scam_obs.seq,
                         "target": getattr(scam_obs, "target", ""),
                         "site": s_site, "drive": s_drive,
                         "points_loaded": loaded},
            "pixl": {"sol": ox["sol"], "pos": ox["pos"],
                     "n_points": ox["n_points"]},
            "rows": rows,
            "note": "SuperCam 侧为逐发检出证据（定性+重复率），PIXL 侧为逐点 XRF "
                    "（~120µm 斑）的扫描中位 wt%；rqc 是处理 sol 产品，采集 sol 可能"
                    "不同；并排显示不合成数值、不做 LIBS 定量。",
        },
    }
