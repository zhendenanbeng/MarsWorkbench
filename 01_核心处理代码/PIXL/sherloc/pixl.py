# -*- coding: utf-8 -*-
"""PIXL 元素化学：发现 → 下载 → 配准 → 逐点化学查询。

PIXL（X 射线荧光，~120 µm 斑）与 SHERLOC 同臂，扫同一研磨斑。本模块把
"ACI 像素坐标"接到"PIXL 的 16 种氧化物 wt%+误差"上，作为矿物判定的
第二证据源（拉曼说硫酸盐？化学说 SO3 是斑中位的几倍？）。

数据链路（06_算法实验与依据/pixl/README.md 已验证，Bellegarde 重投影 2.66px≈18µm）：

  SHERLOC RLS(ACI 像素) → ACI↔MCC 单应(SIFT) → MCC 像素 → rxl(PMC) → rqc 化学

产品要点（坑都是 spike 实测的）：
  * `rxl` 坐标表在 data_processed，`rqc` 成分表在 data_oxides_pmc；
    成分表按**处理 sol** 分目录（Bellegarde 采集 186、成分表在 187），
    所以候选按 sol ± 窗口搜，不能拿采集 sol 硬查。
  * `rxl` 与 `rqc` 的序列号可以不同（699274310005 vs …000），同一扫描
    被拆成多个子序列——**用 PMC 做连接键**，选 PMC 交集最大的配对。
  * MCC 图像本体在 mars2020_imgops:data_mcc_imgops（与 ACI/WATSON 同
    一个 bundle）；PIXL 引用名叫 `...rcm...`，归档叫 `...fdr/ecm...`。
  * Wustl 偶发把错误页 HTML 当下载体返回——下载后验头字节。

用法（给质控台/分布图用的最小接口）：

    chem = pixl.chemistry_at_points(
        aci_gray, sol="0186", placement="N0070000",
        points=[(samp, line), ...], folder=数据目录, log=print)
    # chem["median"]  各氧化物斑中位
    # chem["points"]  与 points 等长：None 或 {"pmc", "dist_mcc_px", "elements"}
"""

from __future__ import annotations

import csv
import json
import os
import re
import urllib.parse
import urllib.request

import numpy as np

try:
    import cv2
    HAS_CV2 = True
except ImportError:                                  # noqa: F401
    HAS_CV2 = False

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(os.path.expanduser("~"), ".sherloc_cache", "pixl")
BASE = "https://pds-geosciences.wustl.edu/m2020/urn-nasa-pds-mars2020_pixl"
IMGOPS_API = "https://pds.nasa.gov/api/search/1/products"
UA = {"User-Agent": "sherloc-workbench/1.0"}

COLL_POS = "data_processed"          # rxl 坐标表
COLL_CHEM = "data_oxides_pmc"        # rqc 成分表
SOL_WINDOW = 15                      # 处理 sol 相对采集 sol 的偏移上界
ELEMENTS = ("Na2O", "MgO", "Al2O3", "SiO2", "P2O5", "SO3", "Cl", "K2O",
            "CaO", "TiO2", "Cr2O3", "MnO", "FeO-T", "NiO", "ZnO", "Br")
MAX_MCC_DIST = 25.0                  # 最近 PMC 的门限（MCC px，1px≈32-49µm）
MIN_INLIERS = 15                     # MCC↔ACI 配准的 RANSAC 内点下限


# ---------------------------------------------------------------------------
# 清单与产品名
# ---------------------------------------------------------------------------

def _ensure_cache():
    os.makedirs(CACHE, exist_ok=True)


def load_inventory(coll):
    """collection 的 LID 清单（缓存 CSV，Wustl 直链稳定）。"""
    _ensure_cache()
    p = os.path.join(CACHE, "inv_%s.csv" % coll)
    if not os.path.exists(p):
        url = "%s/%s/collection_%s_inventory.csv" % (BASE, coll, coll)
        req = urllib.request.Request(url, headers=UA)
        with urllib.request.urlopen(req, timeout=180) as r, \
                open(p + ".part", "wb") as f:
            f.write(r.read())
        os.replace(p + ".part", p)
    out = []
    with open(p, encoding="utf-8", errors="replace") as f:
        for line in f:
            lid = line.strip()
            if "," in lid:
                lid = lid.split(",", 1)[1]
            if "::" in lid:
                lid = lid.split("::")[0]
            if lid:
                out.append(lid)
    return out


def parse_lid(lid):
    """拆 PIXL 产品 LID。

    urn:nasa:pds:mars2020_pixl:data_processed:pe__0186_0683483744_000rxl__00700000699274310005___j07
      → {coll, pre: pe, sol: 0186, sclk: 0683483744, type: rxl,
         pos: 0070000, seq: 699274310005, tail: j07}
    """
    stem = lid.split(":")[-1]
    parts = stem.split("__")
    if len(parts) < 3:
        return None
    pre = parts[0]
    m = re.match(r"(\d{4})_(\d{9,10})_0*(\w+)", parts[1])
    posseq = parts[2]
    if not m or len(posseq) < 20:
        return None
    return {"coll": lid.split(":")[-2], "pre": pre, "lid": lid,
            "sol": m.group(1), "sclk": m.group(2), "type": m.group(3),
            "pos": posseq[:7], "seq": posseq[7:], "tail": parts[-1],
            "stem": stem}


def find_scans(placement, sol, window=SOL_WINDOW, inventories=None):
    """按机械臂定位 + sol 窗口找 rxl / rqc 产品 LID。

    placement 传 SHERLOC 侧的 8 位（如 "N0070000"）或 7 位数字。
    """
    want = placement.strip().upper().lstrip("N")
    inv = inventories or {COLL_POS: load_inventory(COLL_POS),
                          COLL_CHEM: load_inventory(COLL_CHEM)}
    lo = int(sol) - window
    hi = int(sol) + window
    want_type = {COLL_POS: "rxl", COLL_CHEM: "rqc"}
    out = {COLL_POS: [], COLL_CHEM: []}
    for coll, lids in inv.items():
        for lid in lids:
            info = parse_lid(lid)
            if not info or info["pos"] != want:
                continue
            if info["type"] != want_type.get(coll):
                continue
            try:
                s = int(info["sol"])
            except ValueError:
                continue
            if lo <= s <= hi:
                out[coll].append(info)
    for coll in out:
        out[coll].sort(key=lambda i: (i["sol"], i["seq"]))
    return out


def _download(url, path):
    _ensure_cache()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".part"
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=300) as r, \
                    open(tmp, "wb") as f:
                f.write(r.read())
            head = open(tmp, "rb").read(64)
            if b"<html" in head.lower() or b"<!DOCTYPE" in head[:64].upper():
                raise ValueError("wustl 返回了错误页")
            os.replace(tmp, path)
            return path
        except Exception:                            # noqa: BLE001
            if attempt == 2:
                if os.path.exists(tmp):
                    os.remove(tmp)
                raise
    return None


_REFS = {}                                   # lid → {文件名: 直链}


def resolve_lid(lid):
    """PDS API：产品 LID → {文件名: 直链}。

    不能自己拼 wustl 路径：**产品名里的 sol ≠ 存放目录的 sol**（处理
    可以晚几天，rxl/rqc 都有采 186 存 187 的例子），必须让 API 给路径。
    一次解析顺带拿到该产品的全部附属文件（rxl 的 4 张 MCC tif 就在
    它自己的 file_ref 列表里）。
    """
    if lid in _REFS:
        return _REFS[lid]
    q = urllib.parse.quote('lid eq "%s"' % lid, safe="")
    req = urllib.request.Request(IMGOPS_API + "?limit=1&q=" + q, headers=UA)
    with urllib.request.urlopen(req, timeout=120) as r:
        d = json.loads(r.read().decode("utf-8", "replace"))
    refs = {}
    if (d.get("summary") or {}).get("hits"):
        p = d["data"][0]["properties"]
        names = p.get("ops:Data_File_Info.ops:file_name") or []
        urls = p.get("ops:Data_File_Info.ops:file_ref") or []
        for n, u in zip(names, urls):
            refs[n.lower()] = u
    _REFS[lid] = refs
    return refs


def get_product(info):
    """下载（缓存）一个 rxl/rqc 产品的 csv 和 xml，返回 csv 本地路径。

    xml 必须要：MCC 上下文图的文件名引用在 rxl 的 xml 里，csv 没有。
    """
    name = info["stem"] + ".csv"
    p = os.path.join(CACHE, info["coll"], name)
    if os.path.exists(p):
        if not os.path.exists(p[:-4] + ".xml"):
            _ensure_xml(info, p)
        return p
    refs = resolve_lid(info["lid"])
    url = refs.get(name.lower())
    if not url:
        raise ValueError("API 没有这个产品的直链：%s" % name)
    _download(url, p)
    _ensure_xml(info, p)
    return p


def _ensure_xml(info, csv_path):
    """补下 label xml。API 的 file_ref 不含 xml，但从 csv 直链换扩展名即可
    （同目录同名）。失败不致命——MCC tif 名单另有来源。"""
    url = resolve_lid(info["lid"]).get((info["stem"] + ".csv").lower())
    if not url:
        return
    try:
        _download(url[:-4] + ".xml", csv_path[:-4] + ".xml")
    except Exception:                                # noqa: BLE001
        pass


# ---------------------------------------------------------------------------
# rxl / rqc 解析
# ---------------------------------------------------------------------------

def read_rxl(path):
    """解析 rxl：返回 {"images": [图像键], "pts": {pmc: {imgkey: (i, j)}},
    "mcc_tifs": [PIXL 引用的 tif 文件名]}。"""
    rows = list(csv.reader(open(path, encoding="utf-8", errors="replace")))
    hdr = rows[0]
    imgkeys = sorted({h.rsplit("_", 1)[0] for h in hdr if h.endswith("_i")})
    pts = {}
    for r in rows[1:]:
        if not r or not r[0].strip().lstrip("-").isdigit():
            continue
        pmc = int(r[0])
        d = {}
        for k in imgkeys:
            try:
                d[k] = (float(r[hdr.index(k + "_i")]),
                        float(r[hdr.index(k + "_j")]))
            except (ValueError, IndexError):
                pass
        pts[pmc] = d
    xml = path[:-4] + ".xml"
    tifs = []
    if os.path.exists(xml):
        text = open(xml, encoding="utf-8", errors="replace").read()
        tifs = sorted(set(re.findall(r"pcw_\w+\.tif", text)))
    return {"images": imgkeys, "pts": pts, "mcc_tifs": tifs}


def read_rqc(path):
    """解析 rqc：{pmc: {elem: (wt%, err)}}。"""
    rows = list(csv.reader(open(path, encoding="utf-8", errors="replace")))
    hdr = rows[0]
    out = {}
    for r in rows[1:]:
        if not r or not r[0].strip().lstrip("-").isdigit():
            continue
        pmc = int(r[0])
        d = {}
        for e in ELEMENTS:
            try:
                v = float(r[hdr.index(e + "_wt%")])
                er = float(r[hdr.index(e + "_wt%_err")])
            except (ValueError, IndexError):
                continue
            d[e] = (v, er)
        if d:
            out[pmc] = d
    return out


# ---------------------------------------------------------------------------
# MCC 图像：PIXL 引用名 → 直链（就在 rxl 产品自己的 file_ref 列表里）
# ---------------------------------------------------------------------------

def get_mcc(tif_name, refs=None):
    """下载（缓存）一张 MCC 上下文图，返回本地路径。

    tif 是 rxl 产品的附属文件，直链随 resolve_lid(rxl LID) 一起拿到；
    名字里 PIXL 内部叫 `rcm`、imgops 叫 `fdr`，别做字符串变换，直接用
    API 给的文件名对上。
    """
    p = os.path.join(CACHE, "mcc", tif_name.replace(".tif", ".IMG"))
    if os.path.exists(p):
        return p
    url = (refs or {}).get(tif_name.lower())
    if not url:
        raise ValueError("MCC 图像不在产品附属文件里：%s" % tif_name)
    return _download(url, p)


def mcc_gray(path):
    """MCC 图读成灰度 uint8。PIXL bundle 里的 .tif 是真 TIFF（PIL 可读），
    imgops 的 .IMG 是 VICAR（vicar.py 读）——两种来源都支持。"""
    try:
        from PIL import Image
        im = Image.open(path)
        arr = np.asarray(im)
        if arr.ndim == 3:
            arr = arr[..., 0]
        return _stretch(arr.astype(np.float64))
    except Exception:                                # noqa: BLE001
        from .vicar import read_vicar
        img, _ = read_vicar(path)
        return _stretch(img.astype(np.float64))


def _stretch(g, lo_pct=0.5, hi_pct=99.5):
    g = g.astype(np.float32)
    lo, hi = np.percentile(g, (lo_pct, hi_pct))
    return np.clip((g - lo) / max(hi - lo, 1e-9) * 255, 0, 255).astype(np.uint8)


# ---------------------------------------------------------------------------
# MCC ↔ ACI 配准（spike 验证过的 SIFT 管线）
# ---------------------------------------------------------------------------

def _prep(g):
    g = g.astype(np.float32)
    bg = cv2.GaussianBlur(g, (0, 0), sigmaX=max(g.shape) / 12.0)
    g = cv2.normalize(g - bg, None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
    return cv2.createCLAHE(2.0, (8, 8)).apply(g)


def register_mcc_aci(mcc_gray, aci_gray):
    """MCC → ACI 的仿射变换。返回 (M 2x3, {"inliers","scale","err"}) 或 (None, 原因)。"""
    sift = cv2.SIFT_create(20000, contrastThreshold=0.01)
    flann = cv2.FlannBasedMatcher(dict(algorithm=1, trees=5), dict(checks=64))
    k2, d2 = sift.detectAndCompute(_prep(mcc_gray), None)
    k1, d1 = sift.detectAndCompute(_prep(aci_gray), None)
    if d2 is None or d1 is None:
        return None, "特征不足"
    ms = flann.knnMatch(d2, d1, k=2)
    bt = {}
    for m, n in ms:
        if m.distance < 0.8 * n.distance:
            t = m.trainIdx
            if t not in bt or m.distance < bt[t].distance:
                bt[t] = m
    good = list(bt.values())
    if len(good) < 10:
        return None, "好匹配太少（%d）" % len(good)
    src = np.float64([k2[m.queryIdx].pt for m in good]).reshape(-1, 1, 2)
    dst = np.float64([k1[m.trainIdx].pt for m in good]).reshape(-1, 1, 2)
    M, inl = cv2.estimateAffinePartial2D(src, dst, method=cv2.RANSAC,
                                         ransacReprojThreshold=4.0,
                                         maxIters=20000, confidence=0.999)
    if M is None or inl is None:
        return None, "RANSAC 失败"
    mask = inl.ravel().astype(bool)
    if mask.sum() < MIN_INLIERS:
        return None, "内点不足（%d）" % mask.sum()
    proj = cv2.transform(src, M).reshape(-1, 2)
    err = float(np.median(np.hypot(*(proj[mask] - dst.reshape(-1, 2)[mask]).T)))
    return M, {"inliers": int(mask.sum()), "good": len(good),
               "scale": float(np.hypot(M[0, 0], M[1, 0])), "err_px": err}


# ---------------------------------------------------------------------------
# 主入口：给一次 SHERLOC 扫描建化学查询
# ---------------------------------------------------------------------------

def chemistry_at_points(aci_gray, sol, placement, points, folder=None,
                        log=None):
    """把一批 ACI 像素坐标接到 PIXL 化学。返回 dict 或 None（原因在 log）。

    points: [(samp, line), ...]（RLS 的 ACI 像素坐标）
    返回 {"median": {elem: 中位wt%}, "scale": MCC→ACI 尺度比,
          "points": [(pmc|None, dist_mcc_px|None, {elem: (v, er)}) 与 points 等长],
          "inliers": n, "err_px": float, "pmc_used": n}
    """
    def say(msg):
        if log:
            log(msg)

    if not HAS_CV2:
        say("PIXL：未安装 opencv，跳过化学佐证")
        return None
    try:
        scans = find_scans(placement, sol)
    except Exception as exc:                         # noqa: BLE001
        say("PIXL：清单拿不到（%s）" % exc)
        return None
    if not scans[COLL_POS]:
        say("PIXL：定位 %s 在 sol %s ±%d 没有扫描" % (placement, sol, SOL_WINDOW))
        return None
    say("PIXL：定位 %s 候选 rxl %d 个、rqc %d 个"
        % (placement, len(scans[COLL_POS]), len(scans[COLL_CHEM])))

    aci8 = _stretch(np.asarray(aci_gray))
    best = None
    for info in scans[COLL_POS]:
        try:
            refs = resolve_lid(info["lid"])
            rxl = read_rxl(get_product(info))
            # MCC tif 名单：产品附属文件里就有（xml 里也有，但 API 的
            # file_ref 不含 label，别依赖 xml）
            tifs = sorted(f for f in refs if f.endswith(".tif")) \
                or rxl["mcc_tifs"]
            chems = []
            for cinfo in scans[COLL_CHEM]:
                try:
                    chem = read_rqc(get_product(cinfo))
                    chems.append((len(set(chem) & set(rxl["pts"])), cinfo, chem))
                except Exception:                        # noqa: BLE001
                    continue
            if not chems:
                say("PIXL：%s 没有可用的成分表" % info["stem"][-20:])
                continue
            chems.sort(key=lambda t: -t[0])
            overlap, cinfo, chem = chems[0]

            # 依次尝试每张 MCC 图，取配准最好的
            for tif in tifs:
                try:
                    mcc = mcc_gray(get_mcc(tif, refs))
                except Exception as exc:             # noqa: BLE001
                    say("PIXL：MCC %s 拿不到（%s）" % (tif[-26:], exc))
                    continue
                imgkey = _tif_to_key(tif, rxl)
                if imgkey is None:
                    continue
                M, stats = register_mcc_aci(mcc, aci8)
                if M is None:
                    say("PIXL：MCC %s 配不上（%s）" % (tif[-26:], stats))
                    continue
                say("PIXL：%s ↔ ACI 内点 %d、重投影 %.2fpx（%s）"
                    % (tif[-26:], stats["inliers"], stats["err_px"], imgkey))
                if best is None or stats["inliers"] > best["stats"]["inliers"]:
                    best = {"M": M, "stats": stats, "rxl": rxl,
                            "chem": chem, "overlap": overlap,
                            "imgkey": imgkey, "info": info,
                            "mcc_gray": mcc}
        except Exception as exc:                     # noqa: BLE001
            say("PIXL：%s 处理失败（%s）" % (info["stem"][-20:], exc))
    if best is None:
        say("PIXL：没有任何一张 MCC 图能与这张 ACI 配准")
        return None

    M = best["M"]
    Minv = cv2.invertAffineTransform(M.astype(np.float64))
    rxl, chem = best["rxl"], best["chem"]
    imgkey = best["imgkey"]
    pmc_pos = {pmc: d[imgkey] for pmc, d in rxl["pts"].items() if imgkey in d}
    chem = {p: c for p, c in chem.items() if p in pmc_pos}
    mc = np.array(list(pmc_pos.values()), np.float64)
    pmcs = list(pmc_pos.keys())

    # 斑中位（用作" elevated 相对什么"的基准）
    per_elem = {e: [chem[p][e][0] for p in chem if e in chem[p]]
                for e in ELEMENTS}
    median = {e: float(np.median(v)) for e, v in per_elem.items() if v}

    pts_out = []
    for (sx, sy) in points:
        m = Minv[:, :2] @ np.array([float(sx), float(sy)])             + Minv[:, 2]                                  # ACI → MCC
        dd = np.hypot(mc[:, 0] - m[0], mc[:, 1] - m[1])
        k = int(np.argmin(dd))
        if dd[k] > MAX_MCC_DIST or pmcs[k] not in chem:
            pts_out.append((None, None, None))
            continue
        pts_out.append((pmcs[k], float(dd[k]), chem[pmcs[k]]))
    n_hit = sum(1 for p in pts_out if p[0] is not None)
    say("PIXL：化学查询 %d/%d 点命中（PMC 交集 %d，斑中位基准 %d 点）"
        % (n_hit, len(points), best["overlap"], len(chem)))
    return {"median": median, "points": pts_out,
            "inliers": best["stats"]["inliers"],
            "err_px": best["stats"]["err_px"],
            "scale": best["stats"]["scale"], "pmc_used": len(chem),
            # 重投影上下文：同一张 ACI 换一批点时，复用配准与成分表只重算
            # 最近邻查询（run_qc 的 memo 跨观测复用就靠它，不然第二份观测
            # 会按位置套用第一份观测的逐点结果）
            "_ctx": {"Minv": Minv, "mc": mc, "pmcs": pmcs, "chem": chem,
                     "src_points": [tuple(map(float, p)) for p in points]},
            # 热图用的原始素材：MCC 灰度 + 每个 PMC 的 MCC 像素坐标与成分
            "mcc_gray": best.get("mcc_gray"),
            "pmc_xy": {p: pmc_pos[p] for p in chem},
            "chem": chem}


def reproject_points(res, points):
    """对同一张 ACI 的另一批点重算逐点化学。

    复用 chemistry_at_points 存下的配准（Minv）与成分表，只重做最近邻查询。
    返回与 points 等长的 [(pmc|None, dist|None, elems|None), ...]；
    res 没有 _ctx（旧结构）时返回 None，调用方自行退回 res["points"]。
    """
    ctx = (res or {}).get("_ctx")
    if ctx is None:
        return None
    pts = [tuple(map(float, p)) for p in points]
    if pts == ctx["src_points"]:
        return res["points"]
    Minv, mc, pmcs, chem = (ctx["Minv"], ctx["mc"],
                            ctx["pmcs"], ctx["chem"])
    out = []
    for (sx, sy) in pts:
        m = Minv[:, :2] @ np.array([sx, sy]) + Minv[:, 2]      # ACI → MCC
        dd = np.hypot(mc[:, 0] - m[0], mc[:, 1] - m[1])
        k = int(np.argmin(dd))
        if dd[k] > MAX_MCC_DIST or pmcs[k] not in chem:
            out.append((None, None, None))
            continue
        out.append((pmcs[k], float(dd[k]), chem[pmcs[k]]))
    return out


def _tif_to_key(tif_name, rxl):
    """MCC tif 文件名 → rxl 里的列键（"PMC_0005_MCC"）。

    名字规则：定位后 13 位序列号的**末 4 位 = 中心 PMC**
    （…699274310005075j → 0005；…699274312816075j → 2816）。
    再用"该 PMC 自己的坐标应落在图心附近"做校验；对不上时退回
    枚举所有键挑图心最近的。
    """
    m = re.search(r"n\d{7}(\d{13})\d{3}j", tif_name.lower())
    if m:
        key = "PMC_%04d_MCC" % int(m.group(1)[-4:])
        if key in rxl["images"]:
            return key
    # 兜底：每个键自己的 PMC 应在图心 (ns/2, nl/2) 附近
    best, best_d = None, 1e9
    for k in rxl["images"]:
        pmc = int(k.split("_")[1])
        pos = rxl["pts"].get(pmc, {}).get(k)
        if pos is None:
            continue
        d = np.hypot(pos[0] - 752 / 2.0, pos[1] - 580 / 2.0)
        if d < best_d:
            best, best_d = k, d
    return best


# ---------------------------------------------------------------------------
# 化学佐证（质控台逐点表用的判定规则）
# ---------------------------------------------------------------------------

# 类别 key → (元素, 说明)。carbonate 特殊处理（无 C，看 Ca/Mg/Fe 高 + Si 不高）。
_CHEM_RULES = {
    "sulfate": ("SO3", "硫酸盐"),
    "sulfate_or_pyroxene": ("SO3", "硫酸盐?"),
    "perchlorate_chlorate": ("Cl", "高氯酸盐"),
    "sodium_perchlorate": ("Cl", "钠高氯酸盐"),
    "perchlorate_or_phosphate": ("P2O5", "磷酸盐?"),
    "phosphate": ("P2O5", "磷酸盐"),
}


def evidence(class_key, elem_dict, median):
    """对一个点的某个拉曼判定给化学佐证。返回 "SO3 3.2%(斑中位1.1,×2.9)✓" 形式。

    标记语义（2026-09-20 起用显著性，不再只看倍数）：
      ✓ = 显著富集：≥1.5× 且 (v−斑中位)/误差 ≥ 3（rqc 每个元素自带误差列）
      ✗ = 低于斑中位（ratio < 1.0）——硫酸盐/高氯酸盐这类"该元素应当富集"
          的判定，点值贴着斑中位就该质疑
      ～ = 其余（富集但未达 3σ，或略高）

    elem_dict/median 任一为空、或类别没有对应元素规则时返回 ""。
    """
    if not elem_dict or not median:
        return ""
    if class_key == "carbonate":
        return _carbonate_evidence(elem_dict, median)
    if class_key == "perchlorate_or_phosphate":
        # 双元素分支：Cl 和 P2O5 各自做显著性检验，谁显著富集且更强判谁；
        # 都不显著 → 这个"二选一"判定失去化学支撑
        best = None
        for e in ("Cl", "P2O5"):
            if e in elem_dict and median.get(e, 0) > 0:
                v, er = elem_dict[e]
                med = median[e]
                ratio, z = v / med, (v - med) / max(abs(er), 1e-9)
                if ratio >= 1.5 and z >= 3.0 and (best is None or ratio > best[2]):
                    best = (e, v, ratio, z, med)
        if best:
            e, v, ratio, z, med = best
            return "%s %.2f%%(斑中位%.2f,×%.1f,z%+.0f)✓" % (e, v, med, ratio, z)
        return "P/Cl 均不显著富集 ✗"
    rule = _CHEM_RULES.get(class_key)
    if not rule:
        return ""
    elem, _label = rule
    if elem not in elem_dict or elem not in median:
        return ""
    v, er = elem_dict[elem]
    med = median[elem]
    if med <= 0:
        # 斑里几乎没有该元素，但这个点报出了对应矿物：查绝对值
        return "%s %.2f%%(斑中位~0)%s" % (elem, v, "✓" if v >= 1.0 else "？")
    ratio = v / med
    z = (v - med) / max(abs(er), 1e-9)
    if ratio >= 1.5 and z >= 3.0:
        mark = "✓"
    elif ratio < 1.0:
        mark = "✗"
    else:
        mark = "～"
    return "%s %.2f%%(斑中位%.2f,×%.1f,z%+.0f)%s" % (elem, v, med, ratio, z, mark)


def _carbonate_evidence(elem_dict, median):
    """碳酸盐佐证：(Ca,Mg,Fe) 相对斑中位最富集者 + SiO2 抬升否决。"""
    cand = []
    for e in ("CaO", "MgO", "FeO-T"):
        if e in elem_dict and median.get(e, 0) > 0:
            v, er = elem_dict[e]
            cand.append((e, v / median[e], (v - median[e]) / max(abs(er), 1e-9),
                         v, median[e]))
    if not cand:
        return ""
    base, ratio, z, v, med = max(cand, key=lambda t: t[1])
    si = elem_dict.get("SiO2", (0, 0))[0]
    si_med = median.get("SiO2", 0.0)
    si_ok = si_med <= 0 or si <= si_med * 1.5
    if ratio >= 1.5 and z >= 3.0 and si_ok:
        mark = "✓"
    elif ratio < 1.2 and si_med > 0 and si >= si_med:
        mark = "✗"
    else:
        mark = "～"
    return "%s %.1f%%(斑中位%.1f,×%.1f)%s" % (base, v, med, ratio, mark)


def arbitrate(class_keys, elems, median):
    """化学参与判定：对"二选一"模糊类别做拆分仲裁。

    只动两类（它们本来就是拉曼缺次级峰逼出来的二选一）。
    判据用显著性（rqc 自带误差列）：z = (v − 斑中位)/误差。
      * sulfate_or_pyroxene：SO3 低于斑中位（<1.0×）→ 判辉石（否决硫酸盐）；
        ≥1.5× 且 z≥3 → 确认硫酸盐；中间不动。
      * perchlorate_or_phosphate：Cl 显著富集（≥1.5× 且 z≥3）且强于 P →
        高氯酸盐；P2O5 显著富集且强于 Cl → 磷酸盐；否则维持二选一。
    其余类别不仲裁（✓/✗ 只作注释，替用户删矿物是更重的科学决定）。

    返回 (被仲裁的键, 新键列表, note)；无需仲裁返回 (None, None, "")。
    注意：new_keys 只替换被仲裁的那一个键，**类别列表里的其他类别保留**——
    首版把整个列表替换掉，顺手丢掉了同点的碳酸盐（实测踩过）。
    """
    if not elems or not median:
        return None, None, ""

    def rp(e):
        """(倍数, 显著性 z)；元素缺失或斑中位为 0 时 (None, None)。"""
        if e not in elems or median.get(e, 0) <= 0:
            return None, None
        v, er = elems[e]
        med = median[e]
        return v / med, (v - med) / max(abs(er), 1e-9)

    if "sulfate_or_pyroxene" in class_keys:
        r, z = rp("SO3")
        if r is not None:
            if r < 1.0:
                return ("sulfate_or_pyroxene", ["pyroxene"],
                        "化学否决硫酸盐（SO3 ×%.1f）→ 辉石" % r)
            if r >= 1.5 and z >= 3.0:
                return ("sulfate_or_pyroxene", ["sulfate"],
                        "化学确认硫酸盐（SO3 ×%.1f，z%.0f）" % (r, z))
    if "perchlorate_or_phosphate" in class_keys:
        cr, cz = rp("Cl")
        pr, pz = rp("P2O5")
        if cr is not None and pr is not None:
            if cr >= 1.5 and cz >= 3.0 and cr > pr:
                # 键名必须是 minerals 里真实存在的类别（perchlorate_chlorate），
                # 否则点阵/图例/配色全部对不上——首版就写错过一个不存在的键
                return ("perchlorate_or_phosphate", ["perchlorate_chlorate"],
                        "化学判高氯酸盐（Cl ×%.1f，z%.0f）" % (cr, cz))
            if pr >= 1.5 and pz >= 3.0 and pr > cr:
                return ("perchlorate_or_phosphate", ["phosphate"],
                        "化学判磷酸盐（P2O5 ×%.1f，z%.0f）" % (pr, pz))
    return None, None, ""


def save_heatmap(chem_res, element, out_path, title=None):
    """把某元素的 PIXL 浓度场画成 MCC 底图上的热图（PIL+numpy，无新依赖）。

    chem_res 是 chemistry_at_points 的返回值（要有 mcc_gray/pmc_xy/chem）。
    网格化用最近邻溅落 + 归一化卷积（WATSON 空洞填充同款思路），
    色带用 cv2 TURBO。返回输出路径；素材缺失返回 None。
    """
    if not HAS_CV2 or not chem_res             or chem_res.get("mcc_gray") is None:
        return None
    mcc = chem_res["mcc_gray"]
    h, w = mcc.shape[:2]
    xy = chem_res["pmc_xy"]
    chem = chem_res["chem"]

    acc = np.zeros((h, w), np.float64)
    cnt = np.zeros((h, w), np.float64)
    n = 0
    for pmc, (i, j) in xy.items():
        d = chem.get(pmc)
        if not d or element not in d:
            continue
        x0, y0 = int(round(i)), int(round(j))
        if not (0 <= x0 < w and 0 <= y0 < h):
            continue
        acc[y0, x0] += d[element][0]
        cnt[y0, x0] += 1.0
        n += 1
    if cnt.sum() == 0:
        return None
    # 归一化卷积：值和掩模各糊一次再相除，空洞不沾染色
    fill = cv2.GaussianBlur(acc, (0, 0), 3.5)
    wgt = cv2.GaussianBlur(cnt, (0, 0), 3.5)
    valid = wgt > 1e-6
    field = np.where(valid, fill / np.maximum(wgt, 1e-9), 0)

    lo, hi = np.percentile(field[valid], (2, 98))
    norm = np.clip((field - lo) / max(hi - lo, 1e-9) * 255.0, 0, 255)
    colored = cv2.applyColorMap(norm.astype(np.uint8), cv2.COLORMAP_TURBO)
    colored[~valid] = 0

    base = cv2.cvtColor(mcc, cv2.COLOR_GRAY2BGR)
    base = np.clip(base.astype(np.float64) * 1.15, 0, 255)
    alpha = np.where(valid, 0.58, 0.0)[..., None]
    out = (base * (1 - alpha) + colored.astype(np.float64) * alpha).astype(np.uint8)

    # 右侧色带 + 刻度（PIL 画字，中文用系统微软雅黑）
    from PIL import Image, ImageDraw
    BAR = 112
    img = Image.fromarray(cv2.cvtColor(out, cv2.COLOR_BGR2RGB))
    canvas = Image.new("RGB", (w + BAR, h), (250, 250, 248))
    canvas.paste(img, (0, 0))
    draw = ImageDraw.Draw(canvas)
    top, bot, x0b = 60, h - 60, w + 12
    for yy in range(top, bot):
        t = 1.0 - (yy - top) / max(bot - top - 1, 1)
        r_, g_, b_ = _turbo(t)
        draw.line([(x0b, yy), (x0b + 34, yy)], fill=(r_, g_, b_))
    try:
        from make_mineral_map import get_font
        f_big, f_small = get_font(22), get_font(15)
    except Exception:                                # noqa: BLE001
        from PIL import ImageFont
        f_big = f_small = ImageFont.load_default()
    draw.text((x0b, 24), "%s wt%%" % element, fill=(20, 20, 20), font=f_big)
    # TURBO 红端在色带顶部（t=1 对应数据 hi），所以 hi 标顶、lo 标底
    draw.text((x0b, top - 20), "%.1f" % hi, fill=(40, 40, 40), font=f_small)
    draw.text((x0b, bot + 4), "%.1f" % lo, fill=(40, 40, 40), font=f_small)
    if title:
        draw.text((10, 8), title, fill=(255, 255, 255), font=f_big)
    canvas.save(out_path)
    return out_path


def _turbo(t):
    """TURBO 色带的近似锚点插值（避免为了一个色带引入 matplotlib）。"""
    anchors = [(48, 18, 59), (70, 107, 227), (40, 187, 219), (27, 229, 131),
               (110, 247, 38), (215, 197, 8), (251, 113, 5), (222, 38, 8),
               (122, 4, 3)]
    t = min(max(t, 0.0), 1.0) * (len(anchors) - 1)
    k = min(int(t), len(anchors) - 2)
    f = t - k
    a, b = anchors[k], anchors[k + 1]
    return tuple(int(a[i] + (b[i] - a[i]) * f) for i in range(3))
