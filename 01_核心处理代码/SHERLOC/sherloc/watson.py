# -*- coding: utf-8 -*-
"""给 ACI 图找配套的 WATSON 彩色图（查找 + 下载）。

"配套"怎么判：不能按时间猜（sol 141 实测 WATSON 比 ACI 早 12 分钟，
看的是两块不同的地方），也不能按 sol 找（跨 sol 的 Alfalfa 配对照样能
配准成功）。唯一的硬判据是 **几何共视**——两台相机的 CAHVOR 模型
（各自的 .xml 标签里）都在 ROVER_NAV_FRAME 坐标系下，把 ACI 的视场
投到 WATSON 上看落不落得进去。

流程：
  1. `candidates()`  从全量清单里挑出同定位（同 sol + 同 `_Nxxxxxxx`
     姿态段）的 WATSON 彩色产品（SIF_ 前缀 = 3 波段）；
  2. `cosees()`      下载候选的 .xml（~80 KB，比 2 MB 的图便宜得多），
     用 CAHVOR 扫 standoff，看 ACI 四角能否投进 WATSON 画幅；
  3. `ensure()`      通过共视判定的候选下载 .IMG 到缓存目录。

清单 CSV 在 PDS 归档里有全量版（一次下载长期可用）；没有清单时
退化到 PDS 检索 API 逐个查。两台相机的几何函数与 06_算法实验与依据/refine_any.py
保持一致（那边是验证过的实现）。

两个坑（都是实测踩过的）：
  * 图名末尾两位（如 LMJ04 的 04）是产品版本号，不属于 LID；
  * PDS 检索 API 的返回字段有时是列表，要取第一个元素。
"""

from __future__ import annotations

import csv
import io
import json
import os
import re
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))

API = "http://pds.nasa.gov/api/search/1/products"
UA = {"User-Agent": "sherloc-workbench/1.0",
      "Accept": "application/json"}

# 全量清单的直链（PDS 归档，一次下载长期可用）
INV_URL = ("https://pds-imaging.jpl.nasa.gov/archive/m20/cumulative/"
           "mars2020_imgops/")

# 缓存目录：清单 + 相机模型 xml 都放这里
CACHE = os.path.join(os.path.expanduser("~"), ".sherloc_cache", "watson")

# standoff 扫描范围（mm）：覆盖 SHERLOC 的工作距离量级，步长 10
STAND_RANGE = range(60, 401, 10)


# ---------------------------------------------------------------------------
# 产品名解析
# ---------------------------------------------------------------------------

def parse_name(stem):
    """拆图像产品名。返回 dict 或 None。

    SC3_0186_0683479674_054FDR_N0070000SRLC11420_0000LMJ05
     │    │     │       │└─┬┘ │        │          │└ 版本号
     前缀  sol  SCLK    序号 类型 姿态段    序列       编码

    第 4 段 = 3 位序号 + 3 位产品类型（FDR/RAD/ECM…）。
    姿态段 `NxxxxxxxSRLCxxxxx` 是"机械臂定位 + 序列"的标识。
    """
    s = os.path.splitext(os.path.basename(stem))[0]
    m = re.match(r"^(S[CI][F0-9])_(\d{4})_(\d{10})_\d{3}(\w{3})_"
                 r"([A-Z]\d{7}[A-Z]+\d+)_([A-Z0-9]+?)(\d{2})?$",
                 s.upper())
    if not m:
        return None
    return {"prefix": m.group(1), "sol": m.group(2), "sclk": m.group(3),
            "kind": m.group(4), "att": m.group(5), "enc": m.group(6)}


def normalize(stem):
    """归一成 LID 用的主干：去目录、去扩展名、去尾部版本号、转小写。"""
    s = os.path.splitext(os.path.basename(stem))[0]
    return re.sub(r"\d+$", "", s).lower()


# ---------------------------------------------------------------------------
# 清单：lid 列表 -> 按定位分组
# ---------------------------------------------------------------------------

def _inventory_path():
    return os.path.join(CACHE, "watson_inventory.csv")


def load_inventory(folder=None, refresh=False):
    """WATSON 彩色（SIF_）产品的 lid 主干列表。优先用本地缓存。

    查找顺序：缓存 CSV -> spike 的清单（开发机） -> PDS API 全量拉取。
    """
    p = _inventory_path()
    if not refresh and os.path.exists(p):
        return _read_inventory(p)
    if folder and not refresh:
        # 数据文件夹附近可能有 spike 留下的全量清单
        for cand in (os.path.join(folder, "watson_inventory.csv"),
                     os.path.join(os.path.dirname(folder),
                                  "06_算法实验与依据", "watson_inventory.csv")):
            if os.path.exists(cand):
                out = _read_inventory(cand)
                if out:
                    _ensure_cache_dir()
                    with open(p, "w", encoding="utf-8") as f:
                        f.write("\n".join(out) + "\n")
                    return out
    return _fetch_inventory()


def _read_inventory(p):
    """认三种行格式：裸 stem、裸 lidvid、PDS 全量清单的 "P,lidvid"。"""
    out = []
    with open(p, encoding="utf-8", errors="replace") as f:
        for line in f:
            lid = line.strip()
            if "," in lid:
                lid = lid.split(",", 1)[1]
            if "::" in lid:
                lid = lid.split("::")[0]
            parts = lid.split(":")
            stem = parts[-1] if parts else ""
            if stem.startswith("sif_"):
                out.append(stem)
    return out


def _ensure_cache_dir():
    os.makedirs(CACHE, exist_ok=True)


def _fetch_inventory():
    """从 PDS 归档拿全量清单。归档不给目录列表时返回空。"""
    _ensure_cache_dir()
    try:
        req = urllib.request.Request(
            INV_URL + "data_watson_imgops.csv", headers=UA)
        with urllib.request.urlopen(req, timeout=120) as r:
            data = r.read().decode("utf-8", "replace")
        if "sif_" in data:
            with open(_inventory_path(), "w", encoding="utf-8") as f:
                f.write(data)
            return _read_inventory(_inventory_path())
    except Exception:                                # noqa: BLE001
        pass
    return []


# ---------------------------------------------------------------------------
# 候选：同定位的 WATSON 彩色产品
# ---------------------------------------------------------------------------

def candidates(aci_image_name, inventory, same_sol=True, sol_window=30):
    """同定位的 WATSON 彩色 FDR/RAD 产品主干列表。

    "同定位" = 姿态段的 N 后 7 位相同（机械臂定位编号）。
    注意**不能比整段姿态串**：两台相机同定位下跑各自的序列号
    （ACI 是 SRLC10600，WATSON 是 SRLC00451 之类），整段比永远配不上。

    same_sol=True 只看同 sol；False 时放宽到 ±sol_window 个 sol
    （漫游车常先用 WATSON 拍研磨斑、隔几天才做 SHERLOC——sol 865 的
    教训：同 sol 只有别的采集点，跨 sol 才有真共视）。跨 sol 候选量
    大一个量级，靠几何共视判定兜底。

    同定位下 WATSON 常有一整组对焦堆栈（间隔 3 秒十几帧），全算候选，
    交给几何判定去挑。
    """
    info = parse_name(aci_image_name)
    if not info:
        return []
    arm = info["att"][1:8]                     # N 后 7 位 = 机械臂定位
    try:
        sol_n = int(info["sol"])
    except ValueError:
        sol_n = None
    out = []
    for stem in inventory:
        w = parse_name(stem)
        if not w or w["att"][1:8] != arm:
            continue
        if same_sol:
            if w["sol"] != info["sol"]:
                continue
        elif sol_n is not None:
            try:
                if abs(int(w["sol"]) - sol_n) > sol_window:
                    continue
            except ValueError:
                pass
        if w["kind"] not in ("FDR", "RAD"):     # 要彩色成品，不要原始帧
            continue
        out.append(stem)
    return sorted(out)


# ---------------------------------------------------------------------------
# CAHVOR 几何：与 spike 的实现同源
# ---------------------------------------------------------------------------

def _tag_floats(t, tag, n=3):
    i = t.find("<geom:%s>" % tag)
    if i < 0:
        return None
    seg = t[i:t.index("</geom:%s>" % tag, i)]
    return [float(v) for v in
            re.findall(r">\s*(-?[\d.eE+-]+)\s*<", seg)][:n]


def load_cam_xml(text, ns=None, nl=None):
    """从 .xml 文本解析 CAHV。ns/nl 缺省时从标签读。"""
    def _int(name):
        m = re.search(r"<(?:img:|%s)" % r"" + name + r">\s*(\d+)" if False
                      else r"<img:%s>\s*(\d+)" % name, text)
        return int(m.group(1)) if m else None

    C = _tag_floats(text, "Vector_Center")
    A = _tag_floats(text, "Vector_Axis")
    H = _tag_floats(text, "Vector_Horizontal")
    V = _tag_floats(text, "Vector_Vertical")
    if not (C and A and H and V):
        # 占位全 0 的模型（Guillaumes/Dourbes 实测过）：当不可用
        if not any(C or []):
            return None
        return None
    return {"C": C, "A": A, "H": H, "V": V,
            "ns": _int("samples") or ns, "nl": _int("lines") or nl}


def _download(url, path, timeout=240):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".part"
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r, \
            open(tmp, "wb") as f:
        while True:
            buf = r.read(1 << 18)
            if not buf:
                break
            f.write(buf)
    os.replace(tmp, path)


def get_cam_xml(stem, folder=None):
    """拿 WATSON 产品的相机模型 xml（带缓存）。返回文本或 None。"""
    p = os.path.join(CACHE, "xml", stem + ".xml")
    if os.path.exists(p):
        try:
            return open(p, encoding="utf-8", errors="replace").read()
        except Exception:                            # noqa: BLE001
            pass
    url = lookup_url(stem, want_label=True)
    if not url:
        return None
    try:
        _download(url, p)
        return open(p, encoding="utf-8", errors="replace").read()
    except Exception:                                # noqa: BLE001
        return None


def lookup_url(stem, want_label=False):
    """PDS 检索 API：产品主干 -> 数据文件（或标签文件）直链。"""
    coll = "data_watson_imgops"
    q = urllib.parse.quote(
        'lid eq "urn:nasa:pds:mars2020_imgops:%s:%s"' % (coll, stem))
    try:
        req = urllib.request.Request(API + "?limit=1&q=" + q, headers=UA)
        with urllib.request.urlopen(req, timeout=90) as r:
            d = json.loads(r.read().decode("utf-8", "replace"))
    except Exception:                                # noqa: BLE001
        return None
    if not (d.get("summary") or {}).get("hits"):
        return None
    props = (d.get("data") or [{}])[0].get("properties") or {}
    key = ("ops:Label_File_Info.ops:file_ref" if want_label
           else "ops:Data_File_Info.ops:file_ref")
    v = props.get(key) or props.get("ops:Data_File_Info.ops:file_ref")
    if not v:
        return None
    return v[0] if isinstance(v, list) else v


# --- 共视判定 -------------------------------------------------------------

def _ray(cam, x, y):
    return (cross(sub(cam["H"], mul(x, cam["A"])),
                  sub(cam["V"], mul(y, cam["A"]))))


def sub(a, b):
    return [a[i] - b[i] for i in range(3)]


def mul(k, v):
    return [k * x for x in v]


def cross(a, b):
    return [a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]]


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def project(cam, P):
    v = sub(P, cam["C"])
    den = dot(v, cam["A"])
    if abs(den) < 1e-12:
        return None
    return (dot(v, cam["H"]) / den, dot(v, cam["V"]) / den)


def cosees(aci_cam, wat_cam, phys_lo=70, phys_hi=150):
    """ACI 四角能否在**物理工作距离内**完整投进 WATSON 画幅。

    只认 phys_lo–phys_hi mm（SHERLOC 的臂上工作距离量级）内的共视。
    之前扫到 400mm：远处的"虚拟平面"上投影必然缩小，指向相近但
    实际没拍同一处的 WATSON 也会被判成共视（sol 865 实测踩过：
    20 个候选在 320mm 全部"通过"，其实 WATSON 拍的是别的采集点）。

    返回 (bool, standoff)。
    """
    corners = [(0, 0), (aci_cam["ns"] - 1, 0),
               (aci_cam["ns"] - 1, aci_cam["nl"] - 1), (0, aci_cam["nl"] - 1)]
    n = aci_cam["A"]

    feas = []
    for d in STAND_RANGE:
        if not (phys_lo <= d <= phys_hi):
            continue
        dist = d / 1000.0
        P0 = [aci_cam["C"][i] + dist * n[i] for i in range(3)]
        ok = True
        for x, y in corners:
            r = _ray(aci_cam, x, y)
            den = dot(r, n)
            if abs(den) < 1e-12:
                ok = False
                break
            t = dot(sub(P0, aci_cam["C"]), n) / den
            if t <= 0:
                ok = False
                break
            pt = project(wat_cam, [aci_cam["C"][i] + t * r[i]
                                   for i in range(3)])
            if (pt is None or pt[0] < 0 or pt[0] > wat_cam["ns"] - 1
                    or pt[1] < 0 or pt[1] > wat_cam["nl"] - 1):
                ok = False
                break
        if ok:
            feas.append(d)
    if feas:
        return True, float(feas[len(feas) // 2])
    return False, 0.0


# ---------------------------------------------------------------------------
# 主入口
# ---------------------------------------------------------------------------

def find_pair(aci_image_name, folder=None, inventory=None, log=None):
    """给一张 ACI 图找共视的 WATSON 彩色产品。

    返回 {"stem", "url", "standoff"} 或 None。log 是可选的收集函数。
    """
    def say(msg):
        if log:
            log(msg)

    inv = inventory if inventory is not None else load_inventory(folder)
    if not inv:
        say("  拿不到 WATSON 全量清单，跳过")
        return None

    # ACI 自己的相机模型从数据文件夹的 .xml 拿（reader 已有 label 解析，
    # 但 CAHVOR 只在 PDS 的完整标签里；先试数据文件夹，缺了再走 API）
    aci_xml_text = _aci_cam_text(aci_image_name, folder)
    if not aci_xml_text:
        say("  拿不到 ACI 的相机模型")
        return None
    ns = nl = None
    try:
        m = re.search(r"<img:samples>\s*(\d+)", aci_xml_text)
        n2 = re.search(r"<img:lines>\s*(\d+)", aci_xml_text)
        ns, nl = (int(m.group(1)), int(n2.group(1))) if m and n2 else (None, None)
    except Exception:                                # noqa: BLE001
        pass
    aci_cam = load_cam_xml(aci_xml_text, ns, nl)
    if not aci_cam:
        say("  ACI 的 CAHVOR 是占位值")
        return None

    # 两阶段搜索：先同 sol（快），没有再放宽到 ±30 sol。
    # 漫游车常先用 WATSON 拍研磨斑、隔几天才做 SHERLOC——跨 sol 的
    # 共视是真实存在的（Alfalfa 367/370、sol 865 的 857-860）。
    stages = [("同 sol", candidates(aci_image_name, inv, same_sol=True))]
    cands2 = candidates(aci_image_name, inv, same_sol=False, sol_window=30)
    stages.append(("跨 ±30 sol", [c for c in cands2 if c not in set(
        stages[0][1])]))

    tried = 0
    for stage_name, cands in stages:
        if not cands:
            if stage_name.startswith("同"):
                say("  清单里没有同定位的 WATSON 彩色产品，放宽到跨 sol")
            continue
        say("  [%s] 候选 %d 个，逐个做几何共视判定" % (stage_name, len(cands)))
        for stem in cands:
            tried += 1
            text = get_cam_xml(stem, folder)
            if not text:
                continue
            m = re.search(r"<img:samples>\s*(\d+)", text)
            n2 = re.search(r"<img:lines>\s*(\d+)", text)
            cam = load_cam_xml(text,
                               int(m.group(1)) if m else None,
                               int(n2.group(1)) if n2 else None)
            if not cam:
                say("  %s：CAHVOR 占位值，跳过" % stem[-24:])
                continue
            ok, d = cosees(aci_cam, cam)
            if ok:
                say("  %s：共视 ✓（standoff %.0f mm，%s）"
                    % (stem[-24:], d, stage_name))
                return {"stem": stem, "url": lookup_url(stem),
                        "standoff": d}
            say("  %s：不共视" % stem[-24:])
    if tried:
        say("  —— 共 %d 个候选全都不共视：WATSON 在这个定位（含跨 ±30 sol）"
            "拍的是别的采集点，这份 ACI 没有可上色的彩图" % tried)
    return None


def _aci_cam_text(aci_image_name, folder):
    """找 ACI 产品的完整 PDS4 xml 文本。先翻数据文件夹，再走 API。"""
    stem = normalize(aci_image_name)
    if folder:
        for root, _dirs, files in os.walk(folder):
            for f in files:
                if (f.lower().endswith(".xml")
                        and normalize(f) == stem):
                    try:
                        return open(os.path.join(root, f), encoding="utf-8",
                                    errors="replace").read()
                    except Exception:                # noqa: BLE001
                        pass
    url = None
    try:
        q = urllib.parse.quote(
            'lid eq "urn:nasa:pds:mars2020_imgops:data_aci_imgops:%s"' % stem)
        req = urllib.request.Request(API + "?limit=1&q=" + q, headers=UA)
        with urllib.request.urlopen(req, timeout=90) as r:
            d = json.loads(r.read().decode("utf-8", "replace"))
        if (d.get("summary") or {}).get("hits"):
            props = d["data"][0]["properties"]
            v = props.get("ops:Label_File_Info.ops:file_ref")
            if v:
                url = v[0] if isinstance(v, list) else v
    except Exception:                                # noqa: BLE001
        return None
    if not url:
        return None
    try:
        req = urllib.request.Request(url, headers=UA)
        with urllib.request.urlopen(req, timeout=120) as r:
            return r.read().decode("utf-8", "replace")
    except Exception:                                # noqa: BLE001
        return None
