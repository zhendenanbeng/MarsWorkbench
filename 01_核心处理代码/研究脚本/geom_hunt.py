# -*- coding: utf-8 -*-
"""用相机模型判断"哪张 WATSON 真的和某张 ACI 看同一个地方"。

背景（2026-09-15 的教训）：上一轮想靠图像特征（SIFT / 互相关）配准，卡住了。
原因不是算法不好，而是**候选帧根本不对**——sol 141 上那些 WATSON 帧和 ACI
光心相距 35 cm、看的是两个地方。用两台相机各自的 CAHVOR 模型一算就知道，
根本不用做图像匹配。

这一轮换路：论文图注里直接给了每张配套 WATSON 的图像 ID（Fig.5a/6a/7a/8a）。
把那些产品取回来，用相机模型验证是否共视，通过的才值得下载全图做上色。

共视判据（纯几何，不含任何拟合参数）：
    两条视轴 C1+t·A1 与 C2+s·A2 的最近距离。两台相机若在瞄同一个目标，
    两条视轴应在目标面上相交，最近距离只有几毫米到几厘米；差几十厘米就是各看各的。

用法：python geom_hunt.py
输出：out/geom_pairs.txt，XML 标签缓存在 out/xml/
"""

from __future__ import annotations

import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
XMLDIR = os.path.join(OUT, "xml")
os.makedirs(XMLDIR, exist_ok=True)

API = "https://pds.nasa.gov/api/search/1/products"
UA = {"User-Agent": "Mozilla/5.0"}

# 论文图注里给出的配套图像（Corpolongo et al. 2023）
# (靶点, WATSON 图像 ID, 出处, ACI 图像 ID, 出处)
PAIRS = [
    ("Guillaumes", "SI1_0160_0681181210_679FDR_N0060000SRLC00003_000095J01", "Fig.5a",
     "SC3_0162_0681352993_195FDR_N0060000SRLC11420_0000LMJ01", "Fig.5b"),
    ("Bellegarde", "SIF_0185_0683368184_652FDR_N0070000SRLC00720_0000LMJ01", "Fig.7a",
     "SC3_0186_0683479674_054FDR_N0070000SRLC11420_0000LMJ02", "Fig.7b"),
    ("Dourbes", "SI1_0257_0689775622_738FDR_N0080000SRLC00006_000095J01", "Fig.2a",
     "SC3_0257_0689786306_125FDR_N0080000SRLC11421_0000LMJ01", "Fig.2b"),
    ("Montpezat", "SIF_0346_0697661461_972FDR_N0092982SRLC01034_0000LMJ01", "Fig.6a",
     "SC3_0349_0697954126_179FDR_N0092982SRLC10600_0000LMJ01", "Fig.6b"),
    ("Alfalfa", "SIF_0367_0699544277_375FDR_N0110108SRLC08029_0000LMJ01", "Fig.8a",
     "SC3_0370_0699816293_742FDR_N0110108SRLC10600_0000LMJ01", "Fig.8b"),
    # 这两张 ACI 的图注没给 WATSON，但既然能取到 ACI 就一并算着
    ("Garde", None, "",
     "SC3_0208_0685432316_968FDR_N0071836SRLC11370_0000LMJ01", "Fig.3b"),
    ("Quartier", None, "",
     "SC3_0293_0692982585_671FDR_N0090000SRLC10600_0000LMJ01", "Fig.4b"),
]


def lid_for(img_id):
    """图像 ID → PDS LID。

    要剥掉末尾两位版本号：论文里写的是 `..._0000LMJ01`，其中 `01` 是产品版本，
    LID 只到 `..._0000lmj`。带上版本号查会一个都命中不了（实测踩过）。
    集合按前缀分流：SI* 是 WATSON，SC* 是 ACI。
    """
    stem = img_id[:-2] if img_id[-2:].isdigit() else img_id
    coll = "data_watson_imgops" if img_id.startswith("SI") else "data_aci_imgops"
    return "urn:nasa:pds:mars2020_imgops:%s:%s" % (coll, stem.lower())


def api_lookup(lid, tries=3):
    """按精确 LID 查产品，返回 label 文件地址。结果缓存到 out/。"""
    cache = os.path.join(OUT, "lookup_%s.json" % re.sub(r"\W+", "_", lid)[-70:])
    if os.path.exists(cache):
        return json.load(open(cache, encoding="utf-8"))
    url = API + "?limit=3&q=" + urllib.parse.quote('lid eq "%s"' % lid)
    last = None
    j = None
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers=dict(UA, Accept="application/json"))
            with urllib.request.urlopen(req, timeout=90) as fh:
                j = json.loads(fh.read().decode("utf-8"))
            break
        except Exception as exc:                     # noqa: BLE001
            last = exc
            time.sleep(3 + 3 * k)
    if j is None:
        raise RuntimeError("查不到 %s：%s" % (lid, last))
    if not j.get("summary", {}).get("hits"):
        raise RuntimeError("LID 无命中：%s" % lid)
    p = j["data"][0]["properties"]

    def one(key):
        """API 有的字段返回列表，有的返回字符串，统一取成字符串。"""
        v = p.get(key)
        if isinstance(v, list):
            return v[0] if v else None
        return v

    out = {"name": one("ops:Data_File_Info.ops:file_name"),
           "data_ref": one("ops:Data_File_Info.ops:file_ref"),
           "label_ref": one("ops:Label_File_Info.ops:file_ref"),
           "lidvid": one("lidvid")}
    json.dump(out, open(cache, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    return out


def fetch(url, dst):
    """下载到 dst，已有就跳过。"""
    if os.path.exists(dst) and os.path.getsize(dst) > 1000:
        return dst
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=120) as fh:
        blob = fh.read()
    with open(dst, "wb") as fh:
        fh.write(blob)
    return dst


def _vec(text, tag):
    key = "<geom:%s>" % tag
    i = text.find(key)
    if i < 0:
        raise ValueError("XML 里没有 %s" % tag)
    seg = text[i:text.index("</geom:%s>" % tag, i)]
    vals = [float(v) for v in re.findall(r">\s*(-?[\d.eE+-]+)\s*<", seg)]
    return np.array(vals[:3])


def parse_cam(path):
    """从 PDS4 XML 里取 CAHVOR 的 C/A/H/V，外加观测时刻和位姿序号。"""
    t = open(path, encoding="utf-8", errors="replace").read()
    cam = dict()
    cam["C"] = _vec(t, "Vector_Center")
    cam["A"] = _vec(t, "Vector_Axis")
    cam["H"] = _vec(t, "Vector_Horizontal")
    cam["V"] = _vec(t, "Vector_Vertical")
    m = re.search(r"<pds:start_date_time>([^<]+)", t)
    cam["time"] = m.group(1) if m else ""
    m = re.search(r"<geom:index_id>POSE</geom:index_id>\s*<geom:index_value_number>(\d+)", t)
    cam["pose"] = m.group(1) if m else "?"
    # SITE / DRIVE 决定 ROVER_NAV_FRAME 的锚点：跨 sol 比较时，
    # 这两个数不一样就说明车开走了，两个坐标系没有共同的绝对位置，
    # 光凭相机模型算不出相对关系（那需要 SPICE 定位）。
    frame = dict()
    for key in ("SITE", "DRIVE", "POSE"):
        m = re.search(r"<geom:index_id>%s</geom:index_id>\s*"
                      r"<geom:index_value_number>(\d+)" % key, t)
        frame[key] = m.group(1) if m else "?"
    cam["frame"] = "%s/%s/%s" % (frame["SITE"], frame["DRIVE"], frame["POSE"])
    return cam


def boresight_gap(c1, c2):
    """两条视轴的最短距离，以及最近点在各自轴上离光心的距离。"""
    a1, a2 = c1["A"], c2["A"]
    w0 = c1["C"] - c2["C"]
    b = float(a1 @ a2)
    d = float(a1 @ w0)
    e = float(a2 @ w0)
    den = 1.0 - b * b
    if den < 1e-9:
        return float(np.linalg.norm(w0)), float("nan"), float("nan")
    t = (b * e - d) / den
    s = (e - b * d) / den
    p1 = c1["C"] + t * a1
    p2 = c2["C"] + s * a2
    return float(np.linalg.norm(p1 - p2)), float(t), float(s)


def image_size(path):
    """从 XML 里取像幅（列, 行）。imgops 标签写作 <img:samples> / <img:lines>。"""
    t = open(path, encoding="utf-8", errors="replace").read()
    # 分开取，不能靠文档顺序——标签里 <lines> 在 <samples> 前面
    ns = re.search(r"<img:samples>\s*(\d+)", t)
    nl = re.search(r"<img:lines>\s*(\d+)", t)
    if not ns or not nl:
        return (None, None)
    return (int(ns.group(1)), int(nl.group(1)))


def ray(cam, x, y):
    """像素 (x,y) → 从光心出发的射线方向（未归一化）。"""
    return np.cross(cam["H"] - x * cam["A"], cam["V"] - y * cam["A"])


def project(cam, P):
    """三维点 → 像素坐标。"""
    v = np.asarray(P, dtype=float) - cam["C"]
    den = float(v @ cam["A"])
    if abs(den) < 1e-9:
        return None
    return (float(v @ cam["H"]) / den, float(v @ cam["V"]) / den)


def footprint_in(cam_from, size, cam_to):
    """把 cam_from 的整幅像（四角）投到 cam_to 上。

    返回 (四角像素, 是否全在画面内, 覆盖宽度像素, 覆盖高度像素)。
    目标面取「过 cam_from 光心前方 d、法向等于其光轴」的平面——即假定
    岩面正对这台相机。这是唯一需要假定的东西，d 扫一遍即可。
    """
    ns, nl = size
    out = []
    for x, y in ((0, 0), (ns - 1, 0), (ns - 1, nl - 1), (0, nl - 1)):
        out.append(ray(cam_from, x, y))
    return out


def corner_projection(cam_from, size_from, cam_to, dist_m, out_size):
    """计算 cam_from 的四角在距离 dist_m 处打到平面后、投到 cam_to 的像素。"""
    ns, nl = size_from
    P0 = cam_from["C"] + dist_m * cam_from["A"]
    n = cam_from["A"]
    pts = []
    for x, y in ((0, 0), (ns - 1, 0), (ns - 1, nl - 1), (0, nl - 1)):
        d = ray(cam_from, x, y)
        den = float(d @ n)
        if abs(den) < 1e-12:
            return None
        t = float((P0 - cam_from["C"]) @ n) / den
        if t <= 0:
            return None
        p = project(cam_to, cam_from["C"] + t * d)
        if p is None:
            return None
        pts.append(p)
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    wn, wl = out_size
    inside = min(xs) >= 0 and max(xs) < wn and min(ys) >= 0 and max(ys) < wl
    return {"corners": pts, "inside": inside,
            "w": max(xs) - min(xs), "h": max(ys) - min(ys),
            "cx": 0.5 * (min(xs) + max(xs)), "cy": 0.5 * (min(ys) + max(ys))}


def load_cam(img_id, cams):
    """取一张图的相机模型（带缓存）。"""
    if img_id in cams:
        return cams[img_id]
    info = api_lookup(lid_for(img_id))
    xml = fetch(info["label_ref"],
                os.path.join(XMLDIR, os.path.basename(info["label_ref"])))
    cam = parse_cam(xml)
    cam["name"] = info["name"]
    cam["ref"] = info["data_ref"]
    cam["xml"] = xml
    cam["size"] = image_size(xml)
    cams[img_id] = cam
    return cam


def main():
    cams = dict()
    for target, wid, _ws, aid, _as in PAIRS:
        for img_id in (wid, aid):
            if not img_id:
                continue
            try:
                cam = load_cam(img_id, cams)
                print("取到 %-12s %-58s SITE/DRIVE/POSE=%-12s C=(%6.3f,%6.3f,%6.3f)"
                      % (target, cam["name"][:58], cam["frame"], *cam["C"]))
            except Exception as exc:                 # noqa: BLE001
                print("!! %-12s %s 取不到：%s" % (target, img_id, exc))

    # 逐对做投影测试：把 ACI 整幅像投到 WATSON 上，看落不落在画面里。
    # 目标面离 ACI 的距离 d 是唯一要扫的量（岩面正对 ACI 是唯一假定）。
    D_MM = [40, 60, 80, 100, 120, 140, 160, 180, 200, 240, 300]
    rows = []
    print("\n" + "-" * 100)
    for target, wid, _ws, aid, _as in PAIRS:
        if not aid or aid not in cams:
            continue
        aci = cams[aid]
        w_size = cams[wid]["size"] if wid in cams else None
        if not wid or wid not in cams or not w_size or not aci["size"][0]:
            line = "%-11s 只有 ACI（%s %s），没有可配的 WATSON" % (target, aci["name"][:40], aci["size"])
            print(line)
            rows.append(line)
            continue
        wat = cams[wid]
        gap, t, s = boresight_gap(wat, aci)
        dsep = float(np.linalg.norm(wat["C"] - aci["C"]))
        print("%-11s 光心相距 %5.1f cm；视轴最近 %5.1f cm（在 WATSON %.1f cm / ACI %.1f cm 处）"
              % (target, dsep * 100, gap * 100, t * 100, s * 100))
        hits = []
        for d in D_MM:
            r = corner_projection(aci, aci["size"], wat, d / 1000.0, w_size)
            if r is None:
                continue
            tag = "画面内" if r["inside"] else "部分出界"
            print("     ACI standoff %3d mm → 投到 WATSON (%5.0f,%5.0f)，占 %4.0f×%4.0f px，%s"
                  % (d, r["cx"], r["cy"], r["w"], r["h"], tag))
            if r["inside"]:
                hits.append((d, r))
        if hits:
            best = max(hits, key=lambda kv: kv[1]["w"] * kv[1]["h"])
            line = ("%-11s ** 可用 **  ACI standoff %d mm 时投到 WATSON 的 %dx%d px 区域"
                    "（ACI 那 1.66 cm 覆盖约 %.0f px）"
                    % (target, best[0], int(best[1]["w"]), int(best[1]["h"]), best[1]["w"]))
            print("  ==> " + line)
            rows.append(line)
        else:
            line = "%-11s 所有 standoff 下 ACI 都投不到 WATSON 画面里 —— 不共视" % target
            print("  ==> " + line)
            rows.append(line)
    with open(os.path.join(OUT, "geom_pairs.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(rows) + "\n")
    print("\n视轴间距 = 两条视轴的最短距离；WATSON处/ACI处 = 最近点离各自光心的距离")
    print("结果写到 out/geom_pairs.txt")


if __name__ == "__main__":
    sys.exit(main())
