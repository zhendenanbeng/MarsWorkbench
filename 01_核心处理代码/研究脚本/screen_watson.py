# -*- coding: utf-8 -*-
"""不再相信"几何配对"，改用**内容**配对——分两步，第一步只花 80 KB/张。

第一步（只下 XML，不下图）：枚举某几个 sol 上所有 WATSON 产品，从标签里取
CAHVOR，判断"这张 WATSON 当时是不是正对着 ACI 那个靶点"，顺便算出
**它当时的像素尺度（µm/px）**：

    T = C_aci + d_aci · A_aci        （靶点，d_aci 由 RLS 锚定 ≈108 mm）
    指向角 = angle(A_wat, T − C_wat)  （< 0.5° 才算瞄着同一个点）
    standoff = |T − C_wat|  →  像素尺度 = standoff / |H_wat|  (µm/px)

这一步能回答三件事：
  * 到底有哪些 WATSON 帧瞄过这个点（现在用的那张是不是唯一/最好的）
  * 有没有更近的帧（更细的采样，直接解决 2.83 倍的分辨率差）
  * 每一张的预测尺度，用来判断"跟 ACI 的 6.66 µm/px 差多少"

第二步（`--fetch N`）：把指向角最小的 N 张图下下来（2–11 MB），用
带极性的粗尺度全搜（`sign_test.py` 里那套：位置 × 尺度 × 转角，max 和 min
都要看）验一验内容。真对上时 |NCC| > 0.4；< 0.15 就是不同的地方。

用法：
    python screen_watson.py                 # 只做第一步（需要网络）
    python screen_watson.py --sols 0185 0186 0187
    python screen_watson.py --fetch 3       # 再做第二步
"""
from __future__ import annotations

import argparse
import os
import re
import sys
import urllib.request

import numpy as np
import cv2

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from colorize_pair import (ACI_IMG, ACI_XML, OUT, load_cam, read_bsq, to8,  # noqa: E402
                           homography)
from sign_test import bandpass, zs, search                                 # noqa: E402

INV = os.path.join(HERE, "watson_inventory.csv")
XMLDIR = os.path.join(OUT, "xml")
BASE = ("https://pds-imaging.jpl.nasa.gov/archive/m20/cumulative/"
        "mars2020_imgops/data_watson_imgops/sol/%s/ids/fdr/shrlc/%s")
UA = {"User-Agent": "Mozilla/5.0"}
ACI_UM = 6.66
D_ACI_MM = 108.0        # RLS 锚定的 ACI 工作距离（geom 给的是 129，见 README）


def urls_from_lidvid(lidvid):
    """LIDVID → (标签地址, 图像地址)，**不需要检索 API**。

    文件名 = LID 名字大写 + 版本号（两位）；目录里的 sol 是 5 位补零，
    而文件名里是 4 位（`SIF_0185_...` ↔ `sol/00185/`）——这两个都好踩。
        urn:...:sif_0185_0683368184_652fdr_n0070000srlc00720_0000lmj::4.0
        → .../sol/00185/ids/fdr/shrlc/SIF_0185_0683368184_652FDR_..._0000LMJ04.xml
    """
    lid, _, ver = lidvid.partition("::")
    stem = lid.rsplit(":", 1)[-1]
    m = re.match(r"[a-z0-9]{3}_(\d{4})_(\d{10})_", stem)
    if not m:
        return None, None
    sol = m.group(1).zfill(5)
    name = stem.upper() + ver.split(".")[0].zfill(2)
    return BASE % (sol, name + ".xml"), BASE % (sol, name + ".IMG")


def read_cam_from(lidvids):
    """取某个"帧"（同一时刻的若干产品变体）的相机模型：逐个变体试，成一个就行。

    坑：文件名末尾两位是产品版本，**不一定**等于 LIDVID 的版本号，
    所以不能只拼一个地址——同一帧的 ecm/edr/sif 变体都试一遍最省事。
    """
    if isinstance(lidvids, str):
        lidvids = [lidvids]
    for lidvid in lidvids:
        xml_url, _ = urls_from_lidvid(lidvid)
        if xml_url is None:
            continue
        dst = os.path.join(XMLDIR, os.path.basename(xml_url))
        if not (os.path.exists(dst) and os.path.getsize(dst) > 5000):
            try:
                req = urllib.request.Request(xml_url, headers=UA)
                with urllib.request.urlopen(req, timeout=60) as fh:
                    blob = fh.read()
                if len(blob) < 5000:
                    continue
                open(dst, "wb").write(blob)
            except Exception:                                    # noqa: BLE001
                continue
        cam = load_cam(dst)
        if cam["H"] is not None and np.linalg.norm(cam["H"]) > 1:
            cam["lidvid"] = lidvid
            return cam
    return None


def collect(sols):
    """从全量清单里挑出这些 sol 的 WATSON 帧（同一时刻算一帧），按优先级排列。"""
    frames = {}
    for line in open(INV, encoding="utf-8"):
        if not line.startswith("P,"):
            continue
        lidvid = line.strip().split(",", 1)[1]
        if "data_watson_imgops" not in lidvid:
            continue
        stem = lidvid.partition("::")[0].rsplit(":", 1)[-1]
        if stem[4:8] not in sols:
            continue
        m = re.match(r"[a-z0-9]{3}_\d{4}_(\d{10})_", stem)
        if not m:
            continue
        # 越靠前越优先：彩图 SIF > SI1 > 其它；FDR/EDR 优先于 ECM；
        # 跳过 "t..."（大概率是缩略图/工程件）和 eby（工程相机）
        prio = 0 if stem.startswith("sif") else (1 if stem.startswith("si1") else 9)
        if "eby" in stem or "_t0" in stem:
            prio += 5
        if "ecm" in stem:
            prio += 2
        frames.setdefault(m.group(1), []).append((prio, stem, lidvid))
    out = []
    for _ts, vs in frames.items():
        vs.sort()
        out.append((vs[0][1], [v[2] for v in vs]))
    out.sort()
    return out


def axis_gap(c1, c2):
    a1, a2 = c1["A"], c2["A"]
    w0 = c1["C"] - c2["C"]
    b = float(a1 @ a2)
    d, e = float(a1 @ w0), float(a2 @ w0)
    den = 1.0 - b * b
    if den < 1e-9:
        return float(np.linalg.norm(w0)), float("nan")
    t = (b * e - d) / den
    s = (e - b * d) / den
    return float(np.linalg.norm(c1["C"] + t * a1 - (c2["C"] + s * a2))), float(t)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sols", nargs="*", default=["0185", "0186", "0187"])
    ap.add_argument("--fetch", type=int, default=0,
                    help="再下前 N 张图做内容验证（0 = 只筛 XML）")
    args = ap.parse_args()

    aci = load_cam(ACI_XML)
    print("ACI %s  C=(%.4f,%.4f,%.4f)  |H|=%.0f px（IFOV %.2f µrad/px）"
          % (os.path.basename(ACI_XML)[:46], *aci["C"], np.linalg.norm(aci["H"]),
             1e6 / np.linalg.norm(aci["H"])))
    T = aci["C"] + (D_ACI_MM / 1000.0) * aci["A"]
    print("靶点 T = C_aci + %.0f mm·A_aci = (%.4f, %.4f, %.4f)\n"
          % (D_ACI_MM, *T))

    cands = collect(args.sols)
    print("清单里 sol %s 的 WATSON 产品去重后 %d 个，开始取标签……"
          % (",".join(args.sols), len(cands)))
    rows = []
    for k, (stem, lidvids) in enumerate(cands, 1):
        cam = read_cam_from(lidvids)
        if cam is None:
            if k % 25 == 0:
                print("   …%d/%d（标签取不到）" % (k, len(cands)), flush=True)
            continue
        v = T - cam["C"]
        dist = float(np.linalg.norm(v))
        aim = float(np.degrees(np.arccos(np.clip(cam["A"] @ (v / dist), -1, 1))))
        ifov = 1e6 / float(np.linalg.norm(cam["H"]))          # µrad/px
        um_px = ifov * 1e-6 * dist * 1e6                      # µm/px
        gap, _t = axis_gap(cam, aci)
        rows.append((aim, stem, cam["lidvid"], dist * 1000, um_px, gap * 1000,
                     float(np.linalg.norm(cam["C"] - aci["C"])) * 1000,
                     dist * 1000 * np.tan(np.radians(aim))))
        if k % 25 == 0:
            print("   …%d/%d" % (k, len(cands)), flush=True)

    rows.sort()
    print("\n按『指向角』排序（越小越说明当时正对着同一个点）：")
    print("  指向角°  standoff(mm)  µm/px   偏移(mm)  轴间距(mm)  |C差|(mm)  产品")
    keep = []
    for aim, stem, lidvid, dmm, umpx, gapmm, dcmm, offmm in rows[:45]:
        flag = ""
        # 注意：WATSON 是 32° 广角机，"对着靶点"不需要 0.5° 那么严；
        # 真正有意义的是"靶点偏离它光轴多少毫米"（偏移列）。
        if offmm < 5.5:
            flag = "  ← 靶点还在它视场内"
            if umpx < 0.8 * ACI_UM:
                flag += "，且比现在这张细"
            keep.append((aim, stem, lidvid))
        print("  %6.2f %12.1f %8.2f %9.1f %11.1f %10.1f  %s%s"
              % (aim, dmm, umpx, offmm, gapmm, dcmm, stem, flag))
    print("\n（ACI 自己 %.2f µm/px；现在用的那张 WATSON 是 %.2f µm/px）"
          % (ACI_UM, ACI_UM / 0.353))
    with open(os.path.join(OUT, "screen_watson.txt"), "w", encoding="utf-8") as fh:
        fh.write("aim_deg standoff_mm um_per_px offset_mm axis_gap_mm dC_mm product\n")
        for aim, stem, _lidvid, dmm, umpx, gapmm, dcmm, offmm in rows:
            fh.write("%.2f %.1f %.2f %.1f %.2f %.1f %s\n"
                     % (aim, dmm, umpx, offmm, gapmm, dcmm, stem))
    print("全表写到 out/screen_watson.txt（%d 行）" % len(rows))

    if args.fetch and keep:
        n = min(args.fetch, len(keep))
        print("\n下载指向角最小的 %d 张做内容验证……" % n)
        ar = read_bsq(ACI_IMG, aci)
        aci8 = to8(ar[0])
        A = np.ascontiguousarray(zs(bandpass(
            aci8, 0.06 * (1000.0 / ACI_UM) / 2.355,
            0.5 * (1000.0 / ACI_UM) / 2.355)).astype(np.float32))
        for aim, stem, lidvid in keep[:n]:
            _x, img_url = urls_from_lidvid(lidvid)
            dst = os.path.join(OUT, "img", os.path.basename(img_url))
            if not os.path.exists(dst):
                try:
                    req = urllib.request.Request(img_url, headers=UA)
                    with urllib.request.urlopen(req, timeout=300) as fh:
                        open(dst, "wb").write(fh.read())
                except Exception as exc:                          # noqa: BLE001
                    print("  %s 下载失败：%s" % (stem, exc))
                    continue
            cam = read_cam_from(lidvid)
            wr = read_bsq(dst, cam)
            b = [i - 1 for i in (cam["rgb"] if cam["rgb"][0] else [1, 1, 1])]
            wg = to8(wr[0]) if cam["bands"] == 1 else cv2.cvtColor(
                np.stack([to8(wr[b[i]]) for i in range(3)], axis=-1),
                cv2.COLOR_RGB2GRAY)
            if wg.shape != (cam["nl"], cam["ns"]):
                wg = wg[:cam["nl"], :cam["ns"]]
            W = np.ascontiguousarray(zs(bandpass(
                wg, 0.06 * (1000.0 / (ACI_UM / 0.353)) / 2.355,
                0.5 * (1000.0 / (ACI_UM / 0.353)) / 2.355)).astype(np.float32))
            bmax = bmin = 0.0
            for scale in (0.30, 0.353, 0.44):
                for ang in range(0, 360, 10):
                    mx, _ = search(W, A, ang, scale, "max")
                    mn, _ = search(W, A, ang, scale, "min")
                    if np.isfinite(mx):
                        bmax = max(bmax, mx)
                    if np.isfinite(mn):
                        bmin = min(bmin, mn)
            verdict = ("** 可能就是对的那张 **" if max(bmax, -bmin) > 0.4
                       else "内容上没有共同结构")
            print("  %s：max %+.3f  min %+.3f  → %s"
                  % (stem, bmax, bmin, verdict))
    return 0


if __name__ == "__main__":
    sys.exit(main())
