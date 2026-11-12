# -*- coding: utf-8 -*-
"""两件判决实验（都只用 WATSON 自己的数据）：

A. **同一帧的不同产品变体，位姿是否逐位一致？**
   清单里同一时刻有十几个变体（ecm/edr/fdr × lmj/luj/085j/300j…）。
   如果同一帧的两种产品给出的 C/A 不同，那"标签位姿不代表这张图"就直接成立。
   顺手统计每个变体的像幅（480×480 子帧 vs 1648×1200 全画幅）和不需下载就能
   拿到的文件大小（HEAD 请求），好挑出"全画幅 + 小文件"来做 B。

B. **同一落位、两张全画幅 WATSON 互验**（WATSON 侧的阳性对照）。
   同落位时模型预测的配准就是"几乎零位移"，所以直接比 NCC 即可：
   若产品与标签一致，应当 0.9 上下；若只有 0.1–0.2，则产品/标签层面就有问题。
   作为对照，再拿**不同落位**的一对算一遍（那个本来就该低）。

用法：python variant_check.py [--sol 0186] [--max-mb 3]
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

from colorize_pair import OUT, read_bsq, to8
from screen_watson import INV, urls_from_lidvid

UA = {"User-Agent": "Mozilla/5.0"}
XDIR = os.path.join(OUT, "xml")


def vec(t, tag):
    i = t.find("<geom:%s>" % tag)
    if i < 0:
        return None
    seg = t[i:t.index("</geom:%s>" % tag, i)]
    v = [float(x) for x in re.findall(r">\s*(-?[\d.eE+-]+)\s*<", seg)][:3]
    return np.array(v) if len(v) == 3 else None


def num(t, tag, cast=float):
    m = re.search(r"<%s>\s*(-?[\d.eE+-]+)" % tag, t)
    return cast(m.group(1)) if m else None


def load_xml(url):
    dst = os.path.join(XDIR, os.path.basename(url))
    if not (os.path.exists(dst) and os.path.getsize(dst) > 5000):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA),
                                        timeout=60) as fh:
                blob = fh.read()
            if len(blob) < 5000:
                return None
            open(dst, "wb").write(blob)
        except Exception:                                          # noqa: BLE001
            return None
    t = open(dst, encoding="utf-8", errors="replace").read()
    out = dict(C=vec(t, "Vector_Center"), A=vec(t, "Vector_Axis"),
               H=vec(t, "Vector_Horizontal"), V=vec(t, "Vector_Vertical"))
    if out["C"] is None:
        return None
    out["ns"] = num(t, "img:samples", int)
    out["nl"] = num(t, "img:lines", int)
    out["fl"] = num(t, "img:first_line", int)
    out["fs"] = num(t, "img:first_sample", int)
    out["exp"] = num(t, "img:exposure_duration")
    out["bands"] = num(t, "img:bands", int) or 3
    rgb = []
    for c in ("red", "green", "blue"):
        m = re.search(r"<disp:%s_channel_band>\s*(\d+)" % c, t)
        rgb.append(int(m.group(1)) if m else None)
    out["rgb"] = rgb if rgb[0] else [1, 2, 3]
    out["file"] = os.path.basename(dst)
    return out


def head_size(url):
    try:
        req = urllib.request.Request(url, headers=UA, method="HEAD")
        with urllib.request.urlopen(req, timeout=60) as fh:
            n = fh.headers.get("Content-Length")
            return int(n) if n else None
    except Exception:                                              # noqa: BLE001
        return None


def frames(sol):
    """按时刻分组，返回 [(ts, [(lidvid, stem), …]), …]"""
    g = {}
    for line in open(INV, encoding="utf-8"):
        if not line.startswith("P,"):
            continue
        lv = line.strip().split(",", 1)[1]
        if "data_watson_imgops" not in lv:
            continue
        stem = lv.partition("::")[0].rsplit(":", 1)[-1]
        if stem[4:8] != sol:
            continue
        m = re.match(r"[a-z0-9]{3}_\d{4}_(\d{10})_", stem)
        if m:
            g.setdefault(m.group(1), []).append((lv, stem))
    return sorted(g.items())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sol", default="0186")
    ap.add_argument("--max-mb", type=float, default=3.0)
    args = ap.parse_args()

    allf = frames(args.sol)
    # 只看 SHERLOC 相关的、且产品数不太夸张的时刻
    allf = [(ts, v) for ts, v in allf if "srlc" in v[0][1]]
    print("sol %s：清单里 WATSON 帧 %d 个（去重后按时刻）" % (args.sol, len(allf)))

    # ---- A：同一帧不同变体的位姿一致性 ----
    print("\nA. 同一帧不同产品变体的位姿比较（只报不一致的）")
    bad = 0
    full = []      # (ts, stem, lidvid, bytes)
    for ts, variants in allf:
        cams = []
        for lidvid, stem in variants[:6]:
            x, i = urls_from_lidvid(lidvid)
            c = load_xml(x)
            if c is None:
                continue
            cams.append((stem, lidvid, c))
        if len(cams) < 2:
            continue
        base = cams[0][2]
        for stem, lidvid, c in cams[1:]:
            dC = float(np.linalg.norm(c["C"] - base["C"])) * 1000
            dA = float(np.degrees(np.arccos(np.clip(c["A"] @ base["A"], -1, 1))))
            dH = float(np.linalg.norm(c["H"] - base["H"]) / np.linalg.norm(base["H"]))
            if dC > 0.05 or dA > 0.01 or dH > 1e-4:
                bad += 1
                print("  !! %s vs %s：ΔC %.3f mm  ΔA %.4f°  Δ|H| %.2e"
                      % (cams[0][0][:34], stem[:34], dC, dA, dH))
        for stem, lidvid, c in cams:
            if c["ns"] == 1648 and c["nl"] == 1200:
                _x, iu = urls_from_lidvid(lidvid)
                b = head_size(iu)
                if b and b <= args.max_mb * 1e6:
                    full.append((ts, stem, lidvid, b))
    print("  → 共检查 %d 帧，位姿不一致的 %d 处" % (len(allf), bad))

    # ---- B：同落位两张全画幅互验 ----
    print("\nB. 可用的“全画幅 + 小文件”候选（按落位分组：看 C 是否相同）")
    seen = {}
    picks = []
    for ts, stem, lidvid, b in full:
        x, _ = urls_from_lidvid(lidvid)
        c = load_xml(x)
        key = tuple(np.round(c["C"], 5))
        seen.setdefault(key, []).append((ts, stem, lidvid, b))
    for key, v in sorted(seen.items(), key=lambda kv: -len(kv[1])):
        print("  落位 C=%s：%d 张" % (str(key), len(v)))
        for ts, stem, lidvid, b in v[:4]:
            print("      %s  %.1f MB  %s" % (ts, b / 1e6, stem[:56]))
        if len(v) >= 2:
            picks.append(v[:2])
    if not picks:
        print("\n  没有找到『同一落位两张全画幅』的组合——说明这一帧的全画幅产品不多，"
              "\n  可以放宽 --max-mb 或换一个 sol/活动再试。")
        return 0

    def gray(lidvid):
        x, iu = urls_from_lidvid(lidvid)
        c = load_xml(x)
        dst = os.path.join(OUT, "img", os.path.basename(iu))
        if not os.path.exists(dst):
            print("    下载 %s（%.1f MB）…" % (os.path.basename(iu),
                                              (head_size(iu) or 0) / 1e6), flush=True)
            with urllib.request.urlopen(urllib.request.Request(iu, headers=UA),
                                        timeout=900) as fh:
                open(dst, "wb").write(fh.read())
        r = read_bsq(dst, c)
        b = [k - 1 for k in c["rgb"]]
        g = cv2.cvtColor(np.stack([to8(r[b[k]]) for k in range(3)], axis=-1),
                         cv2.COLOR_RGB2GRAY)[:c["nl"], :c["ns"]].astype(np.float32)
        return g, c

    def ncc(x, y):
        x = x - x.mean()
        y = y - y.mean()
        return float((x * y).sum() / np.sqrt((x ** 2).sum() * (y ** 2).sum()))

    def coarse(g):
        return cv2.GaussianBlur(g, (0, 0), 1.0) - cv2.GaussianBlur(g, (0, 0), 12.0)

    for ts_a, stem_a, lid_a, _b in [p[0] for p in picks[:2]]:
        ts_b, stem_b, lid_b, _b2 = [p[1] for p in picks[:2]][0]
        print("\n同落位互验：%s × %s" % (ts_a, ts_b))
        ga, ca = gray(lid_a)
        gb, cb = gray(lid_b)
        print("  %dx%d 与 %dx%d；C 相差 %.4f mm"
              % (ca["ns"], ca["nl"], cb["ns"], cb["nl"],
                 np.linalg.norm(ca["C"] - cb["C"]) * 1000))
        xa, xb = coarse(ga), coarse(gb)
        print("  零位移 NCC = %.3f" % ncc(xa[80:-80, 80:-80], xb[80:-80, 80:-80]))
        best = (-9, None)
        for dy in range(-40, 41, 2):
            for dx in range(-40, 41, 2):
                v = ncc(xa[120:-120, 120:-120],
                        np.roll(np.roll(xb, dy, 0), dx, 1)[120:-120, 120:-120])
                if v > best[0]:
                    best = (v, (dx, dy))
        print("  ±40 px 平移搜索最好 NCC = %.3f 在 %s" % best)
        break
    return 0


if __name__ == "__main__":
    sys.exit(main())
