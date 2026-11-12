# -*- coding: utf-8 -*-
"""按内容验候选：把 ACI 的结构图拿到每张候选 WATSON 上全搜一遍。

和 screen_watson.py 配套：那边按"指向角 + 预测像素尺度"筛出候选，这边下图验内容。
关键区别（也是必须按候选分别算的原因）：
  * 每张候选的**像素尺度不一样**（23 mm 处 7.8 µm/px 到 224 mm 处 77 µm/px），
    所以几何预测的倍率 scale = ACI 的 µm/px ÷ 候选的 µm/px 每张都不同；
  * 带通也必须按各自的 µm/px 折算成像素 σ，否则比的根本不是同一个物理波段；
  * 阳性对照（search_control.py）说明：真对上时 |NCC| ≈ 0.95，无关纹理 ≈ 0.1。

用法：python verify_candidates.py [--top 4] [--max-um 30]
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

from colorize_pair import ACI_IMG, ACI_XML, OUT, load_cam, read_bsq, to8   # noqa: E402
from screen_watson import urls_from_lidvid, read_cam_from                   # noqa: E402
from sign_test import bandpass, zs, search                                  # noqa: E402

ACI_UM = 6.66
UA = {"User-Agent": "Mozilla/5.0"}


def read_gray(path, cam):
    """读图 → 8 位灰度。3 波段按标签里的 R/G/B 通道组合。"""
    r = read_bsq(path, cam)
    if cam["bands"] == 1 or not cam["rgb"] or cam["rgb"][0] is None:
        g = to8(r[0])
    else:
        b = [i - 1 for i in cam["rgb"]]
        g = cv2.cvtColor(np.stack([to8(r[b[i]]) for i in range(3)], axis=-1),
                         cv2.COLOR_RGB2GRAY)
    return g[:cam["nl"], :cam["ns"]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", type=int, default=4)
    ap.add_argument("--max-um", type=float, default=40.0,
                    help="只看预测尺度比这更细的候选（太粗的没意义）")
    args = ap.parse_args()

    tab = os.path.join(OUT, "screen_watson.txt")
    rows = []
    for line in open(tab, encoding="utf-8").read().splitlines()[1:]:
        p = line.split()
        if len(p) >= 7:
            rows.append((float(p[0]), p[6], float(p[1]), float(p[2]),
                         float(p[3]), float(p[5])))   # aim, stem, dist, um, off, dC
    rows.sort()
    picked = []
    for aim, stem, dist, um, off, dc in rows:
        if um > args.max_um:
            continue
        if any(abs(um - u) / u < 0.05 for _a, _s, _d, u, _o, _c in picked):
            continue
        picked.append((aim, stem, dist, um, off, dc))
        if len(picked) >= args.top:
            break

    aci = load_cam(ACI_XML)
    aci8 = to8(read_bsq(ACI_IMG, aci)[0])
    aci_ppm = 1000.0 / ACI_UM
    A = np.ascontiguousarray(zs(bandpass(aci8, 0.06 * aci_ppm / 2.355,
                                         0.5 * aci_ppm / 2.355)).astype(np.float32))

    print("候选（按指向角排序，已按预测尺度去重）：")
    for aim, stem, dist, um, off, dc in picked:
        print("  %-52s 指向角 %.2f°  靶点偏移 %.1f mm  standoff %.1f mm  %.2f µm/px"
              % (stem, aim, dist * np.tan(np.radians(aim)), dist, um))
    print("\nACI 自己 %.2f µm/px。开始逐张下载 + 全搜（位置 × 3 档尺度 × 36 转角 × 正负）\n"
          % ACI_UM)

    results = []
    for aim, stem, dist, um, off, dc in picked:
        # 该候选的 LIDVID：从清单里找同名帧，取任意能下成功的变体
        from screen_watson import collect
        lv = None
        for s, variants in collect([stem[4:8]]):
            if s == stem:
                lv = variants
                break
        if lv is None:
            print("  %s：清单里找不到" % stem)
            continue
        cam = read_cam_from(lv)
        if cam is None:
            print("  %s：标签取不到" % stem)
            continue
        _x, img_url = urls_from_lidvid(cam["lidvid"])
        dst = os.path.join(OUT, "img", os.path.basename(img_url))
        if not os.path.exists(dst):
            try:
                print("  下载 %s（%s）…" % (os.path.basename(img_url), stem), flush=True)
                with urllib.request.urlopen(urllib.request.Request(img_url, headers=UA),
                                            timeout=600) as fh:
                    open(dst, "wb").write(fh.read())
            except Exception as exc:                              # noqa: BLE001
                print("  %s 下载失败：%s" % (stem, exc))
                continue
        g = read_gray(dst, cam)
        ppm = 1000.0 / um
        ns, nl = cam["ns"], cam["nl"]
        field_mm = ns * um / 1000.0
        # 候选帧可能比 ACI 的视场小得多（这批 SIF 只有 480×480）：
        # 模板要按比例裁到能塞进候选帧，否则 matchTemplate 直接报错/被跳过。
        pred = ACI_UM / um
        max_t = int(min(ns, nl) * 0.95 / pred)
        if max_t < 64:
            print("  %s：预测倍率 %.2f，模板只剩 %d px，跳过" % (stem, pred, max_t))
            continue
        if max_t < aci8.shape[0] or max_t < aci8.shape[1]:
            cy, cx = aci8.shape[0] // 2, aci8.shape[1] // 2
            h = min(max_t, aci8.shape[0])
            w = min(max_t, aci8.shape[1])
            sub = aci8[cy - h // 2:cy + h // 2, cx - w // 2:cx + w // 2]
            print("     （ACI 视场 %.1f mm > 候选帧 %.1f mm，改用 ACI 中心 %dx%d px 那块）"
                  % (aci8.shape[1] * ACI_UM / 1000, field_mm, sub.shape[1], sub.shape[0]))
            Ause = np.ascontiguousarray(zs(bandpass(
                sub, 0.06 * aci_ppm / 2.355, 0.5 * aci_ppm / 2.355)).astype(np.float32))
        else:
            Ause = A
        # 带通上限不能超过候选帧本身能容纳的尺度（否则等于全帧低频）
        hi_mm = min(0.5, min(ns, nl) * um / 1000.0 / 8.0)
        W = np.ascontiguousarray(zs(bandpass(
            g, 0.06 * ppm / 2.355, hi_mm * ppm / 2.355)).astype(np.float32))
        # 阴影/暗区置零（尺度随候选的 µm/px 变）
        k = max(3, int(round(min(field_mm, 8.0) * ppm)) | 1)
        loc = cv2.blur(g.astype(np.float32), (k, k))
        W[loc < 0.55 * float(np.median(loc))] = 0.0
        print("  %s：%dx%d px（视场 %.1f mm），带通上限 %.2f mm，预测倍率 %.3f"
              % (stem, ns, nl, field_mm, hi_mm, pred), flush=True)
        bmax, bmin, amax, amin = -9.0, 9.0, None, None
        for s in (pred * 0.75, pred * 0.9, pred, pred * 1.24, pred * 1.5):
            if s <= 0.05 or s > 10:
                continue
            if max(Ause.shape) * s > min(ns, nl):
                continue
            for ang in range(0, 360, 10):
                mx, lm = search(W, Ause, ang, s, "max")
                mn, ln = search(W, Ause, ang, s, "min")
                if np.isfinite(mx) and mx > bmax:
                    bmax, amax = mx, (s, ang, lm)
                if np.isfinite(mn) and mn < bmin:
                    bmin, amin = mn, (s, ang, ln)
        if amax is None:
            print("     没有一组能算（模板仍太大）")
            continue
        score = max(bmax, -bmin)
        verdict = ("** 内容对上了（这才是那把钥匙）**" if score > 0.4
                   else ("有点意思，值得再看" if score > 0.2 else "内容上没有共同结构"))
        print("     max %+.3f%s   min %+.3f%s   → %s"
              % (bmax, (" @尺度%.3f/%d°%s" % amax) if amax else "",
                 bmin, (" @尺度%.3f/%d°%s" % amin) if amin else "", verdict),
              flush=True)
        results.append((score, stem, um, bmax, bmin))

    print("\n汇总（按 |相关| 排序）：")
    for score, stem, um, bmax, bmin in sorted(results, reverse=True):
        print("  |%.3f|  %.2f µm/px  max %+.3f  min %+.3f  %s"
              % (score, um, bmax, bmin, stem))
    return 0


if __name__ == "__main__":
    sys.exit(main())
