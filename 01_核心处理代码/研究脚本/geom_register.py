# -*- coding: utf-8 -*-
"""用两台相机各自的 CAHVOR 模型做几何配准（不靠图像特征）。

为什么这条路比 SIFT / 互相关靠谱：

* ACI 和 WATSON 的 label 里都给了 **CAHVOR 相机模型**，坐标系都是
  ROVER_NAV_FRAME，单位是米。也就是说两台相机在同一时刻各自相对车体的
  位姿是已知的——**连转塔姿态不同这件事都已经被模型吸收掉了**。
* 于是 ACI 图像里任何一个像素，都能算出一条三维射线；射线打到目标面上
  得到一个三维点；再投影进 WATSON 的模型，就得到它在 WATSON 里的位置。
  整条链路没有需要拟合的自由参数——除了**目标面离 ACI 多远**（standoff）。

所以搜索空间从"尺度 × 旋转 × 平移"（上次卡住的地方）缩成**一个标量**：
ACI 的 standoff。就扫它。

CAHVOR 的线性部分（忽略 O/R 径向项，M2020 这两台相机的畸变很小）：
    x = (v·H)/(v·A)      y = (v·V)/(v·A)        v = P - C
反过来，已知像素 (x,y) 求射线方向：
    dir = (H - x·A) × (V - y·A)

用法：python geom_register.py
输出：geom_rank.txt（每张帧 × 每个 standoff 的 NCC），best_warp.png
"""

from __future__ import annotations

import glob
import os
import sys

import numpy as np
import cv2

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "01_软件工作台"))
from sherloc.vicar import read_vicar          # noqa: E402

ACI_DIR = os.path.join(os.path.dirname(HERE), "02_练习数据与验证基准", "images", "aci")
ACI_NAME = "SC3_0141_0679457520_796ECM_N0051812SRLC10203_0000LUJ04.IMG"

# 扫 ACI 到目标面的距离。ACI 的等效焦距约 16250 px，
# 所以 10.07 µm/px 对应 standoff ≈ 16.4 cm。范围给宽一点。
STANDS_MM = [40, 60, 80, 100, 120, 140, 150, 160, 164, 170, 180, 200, 230, 260]


def parse_cahvor(path):
    """从 ODL 标签里取出 CAHVOR 的 C/A/H/V 和像幅。"""
    raw = open(path, "rb").read()
    cut = raw.find(b"END\r\n")
    cut = cut + 5 if cut >= 0 else raw.find(b"END\n") + 4
    text = raw[:cut].decode("ascii", "replace")

    def vec(n):
        key = "MODEL_COMPONENT_%d" % n
        i = text.find(key)
        if i < 0:
            raise ValueError("标签里没有 " + key)
        s = text[text.index("(", i) + 1: text.index(")", i)]
        return np.array([float(v) for v in s.split(",")])

    lab = {}
    for line in text.splitlines():
        s = line.strip()
        if "=" not in s or s.startswith("/*"):
            continue
        k, v = s.split("=", 1)
        lab[k.strip()] = v.strip()
    # 像幅在 ODL 的 SUBFRAME_REQUEST_PARMS 里
    nl = int(str(lab["LINES"]).split("<")[0])
    ns = int(str(lab["LINE_SAMPLES"]).split("<")[0])
    return {"C": vec(1), "A": vec(2), "H": vec(3), "V": vec(4),
            "O": vec(5), "R": vec(6), "nl": nl, "ns": ns}


def project(cam, P):
    """三维点 → 像素。P 可以是 Nx3。"""
    P = np.atleast_2d(P)
    v = P - cam["C"]
    den = v @ cam["A"]
    x = (v @ cam["H"]) / den
    y = (v @ cam["V"]) / den
    return np.stack([x, y], axis=1)


def ray(cam, x, y):
    """像素 → 射线方向（未归一化）。"""
    return np.cross(cam["H"] - x * cam["A"], cam["V"] - y * cam["A"])


def aci_to_watson_homography(aci, wat, stand_mm):
    """给定 ACI 的 standoff，算 ACI 像素 → WATSON 像素的单应。"""
    d = stand_mm / 1000.0
    P0 = aci["C"] + d * aci["A"]          # 目标面过这点，法向取 ACI 的光轴
    n = aci["A"]
    src = np.array([[0, 0], [aci["ns"] - 1, 0],
                    [aci["ns"] - 1, aci["nl"] - 1], [0, aci["nl"] - 1]],
                   dtype=np.float64)
    dst = []
    for x, y in src:
        dirv = ray(aci, x, y)
        t = np.dot(P0 - aci["C"], n) / np.dot(dirv, n)
        if t <= 0:
            return None, None
        dst.append(project(wat, aci["C"] + t * dirv)[0])
    dst = np.array(dst, dtype=np.float64)
    Hm = cv2.getPerspectiveTransform(src.astype(np.float32), dst.astype(np.float32))
    return Hm, dst


def feature(img, blur=None):
    """去低频 + CLAHE + 梯度幅值：跨照明也能比的"结构图"。"""
    g = img.astype(np.float32)
    bg = cv2.GaussianBlur(g, (0, 0), sigmaX=max(g.shape) / 12.0)
    g = cv2.normalize(g - bg, None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
    g = cv2.createCLAHE(2.0, (8, 8)).apply(g)
    gx = cv2.Sobel(g, cv2.CV_32F, 1, 0, ksize=3)
    gy = cv2.Sobel(g, cv2.CV_32F, 0, 1, ksize=3)
    f = cv2.magnitude(gx, gy)
    return cv2.GaussianBlur(f, (0, 0), 0.8) if blur is None else cv2.GaussianBlur(f, (0, 0), blur)


def ncc_masked(a, b, mask):
    """只在 mask 内算归一化互相关。"""
    if mask.sum() < 500:
        return float("nan")
    av, bv = a[mask], b[mask]
    av = av - av.mean()
    bv = bv - bv.mean()
    den = np.sqrt((av ** 2).sum() * (bv ** 2).sum())
    return float((av * bv).sum() / den) if den > 0 else float("nan")


def main():
    aci_p = os.path.join(ACI_DIR, ACI_NAME)
    aci_cam = parse_cahvor(aci_p)
    aci_img, _ = read_vicar(aci_p)
    aci_f = feature(aci_img)

    f_aci = np.sqrt(np.linalg.norm(aci_cam["H"]) ** 2 - np.dot(aci_cam["A"], aci_cam["H"]) ** 2)
    print("ACI  : 像幅 %dx%d  等效焦距 %.0f px" % (aci_cam["ns"], aci_cam["nl"], f_aci))
    print("       光心 %s m" % np.array2string(aci_cam["C"], precision=4))

    frames = []
    for p in sorted(glob.glob(os.path.join(HERE, "ecm", "*.IMG"))):
        cam = parse_cahvor(p)
        img, _ = read_vicar(p)
        f_w = np.sqrt(np.linalg.norm(cam["H"]) ** 2 - np.dot(cam["A"], cam["H"]) ** 2)
        frames.append((os.path.basename(p), cam, img, feature(img), f_w))
    print("WATSON: %d 帧  等效焦距 %.0f px（第一帧）" % (len(frames), frames[0][4]))
    # 两台相机光心间距：转塔姿态不同，所以这个数每帧都不一样
    print("       光心到 ACI 光心的距离：%.3f – %.3f m"
          % (min(np.linalg.norm(fr[1]["C"] - aci_cam["C"]) for fr in frames),
             max(np.linalg.norm(fr[1]["C"] - aci_cam["C"]) for fr in frames)))

    lines = []
    best = None
    for name, cam, img, feat, f_w in frames:
        wn, wl = cam["ns"], cam["nl"]
        for d in STANDS_MM:
            Hm, dst = aci_to_watson_homography(aci_cam, cam, d)
            if Hm is None:
                continue
            if dst[:, 0].max() < 0 or dst[:, 0].min() > wn or \
               dst[:, 1].max() < 0 or dst[:, 1].min() > wl:
                lines.append("%-58s d=%3d  ACI 完全落在 WATSON 视野外" % (name[:58], d))
                continue
            warp = cv2.warpPerspective(feat, Hm, (aci_cam["ns"], aci_cam["nl"]),
                                       flags=cv2.INTER_LINEAR,
                                       borderMode=cv2.BORDER_CONSTANT, borderValue=0)
            valid = cv2.warpPerspective(np.full((wl, wn), 255, np.uint8), Hm,
                                        (aci_cam["ns"], aci_cam["nl"]),
                                        flags=cv2.INTER_NEAREST) > 127
            # 边界一圈不算，避免黑边把相关拉高
            inner = np.zeros_like(valid)
            inner[40:-40, 40:-40] = True
            mask = valid & inner
            score = ncc_masked(aci_f, warp, mask)
            lines.append("%-58s d=%3d  NCC=%s  覆盖 %.0f%%"
                         % (name[:58], d,
                            ("%.3f" % score) if np.isfinite(score) else "  --  ",
                           100.0 * mask.sum() / inner.sum()))
            if np.isfinite(score) and (best is None or score > best[0]):
                best = (score, name, d, Hm)

    with open(os.path.join(HERE, "geom_rank.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print("\n前 15 名：")
    scored = [l for l in lines if "NCC=" in l and "--" not in l]
    for l in sorted(scored, key=lambda s: -float(s.split("NCC=")[1].split()[0]))[:15]:
        print("  " + l)

    if best:
        score, name, d, Hm = best
        img = dict((fr[0], fr[2]) for fr in frames)[name]
        warp = cv2.warpPerspective(img, Hm, (aci_cam["ns"], aci_cam["nl"]),
                                   flags=cv2.INTER_LINEAR)
        panel = np.hstack([aci_img, warp])
        cv2.imwrite(os.path.join(HERE, "best_warp.png"), panel)
        print("\n最佳：%s  standoff=%d mm  NCC=%.3f  → best_warp.png"
              % (name, d, score))


if __name__ == "__main__":
    main()
