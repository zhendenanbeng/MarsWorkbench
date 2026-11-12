# -*- coding: utf-8 -*-
"""从 XML 的 CAHVOR 里直接读出两台相机的 IFOV / FOV，并把"采样率之比"重算一遍。

目的：核对 README 里那条"几何预测 ACI 1 px = WATSON 0.35 px"。
CAHV 的约定（Gennery）：x = (H·(P-C))/(A·(P-C))，y = (V·(P-C))/(A·(P-C))。
所以 |H| 就是"像素焦距"，IFOV ≈ 1/|H| rad/px。
"""
from __future__ import annotations

import glob
import os
import re

import numpy as np
import cv2

HERE = os.path.dirname(os.path.abspath(__file__))
XMLDIR = os.path.join(HERE, "out", "xml")


def xml_text(path):
    return open(path, encoding="utf-8", errors="replace").read()


def tag_int(t, name, default=None):
    m = re.search(r"<%s>\s*(-?\d+)" % name.replace(":", r"\:"), t)
    return int(m.group(1)) if m else default


def tag_floats(t, tag, n=3):
    i = t.find("<geom:%s>" % tag)
    if i < 0:
        return None
    seg = t[i:t.index("</geom:%s>" % tag, i)]
    v = [float(x) for x in re.findall(r">\s*(-?[\d.eE+-]+)\s*<", seg)][:n]
    return np.array(v) if len(v) == 3 else None


def load(path):
    t = xml_text(path)
    cam = {}
    for k, tag in (("C", "Vector_Center"), ("A", "Vector_Axis"),
                   ("H", "Vector_Horizontal"), ("V", "Vector_Vertical")):
        cam[k] = tag_floats(t, tag)
    cam["ns"] = tag_int(t, "img:samples")
    cam["nl"] = tag_int(t, "img:lines")
    cam["units"] = (re.search(r"<geom:unit>\s*([^<]+)", t) or [None, "?"])[1] \
        if re.search(r"<geom:unit>\s*([^<]+)", t) else "?"
    cam["name"] = os.path.basename(path)
    return cam


def ifov(cam):
    """每像素的角分辨率（µrad），H/V 两个方向都算，顺便给出全画幅视场（度）。"""
    fa = 1.0 / np.linalg.norm(cam["H"])
    fb = 1.0 / np.linalg.norm(cam["V"])
    return (fa * 1e6, fb * 1e6,
            np.degrees(fa * cam["ns"]), np.degrees(fb * cam["nl"]))


def report(cam, label):
    ix, iy, fx, fy = ifov(cam)
    print("%-8s %-52s %4dx%-4d  |H|=%.1f px  |V|=%.1f px" %
          (label, cam["name"][:52], cam["ns"], cam["nl"],
           np.linalg.norm(cam["H"]), np.linalg.norm(cam["V"])))
    print("         IFOV %.2f x %.2f µrad/px   FOV %.3f° x %.3f°   |A|=%.4f  |C|=%.1f %s"
          % (ix, iy, fx, fy, np.linalg.norm(cam["A"]), np.linalg.norm(cam["C"]),
             cam["units"]))
    return ix, iy, fx, fy


def main():
    for p in sorted(glob.glob(os.path.join(XMLDIR, "*.xml"))):
        c = load(p)
        if c["H"] is None:
            print("!! %s 没有 CAHVOR" % os.path.basename(p))
            continue
        report(c, "WATSON" if "watson" in p.lower() or "SIF" in p or "SI1" in os.path.basename(p)
               else "ACI")

    print("\n=== Bellegarde 这一对的几何 ===")
    wat = load(os.path.join(XMLDIR, "SIF_0185_0683368184_652FDR_N0070000SRLC00720_0000LMJ04.xml"))
    aci = load(os.path.join(XMLDIR, "SC3_0186_0683479674_054FDR_N0070000SRLC11420_0000LMJ05.xml"))
    wix, wiy, wfx, wfy = ifov(wat)
    aix, aiy, afx, afy = ifov(aci)
    print("光心相距 %.2f mm" % (np.linalg.norm(wat["C"] - aci["C"]) * 1000))
    print("WATSON FOV %.3f° x %.3f°   ACI FOV %.3f° x %.3f°"
          % (wfx, wfy, afx, afy))

    # 视轴夹角与最近点
    A1, A2 = wat["A"], aci["A"]
    C1, C2 = wat["C"], aci["C"]
    print("视轴夹角 %.3f°" % np.degrees(np.arccos(np.clip(A1 @ A2, -1, 1))))

    # 采样率之比：同一块岩面上，ACI 一个像素的footprint 对应 WATSON 几个像素
    print("\nstandoff(ACI mm)  WATSON standoff(mm)  WATSON px per ACI px   "
          "ACI µm/px   WATSON µm/px")
    for d_aci_mm in (60, 80, 100, 110, 129, 150):
        d_wat_mm = d_aci_mm - np.linalg.norm(wat["C"] - aci["C"]) * 1000  # 共轴，岩面在两者之间
        # 通用写法：沿 ACI 轴从 ACI 光心走 d_aci_mm 得到岩面点，再算它到 WATSON 的距离
        P = aci["C"] + (d_aci_mm / 1000.0) * aci["A"]
        d_w = float((P - wat["C"]) @ wat["A"])
        ratio = (aix * d_aci_mm) / (wix * d_w)
        print("%8.0f %18.1f %21.3f %13.2f %13.2f"
              % (d_aci_mm, d_w, ratio,
                 aix * 1e-6 * d_aci_mm * 1000, wix * 1e-6 * d_w * 1000))

    # 和单应本身的局部尺度对一下（避免上面公式写错）
    def ray(cam, x, y):
        return np.cross(cam["H"] - x * cam["A"], cam["V"] - y * cam["A"])

    def project(cam, P):
        v = np.asarray(P, float) - cam["C"]
        return np.array([float(v @ cam["H"]) / float(v @ cam["A"]),
                         float(v @ cam["V"]) / float(v @ cam["A"])])

    def homography(stand_mm):
        d = stand_mm / 1000.0
        P0 = aci["C"] + d * aci["A"]
        n = aci["A"]
        src = np.array([[0, 0], [aci["ns"] - 1, 0], [aci["ns"] - 1, aci["nl"] - 1],
                        [0, aci["nl"] - 1]], np.float32)
        dst = []
        for x, y in src:
            dv = ray(aci, x, y)
            t = float((P0 - aci["C"]) @ n) / float(dv @ n)
            dst.append(project(wat, aci["C"] + t * dv))
        dst = np.array(dst, np.float32)
        return cv2.getPerspectiveTransform(src, dst)

    print("\n单应在中点的局部尺度（WATSON px / ACI px）：")
    for d in (60, 80, 100, 110, 129, 150):
        Hm = homography(d)
        cx, cy = (aci["ns"] - 1) / 2.0, (aci["nl"] - 1) / 2.0
        p0 = np.array([cx, cy, 1.0])
        q0 = Hm @ p0
        q0 = q0[:2] / q0[2]
        s = []
        for dx, dy in ((10, 0), (0, 10)):
            q = Hm @ np.array([cx + dx, cy + dy, 1.0])
            q = q[:2] / q[2]
            s.append(np.linalg.norm(q - q0) / np.hypot(dx, dy))
        print("  %4d mm   水平 %.3f   垂直 %.3f" % (d, s[0], s[1]))


if __name__ == "__main__":
    main()
