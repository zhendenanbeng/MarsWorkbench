# -*- coding: utf-8 -*-
"""把 WATSON 的彩色搬到 ACI 灰图上（Bellegarde，sol 185 / 186）。

这是上一轮"spike 实验"的收口脚本。走到这里经历了三次方向修正：

1. 靠文件时刻猜配套（错：sol 141 的 WATSON 比 ACI 早 12 分钟，看的是别处）
2. 靠图像特征配准（错：候选帧本身就不对，NCC 峰值顶在搜索边界上）
3. **靠两台相机各自的 CAHVOR 模型做几何**（对：坐标系都是 ROVER_NAV_FRAME，
   转塔姿态已经被模型吸收，剩下只有一个未知量——岩面离 ACI 多远）

本脚本：扫 ACI 的 standoff → 每个值算一遍单应 → 用灰度/A CI 的结构相似度
选最优 → 用 LAB 迁移颜色。颜色只取 WATSON 的 a/b（色度），亮度仍用 ACI
自己的，这样纹理一点不丢。

用法：python colorize_pair.py
输出：out/bellegarde_compare.png 等
"""

from __future__ import annotations

import os
import re
import sys

import numpy as np
import cv2

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
XMLDIR = os.path.join(OUT, "xml")
IMGDIR = os.path.join(OUT, "img")

WAT_XML = os.path.join(XMLDIR, "SIF_0185_0683368184_652FDR_N0070000SRLC00720_0000LMJ04.xml")
ACI_XML = os.path.join(XMLDIR, "SC3_0186_0683479674_054FDR_N0070000SRLC11420_0000LMJ05.xml")
WAT_IMG = os.path.join(IMGDIR, "SIF_0185_0683368184_652FDR_N0070000SRLC00720_0000LMJ04.IMG")
ACI_IMG = os.path.join(IMGDIR, "SC3_0186_0683479674_054FDR_N0070000SRLC11420_0000LMJ05.IMG")

# 扫岩面到 ACI 的距离。几何测试给出的可行区间是 100–300 mm。
STANDS_MM = [50, 60, 65, 70, 75, 80, 90, 100, 110, 120, 130, 140, 150, 160,
             180, 200, 220, 250, 300]

# 出图用的 standoff：取**几何解**，不是 NCC 最优。
# 理由：两条视轴相交在 ACI 光心前 12.9 cm 处，那就是目标所在的距离；
# 而 NCC / 互信息在这对图上是平的（见 ncc_search.py），没有分辨力，
# 用它挑 standoff 等于拿噪声做决定。
GEOM_STAND_MM = 129


def xml_text(path):
    return open(path, encoding="utf-8", errors="replace").read()


def tag_int(t, name, default=None):
    m = re.search(r"<%s>\s*(-?\d+)" % name.replace(":", r"\:"), t)
    return int(m.group(1)) if m else default


def tag_floats(t, tag, n=3):
    i = t.find("<geom:%s>" % tag)
    seg = t[i:t.index("</geom:%s>" % tag, i)]
    return np.array([float(v) for v in re.findall(r">\s*(-?[\d.eE+-]+)\s*<", seg)][:n])


def load_cam(path):
    t = xml_text(path)
    cam = {"C": tag_floats(t, "Vector_Center"), "A": tag_floats(t, "Vector_Axis"),
           "H": tag_floats(t, "Vector_Horizontal"), "V": tag_floats(t, "Vector_Vertical")}
    cam["ns"] = tag_int(t, "img:samples")
    cam["nl"] = tag_int(t, "img:lines")
    # 波段数写在 <Axis_Array><axis_name>Band</axis_name><elements>3 里；
    # 单波段产品（ACI）没有这个轴，按 1 算。
    m = re.search(r"<axis_name>Band</axis_name>\s*<elements>\s*(\d+)", t)
    cam["bands"] = int(m.group(1)) if m else 1
    # 显示方向：原始栅格未必是正的，标签里给了水平/垂直轴和方向
    cam["hax"] = re.search(r"<disp:horizontal_display_axis>\s*([^<]+)", t)
    cam["vax"] = re.search(r"<disp:vertical_display_axis>\s*([^<]+)", t)
    cam["hdir"] = re.search(r"<disp:horizontal_display_direction>\s*([^<]+)", t)
    cam["vdir"] = re.search(r"<disp:vertical_display_direction>\s*([^<]+)", t)
    cam["hax"] = cam["hax"].group(1).strip() if cam["hax"] else "?"
    cam["vax"] = cam["vax"].group(1).strip() if cam["vax"] else "?"
    cam["hdir"] = cam["hdir"].group(1).strip() if cam["hdir"] else "?"
    cam["vdir"] = cam["vdir"].group(1).strip() if cam["vdir"] else "?"
    rgb = []
    for c in ("red", "green", "blue"):
        m = re.search(r"<disp:%s_channel_band>\s*(\d+)" % c, t)
        rgb.append(int(m.group(1)) if m else None)
    cam["rgb"] = rgb
    return cam


def read_bsq(path, cam):
    """读 16 位有符号大端 BSQ 栅格。栅格在文件末尾（imgops 的排布）。"""
    nb = cam["bands"]
    need = cam["ns"] * cam["nl"] * nb * 2
    raw = open(path, "rb").read()
    buf = np.frombuffer(raw[-need:], dtype=">i2").astype(np.float32)
    return buf.reshape(nb, cam["nl"], cam["ns"])


def to8(a, lo_pct=0.5, hi_pct=99.5):
    lo, hi = np.percentile(a, (lo_pct, hi_pct))
    return np.clip((a - lo) / max(hi - lo, 1e-9) * 255.0, 0, 255).astype(np.uint8)


def feature(img):
    """去低频 + CLAHE + 梯度幅值：跨照明也能比的"结构图"。"""
    g = img.astype(np.float32)
    bg = cv2.GaussianBlur(g, (0, 0), sigmaX=max(g.shape) / 12.0)
    g = cv2.normalize(g - bg, None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
    g = cv2.createCLAHE(2.0, (8, 8)).apply(g)
    gx = cv2.Sobel(g, cv2.CV_32F, 1, 0, ksize=3)
    gy = cv2.Sobel(g, cv2.CV_32F, 0, 1, ksize=3)
    return cv2.GaussianBlur(cv2.magnitude(gx, gy), (0, 0), 0.8)


def ray(cam, x, y):
    return np.cross(cam["H"] - x * cam["A"], cam["V"] - y * cam["A"])


def project(cam, P):
    v = np.asarray(P, float) - cam["C"]
    den = float(v @ cam["A"])
    return (float(v @ cam["H"]) / den, float(v @ cam["V"]) / den)


def homography(aci, wat, stand_mm):
    """ACI 像素 → WATSON 像素的单应：把 ACI 四角打到岩面上再投进 WATSON。"""
    d = stand_mm / 1000.0
    P0 = aci["C"] + d * aci["A"]
    n = aci["A"]
    src = np.array([[0, 0], [aci["ns"] - 1, 0],
                    [aci["ns"] - 1, aci["nl"] - 1], [0, aci["nl"] - 1]], np.float32)
    dst = []
    for x, y in src:
        dirv = ray(aci, x, y)
        t = float((P0 - aci["C"]) @ n) / float(dirv @ n)
        dst.append(project(wat, aci["C"] + t * dirv))
    dst = np.array(dst, np.float32)
    return cv2.getPerspectiveTransform(src, dst), dst


def ncc(a, b, mask):
    if mask.sum() < 2000:
        return float("nan")
    av = a[mask] - a[mask].mean()
    bv = b[mask] - b[mask].mean()
    den = np.sqrt((av ** 2).sum() * (bv ** 2).sum())
    return float((av * bv).sum() / den) if den > 0 else float("nan")


# warpPerspective 的坑（实测踩过一次）：不写 WARP_INVERSE_MAP 时，OpenCV 把
# 传入的矩阵当成「源→目标」并**内部求逆**再采样。我们的单应是 ACI→WATSON，
# 而源图是 WATSON、目标网格是 ACI，所以必须加 WARP_INVERSE_MAP，
# 否则等于拿反了的映射去采样——出图看起来"投进去了"，其实完全错位。
WARP_FLAGS = cv2.INTER_LINEAR | cv2.WARP_INVERSE_MAP


def sample_watson(img, Hm, size):
    """把 WATSON 按 ACI→WATSON 的单应采样到 ACI 网格上。"""
    return cv2.warpPerspective(img, Hm, size, flags=WARP_FLAGS)


def main():
    wat = load_cam(WAT_XML)
    aci = load_cam(ACI_XML)
    print("WATSON %dx%d %d 波段 RGB 通道=%s 显示轴 %s(%s) / %s(%s)"
          % (wat["ns"], wat["nl"], wat["bands"], wat["rgb"], wat["hax"], wat["hdir"],
             wat["vax"], wat["vdir"]))
    print("ACI    %dx%d %d 波段 显示轴 %s(%s) / %s(%s)"
          % (aci["ns"], aci["nl"], aci["bands"], aci["hax"], aci["hdir"],
             aci["vax"], aci["vdir"]))

    wr = read_bsq(WAT_IMG, wat)
    ar = read_bsq(ACI_IMG, aci)
    print("栅格：WATSON %s  ACI %s" % (wr.shape, ar.shape))

    # 组装 WATSON 彩色（BSQ → HWC），通道顺序照标签
    b = [i - 1 for i in wat["rgb"]]
    wcol = np.stack([wr[b[0]], wr[b[1]], wr[b[2]]], axis=-1)      # R,G,B
    wcol8 = np.stack([to8(wcol[..., i]) for i in range(3)], axis=-1)
    aci8 = to8(ar[0])

    wf = feature(cv2.cvtColor(wcol8, cv2.COLOR_RGB2GRAY))
    af = feature(aci8)

    print("\n--- 扫 ACI standoff，看哪一档的结构最吻合 ---")
    # 两台相机的有效分辨率不一样：ACI 是 6–10 µm/px，WATSON 投过来约 11 µm/px。
    # 比对前把细的那一边糊到粗的那一边，否则 NCC 被"高频不对齐"压低。
    best = None
    null_samples = {}

    def coverage(d):
        """这一档 standoff 下，ACI 网格里能采到 WATSON 的比例（%）。"""
        Hm, _ = homography(aci, wat, d)
        m = sample_watson(np.full((wat["nl"], wat["ns"]), 255, np.uint8), Hm,
                          (aci["ns"], aci["nl"])) > 127
        return 100.0 * m.sum() / m.size

    print("   覆盖率%%   %s" % " ".join("%6.0f" % coverage(d) for d in STANDS_MM))
    for sig in (0.0, 1.0, 2.0, 3.0):
        af_b = af if sig == 0 else cv2.GaussianBlur(af, (0, 0), sig)
        for d in STANDS_MM:
            Hm, dst = homography(aci, wat, d)
            # 不要求四角全在画面内：只要求还有一块重叠可比。
            # （小 standoff 时 ACI 的视场会大于 WATSON 的画幅，但重叠部分一样能比。）
            warp = sample_watson(wf, Hm, (aci["ns"], aci["nl"]))
            m = sample_watson(np.full((wat["nl"], wat["ns"]), 255, np.uint8), Hm,
                              (aci["ns"], aci["nl"])) > 127
            m[:80, :] = False
            m[-80:, :] = False
            m[:, :80] = False
            m[:, -80:] = False
            s = ncc(af_b, warp, m)
            if np.isfinite(s) and (best is None or s > best[0]):
                best = (s, d, Hm, sig)
            null_samples.setdefault(d, []).append(s)
    for sig in (0.0, 1.0, 2.0, 3.0):
        row = []
        for d in STANDS_MM:
            Hm, dst = homography(aci, wat, d)
            af_b = af if sig == 0 else cv2.GaussianBlur(af, (0, 0), sig)
            warp = sample_watson(wf, Hm, (aci["ns"], aci["nl"]))
            m = sample_watson(np.full((wat["nl"], wat["ns"]), 255, np.uint8), Hm,
                              (aci["ns"], aci["nl"])) > 127
            m[:80, :] = False
            m[-80:, :] = False
            m[:, :80] = False
            m[:, -80:] = False
            row.append("%6.3f" % ncc(af_b, warp, m))
        print("   blur σ=%.0f  %s" % (sig, " ".join(row)))
    print("   standoff(mm) " + " ".join("%6d" % d for d in STANDS_MM))

    # 零假设对照：把同样的单应加一个平移扰动，看 NCC 掉不掉。
    # 真配准的话，扰动几十像素 NCC 应显著下降；否则说明这点相关是"通用相似"。
    if best is not None:
        _s, d0, H0, sig0 = best
        af_b = af if sig0 == 0 else cv2.GaussianBlur(af, (0, 0), sig0)
        m = np.zeros_like(aci8, bool)
        m[80:-80, 80:-80] = True
        nulls = []
        for dx, dy in ((60, 0), (-60, 0), (0, 60), (0, -60),
                       (120, 120), (-120, -120), (200, 0), (0, 200)):
            T = np.array([[1, 0, dx], [0, 1, dy], [0, 0, 1]], np.float64)
            w2 = sample_watson(wf, H0 @ np.linalg.inv(T), (aci["ns"], aci["nl"]))
            nulls.append(ncc(af_b, w2, m))
        print("\n零假设对照（把配准平移几十像素后的 NCC）：%s"
              % " ".join("%.3f" % v for v in nulls))
        print("   最好 %.3f ／ 扰动后中位 %.3f ／ 干净处 %.3f"
              % (_s, float(np.median(nulls)), min(nulls)))
    if best is None:
        print("没有得到有效的 standoff")
        return 1
    score, d_best, Hm, sig_best = best
    print("==> 最佳 standoff %d mm（ACI 比对前 blur σ=%.0f），NCC %.4f"
          % (d_best, sig_best, score))
    print("    注意：NCC 在这对图上没有分辨力，出图改用几何解 %d mm" % GEOM_STAND_MM)
    Hm, dst = homography(aci, wat, GEOM_STAND_MM)
    d_best = GEOM_STAND_MM

    warped = sample_watson(wcol8, Hm, (aci["ns"], aci["nl"]))
    valid = cv2.warpPerspective(np.full((wat["nl"], wat["ns"]), 255, np.uint8), Hm,
                                (aci["ns"], aci["nl"]),
                                flags=cv2.INTER_NEAREST | cv2.WARP_INVERSE_MAP) > 127

    # LAB 迁移：亮度用 ACI，色度用 WATSON
    lab_w = cv2.cvtColor(warped, cv2.COLOR_RGB2LAB).astype(np.float32)
    lab_a = cv2.cvtColor(cv2.cvtColor(aci8, cv2.COLOR_GRAY2RGB),
                         cv2.COLOR_RGB2LAB).astype(np.float32)
    out = lab_a.copy()
    out[..., 1] = lab_w[..., 1]
    out[..., 2] = lab_w[..., 2]
    color = cv2.cvtColor(np.clip(out, 0, 255).astype(np.uint8), cv2.COLOR_LAB2RGB)
    color[~valid] = 0

    def panel(img, text):
        bar = np.full((40, img.shape[1], 3), 255, np.uint8)
        cv2.putText(bar, text, (10, 28), cv2.FONT_HERSHEY_SIMPLEX, 0.7, (0, 0, 0), 2,
                    cv2.LINE_AA)
        return np.vstack([bar, img])

    h = 760
    def rs(im):
        s = h / im.shape[0]
        return cv2.resize(im, (int(im.shape[1] * s), h))
    comp = np.hstack([
        panel(rs(cv2.cvtColor(aci8, cv2.COLOR_GRAY2RGB)), "ACI gray (no color)"),
        panel(rs(warped), "WATSON warped (%d mm, NCC %.3f)" % (d_best, score)),
        panel(rs(color), "colorized: L=ACI, ab=WATSON"),
    ])
    cv2.imwrite(os.path.join(OUT, "bellegarde_compare.png"), cv2.cvtColor(comp, cv2.COLOR_RGB2BGR))
    cv2.imwrite(os.path.join(OUT, "bellegarde_color.png"), cv2.cvtColor(color, cv2.COLOR_RGB2BGR))
    cv2.imwrite(os.path.join(OUT, "bellegarde_aci.png"), aci8)
    cv2.imwrite(os.path.join(OUT, "bellegarde_watson.png"), cv2.cvtColor(warped, cv2.COLOR_RGB2BGR))
    # 再出一张"standoff 联络表"：把 WATSON 分别在几个距离上投进 ACI 网格，
    # 让肉眼判断哪一个才真的对得上（NCC 做不到这件事）。
    tiles = []
    for d in (60, 70, 80, 90, 100, 115, 129):
        H2, _ = homography(aci, wat, d)
        wv = sample_watson(wcol8, H2, (aci["ns"], aci["nl"]))
        tiles.append(panel(rs(wv), "%d mm" % d))
    cv2.imwrite(os.path.join(OUT, "bellegarde_standoffs.png"),
                cv2.cvtColor(np.hstack(tiles), cv2.COLOR_RGB2BGR))
    cv2.imwrite(os.path.join(OUT, "bellegarde_watson_full.png"),
                cv2.cvtColor(wcol8, cv2.COLOR_RGB2BGR))
    print("写到 out/bellegarde_compare.png 等 4 张")
    return 0


if __name__ == "__main__":
    sys.exit(main())
