# -*- coding: utf-8 -*-
"""ACI 上色可行性试验（一次性脚本，不进项目）。

思路：
  1. 读 WATSON 彩色 RAD（3 波段 BSQ 16 位）和本地 ACI 灰度 ECM；
  2. SIFT 特征 + FLANN 匹配 + RANSAC 估计「相似变换」（平移/旋转/等比缩放），
     把 WATSON 重采样到 ACI 网格——两台相机都垂直朝下装在转台上，
     相似变换比单应更符合物理，也能抑制错误匹配；
  3. 颜色迁移：两边转 LAB，取 WATSON 的 a/b 色度（按分辨率差做模糊），
     亮度 L 用 ACI 自己的——纹理全保真，颜色来自 WATSON。

用法：python colorize_spike.py
输出：aci_gray.png / watson.png / colorized.png / compare.png
"""

import sys
from pathlib import Path

import numpy as np
import cv2

# 旧开发树位置：重做实验时改成你的本地路径
DEV_ROOT = Path.home() / "Desktop" / "Sherlock"
sys.path.insert(0, str(DEV_ROOT / "01_软件工作台"))
from sherloc.vicar import _parse_label  # noqa: E402

HERE = str(DEV_ROOT / "06_算法实验与依据")
ACI = str(DEV_ROOT / "02_练习数据与验证基准" / "images" / "aci"
          / "SC3_0141_0679457520_796ECM_N0051812SRLC10203_0000LUJ04.IMG")
WAT = HERE + r"\wat_rad.IMG"


def read_img_bsq(path):
    """读 3 波段 BSQ 16 位 IMG，返回 HxWx3 uint8。"""
    raw = open(path, "rb").read()
    m = raw.find(b"END\r\n")
    if m < 0:
        m = raw.find(b"END\n")
    lab = _parse_label(raw[:m].decode("ascii", "replace"))
    lines, samps, bands = lab["LINES"], lab["LINE_SAMPLES"], lab["BANDS"]
    rec = lab["RECORD_BYTES"]
    start = (int(lab["^IMAGE"]) - 1) * rec
    dt = np.dtype(">i2") if "MSB" in lab.get("SAMPLE_TYPE", "") else np.dtype("<i2")
    arr = np.frombuffer(raw, dtype=dt, offset=start,
                        count=lines * samps * bands).reshape(bands, lines, samps)
    arr = arr.astype(np.float64)
    # 百分位拉伸到 0-255，仰赖 RAD 已做过辐射校正
    lo, hi = np.percentile(arr, (0.5, 99.5))
    arr = np.clip((arr - lo) / max(hi - lo, 1e-9) * 255.0, 0, 255)
    return np.transpose(arr, (1, 2, 0)).astype(np.uint8), lab


def main():
    wat, wlab = read_img_bsq(WAT)
    import PIL.Image as PImage
    from sherloc.vicar import load_aci
    aci = np.asarray(load_aci(ACI))
    print("ACI    :", aci.shape, aci.dtype)
    print("WATSON :", wat.shape, wat.dtype,
          "lines/samps/bits =", wlab["LINES"], wlab["LINE_SAMPLES"],
          wlab["SAMPLE_BITS"])

    # ---------- 配准 ----------
    # 两机像素尺度差约 3 倍（WATSON ~31.85 um/px vs ACI 10.1 um/px），
    # 直接跨尺度 SIFT 匹配质量差。先把 WATSON 预放大到接近 ACI 的尺度，
    # 匹配完再把预放大折回变换矩阵里。
    WATSON_UM, ACI_UM = 31.85, 10.1
    f = WATSON_UM / ACI_UM
    wat_up = cv2.resize(wat, None, fx=f, fy=f, interpolation=cv2.INTER_CUBIC)

    clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8, 8))
    g_aci = clahe.apply(aci)
    g_wat = clahe.apply(cv2.cvtColor(wat_up, cv2.COLOR_BGR2GRAY))

    sift = cv2.SIFT_create(nfeatures=20000, contrastThreshold=0.01)
    k1, d1 = sift.detectAndCompute(g_aci, None)
    k2, d2 = sift.detectAndCompute(g_wat, None)
    print("keypoints: aci %d, watson_up %d" % (len(k1), len(k2)))

    flann = cv2.FlannBasedMatcher(dict(algorithm=1, trees=5), dict(checks=64))
    matches = flann.knnMatch(d2, d1, k=2)          # watson -> aci
    good = [m for m, n in matches if m.distance < 0.8 * n.distance]
    print("good matches:", len(good))
    if len(good) < 20:
        print("!! 匹配太少，放弃")
        return 1

    src = np.float64([k2[m.queryIdx].pt for m in good]).reshape(-1, 1, 2)
    dst = np.float64([k1[m.trainIdx].pt for m in good]).reshape(-1, 1, 2)
    M, inl = cv2.estimateAffinePartial2D(
        src, dst, method=cv2.RANSAC, ransacReprojThreshold=3.0,
        maxIters=20000, confidence=0.999)
    if M is None or inl is None:
        print("!! RANSAC 失败")
        return 1
    n_inl = int(inl.sum())
    scale = np.hypot(M[0, 0], M[1, 0])
    mask = inl.ravel().astype(bool)
    proj = cv2.transform(src, M).reshape(-1, 2)[mask]
    ref = dst.reshape(-1, 2)[mask]
    err = np.hypot(*(proj - ref).T)
    print("inliers %d/%d  scale=%.4f (expect ~%.2f)  median err=%.2f px"
          % (n_inl, len(good), scale, f, float(np.median(err))))
    if n_inl < 15 or not (0.5 < scale < 2.0):
        print("!! 配准置信度不足，放弃")
        return 1
    # 折回预放大：upsampled 坐标 = 原图坐标 * f
    M = M.copy()
    M[:, :2] *= f
    scale *= f

    h, w = aci.shape
    warped = cv2.warpAffine(wat, M, (w, h), flags=cv2.INTER_LINEAR,
                            borderMode=cv2.BORDER_CONSTANT,
                            borderValue=(0, 0, 0))
    valid = cv2.warpAffine(np.full(wat.shape[:2], 255, np.uint8),
                           M, (w, h), flags=cv2.INTER_NEAREST) > 127

    # ---------- 颜色迁移（LAB）----------
    wat_f = cv2.cvtColor(warped, cv2.COLOR_BGR2LAB).astype(np.float32)
    aci_l = cv2.cvtColor(aci, cv2.COLOR_GRAY2BGR)
    aci_lab = cv2.cvtColor(aci_l, cv2.COLOR_BGR2LAB).astype(np.float32)

    # WATSON 像素比 ACI 粗，色度先模糊，避免色块锯齿
    k = max(int(round(scale * 3)) | 1, 3)
    a = cv2.GaussianBlur(wat_f[..., 1], (k, k), 0)
    b = cv2.GaussianBlur(wat_f[..., 2], (k, k), 0)

    out = aci_lab.copy()
    out[..., 1] = a
    out[..., 2] = b
    out = np.clip(out, 0, 255).astype(np.uint8)
    colorized = cv2.cvtColor(out, cv2.COLOR_LAB2BGR)
    colorized[~valid] = 0

    # ---------- 输出 ----------
    PImage.fromarray(aci).save(HERE + r"\aci_gray.png")
    PImage.fromarray(cv2.cvtColor(wat, cv2.COLOR_BGR2RGB)).save(HERE + r"\watson.png")
    PImage.fromarray(cv2.cvtColor(colorized, cv2.COLOR_BGR2RGB)).save(
        HERE + r"\colorized.png")

    def label(img, text):
        bar = np.full((36, img.shape[1], 3), 255, np.uint8)
        cv2.putText(bar, text, (10, 26), cv2.FONT_HERSHEY_SIMPLEX, 0.8,
                    (0, 0, 0), 2, cv2.LINE_AA)
        return np.vstack([bar, img])

    cmp_h = 700
    def rs(im):
        s = cmp_h / im.shape[0]
        return cv2.resize(im, (int(im.shape[1] * s), cmp_h))
    panel = np.hstack([
        label(rs(cv2.cvtColor(aci, cv2.COLOR_GRAY2BGR)), "ACI gray"),
        label(rs(warped), "WATSON warped (SIFT+RANSAC, %d inl)" % n_inl),
        label(rs(colorized), "colorized: L=ACI, ab=WATSON"),
    ])
    PImage.fromarray(panel).save(HERE + r"\compare.png")
    print("done -> compare.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
