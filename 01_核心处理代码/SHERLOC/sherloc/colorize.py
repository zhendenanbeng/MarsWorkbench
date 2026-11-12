# -*- coding: utf-8 -*-
"""WATSON 彩色 → ACI 灰度的配准与上色（spike 验证过的算法，产品化封装）。

输入一对图像（ACI 灰度 + WATSON 彩色）和各自的 CAHVOR 相机模型，
输出上色图和一份质量报告。管线：

  几何单应当初值（standoff 由调用方给，通常来自 watson.find_pair）
  → ACI 亮斑簇锚点（程序化检测）
  → 缩放×模糊联合搜索（显示方向一致时不需要翻转）
  → 刚性拟合残差 + 锚点峰噪 + 受光率 + 旋转角四道闸门自动判 GO/NO-GO
  → GO：LAB 迁移（L=ACI 保纹理，ab=WATSON 给颜色），阴影区中性

判 NO-GO 时返回 reason，调用方回退灰度底图即可——用户不需要知道细节。

依赖：opencv（可选依赖，import 失败时上层应隐藏上色功能）。
算法的来龙去脉见 06_算法实验与依据/README.md 第七、八节。
"""

from __future__ import annotations

import numpy as np

try:
    import cv2
    HAS_CV2 = True
except ImportError:                                  # noqa: F401
    HAS_CV2 = False

# 搜索网格：μ 覆盖 0.40–2.40。曾收窄到 0.60–1.45（按当时的两个 GO 案例），
# sol 865 的宽幅 survey（SRLC10000）实测最优 μ=1.85 顶穿网格——不同类型的
# 扫描（HDR/Survey/宽幅）几何差异比想象大，别按样本收窄。σ 覆盖两机
# 5.4 倍的角分辨率差需要 ≥3。
MUS = np.round(np.arange(0.40, 2.41, 0.08), 2)
SIGMAS = (0.0, 3.0, 5.0)
TPL = 280                 # 锚点模板边长（ACI px）
SEARCH_PX = 300           # 模板匹配邻域
MIN_LIT = 0.10            # 投影区受光率下限
MAX_RESID = 90            # 刚性残差上限（跨 sol keystone 吃不满是物理必然）
MAX_ROT = 3.0             # 刚性旋转角上限（度）——转塔相机取向差只有几度
MIN_PEAKNOISE = 0.20      # 锚点响应面峰噪中位下限（真匹配实测 0.27–0.50）


def _stretch(a, lo_pct=0.5, hi_pct=99.5):
    """百分位拉伸到 0–255 uint8。"""
    a = a.astype(np.float32)
    lo, hi = np.percentile(a, (lo_pct, hi_pct))
    return np.clip((a - lo) / max(hi - lo, 1e-9) * 255.0, 0, 255).astype(np.uint8)


def feature(img):
    """去低频 + CLAHE + 梯度幅值：跨照明也能比的结构图。"""
    g = img.astype(np.float32)
    bg = cv2.GaussianBlur(g, (0, 0), sigmaX=max(g.shape) / 12.0)
    g = cv2.normalize(g - bg, None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
    g = cv2.createCLAHE(2.0, (8, 8)).apply(g)
    gx = cv2.Sobel(g, cv2.CV_32F, 1, 0, ksize=3)
    gy = cv2.Sobel(g, cv2.CV_32F, 0, 1, ksize=3)
    return cv2.GaussianBlur(cv2.magnitude(gx, gy), (0, 0), 0.8)


def detect_anchors(aci8, want=6):
    """程序化找亮斑簇：局部极大 → 网格聚类 → 对比度×数量最高的几簇。"""
    g = cv2.GaussianBlur(aci8.astype(np.float32), (0, 0), 1.2)
    dil = cv2.dilate(g, np.ones((25, 25), np.float32))
    peaks = (g == dil) & (g > g.mean() + 2.2 * g.std())
    ys, xs = np.nonzero(peaks)
    if len(xs) == 0:
        return []
    bins = {}
    for x, y in zip(xs, ys):
        bins.setdefault((x // 100, y // 100), []).append((x, y))
    clusters = []
    for (_bx, _by), pts in bins.items():
        xs_c = np.array([p[0] for p in pts])
        ys_c = np.array([p[1] for p in pts])
        clusters.append((float(xs_c.mean()), float(ys_c.mean()),
                         float(np.mean(g[ys_c, xs_c])) * len(pts)))
    clusters.sort(key=lambda c: -c[2])
    picked = []
    for cx, cy, sc in clusters:
        if all(np.hypot(cx - px, cy - py) > 260 for px, py, _s in picked):
            picked.append((cx, cy, sc))
        if len(picked) >= want:
            break
    return picked


def zoom_matrix(mu, ns, nl):
    c = np.array([ns / 2.0, nl / 2.0])
    return np.float64([[mu, 0, (1 - mu) * c[0]],
                       [0, mu, (1 - mu) * c[1]],
                       [0, 0, 1]])


def warp(img, M, size):
    """把 WATSON 侧图像按 ACI→WATSON 映射采样到 ACI 网格。

    坑（spike 实测踩过）：我们的 M 是 ACI→WATSON 方向、源图是 WATSON，
    必须加 WARP_INVERSE_MAP，否则等于拿反了的映射去采样。
    """
    return cv2.warpPerspective(img, M, size,
                               flags=cv2.INTER_LINEAR | cv2.WARP_INVERSE_MAP)


def ncc(a, b, mask):
    if mask.sum() < 2000:
        return float("nan")
    av = a[mask] - a[mask].mean()
    bv = b[mask] - b[mask].mean()
    den = np.sqrt((av ** 2).sum() * (bv ** 2).sum())
    return float((av * bv).sum() / den) if den > 0 else float("nan")


def rigid_fit(src, dst):
    """Kabsch 刚性拟合（旋转+平移，无缩放——尺度由 μ 搜索负责）。"""
    src_c = src - src.mean(0)
    dst_c = dst - dst.mean(0)
    H = src_c.T @ dst_c
    U, _s, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1] *= -1
        R = Vt.T @ U.T
    t = dst.mean(0) - R @ src.mean(0)
    S = np.vstack([np.hstack([R, t.reshape(2, 1)]), [0, 0, 1.0]])
    ang = float(np.degrees(np.arctan2(R[1, 0], R[0, 0])))
    return S, ang


def rigid_residual(src, dst):
    S, ang = rigid_fit(src, dst)
    proj = src @ S[:2, :2].T + S[:2, 2]
    return float(np.hypot(*(proj - dst).T).mean()), ang


def homography_from_cahvor(aci_cam, wat_cam, stand_mm):
    """ACI 像素 → WATSON 像素的单应（把 ACI 四角打到岩面距离再投进 WATSON）。

    aci_cam/wat_cam 是 dict：C/A/H/V 各 3 元列表 + ns/nl + bands。
    """
    d = stand_mm / 1000.0
    n = np.asarray(aci_cam["A"], float)
    C = np.asarray(aci_cam["C"], float)
    P0 = C + d * n
    src = np.array([[0, 0], [aci_cam["ns"] - 1, 0],
                    [aci_cam["ns"] - 1, aci_cam["nl"] - 1],
                    [0, aci_cam["nl"] - 1]], np.float32)
    dst = []
    for x, y in src:
        dirv = np.cross(np.asarray(aci_cam["H"]) - x * n,
                        np.asarray(aci_cam["V"]) - y * n)
        t = float((P0 - C) @ n) / float(dirv @ n)
        v = (C + t * dirv) - np.asarray(wat_cam["C"], float)
        den = float(v @ np.asarray(wat_cam["A"], float))
        dst.append((float(v @ np.asarray(wat_cam["H"], float)) / den,
                    float(v @ np.asarray(wat_cam["V"], float)) / den))
    dst = np.array(dst, np.float32)
    return cv2.getPerspectiveTransform(src, dst)


def colorize(aci_gray, wat_rgb, aci_cam, wat_cam, stand_mm, log=None):
    """主入口。返回 dict：

    GO:  {"ok": True, "color": HxWx3, "scale": 尺度比, "resid": px,
          "rotation": 度, "n_anchors": n, "peak_noise": 中位峰噪,
          "lit": 受光率, "ncc": 全图相关（参考值）}
    NO-GO: {"ok": False, "reason": 文案}
    """
    def say(msg):
        if log:
            log(msg)

    if not HAS_CV2:
        return {"ok": False, "reason": "未安装 opencv（pip install "
                                       "opencv-python-headless）"}
    ns, nl = aci_gray.shape[1], aci_gray.shape[0]
    wgray = cv2.cvtColor(wat_rgb, cv2.COLOR_RGB2GRAY)
    # 亮度规范化（百分位拉伸）。必须做：不同读图器对 16 位 VICAR 的
    # 归一差异极大（load_aci 中位 124 vs spike to8 的 62，实测差 81），
    # 不规范的话特征图结构跟着变，锚点匹配会整体走样。
    aci8 = _stretch(aci_gray)
    wgray = _stretch(wgray)

    H0 = homography_from_cahvor(aci_cam, wat_cam, stand_mm)
    say("几何初值 standoff %.0f mm" % stand_mm)

    anchors = detect_anchors(aci8)
    if len(anchors) < 3:
        return {"ok": False, "reason": "ACI 上可作锚点的亮斑簇不足"}
    say("锚点 %d 个" % len(anchors))

    af, wf = feature(aci8), feature(wgray)
    boxes = []
    for cx, cy, _s in anchors:
        x0 = int(np.clip(cx - TPL / 2, 0, ns - TPL))
        y0 = int(np.clip(cy - TPL / 2, 0, nl - TPL))
        boxes.append((x0, y0))
    src_pts = np.float32([[a[0], a[1]] for a in anchors])

    Zs = {mu: zoom_matrix(mu, ns, nl) for mu in MUS}
    # 双边低通：模板和投影图都按 σ 模糊，匹配才真正在同一频段。
    # 只糊模板的旧做法在 μ 偏离 1 较多时对不齐（sol 865 实测）。
    afs = {s: (af if s == 0 else cv2.GaussianBlur(af, (0, 0), s))
           for s in SIGMAS}
    tpls = {s: [afs[s][y0:y0 + TPL, x0:x0 + TPL] for x0, y0 in boxes]
            for s in SIGMAS}
    lit_raw0 = (wgray > np.percentile(wgray, 25)).astype(np.uint8) * 255

    combos = []
    for mu in MUS:
        Mz = H0 @ Zs[mu]
        wfl_raw = sample_level(wf, Mz, (ns, nl))
        cover = warp(np.full_like(wgray, 255), Mz, (ns, nl)) > 127
        lit_frac = float((cover & (warp(lit_raw0, Mz, (ns, nl)) > 127)).sum()) \
            / cover.size
        for s in SIGMAS:
            wfl = wfl_raw if s == 0 else cv2.GaussianBlur(wfl_raw, (0, 0), s)
            row = []
            for (x0, y0), tpl in zip(boxes, tpls[s]):
                sx0, sy0 = max(0, x0 - SEARCH_PX), max(0, y0 - SEARCH_PX)
                sx1 = min(ns, x0 + TPL + SEARCH_PX)
                sy1 = min(nl, y0 + TPL + SEARCH_PX)
                if sx1 - sx0 < TPL or sy1 - sy0 < TPL:
                    row.append((0.0, 0, 0, 0.0))
                    continue
                res = cv2.matchTemplate(wfl[sy0:sy1, sx0:sx1], tpl,
                                        cv2.TM_CCOEFF_NORMED)
                _, mx, _, loc = cv2.minMaxLoc(res)
                row.append((float(mx), (loc[0] + sx0) - x0,
                            (loc[1] + sy0) - y0,
                            float(mx - np.median(res))))
            scores = np.array([r[0] for r in row])
            offs = np.array([[r[1], r[2]] for r in row])
            nok = int((scores > 0.18).sum())
            resid, ang = 9999.0, 0.0
            ok = scores > 0.18
            if ok.sum() >= 3:
                resid, ang = rigid_residual(src_pts[ok], src_pts[ok] + offs[ok])
            med = float(np.median(scores)) if nok else 0.0
            pns = np.array([r[3] for r in row])
            pn_med = float(np.median(pns[ok])) if nok else 0.0
            if lit_frac < MIN_LIT:
                nok, med, resid, pn_med = 0, 0.0, 9999.0, 0.0
            crit = med - 0.006 * resid if nok >= 3 else -9
            combos.append((crit, med, resid, nok, float(mu), s, row,
                           offs, scores, lit_frac, pn_med, ang))

    if not combos:
        return {"ok": False, "reason": "搜索无结果"}
    combos.sort(key=lambda c: -c[0])
    (crit, med, resid, nok, mu, sig, row, offs, scores,
     lit_frac, pn_med, ang) = combos[0]
    say("最优：μ=%.2f σ=%.1f（%d 锚点，残差 %.0fpx，旋转 %.1f°）"
        % (mu, sig, nok, resid, ang))

    if crit <= 0:
        if nok < 3:
            return {"ok": False, "reason":
                    "锚点匹配不成（仅 %d 个过线，WATSON 里找不到对应纹理）"
                    % nok}
        return {"ok": False, "reason":
                "锚点得分过低（最佳 %d 锚点、中位 %.2f、残差 %.0fpx，"
                "未达置信标准——两次拍摄可能不是同一小区域）"
                % (nok, med, resid)}
    if resid >= MAX_RESID:
        return {"ok": False, "reason": "锚点残差过大（%.0f px ≥ %d）"
                                       % (resid, MAX_RESID)}
    if abs(ang) > MAX_ROT:
        return {"ok": False, "reason": "需要的旋转 %.1f° 超出物理范围"
                                       % ang}
    if pn_med < MIN_PEAKNOISE:
        return {"ok": False, "reason": "锚点响应不尖（峰噪 %.2f）" % pn_med}

    # 最终映射：几何 × 缩放 × 刚性校正
    okm = scores > 0.18
    S3, _ = rigid_fit(src_pts[okm], src_pts[okm] + offs[okm])
    M = H0 @ Zs[mu] @ S3

    # 全图 NCC（参考指标）：阴影掩模自适应腐蚀
    cover = warp(np.full_like(wgray, 255), M, (ns, nl)) > 127
    mask = None
    for er in (31, 21, 13, 9, 5):
        lit = cv2.erode(lit_raw0, np.ones((er, er), np.uint8))
        m2 = cover & (warp(lit, M, (ns, nl)) > 127)
        m2[40:-40, 40:-40] &= True
        mask = m2
        if int(m2.sum()) >= 300000:
            break
    n_after = ncc(afs[sig], sample_level(wf, M, (ns, nl)), mask)

    # 上色：L=ACI，ab=WATSON（按尺度比定核模糊），阴影/画幅外中性
    warped = warp(wat_rgb, M, (ns, nl))
    colmask = cover & (warp(lit_raw0, M, (ns, nl)) > 127)
    lab_w = cv2.cvtColor(warped, cv2.COLOR_RGB2LAB).astype(np.float32)
    lab_a = cv2.cvtColor(cv2.cvtColor(aci8, cv2.COLOR_GRAY2RGB),
                         cv2.COLOR_RGB2LAB).astype(np.float32)
    r_emp = _edge_ratio(M, ns, nl)
    k = max(int(round(r_emp * 3)) | 1, 3)
    ab = cv2.GaussianBlur(lab_w[..., 1:], (k, k), 0)
    out = lab_a.copy()
    out[..., 1:] = ab
    out[~colmask, 1:] = 128.0
    color = cv2.cvtColor(np.clip(out, 0, 255).astype(np.uint8),
                         cv2.COLOR_LAB2RGB)

    say("GO：全图 NCC %.3f（参考）" % n_after)
    return {"ok": True, "color": color, "scale": r_emp, "resid": resid,
            "rotation": ang, "n_anchors": nok, "peak_noise": pn_med,
            "lit": lit_frac, "ncc": n_after}


def sample_level(wf, M, size):
    return warp(wf, M, size)


def _edge_ratio(M, ns, nl):
    corners = np.float32([[0, 0], [ns - 1, 0], [ns - 1, nl - 1], [0, nl - 1]])
    proj = cv2.perspectiveTransform(corners.reshape(-1, 1, 2),
                                    np.float64(M)).reshape(-1, 2)
    d_w = sum(np.hypot(*(proj[(i + 1) % 4] - proj[i])) for i in range(4))
    d_a = sum(np.hypot(*(corners[(i + 1) % 4] - corners[i])) for i in range(4))
    return d_w / d_a
