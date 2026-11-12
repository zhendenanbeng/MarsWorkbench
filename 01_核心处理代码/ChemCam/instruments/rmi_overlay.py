# -*- coding: utf-8 -*-
"""ChemCam RMI 元素分布图：CAHV 共帧投影（第六期）。

路线（S1 侦察结论）：RDR 的 cr0 RMI label 不带相机模型，但 EDR 包的同 SCLK
label（cr0_<SCLK>edr_<obsid>m1.lbl）逐帧自带 CAHV，且官方 catalog 佐证
「RDR 坐标与输入 EDR 相同」。LIBS 每点自己的 CAHV 在 cl5 RDR label 里。

投影：光斑 3D 点 P = C_i + t·A_i（cl5 相机光心沿光轴），经 watson.project
投进 RMI 相机得像素。t（距离）对共轴近距系统不敏感——代码里显式验证并写日志。

诚实口径：灼蚀斑人工目检通过前，渲染图角标必须带「未做地面核验」；
无比例尺（无距离标定）。
"""

from __future__ import annotations

import os
import re
import appenv

import numpy as np

from . import pds3, pds_listing

EDR_ROOT = ("https://pds-geosciences.wustl.edu/msl/"
            "msl-m-chemcam-libs-2-edr-v1/mslccm_0xxx/data")
LABEL_CACHE = str(appenv.copy_cached('edr_labels'))

RMI_RX = re.compile(r"^cr(\d)_(\d{9})(?:edr|prc)_(f\d{7}ccam\d{5})",
                    re.I)


def parse_cahv_lbl(path, ns=None, nl=None):
    """PDS3 label 的 GEOMETRIC_CAMERA_MODEL_PARMS → watson.project 兼容 dict。

    MODEL_COMPONENT_1..4 = C/A/H/V；ns/nl 从 IMAGE 组读，缺省 1024（RMI 标称）。
    没有相机模型的 label（RDR cr0）返回 None。
    """
    lab = pds3.parse_lbl(path)
    c = lab.get("MODEL_COMPONENT_1")
    if not c or len(c) < 3:
        return None
    try:
        C = [float(v) for v in c[:3]]
        A = [float(v) for v in lab["MODEL_COMPONENT_2"][:3]]
        H = [float(v) for v in lab["MODEL_COMPONENT_3"][:3]]
        V = [float(v) for v in lab["MODEL_COMPONENT_4"][:3]]
    except (KeyError, ValueError, IndexError):
        return None
    ns = ns or int(float(lab.get("LINE_SAMPLES", 1024) or 1024))
    nl = nl or int(float(lab.get("LINES", 1024) or 1024))
    return {"C": C, "A": A, "H": H, "V": V, "ns": ns, "nl": nl,
            "source": os.path.basename(path)}


def edr_label_name(rmi_rdr_stem):
    """RDR RMI 主名（cr0_<sclk>prc_<obsid>l1）→ EDR label 名（…edr_<obsid>m1.lbl）。

    cr 后的通道位与下划线必须原样保留（cr0_486263086edr_…，不是
    cr486263086edr_…）——PDS 文件名逐字符敏感。
    """
    m = RMI_RX.match(rmi_rdr_stem)
    if not m:
        return None
    return "cr%s_%sedr_%sm1.lbl" % (m.group(1), m.group(2), m.group(3))


def ensure_edr_label(rmi_rdr_stem, sol, download=True):
    """按需取 EDR label（~32KB），缓存 .cache/edr_labels/solXXXXXX/。

    返回本地路径；download=False 或下载失败返回 None。
    """
    name = edr_label_name(rmi_rdr_stem)
    if not name:
        return None
    sol_dir = os.path.join(LABEL_CACHE, "sol%05d" % int(sol))
    path = os.path.join(sol_dir, name)
    if os.path.exists(path) and os.path.getsize(path) > 5000:
        return path
    if not download:
        return None
    try:
        os.makedirs(sol_dir, exist_ok=True)
        data = pds_listing.http_get("%s/sol%05d/%s"
                                    % (EDR_ROOT, int(sol), name), timeout=60)
        if len(data) < 5000:
            return None
        tmp = path + ".part"
        with open(tmp, "wb") as f:
            f.write(data)
        os.replace(tmp, path)
        return path
    except Exception:                                    # noqa: BLE001
        return None


def libspot_cam(cl5_lbl_path):
    """cl5 RDR label → 该点 LIBS 相机模型（C/A/H/V）。"""
    return parse_cahv_lbl(cl5_lbl_path)


def load_edr_raster(edr_lbl_path):
    """EDR RMI label → 12-bit 原始栅格（numpy 1024×1024）。

    栅格偏移从该 label 自己的 ^IMAGE 指针记录号解出（逐帧不同：
    sol01000 前帧 73、后帧 70），并对照 label 的 MEDIAN 自证没读错位。
    """
    import re as _re
    raw_lbl = open(edr_lbl_path, encoding="latin-1").read()
    m = _re.search(r"\^IMAGE\s*=\s*\([^,]+,\s*(\d+)\s*\)", raw_lbl)
    if not m:
        raise ValueError("label 无 ^IMAGE 指针：%s" % edr_lbl_path)
    lab = pds3.parse_lbl(edr_lbl_path)
    rb = int(float(lab.get("RECORD_BYTES", 2048)))
    nl = int(float(lab.get("LINES", 1024)))
    ns = int(float(lab.get("LINE_SAMPLES", 1024)))
    img_path = os.path.splitext(edr_lbl_path)[0] + ".img"
    raw = open(img_path, "rb").read()
    off = (int(m.group(1)) - 1) * rb
    raster = np.frombuffer(raw[off:off + nl * ns * 2], dtype=">u2").reshape(nl, ns)
    median_lbl = lab.get("MEDIAN")
    if median_lbl is not None and abs(float(np.median(raster))
                                      - float(median_lbl)) > 1.0:
        raise ValueError("栅格中位数 %s 与 label MEDIAN %s 不符（偏移错读）"
                         % (float(np.median(raster)), median_lbl))
    return raster


def principal_point(cam):
    """CAHV 主点：A·H / A·V（|A|=1 时即帧坐标主点）。

    ChemCam EDR 模型的 h/v 以主点为原点（sol01000 实测 h₀=7646、v₀=-10564，
    而焦距分量 ≈48k px 与 RMI 1.158°/1024px 自洽）——必须减主点、加帧中心
    才是帧像素。锚验证：修正后首点 (510,510)，与 label SPECIAL (line 497,
    sample 532) 差 ~20px，即 LIBS↔RMI 真实视轴失配，非投影错误。
    """
    a = float(np.dot(cam["A"], cam["A"])) or 1.0
    return float(np.dot(cam["A"], cam["H"])) / a, float(np.dot(cam["A"], cam["V"])) / a


def project_spot(lib_cam, rmi_cam, t_m=2.0):
    """P = C + t·A 投进 RMI，返回帧像素 (line, sample)；射线背向相机返回 None。

    t 为 standoff（米）。投影含主点重定心（见 principal_point）。
    """
    from sherloc import watson
    P = [lib_cam["C"][k] + t_m * lib_cam["A"][k] for k in range(3)]
    hv = watson.project(rmi_cam, P)
    if hv is None:
        return None
    h0, v0 = principal_point(rmi_cam)
    sample = hv[0] - h0 + (rmi_cam["ns"] - 1) / 2.0
    line = hv[1] - v0 + (rmi_cam["nl"] - 1) / 2.0
    return (line, sample)


def project_all(points, rmi_cam, t_m=2.0):
    """points: [(点号, lib_cam)] → {点号: (line, sample)}。"""
    out = {}
    for no, cam in points:
        px = project_spot(cam, rmi_cam, t_m)
        if px is not None:
            out[no] = px
    return out


# ---- 稳定性自检（护栏用，全部返回 (ok, detail)） ---------------------------

def check_in_frame(px, ns=1024, nl=1024):
    bad = {no: p for no, p in px.items()
           if not (2 <= p[0] < nl - 2 and 2 <= p[1] < ns - 2)}
    return (not bad, "出帧: %s" % bad if bad else "%d 点全部在帧内" % len(px))


def check_spread(px, min_extent_px=100.0):
    """点集外接框对角线——共轴系统下光栅铺开应达百像素量级（A 轴证据 ~318px）。"""
    if len(px) < 2:
        return (False, "点数不足")
    ls = np.array([p[0] for p in px.values()])
    ss = np.array([p[1] for p in px.values()])
    extent = float(np.hypot(ls.max() - ls.min(), ss.max() - ss.min()))
    return (extent >= min_extent_px, "外接对角 %.0f px（阈值 %.0f）"
            % (extent, min_extent_px))


def check_t_insensitive(points, rmi_cam, ts_m=(1.5, 2.0, 2.5), max_drift_px=30.0):
    """投影对距离 t 的敏感度（共轴系统应远小于光栅铺开）。"""
    sets = [project_all(points, rmi_cam, t) for t in ts_m]
    keys = set(sets[0]) & set(sets[-1])
    if not keys:
        return (False, "无共同点")
    drift = max(np.hypot(sets[0][k][0] - sets[-1][k][0],
                         sets[0][k][1] - sets[-1][k][1]) for k in keys)
    return (drift <= max_drift_px, "t=%.1f→%.1fm 最大漂移 %.2f px"
            % (ts_m[0], ts_m[-1], drift))


def check_pair_consistency(points, cam_a, cam_b, min_corr=0.95):
    """双 RMI 交叉一致：两张 RMI 各自投影，点间距矩阵的相关性。"""
    pa = project_all(points, cam_a)
    pb = project_all(points, cam_b)
    keys = sorted(set(pa) & set(pb))
    if len(keys) < 3:
        return (False, "共同点不足 3")
    ka = np.array([[pa[k][0], pa[k][1]] for k in keys])
    kb = np.array([[pb[k][0], pb[k][1]] for k in keys])
    # 相对几何：两两距离矩阵的相关（尺度可能略不同，形状应一致）
    da, db = [], []
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            da.append(float(np.hypot(*(ka[i] - ka[j]))))
            db.append(float(np.hypot(*(kb[i] - kb[j]))))
    if len(da) < 2 or np.std(da) < 1e-9 or np.std(db) < 1e-9:
        return (False, "点距矩阵无方差")
    corr = float(np.corrcoef(da, db)[0, 1])
    return (corr >= min_corr, "点距相关性 r=%.3f（阈值 %.2f）" % (corr, min_corr))


# ---- 叠加渲染 -------------------------------------------------------------

MAP_NOTE = ("CAHV 来自 EDR label，主点重定心；未做地面核验，无比例尺；"
            "元素标签为快速判定口径（非定量）")


def _stretch_pil(array, lo=0.3, hi=99.7):
    """numpy 栅格 → 百分位拉伸的 RGB PIL 图（16-bit 直转 L 会全白，见 rmi.preview）。"""
    a = np.asarray(array, dtype=float)
    low, high = np.percentile(a[np.isfinite(a)], [lo, hi])
    disp = np.zeros(a.shape, dtype=np.uint8)
    if high > low:
        disp = np.clip((a - low) / (high - low) * 255, 0, 255).astype(np.uint8)
    from PIL import Image
    return Image.fromarray(disp).convert("RGB")


def _compose(base_img, px, labels, cap, out_path, cap2=MAP_NOTE, highlight=None):
    """共享绘制：底图 + 圈 + 点号 + 元素标签 + 角标。

    highlight=元素符号时进入高亮模式：含该元素的点正常橙圈，其余点灰细圈
    不带标签——「哪些点亮」一眼可见，未点亮只是弱化不是删除。
    """
    from PIL import Image, ImageDraw
    from sherloc.render import font

    d = ImageDraw.Draw(base_img)
    items = sorted(px.items(), key=lambda kv: int(kv[0]))
    r = 12.0
    f15 = font(15)
    n_hit = 0
    for no, (line, sample) in items:
        x, y = sample, line
        elem = (labels or {}).get(no) or ""
        hit = (highlight is None) or (highlight in [e.strip() for e in elem.split(",") if e.strip()])
        if highlight is not None and hit:
            n_hit += 1
        color = (255, 106, 0) if hit else (110, 116, 124)
        width = 3 if hit else 2
        d.ellipse([x - r, y - r, x + r, y + r], outline=color, width=width)
        d.ellipse([x - 2, y - 2, x + 2, y + 2], fill=color)
        if highlight is not None and not hit:
            continue                                   # 未命中点不画标签，只留位置
        text = str(no)
        if elem:
            text += "·" + elem
        # 标签画在圆圈正下方居中：光栅列距 ~110px，右伸会挤到邻列圆圈
        w = d.textlength(text, font=f15)
        d.text((x - w / 2, y + r + 4), text,
               fill=(255, 220, 80) if hit else (150, 158, 168), font=f15)
    if highlight is not None:
        cap += " · 高亮 %s：%d/%d 点检出" % (highlight, n_hit, len(items))
    bar = 34
    img2 = Image.new("RGB", (1024, 1024 + bar * 2), (12, 14, 16))
    img2.paste(base_img, (0, bar))
    d2 = ImageDraw.Draw(img2)
    d2.text((10, 7), cap, fill=(235, 238, 242), font=font(15))
    d2.text((10, 1024 + bar + 5), cap2, fill=(150, 158, 168), font=font(14))
    if isinstance(out_path, str):
        img2.save(out_path)
    else:                                  # BytesIO 等内存目标（高亮重绘）
        img2.save(out_path, format="PNG")
    return out_path


def base_image(edr_lbl_path, fallback_tif):
    """叠加底图：优先 EDR 12-bit 原片（百分位拉伸），退回本地 PRC。

    返回 (PIL RGB, basis)；两者都缺返回 (None, None)。
    """
    # download=False 时 sol 不参与（只查本地缓存存在性）
    img_path = ensure_edr_image(edr_lbl_path, 0, download=False)
    if img_path:
        return _stretch_pil(load_edr_raster(edr_lbl_path)), "EDR 12-bit 原片（百分位拉伸）"
    if fallback_tif and os.path.exists(fallback_tif):
        from sherloc.vicar import load_aci
        return (_stretch_pil(load_aci(fallback_tif).resize((1024, 1024)), 1, 99),
                "本地 PRC 拉伸片（EDR 原片未取到）")
    return None, None


def render_overlay(rmi_tif_path, px, out_path, t_m=2.0, notes=None,
                   title="ChemCam RMI · CAHV 共帧投影"):
    """RMI 灰度图 + 光斑投影点。px: {点号: (line, sample)}。

    无比例尺（无距离标定）；角标固定带「未做地面核验」——人工目检通过前
    这是硬口径。notes 附加到图注。
    """
    from sherloc.vicar import load_aci
    base = _stretch_pil(load_aci(rmi_tif_path).resize((1024, 1024)))
    cap = "%s · %d 点 · standoff %.1fm（t±0.5m 定位不确定 ~±7px）" % (
        title, len(px), t_m)
    if notes:
        cap += " · " + notes
    return _compose(base, px, {}, cap, out_path)


def ensure_edr_image(edr_lbl_path, sol, download=True, min_bytes=2_200_000):
    """按需取 EDR 原片（~2.2MB），与 label 同目录缓存。失败返回 None。"""
    img_path = os.path.splitext(edr_lbl_path)[0] + ".img"
    if os.path.exists(img_path) and os.path.getsize(img_path) >= min_bytes:
        return img_path
    if not download:
        return None
    try:
        data = pds_listing.http_get("%s/sol%05d/%s"
                                    % (EDR_ROOT, int(sol), os.path.basename(img_path)),
                                    timeout=180)
        if len(data) < min_bytes:
            return None
        tmp = img_path + ".part"
        with open(tmp, "wb") as f:
            f.write(data)
        os.replace(tmp, img_path)
        return img_path
    except Exception:                                    # noqa: BLE001
        return None


def _element_labels(backend, pts, warnings):
    """逐点快速元素标签（annotate 口径）。失败不致命，只降级成无标签。"""
    labels = {}
    for no, _cam, pt in pts:
        try:
            d = backend.load_point(pt, None, {})
            _marks, elements = backend.annotate_data(d)
            names = [str(e).split()[0] for e in (elements or [])][:3]
            labels[no] = ",".join(names)
        except Exception as exc:                          # noqa: BLE001
            warnings.append("第%d点元素标签失败：%s" % (no, str(exc)[:60]))
            labels[no] = ""
    return labels


def build_element_map(backend, obs, out_dir):
    """观测级元素分布图（产品入口）：RMI + CAHV 投影点 + 逐点元素标签。

    降级路径全部给明确 message，绝不猜坐标：
      无 RMI / 无 cl5 相机模型 / EDR label 取不到 → {"error": …} 或逐帧 warning；
      EDR 原片取不到 → 退回本地 PRC 拉伸片并在 basis 里注明。
    返回 {"maps": [{file, sclk, basis, points}…], "warnings": […],
    "points": n, "note": 口径说明} 或 {"error": …}。
    """
    warnings = []
    pts = []
    for pt in obs.points:
        prod = "rdr" if "rdr" in pt.files else ("ccs" if "ccs" in pt.files else None)
        if not prod:
            continue
        lbl = os.path.splitext(pt.files[prod])[0] + ".lbl"
        cam = libspot_cam(lbl) if os.path.exists(lbl) else None
        if cam is None:
            warnings.append("第%d点缺 cl5 label 相机模型，未投影" % pt.point_no)
        else:
            pts.append((pt.point_no, cam, pt))
    if not pts:
        return {"error": "本观测没有可用的 cl5 相机模型（缺 label），"
                         "无法投影；不猜坐标。请确认 spectra 的 .lbl 与数据同目录。"}
    images = getattr(obs, "rmi_images", [])
    if not images:
        return {"error": "当前文件夹没有本观测的 RMI（CR0 PRC）图像；"
                         "可在下载面板补齐后再生成。"}

    labels = _element_labels(backend, pts, warnings)
    os.makedirs(out_dir, exist_ok=True)
    target = getattr(obs, "target", "")
    maps = []
    all_elements = set()
    for row in images:
        sclk, sol = row["sclk"], int(row["sol"])
        edr_lbl = ensure_edr_label(row["name"], sol)
        rmi_cam = parse_cahv_lbl(edr_lbl) if edr_lbl else None
        if rmi_cam is None:
            warnings.append("SCLK %s：EDR label 未取到或无相机模型，该帧跳过" % sclk)
            continue
        px = project_all([(no, cam) for no, cam, _ in pts], rmi_cam, 2.0)
        inpx = {no: q for no, q in px.items()
                if 2 <= q[0] < rmi_cam["nl"] - 2 and 2 <= q[1] < rmi_cam["ns"] - 2}
        if len(inpx) < len(px):
            warnings.append("SCLK %s：%d 点投影出帧，未绘制" % (sclk, len(px) - len(inpx)))
        base, basis = base_image(edr_lbl, row["path"])
        if base is None:
            warnings.append("SCLK %s：既无 EDR 原片也无本地 PRC，该帧跳过" % sclk)
            continue
        stem = "map_sol%05d_%s_%s" % (sol, row["seq"], sclk)
        cap = "ChemCam 元素分布图 · %d 点 · SCLK %s%s" % (
            len(inpx), sclk, (" · " + target) if target else "")
        _compose(base, inpx, labels, cap, os.path.join(out_dir, stem + ".png"))

        # 逐点坐标（复算/高亮模式的数据源）：JSON 给面板，CSV 给 Excel
        sclk_of = {pt.sclk: pt.point_no for _no, _c, pt in pts}
        pts_rows = []
        for no, (line, sample) in sorted(inpx.items()):
            elem = [e.strip() for e in (labels.get(no) or "").split(",") if e.strip()]
            all_elements.update(elem)
            pts_rows.append({"point": no, "sclk": next((s for s, n in sclk_of.items() if n == no), ""),
                             "line": round(float(line), 1), "sample": round(float(sample), 1),
                             "elements": elem})
        import json as _json
        with open(os.path.join(out_dir, stem + "_points.json"), "w", encoding="utf-8") as f:
            _json.dump({"frame": sclk, "basis": basis, "target": target,
                        "note": MAP_NOTE, "points": pts_rows}, f,
                       ensure_ascii=False, indent=1)
        with open(os.path.join(out_dir, stem + "_points.csv"), "w",
                  encoding="utf-8-sig", newline="") as f:
            f.write("点号,SCLK,line,sample,元素(快速判定)\n")
            for pr in pts_rows:
                f.write("%d,%s,%.1f,%.1f,%s\n" % (pr["point"], pr["sclk"],
                                                   pr["line"], pr["sample"],
                                                   "/".join(pr["elements"])))
        maps.append({"file": stem + ".png", "csv": stem + "_points.csv",
                     "sclk": sclk, "basis": basis, "points": len(inpx),
                     "point_rows": pts_rows})
    if not maps:
        return {"error": "所有帧都未能生成：" + "；".join(warnings[:3])}
    return {"maps": maps, "warnings": warnings, "points": len(pts),
            "elements": sorted(all_elements),
            "note": MAP_NOTE + "。判定≠定量；元素标签为检出名称，不含含量。"}


def highlight_png(edr_lbl_path, fallback_tif, points_json_path, element, cap_prefix=""):
    """高亮模式重绘（内存 PNG bytes，不落盘）：含该元素的点亮、其余灰。

    points_json_path 是 build_element_map 落盘的 *_points.json。材料缺失
    抛 ValueError（调用方给明确提示，不猜）。
    """
    import io
    import json as _json

    rmi_cam = parse_cahv_lbl(edr_lbl_path)
    if rmi_cam is None:
        raise ValueError("EDR label 无相机模型，不能重绘")
    base, basis = base_image(edr_lbl_path, fallback_tif)
    if base is None:
        raise ValueError("底图素材缺失（EDR 原片与本地 PRC 都没有）")
    data = None
    try:
        data = _json.load(open(points_json_path, encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ValueError("逐点坐标表不可读（%s）：%s" % (points_json_path, exc))
    px = {int(p["point"]): (float(p["line"]), float(p["sample"]))
          for p in data.get("points", [])}
    labels = {int(p["point"]): ",".join(p.get("elements", []))
              for p in data.get("points", [])}
    cap = "%s%d 点 · SCLK %s" % (cap_prefix or "ChemCam 元素分布图 · ", len(px),
                                 data.get("frame", ""))
    if data.get("target"):
        cap += " · " + data["target"]
    buf = io.BytesIO()
    _compose(base, px, labels, cap, buf, highlight=element)
    return buf.getvalue()


def burn_soft_signal(rmi_tif_path, px, disc_px=5, trials=300, seed=7):
    """软信号：投影点邻域亮度 vs 随机基线（灼蚀斑偏暗）。只报告不判死。

    返回 dict(mean_at_points, mean_random, darker_fraction)。
    """
    from sherloc.vicar import load_aci
    rng = np.random.default_rng(seed)
    img = np.asarray(load_aci(rmi_tif_path).convert("L").resize((1024, 1024)),
                     dtype=float)
    yy, xx = np.mgrid[0:1024, 0:1024]

    def disc_mean(line, sample):
        m = (xx - sample) ** 2 + (yy - line) ** 2 <= disc_px ** 2
        return float(img[m].mean()) if m.any() else float("nan")

    at_pts = [disc_mean(*p) for p in px.values()]
    rand = []
    for _ in range(trials):
        rand.append(disc_mean(rng.uniform(80, 944), rng.uniform(80, 944)))
    at_pts = [v for v in at_pts if np.isfinite(v)]
    rand = [v for v in rand if np.isfinite(v)]
    return {"mean_at_points": round(float(np.mean(at_pts)), 1),
            "mean_random": round(float(np.mean(rand)), 1),
            "darker_fraction": round(float(np.mean(np.asarray(rand)[None, :]
                                                    < np.asarray(at_pts)[:, None])),
                                     3)}
