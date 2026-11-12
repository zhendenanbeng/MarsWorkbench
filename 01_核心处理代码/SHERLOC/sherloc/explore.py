# -*- coding: utf-8 -*-
"""SHERLOC 连续特征与探索性分型。只生成派生产物，不修改输入数据。"""
from __future__ import annotations

import csv
import hashlib
import html
import json
import math
import os
from pathlib import Path
from datetime import datetime, timezone
import uuid

import numpy as np
from PIL import Image, ImageDraw

from . import minerals, preprocess, qc, render

DEFAULT_BANDS = ((1060.0, 1110.0), (990.0, 1030.0))
STEPS = {"deadpixel": True, "dark": True, "laser": True,
         "despike": False, "baseline": False}
PALETTE = ("#2563eb", "#ea580c", "#059669", "#9333ea", "#db2777", "#0891b2", "#ca8a04", "#475569")


class CancelledError(Exception):
    pass


def checkpoint(cancel=None):
    if cancel and cancel():
        raise CancelledError("分析已取消，已有产物保留")


def integration_weights(x, low, high):
    """对分段线性光谱精确积分；端点插值的权重也用于噪声传播。"""
    x = np.asarray(x, dtype=float)
    if x.ndim != 1 or len(x) < 2 or not np.isfinite(x).all() or not np.all(np.diff(x) > 0):
        raise ValueError("波数轴必须有限、严格递增且至少有两个通道")
    if not np.isfinite([low, high]).all() or not x[0] <= low < high <= x[-1]:
        raise ValueError("积分窗口必须完整落在有效波数范围内，且下限小于上限")
    nodes = np.r_[low, x[(x > low) & (x < high)], high]
    widths = np.diff(nodes)
    node_weights = np.r_[widths, 0] / 2 + np.r_[0, widths] / 2
    left = np.clip(np.searchsorted(x, nodes, side="right") - 1, 0, len(x) - 2)
    frac = (nodes - x[left]) / (x[left + 1] - x[left])
    weights = np.zeros(len(x))
    np.add.at(weights, left, node_weights * (1 - frac))
    np.add.at(weights, left + 1, node_weights * frac)
    return weights


def band_features(x, spectra, sigma, band, min_snr=3.0):
    weights = integration_weights(x, *band)
    spectra = np.asarray(spectra, dtype=float)
    sigma = np.asarray(sigma, dtype=float)
    if spectra.ndim != 2 or spectra.shape[1] != len(x) or sigma.shape != (len(spectra),):
        raise ValueError("光谱矩阵、波数轴和逐点噪声的维度不一致")
    area = spectra @ weights
    uncertainty = sigma * np.linalg.norm(weights)
    mask = (x >= band[0]) & (x <= band[1])
    # 窄窗无内部采样点时，积分仍成立，但不虚构一个采样峰位。
    peak = np.max(spectra[:, mask], axis=1) if mask.any() else np.full(len(spectra), np.nan)
    snr = np.divide(area, uncertainty, out=np.full_like(area, np.nan), where=uncertainty > 0)
    peak_snr = np.divide(peak, sigma, out=np.full_like(peak, np.nan), where=sigma > 0)
    return {"area": area, "area_sigma": uncertainty, "area_snr": snr,
            "peak_snr": peak_snr, "valid": np.isfinite(snr) & (np.abs(snr) >= min_snr),
            "weights": weights}


def ratio_features(a, b, sigma):
    valid = a["valid"] & b["valid"] & (a["area"] > 0) & (b["area"] > 0)
    values = np.full(len(valid), np.nan)
    uncertainty = np.full(len(valid), np.nan)
    values[valid] = a["area"][valid] / b["area"][valid]
    # 同一条谱上两个窗口可能重叠，保留其共享通道噪声的协方差。
    for i in np.flatnonzero(valid):
        gradient = (a["weights"] / b["area"][i]
                    - a["area"][i] * b["weights"] / b["area"][i] ** 2)
        uncertainty[i] = sigma[i] * np.linalg.norm(gradient)
    return {"value": values, "sigma": uncertainty, "valid": valid}


def spectral_types(x, spectra, candidates, k=3, cancel=None):
    """逐谱 L2 归一化 → PCA(最多10维、目标95%) → 确定性多起点 k-means。"""
    if not isinstance(k, (int, np.integer)) or isinstance(k, bool) or not 2 <= k <= 8:
        raise ValueError("分型组数必须在 2 到 8 之间")
    spectra = np.asarray(spectra, dtype=float)
    norm = np.linalg.norm(spectra, axis=1)
    selected = np.flatnonzero(np.asarray(candidates) & np.isfinite(spectra).all(axis=1) & np.isfinite(norm) & (norm > 0))
    if len(selected) < 2:
        return {"status": "insufficient", "reason": "符合筛选条件的有效光谱不足2条", "selected": selected,
                "labels": np.empty(0, dtype=int), "groups": []}
    y = spectra[selected] / norm[selected, None]
    centered = y - y.mean(axis=0)
    checkpoint(cancel)
    _, singular, loadings = np.linalg.svd(centered, full_matrices=False)
    energy = singular ** 2
    if energy.sum() <= np.finfo(float).eps:
        return {"status": "single_pattern", "reason": "归一化后的光谱形状相同，未强行分组", "selected": selected,
                "labels": np.zeros(len(selected), dtype=int), "effective_k": 1,
                "groups": [{"id":1,"n":len(selected),"points":selected.tolist(),"representative":int(selected[0]),
                            "mean":spectra[selected].mean(axis=0),"mean_normalized":y.mean(axis=0),"small_group":len(selected)<3}]}
    explained = energy / energy.sum()
    dims = max(1, min(10, len(selected) - 1, int(np.searchsorted(np.cumsum(explained), .95)) + 1))
    loadings = loadings[:dims].copy()
    for row in loadings:
        if row[np.argmax(np.abs(row))] < 0:
            row *= -1
    scores = centered @ loadings.T
    distinct = len(np.unique(np.round(scores, 12), axis=0))
    effective_k = min(k, distinct, len(selected))
    best = None
    # 固定起点规则，无随机种子依赖；重复运行可比较。
    starts = [int(np.argmax(np.sum(scores ** 2, axis=1))), int(np.argmin(scores[:, 0])), int(np.argmax(scores[:, 0]))]
    for first in starts:
        centers = [scores[first].copy()]
        while len(centers) < effective_k:
            distances = np.min(np.sum((scores[:, None] - np.asarray(centers)) ** 2, axis=2), axis=1)
            centers.append(scores[int(np.argmax(distances))].copy())
        centers = np.asarray(centers)
        previous = None
        for _ in range(100):
            checkpoint(cancel)
            distances = np.sum((scores[:, None] - centers) ** 2, axis=2)
            labels = np.argmin(distances, axis=1)
            if previous is not None and np.array_equal(labels, previous):
                break
            previous = labels.copy()
            for g in range(effective_k):
                members = scores[labels == g]
                if len(members):
                    centers[g] = members.mean(axis=0)
        inertia = float(np.sum((scores - centers[labels]) ** 2))
        if best is None or inertia < best[0]:
            best = (inertia, labels.copy())
    labels = best[1]
    # 按组平均谱的主峰位置排序，组号仅在本次扫描/本次模型内有效。
    order = sorted(np.unique(labels), key=lambda g: (float(x[np.argmax(y[labels == g].mean(axis=0))]), int(g)))
    remap = {old: new for new, old in enumerate(order)}
    labels = np.array([remap[g] for g in labels])
    groups, distance, outlier = [], np.zeros(len(selected)), np.zeros(len(selected), dtype=bool)
    for g in range(len(order)):
        local = np.flatnonzero(labels == g)
        centroid = scores[local].mean(axis=0)
        d = np.linalg.norm(scores[local] - centroid, axis=1)
        distance[local] = d
        median = float(np.median(d))
        mad = float(np.median(np.abs(d - median))) * 1.4826
        if len(local) >= 5:
            outlier[local] = d > median + 3 * max(mad, 1e-12)
        groups.append({"id": g + 1, "n": len(local), "points": selected[local].tolist(),
                       "small_group": len(local) < 3,
                       "representative": int(selected[local[np.argmin(d)]]),
                       "mean": spectra[selected[local]].mean(axis=0),
                       "mean_normalized": y[local].mean(axis=0)})
    return {"status": "ok", "selected": selected, "labels": labels, "groups": groups,
            "scores": scores, "loadings": loadings, "explained_variance": explained[:dims],
            "retained_variance": float(explained[:dims].sum()), "effective_k": len(groups),
            "requested_k": k, "distance": distance, "outlier": outlier, "inertia": best[0]}


def prepare(scan, cancel=None):
    if not scan.n_points:
        raise ValueError("所选扫描没有可分析的光谱点")
    if (scan.n_points or 0) > 10000:
        raise ValueError("一次分析最多10000点，请先选择单次扫描")
    data = preprocess.run_pipeline(scan, STEPS)
    if data["shift"] is None or data["r1"] is None:
        raise ValueError("分析需要带波长标定的 Region 1 光谱；原始通道数据请先完成标定")
    wide = (data["shift"] >= 800) & (data["shift"] <= 4000)
    x = np.asarray(data["shift"])[wide]
    spectra = np.asarray(data["r1"], dtype=float)[:, wide].copy()
    if len(x) < 15 or not np.all(np.diff(x) > 0):
        raise ValueError("有效波数轴不足或包含重复波数")
    noise = (x >= 2500) & (x <= 3000)
    if noise.sum() < 3:
        raise ValueError("缺少2500–3000 cm⁻¹噪声窗口，无法进行可靠的信噪比筛选")
    finite = np.isfinite(spectra).all(axis=1)
    residual = np.full_like(spectra, np.nan)
    sigma = np.full(len(spectra), np.nan)
    for i in np.flatnonzero(finite):
        checkpoint(cancel)
        residual[i] = spectra[i] - qc.poly_baseline(x, spectra[i])
        sigma[i] = float(np.sqrt(np.mean(residual[i, noise] ** 2)))
    if finite.any():
        sigma[finite], _ = minerals.robust_sigma(sigma[finite], minerals.SIGMA_OUTLIER_FACTOR)
    # 分型的“未归类”和“合格波段”沿用质控结果，不另造一套矿物规则。
    wl, r1, r2, r3 = scan.read_regions()
    unknown, qualified = np.zeros(len(spectra), bool), np.zeros(len(spectra), bool)
    names = [""] * len(spectra)
    valid_rows = np.flatnonzero(finite)
    if len(valid_rows):
        checkpoint(cancel)
        quality = qc.analyse_spectra(wl, r1[valid_rows],
                                    None if r2 is None else r2[valid_rows],
                                    None if r3 is None else r3[valid_rows])
        if len(quality["point_detail"]) != len(valid_rows):
            raise ValueError("质控未返回完整逐点结果，无法确定分型范围")
        for i, detail in zip(valid_rows, quality["point_detail"]):
            named = sorted({c["class"] for c in detail["classes"] if c["key"] != "unassigned"})
            names[i] = "；".join(named)
            unknown[i] = not named
            qualified[i] = detail["n_bands"] > 0
    return {"x": x, "residual": residual, "sigma": sigma, "finite": finite,
            "unknown": unknown, "qualified": qualified, "known_classes": names,
            "log": data["log"] + [{"step": "analysis_baseline", "kind": "done",
                "text": "Region 1：800–4000 cm⁻¹，质控五次多项式基线；2500–3000估计RMS噪声"}]}


def clean_json(value):
    if isinstance(value, dict):
        return {str(k): clean_json(v) for k, v in value.items() if k != "weights"}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [clean_json(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, (float, np.floating)):
        return float(value) if math.isfinite(value) else None
    return value


def write_csv(path, rows, fields):
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def draw_spatial(path, positions, values, title, categorical=False):
    """按 ACI 图像分别出图，坐标系不混合；没有底图仍能保留真实相对位置。"""
    from make_mineral_map import median_spacing, draw_scale_bar
    img = Image.new("RGB", (1100, 900), "white")
    d = ImageDraw.Draw(img)
    font, small = render.font(20), render.font(15)
    points = [{"samp": p["samp"], "line": p["line"]} for p in positions]
    xs, ys = np.array([p["samp"] for p in points]), np.array([p["line"] for p in points])
    dx, dy = max(float(np.ptp(xs)), 1.), max(float(np.ptp(ys)), 1.)
    scale = min(950/dx, 620/dy)
    ox, oy = 550-float((xs.min()+xs.max())/2)*scale, 400-float((ys.min()+ys.max())/2)*scale
    to_px = lambda p: (ox+p["samp"]*scale, oy+p["line"]*scale)
    vals = np.array([values[int(p["number"])] for p in positions], float)
    finite = vals[np.isfinite(vals)]
    low, high = (float(finite.min()), float(finite.max())) if len(finite) else (0., 1.)
    sp = median_spacing([to_px(p) for p in positions])
    radius = max(3, min(16, sp * .35)) if sp else 12
    def color(v):
        if not np.isfinite(v):
            return (190, 190, 190)
        if categorical:
            return PALETTE[(int(v) - 1) % len(PALETTE)]
        if low < 0 < high:
            t = v / max(abs(low), abs(high))
            return (int(250 - 190 * max(0, -t)), int(250 - 170 * abs(t)), int(250 - 190 * max(0, t)))
        t = (v - low) / (high - low) if high > low else .5
        return (int(35 + 210 * t), int(110 - 55 * t), int(210 - 170 * t))
    for p in positions:
        px, py = to_px(p)
        d.ellipse((px-radius, py-radius, px+radius, py+radius),
                  fill=color(values[int(p["number"])]), outline=(255,255,255), width=1)
    d.rectangle((0, 0, 1100, 75), fill="white")
    d.text((20, 12), title, fill=(20,20,20), font=font)
    d.text((20, 43), "RLS相对坐标；未做空间插值；" + str(len(positions)) + "点", fill=(90,90,90), font=small)
    d.rectangle((0, 765, 1100, 900), fill="white")
    if categorical:
        for j, g in enumerate(sorted({int(v) for v in finite})):
            col, row = j % 4, j // 4
            x0, y0 = 25 + col*220, 785 + row*30
            d.ellipse((x0,y0,x0+14,y0+14), fill=color(g))
            d.text((x0+20,y0-3), "类型 %d（%d点）" % (g, np.sum(vals == g)), fill=(30,30,30), font=small)
    elif len(finite):
        for j in range(350):
            d.line((25+j, 790, 25+j, 810), fill=color(low+(high-low)*j/349))
        d.text((25, 817), "%.4g" % low, fill=(30,30,30), font=small)
        d.text((300, 817), "%.4g" % high, fill=(30,30,30), font=small)
        d.text((425, 795), "输入强度 × cm-1" if "积分面积" in title else "无量纲面积比", fill=(60,60,60), font=small)
    else:
        d.text((25, 795), "没有达到门槛的有效数据，不显示数值色标", fill=(90,90,90), font=small)
    d.text((25, 860), "灰色：低信噪比、无效值或未参与分型；非矿物丰度图", fill=(90,90,90), font=small)
    draw_scale_bar(d, 850, 800, scale * 1000 / 10.1, small)
    img.save(path)


def draw_pca(path, model):
    img = Image.new("RGB", (1000, 760), "white")
    d = ImageDraw.Draw(img)
    scores, labels = model["scores"], model["labels"]
    xy = np.column_stack((scores[:, 0], scores[:, 1] if scores.shape[1] > 1 else np.zeros(len(scores))))
    low, high = xy.min(axis=0), xy.max(axis=0)
    span = np.maximum(high-low, 1e-9)
    d.rectangle((85, 80, 940, 630), outline=(120,120,120))
    for i, (v, g) in enumerate(zip(xy, labels)):
        px, py = 110+(v[0]-low[0])/span[0]*800, 605-(v[1]-low[1])/span[1]*500
        d.ellipse((px-5,py-5,px+5,py+5), fill=PALETTE[g], outline="black" if model["outlier"][i] else PALETTE[g])
    explained = model["explained_variance"]
    d.text((30,20), "PCA光谱形状分型（非矿物鉴定）", fill=(20,20,20), font=render.font(25))
    d.text((30,55), "图中显示前2个PC；分组使用%d个PC。少于3点的组仅作异常候选。"%scores.shape[1],fill=(90,90,90),font=render.font(15))
    d.text((350,645), "PC1（%.1f%%）" % (100*explained[0]), fill="black", font=render.font(18))
    d.text((15,85), "PC2 %.1f%%" % (100*explained[1] if len(explained)>1 else 0), fill="black", font=render.font(15))
    for j, g in enumerate(model["groups"]):
        d.text((30+(j%4)*235, 690+(j//4)*25), "类型%d：%d点" % (g["id"],g["n"]), fill=PALETTE[j], font=render.font(17))
    img.save(path)


def validate_parameters(bands, min_snr, k, scope):
    if scope not in ("unknown", "all") or not math.isfinite(min_snr) or not 0 < min_snr <= 100:
        raise ValueError("分型范围必须为unknown/all，信噪比门槛须大于0且不超过100")
    if len(bands) != 2 or any(len(b) != 2 or not 800 <= b[0] < b[1] <= 1800 for b in bands):
        raise ValueError("A、B两个积分窗口必须位于800–1800 cm⁻¹内")
    if not isinstance(k, (int, np.integer)) or isinstance(k, bool) or not 2 <= k <= 8:
        raise ValueError("分型组数须在2–8之间")


def run(scan, out_root, bands=DEFAULT_BANDS, min_snr=3., k=3, scope="unknown", progress=None, cancel=None):
    validate_parameters(bands, min_snr, k, scope)
    name = "sol%s_%s_%s_%s" % (scan.sol, scan.sequence_id,
        datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S"), uuid.uuid4().hex[:6])
    out = Path(out_root).resolve() / name
    out.mkdir(parents=True, exist_ok=False)
    def step(done, message):
        checkpoint(cancel)
        if progress:
            progress(done, 6, message)
    try:
        step(0, "读取光谱、扣除基线并计算质控标签")
        data = prepare(scan, cancel)
        x, spectra, sigma = data["x"], data["residual"], data["sigma"]
        a, b = [band_features(x, spectra, sigma, band, min_snr) for band in bands]
        ratio = ratio_features(a, b, sigma)
        step(1, "积分面积、噪声传播与比值计算完成")
        fingerprint = (x >= 800) & (x <= 1800)
        candidates = data["finite"] & data["qualified"] & np.isfinite(sigma) & (sigma > 0)
        peak_snr = np.divide(np.max(spectra[:, fingerprint], axis=1), sigma,
                             out=np.zeros(n := len(spectra)), where=sigma > 0)
        candidates &= np.nan_to_num(peak_snr, nan=0., posinf=0.) >= min_snr
        if scope == "unknown":
            candidates &= data["unknown"]
        model = spectral_types(x[fingerprint], spectra[:, fingerprint], candidates, k, cancel)
        step(2, "光谱形状分型完成，正在保存逐点数据")
        n = len(spectra)
        labels, distances, outliers = np.full(n, np.nan), np.full(n, np.nan), np.zeros(n, bool)
        selected = model["selected"]
        if model["status"] in ("ok", "single_pattern"):
            labels[selected] = model["labels"] + 1
        if model["status"] == "ok":
            distances[selected], outliers[selected] = model["distance"], model["outlier"]
        positions = scan.read_shot_positions() or []
        pos = {int(p["number"]): p for p in positions if 0 <= int(p["number"]) < n}
        rows = []
        for i in range(n):
            if not data["finite"][i]: reason = "nonfinite_spectrum"
            elif not data["qualified"][i]: reason = "no_qualified_band"
            elif scope == "unknown" and not data["unknown"][i]: reason = "already_classified"
            elif not candidates[i]: reason = "low_snr_or_invalid_noise"
            elif model["status"] == "insufficient": reason = "insufficient_candidates"
            else: reason = ""
            row = {"point_index_0based": i, "point_number_1based": i+1,
                "known_classes": data["known_classes"][i], "noise_sigma": sigma[i],
                "aci_image": pos.get(i,{}).get("image", ""), "samp": pos.get(i,{}).get("samp"),
                "line": pos.get(i,{}).get("line"), "type": labels[i], "excluded_reason": reason,
                "distance_to_centroid": distances[i], "exploratory_outlier": bool(outliers[i]),
                "A_area": a["area"][i], "A_area_sigma": a["area_sigma"][i], "A_area_snr": a["area_snr"][i],
                "A_valid": a["valid"][i], "B_area": b["area"][i], "B_area_sigma": b["area_sigma"][i],
                "B_area_snr": b["area_snr"][i], "B_valid": b["valid"][i],
                "A_div_B": ratio["value"][i], "ratio_sigma": ratio["sigma"][i],
                "ratio_reason": "" if ratio["valid"][i] else "low_snr_nonpositive_or_invalid_area"}
            rows.append(clean_json(row))
        write_csv(out / "points.csv", rows, list(rows[0]))
        if model["status"] in ("ok", "single_pattern"):
            prototypes = []
            for j, wave in enumerate(x[fingerprint]):
                prototypes.append({"shift_cm_1": float(wave), **{
                    "type%d_mean" % g["id"]: float(g["mean"][j]) for g in model["groups"]}, **{
                    "type%d_shape" % g["id"]: float(g["mean_normalized"][j]) for g in model["groups"]}})
            write_csv(out / "type_spectra.csv", prototypes, list(prototypes[0]))
            render.render_png(str(out / "type_spectra.png"),
                [("类型%d (%d点)" % (g["id"],g["n"]),g["mean_normalized"],PALETTE[j])
                 for j,g in enumerate(model["groups"])], x[fingerprint], "拉曼位移 (cm-1)",
                "光谱类型的平均归一化形状", "组号不是矿物名称；分型仅适用于本次扫描", ylabel="L2归一化强度")
        if model["status"] == "ok":
            draw_pca(out / "pca.png", model)
            pca_rows = [{"point_index_0based": int(i), "type": int(model["labels"][j])+1,
                         **{"PC%d"%(c+1): float(v) for c,v in enumerate(model["scores"][j])}}
                        for j,i in enumerate(selected)]
            write_csv(out / "pca_scores.csv", pca_rows, list(pca_rows[0]))
        step(3, "正在生成空间分布图")
        maps = []
        by_image = {}
        for p in pos.values():
            if p.get("samp") is not None and p.get("line") is not None and np.isfinite([p["samp"],p["line"]]).all():
                by_image.setdefault(p["image"], []).append(p)
        fields = [("A_area", np.where(a["valid"],a["area"],np.nan), "A积分面积 %.1f–%.1f cm-1" % tuple(bands[0]),False),
                  ("B_area", np.where(b["valid"],b["area"],np.nan), "B积分面积 %.1f–%.1f cm-1" % tuple(bands[1]),False),
                  ("A_div_B",ratio["value"],"A/B积分面积比",False),
                  ("types",labels,"光谱类型（非矿物名称）",True)]
        for image_i,(image,points) in enumerate(sorted(by_image.items())):
            for field,values,title,cat in fields:
                checkpoint(cancel)
                if field == "types" and model["status"] == "insufficient": continue
                filename = "%s_image%d.png" % (field,image_i+1)
                draw_spatial(out/filename, points, values, title,cat)
                maps.append({"file":filename,"image":image,"quantity":field})
        step(4, "记录分析参数、输入校验值与结果摘要")
        sources = []
        from .label import label_path_for
        seen_sources = set()
        for product in scan.products:
            for source in (product.data_path, label_path_for(product.data_path)):
                if not source: continue
                p = Path(source).resolve()
                if p.is_file() and str(p) not in seen_sources:
                    seen_sources.add(str(p))
                    sources.append({"path":str(p),"sha256":file_hash(p)})
        public_model = {key:value for key,value in model.items() if key not in ("scores","loadings")}
        manifest = {"schema":1,"sol":scan.sol,"sequence_id":scan.sequence_id,"sclks":scan.all_sclks,
            "created_utc":datetime.now(timezone.utc).isoformat(),"points":n,
            "parameters":{"bands_cm_1":bands,"min_snr":min_snr,"scope":scope,"k":k,"preprocessing":STEPS},
            "methods":{"area":"signed piecewise-linear integral of baseline-subtracted R1",
                "area_unit":"input intensity × cm^-1 (not abundance)",
                "uncertainty":"independent channel noise; shared-window covariance included for ratio; baseline/model uncertainty excluded",
                "clustering":"L2 shape normalization; centered PCA <=10 components / target 95%; deterministic multi-start k-means",
                "candidate_gate":"existing QC qualified bands and peak SNR; unknown excludes any named QC class",
                "baseline":{"window_cm_1":[800,4000],"degree":5,"iterations":10,"clip":1.5},
                "spatial":"separate RLS/ACI coordinate frames; no interpolation; ACI scale convention 10.1 um/px"},
            "software":{"numpy":np.__version__,"code_sha256":{
                Path(p).name:file_hash(Path(p)) for p in (__file__,qc.__file__,preprocess.__file__,minerals.__file__)}},
            "preprocessing_log":data["log"],"input_files":sources,"classification":public_model,"maps":maps,
            "valid_feature_points":{"A":int(a["valid"].sum()),"B":int(b["valid"].sum()),"ratio":int(ratio["valid"].sum())},
            "warnings":["光谱类型不是矿物类别；组号不能直接跨扫描比较。",
                "特征积分及面积比不是矿物丰度。","无RLS时只输出表格和光谱分型，不虚构空间坐标。",
                "单像素伪影可能影响积分，请结合质控和原谱复核。",
                "少于3点的分组仅作为异常候选，不能视为稳定光谱类型。"]}
        (out/"analysis.json").write_text(json.dumps(clean_json(manifest),ensure_ascii=False,indent=2,allow_nan=False),encoding="utf-8")
        step(5, "生成可直接打开的分析报告")
        images = sorted(out.glob("*.png"))
        body = '<!doctype html><html lang="zh"><meta charset="utf-8"><title>连续特征与光谱分型</title><style>body{font:16px/1.6 sans-serif;max-width:1100px;margin:30px auto;padding:0 20px}img{max-width:100%}table{border-collapse:collapse}td,th{border:1px solid #ccc;padding:7px}</style>'
        body += '<h1>Sol %s · %s 连续特征与光谱分型</h1>' % (html.escape(scan.sol),html.escape(scan.sequence_id))
        scope_text = "未归类且有合格波段的点" if scope == "unknown" else "全部有合格波段的点"
        status_text = {"ok":"完成分型","insufficient":"候选不足，未分型","single_pattern":"仅一种谱形"}[model["status"]]
        body += '<p>%d点；A窗口%s，B窗口%s cm⁻¹；积分SNR门槛%g；分型范围%s。</p>' % (n,bands[0],bands[1],min_snr,scope_text)
        body += '<p>面积单位：输入强度×cm⁻¹。面积比不是含量比；分型组号不是矿物鉴定。</p>'
        body += '<p>有效面积点：A %d、B %d；有效比值 %d。分型状态：%s；候选%d点。</p>' % (a["valid"].sum(),b["valid"].sum(),ratio["valid"].sum(),status_text,len(selected))
        if model.get("reason"): body += '<p>'+html.escape(model["reason"])+'</p>'
        if not by_image: body += '<p>未找到RLS位置数据，未生成空间地图。</p>'
        body += '<p><a href="points.csv">逐点特征CSV</a> · <a href="analysis.json">参数与输入校验值JSON</a></p>'
        if model["status"] in ("ok", "single_pattern"):
            body += '<p><a href="type_spectra.csv">类型平均谱CSV</a>；少于3点的组仅作异常候选。</p>'
            if model["status"] == "ok":
                body += '<p><a href="pca_scores.csv">PCA得分CSV</a>；保留方差 %.1f%%。</p>' % (100*model["retained_variance"])
            body += '<table><tr><th>光谱类型</th><th>点数</th><th>代表点（从0开始）</th></tr>'
            for g in model["groups"]: body += '<tr><td>%d</td><td>%d</td><td>%d</td></tr>'%(g["id"],g["n"],g["representative"])
            body += '</table>'
        for p in images:
            image_note = next((m["image"] for m in maps if m["file"] == p.name), "")
            body += '<h2>%s</h2><p>%s</p><img src="%s">' % (html.escape(p.stem),html.escape(image_note),p.name)
        (out/"index.html").write_text(body+"</html>",encoding="utf-8")
        step(6, "分析完成")
        return {"out":str(out),"index":str(out/"index.html"),"maps":len(maps),"points":n,
                "types":model.get("effective_k",0),"type_status":model["status"],"candidate_points":len(selected)}
    except BaseException as exc:
        exc.output_path = str(out)
        (out/"job.json").write_text(json.dumps({"complete":False,"cancelled":isinstance(exc,CancelledError),"error":str(exc)},ensure_ascii=False),encoding="utf-8")
        raise


def file_hash(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024*1024), b""):
            h.update(chunk)
    return h.hexdigest()
