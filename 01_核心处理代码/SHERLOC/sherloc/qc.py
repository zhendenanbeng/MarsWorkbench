# -*- coding: utf-8 -*-
"""质控指标：回答这次观测里到底有没有可用的拉曼信号。

当前默认的工程判据参考 Corpolongo (2023)，并非所有文献通用的检测标准：

1. 候选按30 cm-1宽度阈值筛选；当前候选筛选保留离散格点近似，另输出插值
   宽度和截断诊断。阈值不是断言所有更窄的峰都是宇宙射线。
2. 默认噪声窗口2500-3000 cm-1，可显式配置其他论文使用的窗口。
3. 同时越过3 sigma与宽度门槛只产生候选，不等于矿物学证实。
"""

from __future__ import annotations

import numpy as np

from . import minerals
from .peaks import half_maximum

LASER_NM = 248.6
RESOLUTION_NM = 0.269

# 判据全部对齐 Corpolongo et al. (2023)：
#   噪声 RMSD 取自 2500-3000 cm-1
#   潜在拉曼信号：FWHM >= 30 cm-1（约 3 像素）且强度 >= 3 倍噪声
#
# 注意 30 cm-1 比仪器标称分辨本领（0.269 nm ~ 41 cm-1）更宽松，
# 这是论文采用的行业口径；它靠 3 sigma 阈值和平均谱来补足可靠性。
PIXEL_LO, PIXEL_HI = 50, 2098

# 指纹区窗口的**唯一定义**在 minerals.py（那边是判类内核）。这里只做别名，
# 免得两个文件各写一份、过两天就走散了。
FINGERPRINT = minerals.FINGERPRINT
G_BAND = minerals.G_BAND
WATER = (3200.0, 3700.0)
NOISE_WINDOW = (2500.0, 3000.0)
MIN_FWHM_CM = 30.0
SNR_MIN = 3.0
MIN_WIDTH_CHANNELS = 3          # 兼容旧字段，仅用于展示


def nm_to_raman(nm):
    nm = np.asarray(nm, dtype=float)
    return 1.0e7 / LASER_NM - 1.0e7 / nm


def poly_baseline(x, y, degree=5, n_iter=10, clip=1.5):
    keep = np.ones_like(y, dtype=bool)
    coef = np.polyfit(x[keep], y[keep], degree)
    for _ in range(n_iter):
        resid = y - np.polyval(coef, x)
        sigma = resid[keep].std()
        keep = resid < clip * sigma
        if keep.sum() < degree + 2:
            break
        coef = np.polyfit(x[keep], y[keep], degree)
    return np.polyval(coef, x)


def _verdict(peak, min_ch=None, min_fwhm=MIN_FWHM_CM):
    """按论文口径判定：FWHM >= min_fwhm（默认 30 cm-1）才算潜在拉曼信号。"""
    if peak is None:
        return "none"
    if peak["fwhm"] >= min_fwhm:
        return "band"
    if peak["fwhm"] >= 20.0:
        return "edge"
    return "spike"


VERDICT_TEXT = {
    "band": "宽峰候选",
    "edge": "临界",
    "spike": "窄尖刺",
    "none": "无候选峰",
}


def analyse_spectra(wl_nm, r1, r2=None, r3=None, max_points=None,
                    noise_window=None, min_fwhm=None, robust=True,
                    sigma_factor=minerals.SIGMA_OUTLIER_FACTOR, snr_min=None,
                    raman_region_only=True, prom_sigma=None):
    """对一张光谱图做质控，返回结果字典。

    noise_window / min_fwhm 默认取论文口径（2500-3000、30 cm-1），
    传值可切到操作手册口径（2000-2100、25 cm-1），参数扫描就是这么做的。
    robust=True 表示启用手册 A-5 的 σ 中位数替代规则。
    """
    nw = NOISE_WINDOW if noise_window is None else tuple(noise_window)
    mf = MIN_FWHM_CM if min_fwhm is None else float(min_fwhm)
    snr_gate = SNR_MIN if snr_min is None else float(snr_min)
    # 突出度门槛在调用时才取（这样参数扫描脚本可以改 minerals.PROM_SIGMA
    # 直接生效，不用每处都传参）。
    prom_gate = minerals.PROM_SIGMA if prom_sigma is None else float(prom_sigma)
    wl = np.asarray(wl_nm, dtype=float)
    if (wl.ndim != 1 or len(wl)<3 or not np.isfinite(wl).all() or np.any(wl<=0)
            or len(np.unique(wl))!=len(wl)):
        raise ValueError('波长轴必须是一维、有限、正值且无重复通道')
    if len(nw)!=2 or not np.isfinite(nw).all() or nw[0]>=nw[1]:
        raise ValueError('无效噪声窗口')
    if not np.isfinite([mf,snr_gate,prom_gate,sigma_factor]).all() or mf<=0 or snr_gate<=0 or prom_gate<0 or sigma_factor<=0:
        raise ValueError('无效的质控阈值')
    if max_points is not None and (not isinstance(max_points,(int,np.integer)) or isinstance(max_points,bool) or max_points<=0):
        raise ValueError('max_points必须是正整数')
    r1=np.asarray(r1,dtype=float)
    if r1.ndim!=2 or r1.shape[0]==0 or r1.shape[1]!=len(wl) or not np.isfinite(r1).all():
        raise ValueError('Region1必须是非空、有限且与波长轴等长的二维矩阵')
    for name,extra in [('Region2',r2),('Region3',r3)]:
        if extra is not None and (np.asarray(extra).shape!=r1.shape or not np.isfinite(extra).all()):
            raise ValueError(name+'含非有限数值或与Region1尺寸不一致')
    shift = nm_to_raman(wl)
    order = np.argsort(shift)
    shift = shift[order]
    wl_sorted = wl[order]

    r1 = np.asarray(r1)[:, order]
    r2 = None if r2 is None else np.asarray(r2)[:, order]
    r3 = None if r3 is None else np.asarray(r3)[:, order]

    n_all = r1.shape[0]
    if max_points and n_all > max_points:
        idx = np.linspace(0, n_all - 1, max_points).round().astype(int)
        r1 = r1[idx]
        r2 = None if r2 is None else r2[idx]
        r3 = None if r3 is None else r3[idx]

    usel = slice(PIXEL_LO, PIXEL_HI)
    shift = shift[usel]
    wl_sorted = wl_sorted[usel]
    r1 = r1[:, usel]
    r2 = None if r2 is None else r2[:, usel]
    r3 = None if r3 is None else r3[:, usel]

    # 完整光谱 = 三个 Region 表相加（已用官方 RM* 产品逐点验证，误差 0.000000%），
    # 但那只适用于"重建整条谱"。做拉曼分析要用哪一份，得看信号和噪声各在哪：
    #   R1 覆盖 250-285 nm，即约 225-4124 cm-1 —— 拉曼区 800-1800 cm-1
    #   （253.7-258.6 nm）和噪声区 2500-3000 cm-1（264.7-268.8 nm）都在里面；
    #   R2/R3 在这两处只贡献噪声，不贡献信号。
    # 三表相加会把噪声按平方和抬高（实测约 7 倍）而峰高不变，SNR 就被稀释了。
    # 实测（verify/param_sweep.py）：只用 R1 时 sol 293 得到碳酸盐 48 / 硫酸盐 27
    # （论文 36 / 21），sol 141 归零（论文也是零）；三表相加则只有 3 / 9，且 sol 141
    # 冒出 4 个假类别。所以拉曼分析默认只用 R1。
    combined = r1.copy()
    for extra in (r2, r3):
        if extra is not None:
            combined = combined + extra
    total = r1.copy() if raman_region_only else combined

    # 峰宽必须在整条谱上量，不能被 800-1200 的窗口边界截断
    wide = (shift >= FINGERPRINT[0]) & (shift <= 4000.0)
    xw = shift[wide]
    band = (shift >= FINGERPRINT[0]) & (shift <= FINGERPRINT[1])
    x = shift[band]
    d_chan = float(np.median(np.diff(x))) if x.size > 1 else float("nan")
    mid_lam = float(np.median(wl_sorted[band])) if band.any() else LASER_NM
    res_cm = float(1.0e7 / mid_lam - 1.0e7 / (mid_lam + RESOLUTION_NM))
    min_ch = res_cm / d_chan if d_chan else float("nan")

    quiet_w = (xw >= nw[0]) & (xw <= nw[1])
    if x.size<3 or quiet_w.sum()<3 or xw.size<10:
        raise ValueError('缺少完整指纹区或至少3个噪声窗口通道，不能把无法分析解释成无信号')
    in_win = (xw >= FINGERPRINT[0]) & (xw <= FINGERPRINT[1])

    # 第一遍：算每条谱的基线与原始 σ
    n_spec = total.shape[0]
    resid = np.empty((n_spec, xw.size))
    sub = np.empty((n_spec, x.size))
    sig_raw = np.empty(n_spec)
    for i in range(n_spec):
        resid[i] = total[i, wide] - poly_baseline(xw, total[i, wide])
        sub[i] = total[i, band] - poly_baseline(x, total[i, band])
        sig_raw[i] = (float(np.sqrt(np.mean(resid[i][quiet_w] ** 2)))
                      if quiet_w.any() else float("nan"))

    # 手册 A-5 的 σ 补丁：被荧光抬高的 σ 用扫描中位 σ 顶替
    if robust:
        sig_use, sigma_replaced = minerals.robust_sigma(sig_raw, sigma_factor)
    else:
        sig_use, sigma_replaced = sig_raw, []

    # 第二遍：用最终 σ 找候选峰
    points = []
    for i in range(n_spec):
        ys = resid[i]
        sigma = float(sig_use[i])
        best = None
        n_cand = 0
        for j in range(1, xw.size - 1):
            if not in_win[j]:
                continue
            if not (ys[j] > ys[j - 1] and ys[j] >= ys[j + 1]):
                continue
            if not np.isfinite(sigma) or sigma <= 0 or ys[j] < snr_gate * sigma:
                continue
            n_cand += 1
            half = ys[j] / 2.0
            k = m = j
            while k > 0 and ys[k] > half:k -= 1
            while m < xw.size-1 and ys[m] > half:m += 1
            ch, fw = m-k-1, float(xw[m]-xw[k])
            measured = half_maximum(xw, ys, j, allow_censored=True)
            if best is None or (fw, ys[j]) > (best["fwhm"], best["height"]):
                best = {"pos": float(xw[j]), "height": float(ys[j]),
                        "channels": int(ch), "fwhm": fw,
                        "fwhm_interpolated":measured['fwhm'] if measured else None,
                        "width_censored":measured['censored'] if measured else True,
                        "width_method":"sample_bracket_upper_bound",
                        "snr": float(ys[j] / sigma)}
        points.append({"point": i, "sigma": sigma, "sigma_raw": float(sig_raw[i]),
                       "sigma_replaced": i in sigma_replaced,
                       "n_candidates": n_cand, "peak": best,
                       "verdict": _verdict(best, min_ch, mf)})

    med = np.median(sub, axis=0) if sub.size else np.array([])
    sig_med = float(np.nanmedian([p["sigma"] for p in points])) if points else float("nan")
    se = 1.253 * sig_med / np.sqrt(max(1, len(points)))
    med_peak = float(np.nanmax(med)) if med.size else float("nan")

    fluor_peak = fluor_med = raman = None
    if r2 is not None:
        f_win = (wl_sorted >= 300.0) & (wl_sorted <= 350.0)
        r_win = (wl_sorted >= 253.5) & (wl_sorted <= 276.0)
        if f_win.any() and r_win.any():
            fmed = np.median(combined, axis=0)
            fluor_peak = float(fmed[f_win].max())
            fluor_med = float(np.median(fmed[f_win]))
            raman = float(np.median(np.median(total, axis=0)[r_win]))

    counts = {"band": 0, "edge": 0, "spike": 0, "none": 0}
    for p in points:
        counts[p["verdict"]] += 1

    # 矿物类别判定：按论文口径（逐点检测 -> 按峰位分组 -> 用组内平均谱定类别）
    min_groups, min_classes, n_with_band, point_detail = [], [], None, []
    scan_level = None
    n_split = 0
    try:
        mr = minerals.analyse_scan(xw, total[:, wide], noise_window=nw,
                                   min_fwhm=mf, robust=robust,
                                   sigma_factor=sigma_factor, snr_min=snr_gate,
                                   prom_sigma=prom_gate)
        min_groups = mr["groups"]
        min_classes = [c for g in mr["groups"] for c in g["classes"]]
        n_with_band = mr["n_points_with_band"]
        point_detail = mr["point_detail"]
        n_split = mr.get("n_split", 0)
        scan_level = mr.get("scan_level")
    except Exception as exc:
        raise RuntimeError('矿物分析失败；本次不能输出“未检出”的成功结论：'+str(exc)) from exc

    return {
        "n_points_total": n_all,
        "n_points_analysed": len(points),
        "d_channel_cm": d_chan,
        "resolution_cm": res_cm,
        "min_width_channels": min_ch,
        "min_fwhm_cm": MIN_FWHM_CM,
        "min_fwhm_used": mf,
        "noise_window_used": nw,
        "robust_sigma": bool(robust),
        "sigma_factor": float(sigma_factor),
        "snr_min_used": float(snr_gate),
        "prom_sigma_used": float(prom_gate),
        "raman_region_only": bool(raman_region_only),
        "sigma_replaced": list(sigma_replaced),
        "noise_sigma_median": sig_med,
        "raman_level": float(np.nanmedian(med)) if med.size else float("nan"),
        "raman_p99": float(np.nanpercentile(med, 99)) if med.size else float("nan"),
        "median_peak_height": med_peak,
        "median_peak_sigma": float(med_peak / se) if se else float("nan"),
        "fluorescence_level": fluor_peak,
        "fluorescence_median": fluor_med,
        "raman_median_level": raman,
        "fluor_raman_ratio": (fluor_peak / raman) if (fluor_peak and raman) else None,
        "counts": counts,
        "points": points,
        "sub": sub,
        "median_sub": med,
        "x": x,
        "mineral_groups": min_groups,
        "mineral_classes": min_classes,
        "n_points_with_band": n_with_band,
        "n_points_in_groups": mr['n_points_in_groups'],
        "n_points_with_candidate_bands": mr['n_points_with_candidate_bands'],
        "n_points_with_named_classes": mr['n_points_with_named_classes'],
        "point_detail": point_detail,
        "scan_level": scan_level,
        "n_split": int(n_split),
        "n_review_groups":sum(bool(g.get('requires_review')) for g in min_groups),
    }


def headline(res):
    """一句话结论，带上支撑它的数字。有矿物类别时优先报类别。"""
    groups = res.get("mineral_groups") or []
    named = [g for g in groups if g["classes"] and g["classes"][0]["key"] != "unassigned"]
    if named:
        top = sorted(named, key=lambda g: -g["n_points"])[:4]
        parts = []
        for g in top:
            c = g["classes"][0]
            parts.append("%s（%d 个点，峰位 %.0f cm-1）"
                         % (c["class"], g["n_points"], g["band"]["center"]))
        n_classes=len({c['key'] for g in named for c in g['classes'] if c['key']!='unassigned'})
        review=sum(bool(g.get('requires_review')) for g in named)
        suffix=('其中%d组的峰形或完整峰宽需要人工复核。'%review) if review else ''
        return "检出与 %d 类矿物相符的候选信号：" % n_classes + "、".join(parts) + "。"+suffix

    c = res["counts"]
    n = max(1, res["n_points_analysed"])
    n_band = res.get("n_points_with_band") or 0
    sl = res.get("scan_level")
    if sl and sl.get("bands"):
        b = sl["bands"][0]
        cls = "、".join(x["class"] for x in b["classes"]) or "未归类"
        more = ("；另有 %d 条更弱的候选带" % (len(sl["bands"]) - 1)
                if len(sl["bands"]) > 1 else "")
        return ("只有 %d 个点进入峰组，但全扫描平均谱上有带：最强 %.0f cm-1"
                "（%.1fσ，幅度 %.2f×单点噪声）→ 与%s区间一致%s。"
                "扫描级弱信号候选，不构成矿物判定，不参与与论文的逐点对照。"
                % (n_band, b["center"], b["snr"], b["ratio_to_point_sigma"], cls, more))
    if c["band"] >= max(3, 0.05 * n):
        return ("有 %d 个点出现宽度达标的拉曼带候选（逐点判据），但峰位未落到已知矿物类别区间，"
                "扫描平均谱上也没有够格的候选带，"
                "建议人工复核。" % c["band"])
    if c["band"]:
        return ("只有 %d/%d 个点的峰宽达到判据下限（%d cm-1）；其余 %d 个候选峰太窄，"
                "未通过当前宽度筛选；不能仅凭窄峰认定是宇宙射线。" % (c["band"], n, int(res.get("min_fwhm_cm", 30)), c["spike"]))
    if c["spike"] or c["edge"]:
        return ("没有峰宽达标的拉曼带：%d 个候选峰全部窄于 %d cm-1，"
                "未通过当前宽度筛选；需要排查窄矿物峰或伪影。" % (c["spike"] + c["edge"], int(res.get("min_fwhm_cm", 30))))
    return ("%d-%d cm-1 内没有任何超过 3 sigma 的局部极大值，扫描平均谱上也没有候选带："
            "这次观测没有可用的拉曼信号。"
            % (int(FINGERPRINT[0]), int(FINGERPRINT[1])))


# 逐点表的列。最后两列留给手工对照论文用（论文的逐点标签要人去图上抄）。
POINT_HEADER = [
    "点号", "行", "列", "合格波段数", "最佳峰位_cm-1", "半高宽_cm-1", "信噪比",
    "噪声sigma", "归属聚类_cm-1", "判定类别", "类别缩写", "备注",
    "论文备注", "是否一致", "观测论文对照", "观测论文备注",
]


def point_table(res):
    """逐点表：一行一个点，用来和论文的逐点标签一格一格对照。

    扫描顺序是蛇形：第 1 点在最左上，第 10 点到右上，第 11 点折回下一行右端，
    最后一点在第 10 行最左（与 Corpolongo et al. (2023) 第 2.1 节的描述一致）。
    """
    from . import minerals as _min

    rows = []
    for d in res.get("point_detail") or []:
        i = int(d["point"])
        from instruments.notes import label as paper_label
        paper=res.get('paper_points',{}).get(i) or {};observation=res.get('paper_observation') or {}
        r, c = divmod(i, 10)
        col = c if r % 2 == 0 else 9 - c

        names = [b["class"] for b in d["classes"]]
        shorts = [_min.SHORT.get(b["key"], "?") for b in d["classes"]]
        pos = d["pos"]
        if pos is None:
            verdict, short = "", ""
            note = "无 FWHM>=30 且 >=3sigma 的波段"
        elif not d["classes"]:
            # 这两种「没类别」差别很大，必须分开写：
            # 落在组里但峰位不认识，说明判类规则缺一档；
            # 没成组说明是孤立的单点波段，被 min_points 规则丢了。
            if d["labeled"]:
                verdict, short = "未归类", "?"
                note = ("落在聚类组 %s，但峰位不落在任何已知矿物类别区间"
                        % " ".join("%.0f" % g for g in d["groups"]))
            else:
                verdict, short = "单点未成组", "·"
                note = ("有合格波段（%.1f cm-1），但峰位相同的点太少、未成组，"
                        "按当前规则不计入任何类别" % pos)
        else:
            verdict, short = "/".join(names), "+".join(shorts)
            note = ""
        if pos is not None and any(abs(pos - a) <= 4.0 for a in (970.0, 981.0)):
            note = (note + "；" if note else "") + "峰位在 970/981 假峰带附近"
        if len(d["classes"]) > 1:
            note = (note + "；" if note else "") + "该点同时归属 %d 类" % len(d["classes"])
        if d.get("chem_arbitrated"):
            note = (note + "；" if note else "") + d["chem_arbitrated"]
        if any(c.get('requires_review') for c in d['classes']):
            note = (note + '；' if note else '') + '类别依赖近似峰宽或未通过的峰形拟合，需要复核'

        rows.append([
            i,
            r + 1,
            col + 1,
            d["n_bands"],
            "" if pos is None else "%.1f" % pos,
            "" if d["fwhm"] is None else "%.1f" % d["fwhm"],
            "" if d["snr"] is None else "%.2f" % d["snr"],
            "%.1f" % d["sigma"] if d["sigma"] == d["sigma"] else "",
            " ".join("%.0f" % g for g in d["groups"]),
            verdict,
            short,
            note,
            paper.get('note',''),paper_label(paper),paper_label(observation),observation.get('note',''),
        ])
    return POINT_HEADER, rows


def point_table_text(res):
    """逐点表的 CSV 文本。"""
    import csv
    import io

    header, rows = point_table(res)
    buf = io.StringIO()
    w = csv.writer(buf)
    w.writerow(header)
    w.writerows(rows)
    return buf.getvalue()
