# -*- coding: utf-8 -*-
"""预处理管线：每一步都能单独开关，方便看清每一步到底做了什么。

步骤顺序（与指南里的地面处理流程对应）：

1. deadpixel  去感光死区   SCCD 首尾各 50 个像素不感光，直接裁掉
2. dark       暗背景扣除   ERA - ERB（只有原始产品才需要，RRS 已经扣过）
3. laser      激光能量归一 按光电二极管能量归一到扫描的第一个点
4. despike    去宇宙射线   同一通道跨光谱的稳健统计，超阈值用插值替换
5. baseline   基线扣除     迭代 sigma-clipping 多项式基线

每一步都会在 log 里留一条记录，说明它做了还是跳过了、为什么。
"""

from __future__ import annotations

import numpy as np

DEAD_LO, DEAD_HI = 50, 2098      # SCCD 首尾各 50 个像素不感光
DEFAULT_INTEGRATION_MS = 750.0

STEPS = [
    ("deadpixel", "去感光死区"),
    ("dark", "暗背景扣除"),
    ("laser", "激光能量归一"),
    ("despike", "去宇宙射线"),
    ("baseline", "基线扣除"),
]
STEP_KEYS = [k for k, _ in STEPS]

# 附加开关（不是流程步骤，只在特定情况下才需要打开）。
# despike_force：已处理产品（LASER-NORMALIZED）默认跳过二次去宇宙射线，
# 勾上这个才强制再清一遍。用户选 processed 数据时通常不该勾。
EXTRA_FLAGS = [
    ("despike_force", "强制二次清理"),
]
EXTRA_KEYS = [k for k, _ in EXTRA_FLAGS]


def _match(scan, ptype, sclk):
    """在同一 SCLK 下找另一种产品（例如给 ERA 找 ERB）。"""
    for p in scan.products:
        if p.product_type == ptype and p.pid.sclk == sclk:
            return p
    return None


def load_scan(scan):
    """读入一次扫描的三段光谱，外加配对文件里可能存在的暗背景帧。

    返回 dict：wl_nm / r1 / r2 / r3 / d1 / d2 / d3 / has_dark / wl_source
    """
    wl = None
    wl_source = None
    act = [[], [], []]
    dark = [[], [], []]
    has_dark = False

    # 波长表：优先取同一个扫描里做过波长校正的产品（RRS 等）
    for p in list(scan.parts) + list(scan.products):
        try:
            w = p.wavelength_nm()
        except Exception:                       # noqa: BLE001
            w = None
        if w is not None:
            wl = w
            wl_source = p.name
            break

    # 主动帧
    for p in scan.parts:
        try:
            rs = p.regions()
        except Exception:                       # noqa: BLE001
            continue
        for i in range(3):
            if rs[i] is not None:
                act[i].append(rs[i])
        if p.product_type in ("ERA", "ECA"):
            d = _match(scan, "ERB" if p.product_type == "ERA" else "ECB", p.pid.sclk)
            if d is not None:
                ds = d.regions()
                for i in range(3):
                    if ds[i] is not None:
                        dark[i].append(ds[i])
                has_dark = True

    out = {"wl_nm": wl, "wl_source": wl_source, "has_dark": has_dark,
           "product": scan.parts[0].product_type if scan.parts else "?",
           "n_parts": len(scan.parts)}
    for i in range(3):
        out["r%d" % (i + 1)] = np.vstack(act[i]) if act[i] else None
        out["d%d" % (i + 1)] = np.vstack(dark[i]) if dark[i] else None
    return out


def laser_energy(scan):
    """从 RLI 读激光能量（每个点一个值），被拆分的扫描按段拼接；没有则 None。"""
    chunk = []
    for part in scan.parts:
        rli = None
        for q in scan.products:
            if q.product_type == "RLI" and q.pid.sclk == part.pid.sclk:
                rli = q
                break
        if rli is None:
            for q in scan.products:
                if q.product_type == "RLI":
                    rli = q
                    break
        if rli is None:
            return None
        try:
            t = rli.find_table(r"PHOTODIODE|INTENSITY")
            if t is None:
                return None
            arr = rli.read_table(t).ravel()
        except Exception:                       # noqa: BLE001
            return None
        n = part.n_spectra or arr.size
        chunk.append(arr[:n])
    if not chunk:
        return None
    return np.concatenate(chunk) if len(chunk) > 1 else chunk[0]


def _poly_fit(x, y, degree=5, n_iter=10, clip=1.5):
    """迭代 sigma-clipping 多项式拟合，返回系数。"""
    keep = np.ones_like(y, dtype=bool)
    coef = np.polyfit(x[keep], y[keep], degree)
    for _ in range(n_iter):
        resid = y - np.polyval(coef, x)
        sigma = resid[keep].std()
        keep = resid < clip * sigma
        if keep.sum() < degree + 2:
            break
        coef = np.polyfit(x[keep], y[keep], degree)
    return coef


def run_pipeline(scan, steps, points=None, region=1, despike_sigma=6.0):
    """执行预处理。steps 是开关字典，返回处理前后的数组和一份日志。"""
    steps = dict(steps or {})
    # 「强制二次清理」隐含开启 despike 主开关——否则它单独勾选时完全
    # 无效（despike 分支要求主开关为真才进入，实测踩过：勾了 force、
    # 日志里 despike 仍是"跳过"）。
    if steps.get("despike_force"):
        steps["despike"] = True
    raw = load_scan(scan)
    log = []

    def put(key, text, kind="done"):
        log.append({"step": key, "text": text, "kind": kind})

    first = raw.get("r%d" % region)
    n_max = first.shape[1] if first is not None else 0
    n_pts = first.shape[0] if first is not None else 0

    # 处理标记。**以 label 为准，文件名只当补充**（2026-09-16 修）。
    #
    # 为什么必须这样：文件名第 45-51 位那串标记**不完全可靠**。实测
    # sol 293 的 RRS 文件名是 `...w108cgzj03`——第 51 位是 `z`，
    # 按文件名读就是"没做激光归一"，工具于是又归一了一次；可它的 label 里
    # 三张表的表名都明明白白写着 "LASER-NORMALIZED SPECTRA REGION 1/2/3"，
    # 也就是 PDS 早就归一过了。**重复归一 = 把数据改了两遍。**
    #
    # 实测归纳（全库 36 个光谱产品）：
    #   表名含 "LASER-NORMALIZED" → 已扣暗背景 + 已去宇宙射线 + 已增益校正
    #                              + 已激光归一（就是 RDR 成品）
    #   表名含 "PROCESS DATA"     → 只有车载处理，没归一（对应文件名 b..zpz）
    #   表名 "ACTIVE SPECTRA"     → 原始帧（ECA/ERA 这类，要自己扣暗背景）
    flags = dict(scan.parts[0].pid.flags) if scan.parts else {}
    table_names = " | ".join((t.name or "") for t in
                             (scan.parts[0].tables() if scan.parts else []))
    label_rdr = "LASER-NORMALIZED" in table_names.upper()
    norm_done = label_rdr or flags.get("normalize") == "N"
    cosmic_done = label_rdr or flags.get("cosmic") == "C"
    why = ("label 表名写着 LASER-NORMALIZED" if label_rdr
           else "文件名第 51 位为 N")

    if raw["wl_nm"] is None:
        put("wl", "找不到波长表：这份产品没做过波长校正，只能按通道号作图"
                  "（需要同一次观测里的 RRS 或定标结果提供波长解）", "warn")

    # ---------- 1. 去感光死区 ----------
    sel = slice(None)
    if steps.get("deadpixel", True):
        sel = slice(DEAD_LO, DEAD_HI)
        put("deadpixel", "裁掉首尾各 50 个不感光像素（%d -> %d 个通道，共 %d 条光谱）。"
                         "**这一步不是重复**：PDS 不做它，处理后产品里那 100 个"
                         "通道照样没感光" % (n_max, DEAD_HI - DEAD_LO, n_pts))
    else:
        put("deadpixel", "保留全部 %d 个通道（含首尾不感光的 100 个）" % n_max, "skip")

    R = {}
    for i in (1, 2, 3):
        arr = raw.get("r%d" % i)
        R[i] = None if arr is None else arr[:, sel].astype(float)
    wl = None if raw["wl_nm"] is None else np.asarray(raw["wl_nm"])[sel]

    # ---------- 2. 暗背景扣除 ----------
    if steps.get("dark", True):
        if raw["has_dark"]:
            n_skip = 0
            for i in (1, 2, 3):
                d = raw.get("d%d" % i)
                if R[i] is None or d is None:
                    continue
                if d.shape[0] != R[i].shape[0]:
                    # 拆分扫描只配到一段暗帧时，行数不齐的 Region 直接相减
                    # 会广播崩溃；跳过该 Region 并说明，其余照扣
                    n_skip += 1
                    put("dark", "Region %d 跳过：主动帧 %d 行 vs 暗帧 %d 行对不齐"
                        % (i, R[i].shape[0], d.shape[0]), "warn")
                    continue
                R[i] = R[i] - d[:, sel]
            if n_skip == 0:
                put("dark", "主动帧 − 暗背景帧（ERA − ERB）")
            else:
                put("dark", "暗背景扣除完成（%d 个 Region 因行数不齐跳过）" % n_skip,
                    "warn")
        elif raw["product"] in ("RRS", "RRA", "RCS", "RCA"):
            put("dark", "跳过：%s 是已经扣过暗背景的产品" % raw["product"], "skip")
        else:
            put("dark", "跳过：这份产品没有配对的暗背景帧", "warn")
    else:
        put("dark", "不扣暗背景", "skip")

    # ---------- 3. 激光能量归一 ----------
    if steps.get("laser", True):
        e = laser_energy(scan)
        if norm_done:
            put("laser", "跳过：%s 已经做过激光归一化（%s）——再归一次等于把数据改两遍"
                % (raw["product"], why), "skip")
        elif e is None:
            put("laser", "跳过：找不到 RLI 激光能量表", "warn")
        else:
            n = R[region].shape[0] if R[region] is not None else len(e)
            e = np.asarray(e[:n], dtype=float)
            # 能量为 0/NaN 的点（坏点、丢包）不做归一（增益置 1）。绝不能把
            # NaN 写进矩阵：它会顺着中位数/polyfit 扩散，最后把整条谱
            # 变成 NaN，全程无报错
            bad_e = ~np.isfinite(e) | (e <= 0)
            good_e = ~bad_e
            if good_e.any():
                ref = e[0] if good_e[0] else float(np.mean(e[good_e]))
                gain = np.ones(n)
                gain[good_e] = ref / e[good_e]
                for i in (1, 2, 3):
                    if R[i] is not None:
                        R[i] = R[i] * gain[:R[i].shape[0], None]
                notes = ""
                if bad_e.any():
                    notes += "；%d 个点能量为 0/缺失，未归一" % int(bad_e.sum())
                if e.size < n:
                    notes += "；RLI 只有 %d/%d 个点，缺的未归一" % (e.size, n)
                put("laser", "按光电二极管能量归一到第 0 点（%d 个点，能量 %.1f，"
                             "范围 %.1f-%.1f%s）"
                    % (len(e), ref, float(np.min(e[good_e])),
                       float(np.max(e[good_e])), notes))
                if notes:
                    log[-1]["kind"] = "warn"
            else:
                put("laser", "跳过：RLI 能量全是 0/缺失，无法归一", "warn")
    else:
        put("laser", "不做激光能量归一", "skip")

    # ---------- 4. 去宇宙射线 ----------
    if steps.get("despike", False):
        M = R[region]
        if M is None:
            put("despike", "跳过：没有可用的光谱矩阵", "warn")
        elif cosmic_done and not steps.get("despike_force", False):
            # 已处理产品**默认不再动它**（2026-09-16 改）。以前这里会做"二次清理"，
            # 实测在 sol 293 上改掉了 522 个点——那是 PDS 已经处理过的数据，
            # 用户选 processed 就是想要原件。需要二次清理的话，把这一步的
            # `force` 打开（界面上勾"强制二次清理"）即可。
            put("despike", "跳过：%s 已经去过宇宙射线（%s）。"
                           "如需二次清理，请勾上「强制二次清理」"
                % (raw["product"], why), "skip")
        else:
            med = np.median(M, axis=0)
            mad = np.median(np.abs(M - med), axis=0) * 1.4826
            # 回退值要显式判 NaN：bool(nan) 是 True，`nanmedian(mad) or 1.0`
            # 在最需要回退的全 NaN 场景反而返回 NaN，让这些通道永不剔除
            med_scale = float(np.nanmedian(mad))
            fallback = med_scale if np.isfinite(med_scale) and med_scale > 0 else 1.0
            scale = np.where(np.isfinite(mad) & (mad > 0), mad, fallback)
            bad = (M - med) > despike_sigma * scale
            if bad.any():
                col = np.arange(M.shape[1])
                for i in range(M.shape[0]):
                    idx = np.where(bad[i])[0]
                    if idx.size == 0:
                        continue
                    good = ~bad[i]
                    if good.sum() > 2:
                        M[i, idx] = np.interp(col[idx], col[good], M[i, good])
            put("despike", "按同一通道跨光谱的稳健统计剔除 %d 个尖峰点"
                           "（%d 条光谱 x %d 个通道，阈值 %.0f sigma）"
                % (int(bad.sum()), M.shape[0], M.shape[1], despike_sigma))
    else:
        put("despike", "不去宇宙射线（原始谱里的尖峰原样保留）", "skip")

    # 先算拉曼位移：基线在波数轴上拟合更合理，画图也用这个轴
    shift = None
    if wl is not None:
        shift = 1.0e7 / 248.6 - 1.0e7 / wl
        order = np.argsort(shift)
        shift = shift[order]
        for i in (1, 2, 3):
            if R[i] is not None:
                R[i] = R[i][:, order]
        wl = wl[order]

    # 选点：去宇宙射线要看整幅图，所以放在它之后、基线之前
    if points:
        idx = [int(p) for p in points if 0 <= int(p) < (R[region].shape[0] if R[region] is not None else 0)]
        if idx:
            for i in (1, 2, 3):
                if R[i] is not None:
                    R[i] = R[i][idx]

    # ---------- 5. 基线扣除 ----------
    base = None
    if steps.get("baseline", True) and R[region] is not None:
        M = R[region]
        if shift is not None:
            band = (shift >= 800) & (shift <= 4000)
            if band.sum() > 12:
                base = np.empty_like(M)
                for i in range(M.shape[0]):
                    coef = _poly_fit(shift[band], M[i][band])
                    base[i] = np.polyval(coef, shift)
                    M[i] = M[i] - base[i]
                put("baseline", "在 800-4000 cm-1 上做迭代多项式基线扣除（5 次多项式）。"
                                "**这一步也不是重复**：PDS 从不扣基线，"
                                "荧光背景原样留在数据里——不扣掉，几个波数宽的"
                                "拉曼带根本看不出来")
            else:
                put("baseline", "跳过：有效波段太窄", "warn")
        else:
            put("baseline", "跳过：没有波数轴，无法定位拟合区间", "warn")
    elif not steps.get("baseline", True):
        put("baseline", "不扣基线（含连续背景）", "skip")

    return {
        "wl_nm": wl, "shift": shift, "baseline": base,
        "r1": R[1], "r2": R[2], "r3": R[3],
        "log": log, "product": raw["product"],
        "wl_source": raw["wl_source"], "has_dark": raw["has_dark"],
    }
