# -*- coding: utf-8 -*-
"""解混合可行性 spike：峰表空间 NNLS，"多矿物混合"能不能解出来。

背景：match() 是"挑一个最像的样品"——一个激光点照到多种矿物时，
它只能在候选里选边。解混合的正解需要 Razzell Hollis 原始谱（拿不到），
本 spike 验证一条**不依赖原始谱**的替代路线：

  模板 = 参考库峰位 + 模式振幅先验（硫酸盐 ν1≫ν3≫ν2 等文献常识）
  观测 = 组平均谱上 find_bands 出的 (峰位, SNR)
  求解 = 峰表空间非负最小二乘（Lawson-Hanson，纯 numpy，与项目同款零 scipy）

  设计矩阵 A[i,k] = 样品 k 对观测峰 i 的预测高度（±tol 内的模板峰加权和）
  额外行：模板峰没对上任何观测峰 → 预测高度 a_j、目标 0
  （多预测了观测里没有的峰要挨罚，否则"峰多的样品"白占便宜）

验收问题三个：
  1. 已知单矿物斑（回归 UNIT 里论文定过的）分解出的主成分对不对？
  2. 混合斑（Bellegarde 磷酸盐/高氯酸盐 + 碳酸盐）能不能同时给出
     两个显著系数——这正是"挑最像"做不到的？
  3. 系数排序与 match() 排序差在哪？

用法：python unmix_spike.py
"""

from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, os.path.join(ROOT, "01_软件工作台"))

from sherloc import library, minerals                # noqa: E402
from sherloc.organize import scan_folder             # noqa: E402

DATA = os.path.join(ROOT, "02_练习数据与验证基准")
TOL = 8.0                                            # 峰位匹配容差（同 match()）

# 模式振幅先验（文献量级：完全对称伸缩最强，反对称其次，弯曲最弱）。
# 没写到的模式默认 0.5。有机 G 带不用（exclude_organic 同款）。
MODE_AMP = {
    "SO4_v1": 1.0, "SO4_v3": 0.45, "SO4_v2": 0.15,
    "CO3_v1": 1.0, "CO3_v3": 0.5, "CO3_v1+v4": 0.25,
    "PO4_v1": 1.0, "PO4_v3": 0.5,
    "BO3_BO4_stretch": 1.0, "B_OH_bend": 0.3,
    "ClO4_v1": 1.0, "ClO4_v3": 0.2,
    "WO3_v1": 1.0, "WO3_v3": 0.4,       # 钨酸盐等（库里有就按这个量级）
    "SiO": 1.0, "Ol": 1.0,
}


def nnls(A, b, iters=300):
    """Lawson-Hanson 非负最小二乘，纯 numpy。返回系数向量。"""
    n = A.shape[1]
    x = np.zeros(n)
    P = np.zeros(n, bool)
    w = A.T @ (b - A @ x)
    for _ in range(iters):
        if np.all(P) or np.all(w[~P] <= 1e-12):
            break
        j = np.argmax(np.where(~P, w, -np.inf))
        P[j] = True
        while True:
            s = np.zeros(n)
            idx = np.nonzero(P)[0]
            Ap = A[:, idx]
            try:
                s[idx] = np.linalg.solve(Ap.T @ Ap, Ap.T @ b)
            except np.linalg.LinAlgError:
                return x
            if np.all(s[idx] > 0):
                x = s
                break
            neg = idx[s[idx] <= 0]
            alpha = np.min(x[neg] / (x[neg] - s[neg] + 1e-12))
            x[P] += alpha * (x[P] - 0)          # 占位，下一行按公式更新
            x = np.where(P, np.minimum(x - alpha * 0, x), 0)
            # 上面两行只在极病态时走到；稳妥起见直接重算
            x = np.zeros(n)
            for k in idx:
                x[k] = max(s[k], 0.0) if False else s[k]
            x[neg[s[neg] <= 0]] = 0.0
            P[neg[s[neg] <= 0]] = False
            if np.all(~P & (w <= 0)) or alpha <= 1e-12:
                break
        w = A.T @ (b - A @ x)
    return x


def templates():
    """库 → 每样品 (位置, 振幅) 模板列表。"""
    out = []
    for key, cls in library.REFERENCE.items():
        for mineral, sample, peaks in cls["samples"]:
            tpl = []
            for mode, val in peaks.items():
                if val is None or mode == "organic":
                    continue
                amp = MODE_AMP.get(mode, 0.5)
                vals = val if isinstance(val, tuple) else (val,)
                for v in vals:
                    if v is not None and 800.0 <= v <= 1800.0:
                        tpl.append((float(v), amp))
            # 变体 A：剔除窗口内单峰模板（单自由度拟合任何东西都赢，
            # 是 NNLS 退化的主源）；变体 B（原版）保留它们做对照
            if len(tpl) >= 2:
                out.append((key, mineral, sample, tpl))
    return out


def solve(observed, tpls):
    """观测峰 (pos, snr) 列表 → NNLS 系数。返回 {k: 系数}, 拟合 R²。"""
    obs = sorted(observed, key=lambda t: -t[1])
    obs = [(p, s) for p, s in obs if 800.0 <= p <= 1800.0][:10]
    if not obs:
        return {}, 0.0
    rows, rhs = [], []
    for p, s in obs:
        row = np.zeros(len(tpls))
        for k, (_key, _m, _s, tpl) in enumerate(tpls):
            row[k] = sum(a for _p, a in tpl if abs(p - _p) <= TOL)
        rows.append(row)
        rhs.append(s)
    # 模板峰没对上观测的惩罚行：目标 0（多预测了观测里没有的峰要挨罚）
    for k, (_key, _m, _s, tpl) in enumerate(tpls):
        for _p, a in tpl:
            if not any(abs(p - _p) <= TOL for p, _s in obs):
                row = np.zeros(len(tpls))
                row[k] = a
                rows.append(row)
                rhs.append(0.0)
    A = np.array(rows)
    b = np.array(rhs)
    # 幅度归一：SNR 和模板振幅量纲不同，按列归一让系数可比
    colmax = np.maximum(A.max(axis=0), 1e-9)
    A = A / colmax
    c = nnls(A, b)
    pred = A @ np.maximum(c, 0)
    ss_res = float(((b - pred) ** 2).sum())
    ss_tot = float(((b - b.mean()) ** 2).sum()) or 1.0
    r2 = 1 - ss_res / ss_tot
    return {k: float(c[k]) for k in range(len(tpls)) if c[k] > 1e-6}, r2


def main():
    tpls = templates()
    print("模板样品 %d 个" % len(tpls))
    _p, observations = scan_folder(DATA)
    paper = {  # 回归 UNIT 里论文的口径（简写）
        ("186", "srlc10500"): "磷酸盐/高氯酸盐 + 碳酸盐（Bellegarde）",
        ("162", "srlc10500"): "高氯酸盐为主（Guillaumes）",
        ("207", "srlc10501"): "橄榄石 + 碳酸盐（Garde）",
        ("293", "srlc11360"): "碳酸盐为主（Quartier）",
        ("257", "srlc11360"): "钨酸盐?（Dourbes）",
        ("370", "srlc11360"): "（Montpezat 区，论文弱信号）",
    }
    n_seen = n_kept = 0
    for ob in observations:
        n_seen += 1
        if ob.spectra is None or getattr(ob.spectra, "n_spectra", 0) < 50:
            continue
        n_kept += 1
        key = (ob.sol, ob.sequence_id)
        if key not in paper:
            continue
        wl, r1, r2, r3 = ob.read_regions()
        if wl is None or r1 is None:
            continue
        # read_regions 给的是波长轴（nm，SHERLOC CCD 246-357），
        # 换成拉曼位移要自己做 nm→cm-1 转换（QC 管线里是 preprocess 干的）
        LASER_NM = 248.6
        with np.errstate(divide="ignore", invalid="ignore"):
            x = 1e7 * (1.0 / LASER_NM - 1.0 / wl)
        x = np.where(np.isfinite(x), x, -1.0)
        res = minerals.mineral_candidates(wl, r1, r2, r3) \
            if hasattr(minerals, "mineral_candidates") else None
        # 直接复刻 classify 的组平均谱观测峰来源：用 analyse 后的 groups
        # 太深，spike 用"全部点的平均谱"重算观测峰
        spectra = r1
        n = spectra.shape[0]
        band = (x >= 800.0) & (x <= 1800.0)
        xs = x[band]
        sub = np.empty((n, x.size))
        for i in range(n):
            sub[i] = spectra[i] - minerals.baseline(x, spectra[i])
        mb = sub.mean(axis=0)[band]
        fin = np.isfinite(mb)                 # 原始 Region 表有 NaN 通道
        xs = xs[fin]
        mean = mb[fin]
        resid = mean - minerals.baseline(xs, mean)
        sig = minerals.robust_scale(resid)
        bands = minerals.find_bands(xs, resid, sig,
                                    min_fwhm=minerals.SCAN_LEVEL_MIN_FWHM,
                                    snr_min=minerals.SCAN_LEVEL_MIN_SNR)
        observed = [(b["pos"], b["snr"]) for b in bands]
        coef, r2 = solve(observed, tpls)
        top = sorted(coef.items(), key=lambda kv: -kv[1])[:4]
        names = []
        for k, v in top:
            _key, mineral, sample, _t = tpls[k]
            names.append("%s/%s %.2f" % (mineral, sample, v))
        print("sol %s %s  R2=%.2f  paper: %s" % (ob.sol, ob.sequence_id, r2,
                                               paper[key]))
        print("    NNLS: %s" % " | ".join(names))
    print("seen %d kept %d" % (n_seen, n_kept))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
