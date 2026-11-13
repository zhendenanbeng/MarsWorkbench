# -*- coding: utf-8 -*-
"""ChemCam/SuperCam 通用 LIBS 自研 PLS 定量（Phase 1/2：从谱级原始数据自主算成分）。

自算/对账边界声明（与操作手册口径一致）：
  管线输入 = PDS 谱级原始产品（ChemCam CCS / SuperCam 校准谱 FITS 的发射谱数组）；
  自研科学变换 = 发除尘(base.drop_dust) → 去尖峰(despike_shots) → 基线(poly_baseline)
              → 波长裁剪 → PLS1（本模块 NIPALS 实现）→ 主氧化物 wt%；
  对账基准 = CCS/SCCT 靶文献成分、Wiens 2012 Table 8c 留一 RMSE、Ytsma 2023
           在轨精度、官方 MOC/成分产品——只出现在验证台账，绝不作为输出。

训练口径：仅在轨标定靶（CCS/SCCT）观测谱 × published_calibrations.json 转录的
文献成分。留一靶交叉验证（leave-one-target-out）按靶分组防泄漏——同一靶的
所有观测/发整体进出模型。

诚实边界：官方标定（Clegg 2017 等）使用数百个地面标样 + ICA + 距离校正，
本管线只在轨靶起步，精度预期弱于官方；偏差本身就是交付（对账台账）。
"""
from __future__ import annotations

import numpy as np

OXIDES = ("SiO2", "TiO2", "Al2O3", "FeOT", "MgO", "CaO", "Na2O", "K2O")


class PLS1Model:
    """单氧化物 PLS1（NIPALS）。

    纯 numpy，无第三方依赖。X 是 [样本 x 通道] 光谱矩阵（建议先做
    与仪器管线一致的预处理与裁剪），y 是浓度向量。components 由
    训练/验证循环选择（默认 5，与 Wiens 2012 §Table 8 的量级一致）。
    """

    def __init__(self, components: int = 5):
        self.components = int(components)
        self.x_mean = None
        self.x_std = None
        self.y_mean = None
        self.W = None        # [通道 x 成分] 权重
        self.P = None        # [通道 x 成分] 载荷
        self.q = None        # [成分] y 载荷

    def _center(self, X, mean=None, std=None):
        mean = X.mean(axis=0) if mean is None else mean
        std = X.std(axis=0)
        std = np.where(std > 1e-12, std, 1.0)
        return (X - mean) / std, mean, std

    def fit(self, X, y):
        X = np.asarray(X, dtype=float)
        y = np.asarray(y, dtype=float).ravel()
        if X.shape[0] != y.size or X.shape[0] < self.components + 2:
            raise ValueError("PLS1 训练样本不足或形状不一致")
        Xc, self.x_mean, self.x_std = self._center(X)
        self.y_mean = float(y.mean())
        yc = y - self.y_mean
        E = Xc.copy()
        f = yc.copy()
        W = np.zeros((X.shape[1], self.components))
        P = np.zeros((X.shape[1], self.components))
        q = np.zeros(self.components)
        for a in range(self.components):
            w = E.T @ f
            norm = float(np.linalg.norm(w))
            if norm < 1e-12:
                self.components = a
                W, P, q = W[:, :a], P[:, :a], q[:a]
                break
            w = w / norm
            t = E @ w
            tt = float(t @ t)
            if tt < 1e-12:
                self.components = a
                W, P, q = W[:, :a], P[:, :a], q[:a]
                break
            p = (E.T @ t) / tt
            qa = float(f @ t) / tt
            E = E - np.outer(t, p)
            f = f - qa * t
            W[:, a], P[:, a], q[a] = w, p, qa
        self.W, self.P, self.q = W, P, q
        return self

    def predict(self, X):
        X = np.asarray(X, dtype=float)
        Xc = (X - self.x_mean) / self.x_std
        # ŷ = Xc · W (PᵀW)⁻¹ q：T = XcW，y = Tq，Xc ≈ T Pᵀ ⇒ B = W(PᵀW)⁻¹q
        B = self.W @ np.linalg.inv(self.P.T @ self.W) @ self.q
        return Xc @ B + self.y_mean


def loo_target_cv(spectra_by_target, composition, oxide, components=5, wavelengths=None,
                  lo_hi=None):
    """留一靶交叉验证（leave-one-target-out）。

    spectra_by_target: {靶名: 谱矩阵 [观测 x 通道]}——每靶先对观测求平均谱
                       （发级噪声平均掉，靶间成分差异才是模型要学的）。
    composition: {靶名: {氧化物: wt%}}——文献成分（published_calibrations 转录）。
    返回 {predictions: {靶: 预测值}, rmse, bias, n_targets}；
    靶不足 3 个时抛 ValueError（区间护栏依赖可复现的 CV 数字）。
    """
    names = sorted(spectra_by_target)
    if len(names) < 3:
        raise ValueError("留一靶交叉验证至少需要 3 个靶（现 %d 个）" % len(names))
    X = np.vstack([np.mean(np.asarray(spectra_by_target[n], dtype=float), axis=0)
                   for n in names])
    if lo_hi is not None:
        lo, hi = lo_hi
        X = X[:, lo:hi]
    y = np.array([float(composition[n][oxide]) for n in names])
    preds = {}
    for i, hold in enumerate(names):
        train_idx = [j for j in range(len(names)) if j != i]
        model = PLS1Model(components=components).fit(X[train_idx], y[train_idx])
        preds[hold] = float(model.predict(X[i:i + 1])[0])
    err = np.array([preds[n] - y[j] for j, n in enumerate(names)])
    return {"predictions": preds, "true": dict(zip(names, y.tolist())),
            "rmse": float(np.sqrt(np.mean(err ** 2))),
            "bias": float(np.mean(err)), "n_targets": len(names)}


def summarise_cv(cv_rows):
    """多靶 CV 结果聚合：{氧化物: {rmse, bias, n_targets}}。"""
    return {ox: {"rmse": r["rmse"], "bias": r["bias"], "n_targets": r["n_targets"]}
            for ox, r in cv_rows.items()}


# ---------- 产品化：模型持久化与平均谱预测（Phase 1 工作台入口） ----------

import json as _json
import os as _os
import hashlib as _hashlib
from pathlib import Path as _Path


def save_model(models, path, meta):
    """保存训练好的模型族 {氧化物: PLS1Model} 为 JSON。

    meta 必须携带：成分表 sha256、训练靶清单、各氧化物留一靶 RMSE（作为
    预测精度标注）、components——随模型一起冻结，保证「数字-模型」可追溯。
    """
    doc = {
        "scope": "自研 PLS 成分模型（独立标定，在轨靶最小训练集）。精度以 meta 里的"
                 "留一靶 RMSE 标注，预期弱于官方（Clegg 2017 含 408 地面标样）；"
                 "预测值是候选参考，不是官方成分替代。",
        "meta": meta,
        "models": {ox: {"components": mdl.components,
                        "x_mean": mdl.x_mean.tolist(), "x_std": mdl.x_std.tolist(),
                        "y_mean": mdl.y_mean,
                        "W": mdl.W.tolist(), "P": mdl.P.tolist(), "q": mdl.q.tolist()}
                   for ox, mdl in models.items()},
    }
    path = _Path(path)
    tmp = str(path) + ".part"
    open(tmp, "w", encoding="utf-8").write(_json.dumps(doc, ensure_ascii=False))
    _os.replace(tmp, str(path))
    return path


def load_model(path):
    """加载模型族；文件缺失/损坏/结构非法返回 None（调用方显式降级）。"""
    try:
        doc = _json.loads(_Path(path).read_text(encoding="utf-8"))
        models = {}
        for ox, d in doc["models"].items():
            mdl = PLS1Model(components=d["components"])
            mdl.x_mean = np.asarray(d["x_mean"], dtype=float)
            mdl.x_std = np.asarray(d["x_std"], dtype=float)
            mdl.y_mean = float(d["y_mean"])
            mdl.W = np.asarray(d["W"], dtype=float)
            mdl.P = np.asarray(d["P"], dtype=float)
            mdl.q = np.asarray(d["q"], dtype=float)
            models[ox] = mdl
        return {"models": models, "meta": doc.get("meta", {}), "scope": doc.get("scope", "")}
    except (OSError, ValueError, KeyError, TypeError):
        return None


def predict_mean(mean_spectrum, model_doc):
    """对一条平均谱（预处理后、通道对齐模型）预测全部氧化物。

    返回 {氧化物: {value: wt%, rmse: 该氧化物留一靶 CV RMSE}}。
    模型通道与谱长度不一致时抛 ValueError（由调用方做裁剪对齐）。
    """
    out = {}
    rmse = model_doc["meta"].get("cv_rmse", {})
    x = np.asarray(mean_spectrum, dtype=float).ravel()
    for ox, mdl in model_doc["models"].items():
        if mdl.x_mean.shape != x.shape:
            raise ValueError("模型通道数 %d 与谱通道数 %d 不一致"
                             % (mdl.x_mean.shape[0], x.shape[0]))
        out[ox] = {"value": float(mdl.predict(x[None, :])[0]),
                   "rmse": float(rmse.get(ox)) if rmse.get(ox) is not None else None}
    return out


def file_sha(path):
    h = _hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def auto_components(spectra_by_target, composition, oxide, max_components=4):
    """按留一靶 CV 扫描成分数，返回 (最优 components, 最优 RMSE)。

    范围 1..max_components（留一靶口径下训练集靶数-2 为实际上限）。
    平局取更少成分（更简模型）。
    """
    names = sorted(spectra_by_target)
    best = (1, float("inf"))
    for k in range(1, max_components + 1):
        try:
            r = loo_target_cv(spectra_by_target, composition, oxide,
                              components=k)
        except ValueError:
            break
        if r["rmse"] < best[1]:
            best = (k, r["rmse"])
    return best
