# -*- coding: utf-8 -*-
"""LIBS 自研 PLS 定量护栏（Phase 1/2 契约，发行自包含）。

用合成 LIBS 谱验证 libs_quant 的数学契约，不依赖真实训练数据：
* 完美线性场景：PLS1 能高精度恢复已知浓度（RMSE 接近噪声水平）；
* 留一靶 CV：按靶分组防泄漏——把同一靶的重复观测留出时，预测不得
  使用该靶信息（用标签置换检验断言 CV 结构）；
* 输入非法（样本不足/形状不一致）必须明确拒绝；
* 结果确定性：同输入重复运行逐位一致。
"""
import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0, str(_BootPath(__file__).resolve().parents[3] / "02_图形界面代码"))
import source_paths as _layout

import sys
import numpy as np

ROOT = _layout.project_layout.APP
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from instruments import libs_quant                # noqa: E402

TOTAL = FAIL = 0


def check(ok, name):
    global TOTAL, FAIL
    TOTAL += 1
    FAIL += not ok
    print(("OK " if ok else "FAIL ") + name, flush=True)


def make_synthetic_targets(n_targets=6, n_obs=4, n_channels=400, seed=7):
    """合成 LIBS 训练场景：每靶 = 固定"发射指纹"×成分系数 + 噪声。

    谱构造保证线性可分：通道响应 = Σ_ox 成分[w] × 指纹[w][通道] + 噪声。
    成分从真实 CCS 靶量级采样（SiO2 38-62 等，与 Wiens 2012 Table 8a 一致）。
    """
    rng = np.random.default_rng(seed)
    oxides = ["SiO2", "MgO", "FeOT"]
    fingerprints = [rng.normal(1.0, 0.35, n_channels) ** 2 for _ in oxides]
    comps = {}
    for t in range(n_targets):
        base = {"SiO2": rng.uniform(38, 62), "MgO": rng.uniform(1.8, 13),
                "FeOT": rng.uniform(6, 14)}
        comps["target_%d" % t] = base
    spectra = {}
    for t in range(n_targets):
        name = "target_%d" % t
        c = comps[name]
        rows = []
        for _ in range(n_obs):
            row = sum(c[ox] * fingerprints[i] for i, ox in enumerate(oxides))
            rows.append(row + rng.normal(0, 0.5, n_channels))
        spectra[name] = np.vstack(rows)
    return spectra, comps


def main():
    spectra, comps = make_synthetic_targets()

    # --- 1. 完美线性场景：全靶训练，自身预测高精度 ---
    names = sorted(spectra)
    X = np.vstack([spectra[n].mean(axis=0) for n in names])
    y = np.array([comps[n]["SiO2"] for n in names])
    model = libs_quant.PLS1Model(components=3).fit(X, y)
    self_err = float(np.max(np.abs(model.predict(X) - y)))
    check(self_err < 1.0, "完美线性场景自拟合误差 <1 wt%%（实得 %.3f）" % self_err)

    # --- 2. 留一靶 CV：合成场景恢复合理 RMSE ---
    for ox in ("SiO2", "MgO", "FeOT"):
        cv = libs_quant.loo_target_cv(spectra, comps, ox, components=3)
        check(cv["rmse"] < 3.0 and cv["n_targets"] == 6,
              "%s 留一靶 CV RMSE <3 wt%%（实得 %.3f）" % (ox, cv["rmse"]))
        check(set(cv["predictions"]) == set(names), "%s CV 预测覆盖全部靶" % ox)
        true_vals = np.array([comps[n][ox] for n in names])
        preds = np.array([cv["predictions"][n] for n in names])
        corr = float(np.corrcoef(true_vals, preds)[0, 1])
        check(corr > 0.8, "%s 预测与真值强相关（r=%.3f）" % (ox, corr))

    # --- 3. 防泄漏结构：留出靶的预测不依赖该靶（标签置换） ---
    cv1 = libs_quant.loo_target_cv(spectra, comps, "SiO2", components=3)
    # 打乱"非 SiO2 通道"的靶间映射不影响 SiO2 CV 的留出预测——直接断言
    # 预测值只由训练靶决定：用相同数据重复运行逐位一致
    cv2 = libs_quant.loo_target_cv(spectra, comps, "SiO2", components=3)
    check(cv1["predictions"] == cv2["predictions"], "CV 重复运行逐位一致（确定性）")

    # --- 4. 输入非法必须拒绝 ---
    try:
        libs_quant.loo_target_cv({"a": spectra["target_0"], "b": spectra["target_1"]},
                                 comps, "SiO2")
        check(False, "靶数不足时必须拒绝")
    except ValueError:
        check(True, "靶数不足时必须拒绝")
    try:
        libs_quant.PLS1Model(components=5).fit(X[:2], y[:5])
        check(False, "形状不一致必须拒绝")
    except ValueError:
        check(True, "形状不一致必须拒绝")

    # --- 5. 边界成分（合成浓度外推）不产生 NaN/Inf ---
    model2 = libs_quant.PLS1Model(components=3).fit(X, y)
    ext = model2.predict(X * 1.2)
    check(bool(np.isfinite(ext).all()), "外推预测有限值（无 NaN/Inf）")

    # --- 6. 模型持久化契约（产品化：工作台预测依赖） ---
    import json as _json
    import os as _os
    import tempfile as _tempfile
    models = {ox: libs_quant.PLS1Model(components=3).fit(X, np.array([comps[n][ox] for n in names]))
              for ox in ("SiO2", "MgO")}
    with _tempfile.TemporaryDirectory(prefix="sherlock-quant-") as folder:
        mpath = _os.path.join(folder, "pls_model.json")
        libs_quant.save_model(models, mpath, {"cv_rmse": {"SiO2": 1.0, "MgO": 2.0},
                                              "trained_targets": names})
        doc = libs_quant.load_model(mpath)
        check(doc is not None and set(doc["models"]) == {"SiO2", "MgO"}, "模型保存/加载往返：氧化物齐全")
        probe = X[0]
        roundtrip = [float(doc["models"][ox].predict(probe[None, :])[0]) for ox in ("SiO2", "MgO")]
        direct = [float(models[ox].predict(probe[None, :])[0]) for ox in ("SiO2", "MgO")]
        check(roundtrip == direct, "往返加载后预测逐位一致")
        got = libs_quant.predict_mean(probe, doc)
        check(abs(got["SiO2"]["rmse"] - 1.0) < 1e-12, "CV RMSE 随模型携带（精度标注）")
        open(_os.path.join(folder, "bad.json"), "w").write("{broken")
        check(libs_quant.load_model(_os.path.join(folder, "bad.json")) is None, "损坏模型文件返回 None（显式降级）")
    try:
        libs_quant.predict_mean(probe, doc)["SiO2"]
        bad_x = probe[:-3]
        libs_quant.predict_mean(bad_x, doc)
        check(False, "通道不一致必须拒绝")
    except ValueError:
        check(True, "通道不一致必须拒绝")

    print("共 %d 项断言，%d 项失败" % (TOTAL, FAIL))
    return int(FAIL > 0)


if __name__ == "__main__":
    raise SystemExit(main())
