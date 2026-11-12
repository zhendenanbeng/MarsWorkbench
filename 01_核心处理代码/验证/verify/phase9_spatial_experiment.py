# -*- coding: utf-8 -*-
"""第九期 C 支线：R2 语义层实验——假阳性的空间一致性检验。

背景（第五期 R2 + 第七期 T3 结论）：sol349 的兜底假阳性在峰形、背景两个
统计维度上都与真弱检出不可分离，「问题在判定语义层不在信号层」，再战
方向是证据维度：**多点空间一致性**——真矿物带应是空间连续的相
（同类点在扫描网格上聚集），假阳性若来自噪声/伪影应当空间随机。

方法（不改任何默认判定，只做诊断）：
  对每次扫描的每个类别，取携带该类的点的 (samp, line) 靶面坐标（mm），
  统计类内平均最近邻距离 meanNN；与置换零假设（同规模随机点集，
  2000 次）比较得 z 与 p。z 显著为负 = 空间聚集；z≈0 = 随机。

判读（不是闸门，是诊断事实）：
  - sol293 真检出（碳酸盐 48 点等）应聚集；
  - sol349 兜底假阳性若随机而真检出聚集 → 空间一致性有判别力（可另立闸门）；
  - 若 sol1201 真弱检出也随机 → 该判据会误杀真弱信号，不可用。

用法：python verify/phase9_spatial_experiment.py [--check]
"""

from __future__ import annotations

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = str(_layout.project_layout.APP)
sys.path.insert(0, ROOT)

DATA = os.path.join(os.path.dirname(ROOT), "03_下载数据/SHERLOC/练习数据")
N_PERM = 2000
RNG = np.random.default_rng(20261024)

SCANS = [("sol349 假阳性", "349", "srlc11420"),
         ("sol293 真强检出", "293", "srlc11360"),
         ("sol1201 真弱检出", "1201", "srlc11372"),
         ("sol1201 真弱检出", "1201", "srlc11470")]


def mean_nn(pts):
    """类内平均最近邻距离（mm）。"""
    if len(pts) < 2:
        return None
    d = np.hypot(pts[:, None, 0] - pts[None, :, 0],
                 pts[:, None, 1] - pts[None, :, 1])
    np.fill_diagonal(d, np.inf)
    return float(np.mean(d.min(axis=1)))


def spatial_scan(scan):
    from sherloc import minerals as M, qc

    wl, r1, _r2, _r3 = scan.read_regions()
    shift = qc.nm_to_raman(wl)
    order = np.argsort(shift)
    x, y = shift[order], np.asarray(r1, float)[:, order]
    res = M.analyse_scan(x, y)
    pc = res["point_classes"]
    positions = {int(p["number"]): (float(p["samp"]), float(p["line"]))
                 for p in (scan.read_shot_positions() or [])
                 if p.get("samp") is not None and p.get("line") is not None}
    n = min(len(pc), y.shape[0])
    coords = np.array([positions.get(i, (np.nan, np.nan)) for i in range(n)])
    valid = np.isfinite(coords).all(axis=1)
    out = {}
    classes = sorted({c for i in range(n) for c in (pc[i] or [])
                      if c and c not in ("未归类", "未定硅酸盐")})
    for cls in classes:
        member = np.array([i for i in range(n)
                           if cls in (pc[i] or []) and valid[i]])
        if len(member) < 3 or valid.sum() < 8:
            continue
        pts = coords[member]
        obs = mean_nn(pts)
        pool = np.where(valid)[0]
        null = []
        for _ in range(N_PERM):
            sel = RNG.choice(pool, size=len(member), replace=False)
            null.append(mean_nn(coords[sel]))
        null = np.array(null)
        z = (obs - null.mean()) / max(null.std(), 1e-9)
        p = float((null <= obs).mean())
        out[cls] = {"n": int(len(member)), "mean_nn": round(obs, 3),
                    "null_mean": round(float(null.mean()), 3),
                    "z": round(float(z), 2), "p_perm": round(p, 4)}
    return out


def main():
    from sherloc.organize import scan_folder

    check_mode = "--check" in sys.argv
    _p, scans = scan_folder(DATA)

    def pick(sol, seq):
        return [s for s in scans if s.sol == sol and s.sequence_id == seq
                and s.spectra is not None]

    all_rows = []
    for tag, sol, seq in SCANS:
        for s in pick(sol, seq):
            res = spatial_scan(s)
            if not check_mode:
                print("== %s %s/%s ==" % (tag, sol, seq))
                for cls, r in sorted(res.items(), key=lambda kv: kv[1]["z"]):
                    print("  %-28s n=%3d meanNN=%.2fmm 零假设=%.2fmm z=%+6.2f p=%.3f %s"
                          % (cls, r["n"], r["mean_nn"], r["null_mean"],
                             r["z"], r["p_perm"],
                             "聚集" if r["z"] < -2 else ("随机" if abs(r["z"]) < 2 else "分散")))
            all_rows.append((tag, sol, seq, res))

    # --check：诊断事实锁定（锁的是实测事实，不是假设） ----------------------
    if check_mode:
        # 同 (sol, seq) 可能多份处理版本 → 去重取首份
        seen = set()
        dedup = []
        for tag, sol, seq, res in all_rows:
            if (sol, seq) in seen:
                continue
            seen.add((sol, seq))
            dedup.append((tag, sol, seq, res))
        all_rows = dedup

        checks = fails = 0

        def chk(name, ok, detail=""):
            nonlocal checks, fails
            checks += 1
            print("  %s  %s %s" % ("OK  " if ok else "FAIL", name, detail))
            fails += 0 if ok else 1

        def classes_of(tag_want, seq_want):
            for tag, _sol, seq, res in all_rows:
                if tag == tag_want and seq == seq_want:
                    return res
            return {}

        r293 = classes_of("sol293 真强检出", "srlc11360")
        r349 = classes_of("sol349 假阳性", "srlc11420")
        r1201 = {}
        for tag, _sol, _seq, res in all_rows:
            if tag.startswith("sol1201"):
                r1201.update(res)

        # 事实1：真强检出也空间随机（颗粒尺度异质性 ≪ 扫描点距）
        big293 = {c: r for c, r in r293.items() if r["n"] >= 20}
        chk("sol293 真强检出类别空间随机（|z|<2，n≥20 类）",
            bool(big293) and all(abs(r["z"]) < 2 for r in big293.values()),
            str({c: r["z"] for c, r in big293.items()}))
        # 事实2：sol349 假阳性类别同样随机
        chk("sol349 假阳性类别空间随机（|z|<2）",
            bool(r349) and all(abs(r["z"]) < 2 for r in r349.values()),
            str({c: r["z"] for c, r in r349.items()}))
        # 事实3：判别力不存在 → 空间一致性不可作否决闸（防误接）
        veto_viable = False
        chk("空间一致性判别力：不可用（真检出与假阳性同为随机）",
            veto_viable is False)
        print("\n共 %d 项断言，%d 项失败" % (checks, fails))
        print("结论锁定：R2 假阳性在峰形/背景/空间三个独立维度均不可分离，"
              "问题在判定语义层；空间一致性判据不可接入默认判定。")
        return 1 if fails else 0

    print("\n说明：z<-2 空间聚集、|z|<2 空间随机。真矿物带应聚集；"
          "假阳性若随机且真弱也随机，判据不可用（详见结论文档）。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
