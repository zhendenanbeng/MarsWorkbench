# -*- coding: utf-8 -*-
"""R2 续作（第七期 T3）：背景模型竞争能否否决 sol349 的兜底假阳性。

第五期 R2 结论：非负峰形幅度（G+L 系数和）判别失败——sol349 碳酸盐组
@1083 以 5.46σ 存活，负系数在真/假检出间不判别。本轮换武器：
**窗口多项式背景竞争检验**（BIC）——

  模型A（有峰）：伪 Voigt（名义峰位/名义FWHM，线性系数）+ 二次背景
  模型B（纯背景）：四次多项式（能把 ≤~FWHM 数倍宽的平滑鼓包整个吸掉）

  峰必要 ⇔ BIC_B − BIC_A > 10（强证据）；否则判「背景即可解释」→ 否决。

物理假设：真拉曼带窄（FWHM 30–80 cm⁻¹），四次多项式在 ±(FWHM+40) 窗口内
拟合不出一个窄凸起；而基线残留/荧光边缘等伪影结构平滑宽广，能被吸掉。

闸门（沿第五期 R2，三条同时满足才考虑接产品）：
  1. sol349 srlc11420 的兜底候选 100% 被否决；
  2. sol1201 的真弱检出兜底候选存活（不被误杀）；
  3. 对照组 141/181 判定零变化（兜底候选本就为空）。

用法：python verify/phase7_r2bg_experiment.py
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
BIC_STRONG = 10.0          # BIC 差过 10 = 强证据（常用阈值）


def candidates(x, y, res):
    """R2 同款提取：矿物组（组内点判类非空）+ 名义峰位/FWHM。"""
    from sherloc import minerals as M

    pc = res["point_classes"]
    out = []
    for g in res.get("groups", []):
        pts = [p for p in g["pts"] if p < len(pc)]
        cls = {c for p in pts for c in pc[p]
               if c and c not in ("未归类", "未定硅酸盐")}
        if not cls:
            continue
        band = g["band"]
        out.append({"classes": sorted(cls),
                    "center": float(band["center"]),
                    "fwhm": float(band.get("fwhm", 60.0)),
                    "n_pts": len(pts)})
    return out


def bic(rss, n, k):
    return n * np.log(max(rss, 1e-300) / n) + k * np.log(n)


def bg_competition(xw, yw, center, fwhm):
    """窗口内 A(峰+二次背景) vs B(四次背景) 的 BIC 竞争。

    xw 归一到 [-1,1] 改善多项式条件数。返回 (dBIC=BIC_B−BIC_A, 拟合FWHM口径说明)。
    """
    t = (xw - xw.mean()) / max(float(np.ptp(xw)), 1e-9) * 2.0
    g = np.exp(-0.5 * ((xw - center) / (fwhm / 2.355)) ** 2)
    l = 1.0 / (1.0 + ((xw - center) / (fwhm / 2.0)) ** 2)
    n = len(xw)
    A = np.column_stack((g + l, np.ones_like(t), t, t * t))
    B = np.column_stack((np.ones_like(t), t, t * t, t ** 3, t ** 4))
    ca, *_ = np.linalg.lstsq(A, yw, rcond=None)
    cb, *_ = np.linalg.lstsq(B, yw, rcond=None)
    ra = yw - A @ ca
    rb = yw - B @ cb
    return bic(float(rb @ rb), n, B.shape[1]) - bic(float(ra @ ra), n, A.shape[1])


def run_scan(scan, tag):
    from sherloc import minerals as M, qc

    wl, r1, _r2, _r3 = scan.read_regions()
    shift = qc.nm_to_raman(wl)
    order = np.argsort(shift)
    x, y = shift[order], np.asarray(r1, float)[:, order]
    res = M.analyse_scan(x, y)
    rows = []
    for c in candidates(x, y, res):
        g = next(g for g in res["groups"]
                 if float(g["band"]["center"]) == c["center"])
        pts = [p for p in g["pts"] if p < y.shape[0]]
        mean = y[pts].mean(axis=0)
        resid = mean - M.baseline(x, mean)
        win = (x >= c["center"] - c["fwhm"] - 40.0) & (x <= c["center"] + c["fwhm"] + 40.0)
        dbic = bg_competition(x[win], resid[win], c["center"], c["fwhm"])
        rows.append(dict(tag=tag, sol=scan.sol, seq=scan.sequence_id,
                         classes="/".join(c["classes"]), center=round(c["center"]),
                         n_pts=c["n_pts"], dbic=round(float(dbic), 1),
                         survives=bool(dbic > BIC_STRONG)))
    return rows


def main():
    from sherloc.organize import scan_folder

    check_mode = "--check" in sys.argv
    _, scans = scan_folder(DATA)

    def pick(sol, seq=None):
        out = [s for s in scans if s.sol == sol and s.spectra is not None
               and (seq is None or s.sequence_id == seq)]
        return out

    rows = []
    for s in pick("349", "srlc11420"):
        rows += run_scan(s, "sol349 假阳性")
    for s in pick("1201"):
        rows += run_scan(s, "sol1201 真弱检出")
    for s in pick("141") + pick("181"):
        rows += run_scan(s, "对照组")

    if not check_mode:
        print("%-14s %-10s %-22s %6s %5s %7s %s" % ("组", "sol/seq", "类别", "峰位", "点数", "dBIC", "峰必要?"))
        for r in rows:
            print("%-14s %-10s %-22s %6d %5d %7.1f %s"
                  % (r["tag"], "%s/%s" % (r["sol"], r["seq"][-5:]), r["classes"][:22],
                     r["center"], r["n_pts"], r["dbic"], "存活" if r["survives"] else "否决"))

    s349 = [r for r in rows if r["tag"].startswith("sol349")]
    s1201 = [r for r in rows if r["tag"].startswith("sol1201")]
    ctrl = [r for r in rows if r["tag"] == "对照组"]
    if not check_mode:
        print()
        print("闸门1 sol349 兜底候选 100%% 否决：%s（%d/%d 否决）"
              % ("过" if s349 and all(not r["survives"] for r in s349) else "不过",
                 sum(1 for r in s349 if not r["survives"]), len(s349)))
        print("闸门2 sol1201 真弱检出存活：%s（%d/%d 存活）"
              % ("过" if s1201 and all(r["survives"] for r in s1201) else "不过",
                 sum(1 for r in s1201 if r["survives"]), len(s1201)))
        print("闸门3 对照组零变化：%s（候选数 %d）"
              % ("过" if not ctrl else "不过（出现候选）", len(ctrl)))
    verdict = (s349 and all(not r["survives"] for r in s349)
               and s1201 and all(r["survives"] for r in s1201) and not ctrl)

    if check_mode:
        # 失败事实锁定（防重蹈）：谁把「背景竞争」当否决闸接进默认判定，
        # 这三条会红——它对真弱检出是大规模误杀武器。
        checks = fails = 0

        def chk(name, ok, detail=""):
            nonlocal checks, fails
            checks += 1
            print("  %s  %s %s" % ("OK  " if ok else "FAIL", name, detail))
            fails += 0 if ok else 1

        n_alive = sum(1 for r in s1201 if r["survives"])
        chk("sol1201 真弱检出被大规模误杀（存活 ≤2/%d）" % len(s1201),
            0 < n_alive <= 2, "存活 %d" % n_alive)
        d349 = sorted(r["dbic"] for r in s349)
        d1201 = sorted(r["dbic"] for r in s1201)
        separable = d349[-1] < d1201[0]      # 存在阈值把两组完全分开？
        chk("dBIC 分布无阈值分离（sol349 max %.1f ≥ sol1201 min %.1f）"
            % (d349[-1], d1201[0]), not separable)
        chk("sol349 假阳性有强 BIC 存活者（≥10，非背景伪影）",
            d349[-1] >= BIC_STRONG, "max dBIC %.1f" % d349[-1])
        print("\n共 %d 项断言，%d 项失败" % (checks, fails))
        print("结论锁定：背景竞争检验判别失败，默认判定不动（详见第七期 T3 文档）。")
        return 1 if fails else 0

    print()
    print("结论：%s" % ("三项闸门全过——背景竞争检验有判别力，可考虑接产品（另立护栏）。"
                       if verdict else
                       "闸门未全过——按纪律不接产品，结论留档（默认判定不动）。"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
