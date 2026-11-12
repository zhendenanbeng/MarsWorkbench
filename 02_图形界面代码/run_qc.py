#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""命令行入口：扫描一个 SHERLOC 数据文件夹，生成质控台报告。

用法::

    python run_qc.py <数据文件夹> [-o 输出目录] [--max-points N] [--open]

会对文件夹里每一次观测生成一个 HTML 质控台，并写一个 index.html 汇总。
"""

from __future__ import annotations

import source_paths  # source roots; no-op for bundled imports

import argparse
import csv
import html
import os
import sys
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from sherloc import qc, workbench                       # noqa: E402
from sherloc.organize import scan_folder, summary_rows  # noqa: E402


def _pixl_evidence_column(ob, res, pts, folder, memo, verbose=True):
    """给一次观测的逐点表加"PIXL 化学佐证"列。

    返回 (列名, 与逐点行等长的字符串列表)。任何一步不成（没装 opencv、
    定位没有 PIXL 扫描、配准失败、点没落在化学网格里）都返回
    (列名, 全空列表)——化学是第二证据，不能因为它挂掉拉曼主流程。
    """
    try:
        from sherloc import pixl, watson
    except ImportError:
        return "PIXL 化学佐证", []
    col = [""] * len(res.get("point_detail") or [])
    if not pts:
        return "PIXL 化学佐证", col
    # 点号 → ACI 图像分组；同一图像只做一次配准（memo 跨观测复用）
    by_image = {}
    for p in pts:
        by_image.setdefault(p["image"], []).append(p)
    chem_by_point = {}
    used_medians = []
    n_hit = n_all = 0
    for aci_name, group in by_image.items():
        try:
            from sherloc.watson import normalize as _norm
            key = _norm(aci_name)
        except Exception:                            # noqa: BLE001
            continue
        # 先解析（memo 命中路径也要用 info 拼热图文件名，原来只在未命中
        # 分支里赋值，命中时直接 NameError）
        info = watson.parse_name(aci_name)
        if not info:
            continue
        if key in memo:
            chem = memo[key]
        else:
            aci_path = _find_local_aci(folder, aci_name)
            if not aci_path:
                if verbose:
                    print("  PIXL：本地没有 %s，跳过化学佐证" % aci_name[-24:])
                memo[key] = None
                continue
            try:
                import numpy as np
                from sherloc.vicar import load_aci
                aci_gray = np.asarray(load_aci(aci_path).convert("L"))
                points = [(p["samp"], p["line"]) for p in group]
                chem = pixl.chemistry_at_points(
                    aci_gray, sol=ob.sol, placement=info["att"][:8],
                    points=points, folder=folder)
            except Exception as exc:                 # noqa: BLE001
                if verbose:
                    print("  PIXL：失败（%s: %s）"
                          % (type(exc).__name__, exc))
                chem = None
            memo[key] = chem
        if not chem:
            continue
        used_medians.append(chem["median"])
        # 元素热图（每个 ACI 图只出一次）：SO3/Cl/P2O5 是三类佐证元素
        if chem.get("mcc_gray") is not None and not chem.get("heatmaps"):
            hm = []
            try:
                maps_dir = os.path.join(folder, "_maps")
                os.makedirs(maps_dir, exist_ok=True)
                tag = "sol%s_%s" % (ob.sol, info["att"][:8].lower())
                for elem in ("SO3", "Cl", "P2O5"):
                    outp = os.path.join(
                        maps_dir, "pixl_%s_%s.png" % (tag, elem.lower()))
                    got = pixl.save_heatmap(
                        chem, elem, outp,
                        title="PIXL %s (wt%%) - sol %s %s"
                              % (elem, ob.sol, info["att"][:8]))
                    if got:
                        hm.append(got)
                chem["heatmaps"] = hm
                if hm and verbose:
                    print("  PIXL 元素热图 %d 张 -> _maps/" % len(hm))
            except Exception as exc:                 # noqa: BLE001
                if verbose:
                    print("  PIXL 热图失败（%s: %s）"
                          % (type(exc).__name__, exc))
        # memo 里存的是第一次观测按它自己的点算的逐点结果；换了一批点必须
        # 重算最近邻，否则本次每个点套用的是上一次同序号点的化学
        pts_chem = pixl.reproject_points(
            chem, [(p["samp"], p["line"]) for p in group])
        if pts_chem is None:
            pts_chem = chem["points"]
        for p, (pmc, _dist, elems) in zip(group, pts_chem):
            if pmc is None:
                continue
            num = int(p["number"])
            d = (res.get("point_detail") or [{}])[num] \
                if num < len(res.get("point_detail") or []) else None
            if not d:
                continue
            n_all += 1
            # 化学参与判定：二选一模糊类别（S/Px、P/Cl）按元素丰度拆分。
            # 只动这两类——它们是拉曼缺次级峰逼出来的二选一，化学正是
            # 裁决者；其余类别保持 ✓/✗ 注释，不替用户删矿物。
            from sherloc import minerals as _min
            keys = [c["key"] for c in d["classes"]]
            arb_key, new_keys, note = pixl.arbitrate(keys, elems,
                                                     chem["median"])
            if new_keys:
                d["classes"] = [
                    c if c["key"] != arb_key else
                    {"key": nk, "class": _min.LABELS.get(nk, nk)}
                    for c in d["classes"] for nk in
                    (new_keys if c["key"] == arb_key else [c["key"]])]
                d["chem_arbitrated"] = note
            parts = [pixl.evidence(c["key"], elems, chem["median"])
                     for c in d["classes"]]
            parts = [x for x in parts if x]
            chem_by_point[num] = "；".join(parts)
            if parts:
                n_hit += 1
    if n_all and verbose:
        m = used_medians[0] if used_medians else None
        if m:
            print("  PIXL 化学佐证：%d 点落在化学网格，%d 点有类别可佐证"
                  "（本斑中位 SO3=%.2f Cl=%.2f SiO2=%.1f CaO=%.1f）"
                  % (n_all, n_hit, m.get("SO3", 0), m.get("Cl", 0),
                     m.get("SiO2", 0), m.get("CaO", 0)))
    for i, d in enumerate(res.get("point_detail") or []):
        num = int(d["point"])
        col[i] = chem_by_point.get(num, "")
    return "PIXL 化学佐证", col


def _find_local_aci(folder, aci_name):
    from sherloc.watson import normalize as _norm
    want = _norm(aci_name)
    for root, _dirs, files in os.walk(folder):
        for f in files:
            if _norm(f) == want and f.lower().endswith(
                    (".img", ".png", ".tif", ".jpg")):
                return os.path.join(root, f)
    return None


def analyse_folder(folder, out_dir, max_points=None, verbose=True):
    products, observations = scan_folder(folder)
    pixl_memo = {}                    # ACI 图 → 化学查询结果（跨观测复用）
    os.makedirs(out_dir, exist_ok=True)
    # 清掉上一轮生成的报告（只删本工具自己写的 sol*.html 和 sol*_points.csv）
    removed = 0
    for name in os.listdir(out_dir):
        if name.startswith("sol") and (name.endswith(".html") or name.endswith("_points.csv")):
            try:
                os.remove(os.path.join(out_dir, name))
                removed += 1
            except OSError:
                pass
    if removed and verbose:
        print("清理上一轮报告 %d 个" % removed)
    results = []
    for ob in observations:
        spec = ob.spectra
        if spec is None:
            results.append({"obs": ob, "error": "这次观测里没有光谱产品"})
            continue
        try:
            wl, r1, r2, r3 = ob.read_regions()
            if wl is None:
                raise ValueError("label 里没有波长表（未做波长校正的产品无法换算波数）")
            if r1 is None:
                raise ValueError("label 里没有 Region 1 表")
            res = qc.analyse_spectra(wl, r1, r2, r3, max_points=max_points)
            from instruments import notes as paper_notes
            paper_notes.attach_sherloc(res,folder,ob)
            pts = ob.read_shot_positions()
            out = os.path.join(out_dir, "sol%s_%s_%s.html"
                               % (ob.sol, ob.sequence_id,
                                  "_".join(s.replace(".", "_") for s in ob.all_sclks)))

            # 逐点表：HTML 和 CSV 都写，CSV 用 utf-8-sig，Excel 双击不乱码。
            # 注意顺序：PIXL 佐证要先算——它会对二选一模糊类别做化学仲裁、
            # 改写 point_detail 的类别，逐点表才能反映改判后的结果。
            pixl_header, pixl_col = _pixl_evidence_column(
                ob, res, pts, folder, pixl_memo, verbose=verbose)
            header, rows = qc.point_table(res)
            if any(pixl_col):
                header = header + [pixl_header]
                for row, ev in zip(rows, pixl_col):
                    row.append(ev)
            csv_name = ""
            if rows:
                csv_path = out[:-5] + "_points.csv"
                with open(csv_path, "w", encoding="utf-8-sig", newline="") as f:
                    w = csv.writer(f)
                    w.writerow(header)
                    w.writerows(rows)
                csv_name = os.path.basename(csv_path)

            # 悬停提示里也带上化学佐证（点号 → 佐证文字）
            pixl_ev = {}
            det = res.get("point_detail") or []
            for i, d in enumerate(det):
                if i < len(pixl_col) and pixl_col[i]:
                    pixl_ev[int(d["point"])] = pixl_col[i]
            page = workbench.build_report(ob, res, rls_pts=pts, folder=folder,
                                          points_csv=csv_name,
                                          pixl_evidence=pixl_ev or None)
            with open(out, "w", encoding="utf-8") as f:
                f.write(page)
            results.append({"obs": ob, "res": res, "path": out, "spectra": spec})
            if verbose:
                c = res["counts"]
                print("%-6s %-11s %-22s %-7s %4s 点 | 宽峰 %3d 临界 %3d 尖刺 %3d 无 %3d"
                      % (ob.sol, ob.sequence_id, ob.sclk[:20],
                         ("%d 段" % len(ob.parts)) if ob.is_split else "单段",
                         res["n_points_total"],
                         c["band"], c["edge"], c["spike"], c["none"]))
        except Exception as exc:                        # noqa: BLE001
            results.append({"obs": ob, "error": "%s: %s" % (type(exc).__name__, exc)})
            if verbose:
                print("%-6s %-11s %-16s 失败: %s" % (ob.sol, ob.sequence_id,
                                                     ob.sclk[:14], exc))
                if os.environ.get("SHERLOC_DEBUG"):
                    traceback.print_exc()
    write_index(results, products, folder, out_dir)
    return results


def write_index(results, products, folder, out_dir):
    import appenv
    rows = []
    for r in results:
        ob = r["obs"]
        if "error" in r:
            rows.append("<tr><td>%s</td><td><code>%s</code></td><td>%s</td><td colspan=4 class='bad'>%s</td></tr>"
                        % (ob.sol, html.escape(ob.sequence_id), ob.sclk, html.escape(r["error"])))
            continue
        c = r["res"]["counts"]
        cls = "good" if c["band"] >= 3 else "bad"
        seg = ("%d 段拼接" % len(ob.parts)) if ob.is_split else "单段"
        rows.append(
            "<tr><td>%s</td><td><code>%s</code></td><td>%s</td><td class='num'>%s</td>"
            "<td class='num'>%s</td><td class='num %s'>%s</td><td><a href='%s'>打开质控台</a></td></tr>"
            % (ob.sol, html.escape(ob.sequence_id), ob.sclk + " · " + seg,
               r["res"]["n_points_total"], c["band"] + c["edge"], cls,
               qc.headline(r["res"])[:56] + "…",
               os.path.basename(r["path"])))

    page = """<!doctype html><html lang="zh"><head><meta charset="utf-8">
<title>SHERLOC 质控台 · 汇总</title><style>
body{margin:0;padding:24px;font:14px/1.6 ui-sans-serif,"Microsoft YaHei",sans-serif;
background:#eef0f3;color:#15181c}
h1{font-size:20px;font-weight:500;margin:0 0 4px}
.sub{color:#5f6773;font-size:12px}
table{border-collapse:collapse;margin-top:16px;background:#fff;border-radius:8px;overflow:hidden;width:100%}
th,td{padding:8px 12px;border-bottom:1px solid #edf0f4;text-align:left;font-size:13px;
vertical-align:top;overflow-wrap:anywhere;word-break:break-word}
th{color:#5f6773;font-weight:400;background:#f7f8fa}
.num{font-variant-numeric:tabular-nums}
.good{color:#1a7f4b}.bad{color:#b3352f}
code{background:#f7f8fa;padding:1px 4px;border-radius:3px;font-size:12px}
.wrap{overflow-x:auto;margin-top:16px}
@media (max-width:640px){body{padding:16px 12px}th,td{padding:6px 8px;font-size:12px}
th:nth-child(3),td:nth-child(3){display:none}}
</style></head><body>
<h1>SHERLOC 质控台 · 汇总</h1>
<p class="sub">数据目录 <code>__FOLDER__</code> · __NOBS__ 次观测 · __NFILES__ 个文件</p>
<div class="wrap"><table><thead><tr><th>Sol</th><th>序列号</th><th>SCLK</th><th>点数</th>
<th>宽峰候选</th><th>结论</th><th></th></tr></thead>
<tbody>__ROWS__</tbody></table></div>
<p class="sub" style="margin-top:20px">判据：峰宽须达仪器分辨下限（41.3 cm-1，约 4.3 个通道）且 SNR &gt; 3。
仅"宽峰候选"才可能是真实拉曼带。</p>
</body></html>"""
    page = (page.replace("__FOLDER__", html.escape(folder))
                .replace("__NOBS__", str(len(results)))
                .replace("__NFILES__", str(len(products)))
                .replace("__ROWS__", "\n".join(rows)))
    with open(os.path.join(out_dir, "index.html"), "w", encoding="utf-8") as f:
        f.write(appenv.clean_report(page))


def main(argv=None):
    import runtime_http
    runtime_http.initialize()
    ap = argparse.ArgumentParser(description="SHERLOC 数据质控台")
    ap.add_argument("folder", help="含 PDS4 label 的数据文件夹")
    ap.add_argument("-o", "--out", default=None, help="输出目录（默认 <folder>/_qc）")
    ap.add_argument("--max-points", type=int, default=None, help="每个观测最多分析多少个点")
    ap.add_argument("--open", action="store_true", help="生成后用浏览器打开汇总页")
    args = ap.parse_args(argv)

    folder = os.path.abspath(args.folder)
    out_dir = args.out or os.path.join(folder, "_qc")
    print("扫描:", folder)
    results = analyse_folder(folder, out_dir, max_points=args.max_points)
    print("\n汇总:", os.path.join(out_dir, "index.html"))
    if args.open:
        import webbrowser
        webbrowser.open("file:///" + os.path.join(out_dir, "index.html").replace("\\", "/"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
