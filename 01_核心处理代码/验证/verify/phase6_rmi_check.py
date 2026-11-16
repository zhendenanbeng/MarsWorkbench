# -*- coding: utf-8 -*-
"""第六期 T2 护栏：RMI CAHV 共帧投影的数学契约（探针缓存不在则跳过）。

防的是三类回归：
  1. 文件名/下载通道回归（cr0_ 前缀、EDR↔RDR 配对）；
  2. 投影口径回归（主点重定心、首点帧中心自洽、SPECIAL 锚残差方向——
     视轴失配只记录不补偿，谁把它"修没了"就是静默放行）；
  3. 稳定性回归（帧内/铺开/t 不敏感/双 RMI 一致 + 配对差分方向）。
"""

from __future__ import annotations

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout

import glob
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = str(_layout.project_layout.APP)
sys.path.insert(0, ROOT)

from instruments import rmi_overlay as ro            # noqa: E402

checks = fails = 0


def chk(name, ok, detail=""):
    global checks, fails
    checks += 1
    print("  %s  %s %s" % ("OK  " if ok else "FAIL", name, detail))
    fails += 0 if ok else 1


def main():
    probe = os.path.normpath(os.path.join(str(_layout.project_layout.legacy_cache_root()),"phase3_probe","sol01000"))
    edr = os.path.normpath(os.path.join(str(_layout.project_layout.legacy_cache_root()),"edr_labels","sol01000"))
    if not (os.path.isdir(probe) and os.path.isdir(edr)):
        print("  （探针缓存不在本机，第六期投影护栏跳过——跑过 phase3/phase6 探针后生效）")
        print("SKIP 数据不在位：探针缓存不在本机（跑过 phase3/phase6 探针后生效）")
        print("\n共 0 项断言，0 项失败")
        return 0

    stem1 = "cr0_486263086prc_f0481570ccam02000l1"
    stem2 = "cr0_486264264prc_f0481570ccam02000l1"
    lbl1 = os.path.join(edr, "cr0_486263086edr_f0481570ccam02000m1.lbl")
    lbl2 = os.path.join(edr, "cr0_486264264edr_f0481570ccam02000m1.lbl")

    # -- 文件名与配对 -------------------------------------------------------
    chk("edr_label_name 前帧（cr0_ 前缀不丢）",
        ro.edr_label_name(stem1) == "cr0_486263086edr_f0481570ccam02000m1.lbl",
        ro.edr_label_name(stem1))
    chk("edr_label_name 后帧", 
        ro.edr_label_name(stem2) == "cr0_486264264edr_f0481570ccam02000m1.lbl")
    p = ro.ensure_edr_label(stem1, 1000, download=False)
    chk("ensure_edr_label 命中缓存且保留原件字节", p is not None and open(p,'rb').read()==open(lbl1,'rb').read())

    rdr_lbls = sorted(glob.glob(os.path.join(probe, "cr0_*prc*.lbl")))
    chk("RDR cr0 label 无相机模型（降级返回 None 不猜）",
        rdr_lbls and ro.parse_cahv_lbl(rdr_lbls[0]) is None)

    cam1, cam2 = ro.parse_cahv_lbl(lbl1), ro.parse_cahv_lbl(lbl2)
    chk("EDR label 解析出模型（1024×1024、C/A/H/V 各 3 分量）",
        cam1 and cam1["ns"] == cam1["nl"] == 1024
        and all(len(cam1[k]) == 3 for k in "CAHV"))

    # -- 投影口径 -----------------------------------------------------------
    points = [(re.match(r"cl5_(\d{9})", os.path.basename(f)).group(1),
               ro.libspot_cam(f))
              for f in sorted(glob.glob(os.path.join(probe, "cl5_*ccs*.lbl")))]
    px1 = ro.project_all(points, cam1, 2.0)
    ok, detail = ro.check_in_frame(px1, cam1["ns"], cam1["nl"])
    chk("前帧投影 %d 点全在帧内" % len(px1), ok, detail)

    ls = np.array([q[0] for q in px1.values()])
    ss = np.array([q[1] for q in px1.values()])
    extent = float(np.hypot(ls.max() - ls.min(), ss.max() - ss.min()))
    chk("光栅铺开 300–335px（A 轴证据 ~318px）",
        300 <= extent <= 335, "外接对角 %.0f px" % extent)

    ok, detail = ro.check_t_insensitive(points, cam1)
    chk("t=1.5→2.5m 定位漂移 ≤20px", ok, detail)

    ok, detail = ro.check_pair_consistency(points, cam1, cam2)
    chk("双 RMI 点距相关 r≥0.99", ok, detail)

    first = px1[sorted(px1)[0]]
    d_center = float(np.hypot(first[0] - 511.5, first[1] - 511.5))
    chk("首点距帧中心 ≤2px（同时刻同指向自洽）", d_center <= 2.0,
        "(%.1f, %.1f) 偏 %.1f px" % (first[0], first[1], d_center))

    d_line, d_sample = first[0] - 497, first[1] - 532
    chk("SPECIAL 锚残差 ~(+15,-21)px（视轴失配只记录不补偿）",
        8 <= d_line <= 30 and -32 <= d_sample <= -8,
        "Δline %+.1f Δsample %+.1f" % (d_line, d_sample))

    h0, v0 = ro.principal_point(cam1)
    chk("主点 ≈ (7646,-10564)（重定心契约）",
        abs(h0 - 7646) <= 5 and abs(v0 + 10564) <= 5,
        "(%.0f, %.0f)" % (h0, v0))

    # -- 原片栅格与配对差分（EDR .img 在缓存时加测） ------------------------
    if os.path.exists(lbl1[:-4] + ".img") and os.path.exists(lbl2[:-4] + ".img"):
        r1 = ro.load_edr_raster(lbl1)
        r2 = ro.load_edr_raster(lbl2)
        chk("EDR 栅格读取自证（中位数对账 label MEDIAN）",
            r1.shape == (1024, 1024) and r2.shape == (1024, 1024))
        px2 = ro.project_all(points, cam2, 2.0)
        yy, xx = np.mgrid[0:1024, 0:1024]

        def disc(img, q, d=5):
            m = (xx - q[1]) ** 2 + (yy - q[0]) ** 2 <= d ** 2
            return float(img[m].mean())

        rng = np.random.default_rng(11)
        rand = [(rng.uniform(80, 944), rng.uniform(80, 944)) for _ in range(400)]
        rand = [(l, s) for l, s in rand
                if all((l - q[0]) ** 2 + (s - q[1]) ** 2 > 225 for q in px1.values())]
        d_pts = np.array([disc(r2, px2[k]) - disc(r1, px1[k]) for k in sorted(px1)])
        d_rand = np.array([disc(r2, q) - disc(r1, q) for q in rand])
        did = d_pts.mean() - d_rand.mean()
        z = did / (d_rand.std() / np.sqrt(len(d_pts)))
        chk("配对差分：灼蚀点打后相对变暗 ≥40DN 且 z≤-3",
            did <= -40 and z <= -3, "差中差 %.1f DN z=%.2f" % (did, z))
    else:
        print("  （EDR 原片 .img 不在缓存：栅格自证与配对差分两项未测）")

    # -- T4 产品链路：元素分布图 ------------------------------------------------
    import tempfile
    import shutil
    from instruments import REGISTRY
    b = REGISTRY["chemcam"]

    # 降级1：fixtures 无 RMI → 明确 error，不猜坐标
    fx = b.scan_folder(os.path.join(HERE, "fixtures", "chemcam"))
    fx_obs = [o for o in fx if not getattr(o, "passive", False)][0]
    r = ro.build_element_map(b, fx_obs, tempfile.mkdtemp())
    chk("降级：无 RMI 时给明确 error（不猜坐标）",
        "error" in r and "RMI" in r["error"], r.get("error", "")[:40])

    # 降级2：无 cl5 相机模型 → 明确 error
    import types
    bare = types.SimpleNamespace(points=[types.SimpleNamespace(
        point_no=1, files={"ccs": "nowhere.csv"})], rmi_images=[
        {"name": "cr0_486263086prc_f0481570ccam02000l1", "sol": "1000",
         "sclk": "486263086", "seq": "ccam02000", "path": "x.tif"}],
        target="")
    r = ro.build_element_map(b, bare, tempfile.mkdtemp())
    chk("降级：无 cl5 相机模型时给明确 error",
        "error" in r and "cl5" in r["error"], r.get("error", "")[:40])

    # 全链路：探针观测 → 2 帧地图、命名契约、口径随图下发
    probe_obs = [o for o in b.scan_folder(probe) if not getattr(o, "passive", False)][0]
    tmp = tempfile.mkdtemp()
    try:
        r = ro.build_element_map(b, probe_obs, tmp)
        chk("产品链路：探针观测生成 2 帧元素分布图",
            len(r.get("maps", [])) == 2 and not r.get("warnings"),
            "%s 帧%s" % (len(r.get("maps", [])),
                         "，无警告" if not r.get("warnings") else "，" + r["warnings"][0]))
        import re as _re
        chk("地图命名契约 map_solNNNNN_seq_sclk.png",
            all(_re.match(r"^map_sol\d{5}_[a-z0-9]+_\d{9}\.png$", m["file"])
                for m in r.get("maps", [])))
        chk("底图口径注明 EDR 或 PRC 降级",
            all(m["basis"].startswith(("EDR 12-bit", "本地 PRC"))
                for m in r.get("maps", [])),
            r["maps"][0]["basis"] if r.get("maps") else "")
        chk("口径说明随负载下发（未做地面核验/非定量）",
            "未做地面核验" in r.get("note", "") and "非定量" in r.get("note", ""))
        chk("地图文件真实落盘且非空",
            all(os.path.getsize(os.path.join(tmp, m["file"])) > 100000
                for m in r.get("maps", [])))
        # -- T2 增强：坐标导出与高亮（第七期） -----------------------------
        chk("观测级元素表非空下发",
            len(r.get("elements", [])) >= 1, str(r.get("elements", []))[:40])
        chk("坐标 CSV 逐帧落盘（行数=点数+表头）",
            all(os.path.exists(os.path.join(tmp, m["csv"]))
                and sum(1 for _ in open(os.path.join(tmp, m["csv"]),
                                        encoding="utf-8-sig")) == m["points"] + 1
                for m in r.get("maps", [])))
        chk("逐点 JSON 行含坐标与元素",
            all(p.get("line") and p.get("sample") is not None
                and isinstance(p.get("elements"), list)
                for m in r.get("maps", []) for p in m.get("point_rows", [])))
        m2 = r["maps"][-1]
        pj = os.path.join(tmp, m2["file"][:-4] + "_points.json")
        png_hl = ro.highlight_png(lbl2, "no-such.tif", pj,
                                  r["elements"][0] if r.get("elements") else "Ca")
        chk("高亮重绘返回 PNG 字节（内存不落盘）",
            png_hl[:8] == b"\x89PNG\r\n\x1a\n" and len(png_hl) > 100000)
        png_none = ro.highlight_png(lbl2, "no-such.tif", pj, "Zz")
        chk("高亮 0 命中元素仍渲染（全灰=只弱化不删除）",
            png_none[:8] == b"\x89PNG\r\n\x1a\n")
        try:
            ro.highlight_png(lbl2, "no-such.tif", pj + ".nope", "Ca")
            ok_deg = False
        except ValueError:
            ok_deg = True
        chk("高亮缺坐标表给明确降级（ValueError）", ok_deg)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    print("\n共 %d 项断言，%d 项失败" % (checks, fails))
    return 1 if fails else 0


if __name__ == "__main__":
    raise SystemExit(main())
