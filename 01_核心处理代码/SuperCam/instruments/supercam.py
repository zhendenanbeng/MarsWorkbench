# -*- coding: utf-8 -*-
"""毅力号 SuperCam LIBS 后端。

数据：urn-nasa-pds-mars2020_supercam / data_calibrated_spectra（PDS4）
  文件名：scam_0700_0729077105_539_cp2_scam01700_目标名___01p01.{fits,xml}
          └仪器┘└sol┘└──SCLK──┘└ms┘└类型┘└序列号┘            └点┘└版┘
  cl1 = LIBS 标定光谱（本后端的主产品）；cp2/cp3 = 被动/红外；
  cr0/cr1 = RMI 图像。

FITS 里（按 EXTNAME 找，不靠顺序）：
  WAVELENGTH  7933 行（波长 µm，×1000 转 nm）
  SPECTRA     7933 × 30 列 Spectrum0..29 —— 一列一发激光
  LASERDATA   每发的电压/电流/光能量（写进日志）
  STATISTICS  每通道 Mean/Median/StDev

同一点位出现多个产品版本（01p01 / 01p02）时只保留最高版。

官方成分预测（MOC 表）：
  PDS 的 data_derived_spectra/supercam_libs_moc.csv 给全任务 12,000+ 个 LIBS 点
  的主量氧化物预测（SiO2/TiO2/Al2O3/FeOT/MgO/CaO/Na2O/K2O，各带标准差），
  按 cdr_fname（= 产品主名去掉末尾两位版本号）与本地文件对应。
  这是 SuperCam 团队官方模型的输出，本软件只做查表展示，不做计算——
  界面上必须带着「官方模型预测值」的口径说明。
"""

from __future__ import annotations

import csv
import os
import re

import numpy as np

from . import base, fitsio

MOC_URL = ("https://pds-geosciences.wustl.edu/m2020/"
           "urn-nasa-pds-mars2020_supercam/data_derived_spectra/supercam_libs_moc.csv")
MOC_OXIDES = ("SiO2", "TiO2", "Al2O3", "FeOT", "MgO", "CaO", "Na2O", "K2O")
MOC_META_ROWS = ("RMSEP", "Precision")


def moc_cache_path():
    """官方成分表的缓存位置（软件目录下 .cache/，一次下载全程共用）。"""
    import appenv
    return str(appenv.copy_cached('supercam_moc.csv'))


def _to_f(v):
    try:
        value = float(v)
        return value if np.isfinite(value) else None
    except (TypeError, ValueError):
        return None


def load_moc(path):
    """解析官方成分预测表。

    返回 (dict 产品主名 -> 成分行, dict 元数据行)。成分行 =
    {"oxides": {名称: (值, 标准差)}, "total": 合计, "distance": 距离 mm}。
    元数据 = {"RMSEP": {氧化物: 值}, "Precision": {...}}（官方给的模型精度）。
    """
    with open(path, encoding="utf-8-sig", errors="replace", newline="") as f:
        lines = f.read().splitlines()
    header_idx = None
    for i, line in enumerate(lines):
        if line.startswith("cdr_fname"):
            header_idx = i
            break
    if header_idx is None:
        raise ValueError("成分表里找不到 cdr_fname 表头行：%s" % path)
    header = next(csv.reader([lines[header_idx]]))

    meta = {}
    if header_idx >= 2:
        names = next(csv.reader([lines[0]]))
        for i in range(1, min(header_idx, 9)):
            vals = next(csv.reader([lines[i]]))
            if not vals or vals[0] not in MOC_META_ROWS:
                continue
            stat = {}
            for j in range(1, min(len(vals), len(names))):
                if vals[j] != "" and names[j] in MOC_OXIDES:
                    stat[names[j]] = _to_f(vals[j])
            meta[vals[0]] = stat

    rows = {}
    for line in lines[header_idx + 1:]:
        if not line.strip():
            continue
        vals = next(csv.reader([line]))
        if len(vals) < len(header):
            continue
        rec = dict(zip(header, vals))
        key = moc_key(rec["cdr_fname"])
        ox = {}
        for name in MOC_OXIDES:
            ox[name] = (_to_f(rec.get(name)), _to_f(rec.get(name + "_stdev")))
        rows[key] = {"oxides": ox, "total": _to_f(rec.get("Total")),
                     "distance": _to_f(rec.get("distance_mm"))}
    return rows, meta


def moc_key(stem):
    """本地产品主名（…06p01）→ MOC 表的键（…06p，去掉末尾两位版本号）。"""
    stem = str(stem).replace("\\", "/").rsplit("/", 1)[-1].rsplit(":", 1)[-1].lower()
    if stem.endswith((".fits", ".xml")):
        stem = stem.rsplit(".", 1)[0]
    return stem[:-2] if re.search(r"p\d{2}$", stem) else stem

NAME_RX = re.compile(
    r"^scam_(\d+)_(\d+)_(\d+)_([a-z]+\d+)_([a-z]+\d{5})_(.+?)__(\d+)p(\d+)\.fits$",
    re.I)


def site_drive_from_xml(path):
    """SuperCam PDS4 label 的 RMC → (site, drive)。文件名/顺序无关，缺字段返回 None。

    跨仪器对照（SuperCam×PIXL）的主判据来源；label 缺失时调用方必须降级说明，
    不猜坐标。
    """
    try:
        text = open(path, encoding="utf-8", errors="replace").read()
    except OSError:
        return None
    idx = dict(re.findall(
        r"index_id>([A-Z]+)</\w+:index_id>\s*<\w+:index_value_number>(\d+)", text))
    if "SITE" in idx and "DRIVE" in idx:
        site, drive = int(idx["SITE"]), int(idx["DRIVE"])
        # 健全性钳制：S1 探针见过个别 label 把别的帧索引读进来（site=17040）
        if 0 <= site < 1000 and 0 <= drive < 100000:
            return site, drive
    return None

TYPE_LABELS = {
    "CL1": "CL1 · LIBS 主动激光光谱（标定）",
    "CP2": "CP2 · 被动红外光谱（1.3–2.6 µm）",
    "CP3": "CP3 · 被动可见/红外光谱（合并统计产品）",
    "CR0": "CR0 · RMI 远景图像",
    "CR1": "CR1 · RMI 远景图像",
}
PASSIVE_TYPES = ("cp2", "cp3")          # 被动光谱：无激光发数、无发射线判定


class SuperCamBackend(base.Backend):
    id = "supercam"
    name = "SuperCam LIBS 激光击穿光谱"
    rover = "m2020"
    rover_label = "毅力号"
    subtitle = "下载 · 看谱 · 元素线识别"
    xlabel = "波长 (nm)"
    ylabel = "强度"
    accent = "#1a63c8"
    steps = [("dust", "丢弃前 5 发（除尘）", True),
             ("despike", "去尖峰（跨发 MAD）", True),
             ("baseline", "连续背景扣除", False),
             ("reflectance", "CP2反射吸收参考（I/F）", False)]
    pds_roots = {
        "calibrated": {
            "label": "标定光谱 data_calibrated_spectra（看图用，推荐）",
            "url": "https://pds-geosciences.wustl.edu/m2020/"
                   "urn-nasa-pds-mars2020_supercam/data_calibrated_spectra"},
    }

    # -- 本地扫描 -----------------------------------------------------------

    def scan_folder(self, root):
        obs_map = {}
        for path, fname, dirname in self.walk_files(root, r".*\.fits$"):
            m = NAME_RX.match(fname)
            if not m:
                continue
            ptype = m.group(4).lower()
            if ptype not in ("cl1",) + PASSIVE_TYPES:
                continue
            sol = int(m.group(1))
            sclk = m.group(2)
            seq = m.group(5).lower()
            target = m.group(6).strip("_").replace("_", " ")
            pno = int(m.group(7))
            ver = int(m.group(8))
            key = (sol, seq, ptype)          # 不同产品类型 = 不同观测
            obs = obs_map.get(key)
            if obs is None:
                obs = obs_map[key] = base.Observation(
                    sol, seq, ptype.upper())
                obs.target = target
                xml = os.path.splitext(path)[0] + ".xml"
                sd = site_drive_from_xml(xml) if os.path.exists(xml) else None
                obs.site, obs.drive = (sd if sd else (None, None))
            obs.points.append(base.PointFile(
                os.path.splitext(fname)[0], {ptype: path}, sclk=sclk,
                point_no=pno,
                meta={"version": ver, "target": target, "ptype": ptype}))
        out = []
        for key, obs in sorted(obs_map.items()):
            # 同点号多版本 → 留最高版
            best = {}
            for pt in obs.points:
                cur = best.get(pt.point_no)
                if cur is None or pt.meta.get("version", 0) > cur.meta.get("version", 0):
                    best[pt.point_no] = pt
            obs.points = [best[k] for k in sorted(best)]
            if obs.points:
                out.append(obs)
        return out

    # -- 读点 ---------------------------------------------------------------

    def load_point(self, point, shot, steps_on):
        ptype = point.meta.get("ptype", "cl1")
        path = next(iter(point.files.values()))
        log = [("读取", os.path.basename(path))]
        cols = fitsio.read_bintable(path, "WAVELENGTH")
        wl_key = next((k for k in cols if k.lower().startswith("wavelength")), None)
        if wl_key is None:
            raise ValueError("FITS 里没有 WAVELENGTH 扩展")
        x = np.asarray(cols[wl_key], dtype=float)

        if ptype in PASSIVE_TYPES:
            return self._load_passive(path, ptype, x, log, point, steps_on)

        if np.nanmedian(x) < 10:                 # µm → nm
            x = x * 1000.0
            log[-1] = ("读取", "%s · %d 通道（µm→nm）" % (os.path.basename(path), len(x)))

        spec = fitsio.read_bintable(path, "SPECTRA")
        shot_keys = sorted((k for k in spec if k.lower().startswith("spectrum")
                            and "." not in k),
                           key=lambda k: (len(k), k))
        shots = np.column_stack([np.asarray(spec[k], dtype=float)
                                 for k in shot_keys])
        laser = None
        try:
            laser = fitsio.read_bintable(path, "LASERDATA")
        except Exception:                            # noqa: BLE001
            pass
        if laser:
            for k in laser:
                if "optical power" in k.lower() and "mj" in k.lower():
                    vals = np.asarray(laser[k], dtype=float)
                    vals = vals[np.isfinite(vals)]
                    if vals.size:
                        log.append(("激光", "共 %d 发，平均光能量 %.1f mJ"
                                    % (vals.size, float(np.mean(vals)))))
                    break
        return self._process(x, shots, shot, steps_on, log, point)

    def _load_passive(self, path, ptype, x, log, point, steps_on):
        """被动/红外光谱：单条谱，无发数、无发射线判定。

        CP2 的 SPECTRA 是 6 个物理量列（Spectrum DN / Radiance / …）× 256 行；
        CP3 的 SPECTRA 可能为空，实际数据在 STATISTICS 的 Mean。
        """
        from . import base as _b
        y = None
        ysrc = ""
        try:
            spec = fitsio.read_bintable(path, "SPECTRA")
            n_rows = len(next(iter(spec.values()))) if spec else 0
        except Exception:                             # noqa: BLE001
            spec, n_rows = {}, 0
        if n_rows:
            for cand in ("Spectrum DN", "spectrum", "Radiance"):
                if cand in spec:
                    y = np.asarray(spec[cand], dtype=float)
                    ysrc = cand
                    break
            if y is None and spec:
                y = np.asarray(next(iter(spec.values())), dtype=float)
                ysrc = next(iter(spec))
        if y is None:
            stat = fitsio.read_bintable(path, "STATISTICS")
            if "Mean" in stat:
                y = np.asarray(stat["Mean"], dtype=float)
                ysrc = "STATISTICS.Mean"
        if y is None:
            raise ValueError("这个被动光谱文件里没有可读的谱列")
        reflection=ptype=='cp2' and bool(steps_on.get('reflectance'))
        if reflection:
            lookup={k.lower():k for k in spec}
            ysrc=next((lookup[k] for k in ('i_f_atm','i_f') if k in lookup),'')
            if not ysrc:raise ValueError('原产品缺少I/F列，不能对DN或Radiance作反射吸收评分')
            y=np.asarray(spec[ysrc],float)
        order = np.argsort(x)               # CP2 波长轴是降序的，统一升序方便作图
        x, y = np.asarray(x, dtype=float)[order], y[order]
        in_um = np.nanmedian(x) < 10
        xlabel = "波长 (µm)" if in_um else "波长 (nm)"
        if steps_on.get("baseline") and not reflection:
            y = _b.poly_baseline(x, y, log=log)
        elif steps_on.get('baseline'):log.append(('基线','反射吸收模式保留I/F，跳过全局扣背景；各参考带单独取局部连续谱'))
        log.append(("被动", "%s · %d 通道 · 谱列 %s" % (ptype.upper(), len(x), ysrc)))
        from . import absorption
        references=absorption.payload(x if in_um else x/1000,y if reflection else None,ysrc,1 if in_um else 1000) if ptype=='cp2' else None
        return {"x": x, "y": y, "log": log, "n_shots": 1,
                "ylabel": (ysrc+'（原产品反射数据）') if reflection else ("Spectrum DN（计数）" if ysrc == "Spectrum DN" else ysrc + "（原产品单位）"),
                'absorption':references,
                "passive": True, "xlabel": xlabel, "no_annotate": True,
                "meta": dict(point.meta, sclk=point.sclk, kind=ptype,spectral_column=ysrc)}

    def _process(self, x, shots, shot, steps_on, log, point):
        from . import base as _b
        original_shots = shots.shape[1]
        if steps_on.get("dust"):
            shots = _b.drop_dust(shots, 5, log)
        if steps_on.get("despike"):
            shots = _b.despike_shots(shots, log)
        n_shots = shots.shape[1]
        shot_numbers = list(range(original_shots - n_shots + 1, original_shots + 1))
        if shot is None or shot < 0 or shot >= n_shots:
            y = shots.mean(axis=1)
            tag = "总谱（%d 发平均）" % n_shots
        else:
            y = shots[:, shot]
            tag = "原始第 %d 发" % shot_numbers[shot]
        if steps_on.get("baseline"):
            y = _b.poly_baseline(x, y, log=log)
        log.append(("显示", tag))
        return {"x": x, "y": y, "log": log, "n_shots": n_shots,
                "shot_matrix":np.column_stack([_b.poly_baseline(x,col) for col in shots.T]) if steps_on.get('baseline') else shots.copy(),
                "shot_numbers": shot_numbers,
                "meta": dict(point.meta, sclk=point.sclk)}

    # -- 官方成分预测（MOC 表） ----------------------------------------------

    _moc_table = None            # (rows, meta)，进程内只解析一次

    def moc_table(self, path=None):
        p = path or moc_cache_path()
        if self._moc_table is None or getattr(self, "_moc_path", None) != p:
            if not os.path.exists(p):
                return None, None
            self._moc_table = load_moc(p)
            self._moc_path = p
        return self._moc_table

    def ensure_moc(self, force=False):
        """下载官方成分表到缓存（已存在且不强制则直接用）。返回缓存路径。"""
        p = moc_cache_path()
        if os.path.exists(p) and not force:
            return p
        os.makedirs(os.path.dirname(p), exist_ok=True)
        from . import pds_listing
        data = pds_listing.http_get(MOC_URL, timeout=120)
        tmp = p + ".part"
        with open(tmp, "wb") as f:
            f.write(data)
        try:
            rows, _ = load_moc(tmp)
            if not rows or not any(any(v[0] is not None for v in row["oxides"].values()) for row in rows.values()):
                raise ValueError("下载内容没有有效成分记录")
        except Exception:
            os.remove(tmp)
            raise
        os.replace(tmp, p)
        self._moc_table = None           # 重新解析
        return p

    def moc_payload(self, point_stem):
        """给界面/报告的成分负载数据。表没下载或该点无预测都明确说出来。"""
        rows, meta = self.moc_table()
        if rows is None:
            return {"have_table": False}
        rec = rows.get(moc_key(point_stem))
        if rec is None:
            return {"have_table": True, "matched": False, "table_rows": len(rows)}
        out_rows = [{"oxide": n, "wt": rec["oxides"][n][0], "err": rec["oxides"][n][1]}
                    for n in MOC_OXIDES]
        known = [r for r in out_rows if r["wt"] is not None]
        return {
            "have_table": True, "matched": True, "table_rows": len(rows),
            "rows": out_rows,
            "total": rec["total"], "distance": rec["distance"],
            "total_known": round(sum(r["wt"] for r in known), 2) if known else None,
            "rmsep": meta.get("RMSEP", {}),
            "note": "PDS 官方模型预测值（supercam_libs_moc），非本软件计算；"
                    "FeOT 为全铁。RMSEP 为官方表给出的模型预测误差指标，与逐点标准差不同。",
        }

    # -- 自研 PLS 定量（Phase 2；与 ChemCam quant_payload 同构） -------------

    def quant_payload(self, obs):
        """观测级自研 PLS 成分（SCCT 训练模型；口径与边界同 ChemCam）。"""
        from . import base as _b
        try:
            import project_layout
            model_path = (project_layout.data_root() / "SuperCam" / "缓存"
                          / "pls_model.json")
        except Exception:                                # noqa: BLE001
            return {"available": False, "reason": "数据根不可用"}
        if not model_path.is_file():
            return {"available": False,
                    "reason": "模型未训练：先运行 SuperCam 数据准备脚本 01 --fit-model"}
        from instruments import libs_quant
        doc = libs_quant.load_model(model_path)
        if doc is None:
            return {"available": False, "reason": "模型文件损坏，请重新训练"}
        spectra = []
        used = 0
        for pt in obs.points:
            path = pt.files.get("cl1")
            if not path:
                continue
            try:
                d = self.load_point(pt, None, {"dust": True, "despike": True,
                                               "baseline": True})
            except Exception:                            # noqa: BLE001
                continue
            spectra.append(np.asarray(d["y"], dtype=float))
            used += 1
        if not used:
            return {"available": False, "reason": "本观测没有可读的主动 CL1 点"}
        try:
            values = libs_quant.predict_mean(np.mean(spectra, axis=0), doc)
        except ValueError as exc:
            return {"available": False, "reason": str(exc)}
        official = self._official_moc_for(obs)
        return {"available": True, "values": values, "official": official,
                "n_points_used": used,
                "scope": doc["scope"],
                "trained_targets": doc["meta"].get("trained_targets", [])}

    def _official_moc_for(self, obs):
        """官方 MOC 查本观测的成分预测（对账列）；缓存缺失/未命中返回 None。"""
        try:
            import project_layout
            moc_p = project_layout.legacy_cache("supercam_moc.csv")
            if moc_p is None or not os.path.isfile(str(moc_p)):
                moc_p = (project_layout.data_root() / "SuperCam" / "随包缓存"
                         / "supercam_moc.csv")
            if not moc_p.is_file():
                return None
            want = set()
            for pt in obs.points:
                for f in pt.files.values():
                    # 文件名带版本后缀（_14p01），MOC 行无版本（_14p）——去尾段匹配
                    stem = re.sub(r"_\d+p\d*\.fits$", "", os.path.basename(f).upper(), flags=re.I)
                    want.add(stem)
            if not want:
                return None
            hits = []
            for line in moc_p.read_text(encoding="utf-8", errors="replace").splitlines():
                parts = [p.strip() for p in line.split(",")]
                if len(parts) < 18 or not parts[0].startswith("scam_"):
                    continue
                if re.sub(r"_\d+p\d*$", "", parts[0].upper(), flags=re.I) in want:
                    try:
                        hits.append({ox: [float(parts[1 + 2 * i]), float(parts[2 + 2 * i])]
                                     for i, ox in enumerate(("SiO2", "TiO2", "Al2O3", "FeOT",
                                                             "MgO", "CaO", "Na2O", "K2O"))})
                    except (ValueError, IndexError):
                        continue
            return hits[0] if hits else None
        except Exception:
            import traceback
            traceback.print_exc()
            return None

    # -- 下载面板分组 --------------------------------------------------------

    def parse_listing(self, entries):
        groups = {}
        for e in entries:
            name = e["name"]
            ext = os.path.splitext(name)[1].lower()
            if ext not in (".fits", ".xml"):
                continue
            product = re.match(r"^scam_\d+_\d+_\d+_([a-z]+\d+)_", name, re.I)
            key = product.group(1).upper() if product else "OTHER"
            label = TYPE_LABELS.get(key, key + " · 其它")
            g = groups.setdefault(key, {"key": key, "label": label,
                                        "files": [], "count": 0, "bytes": 0})
            g["files"].append(e)
            g["count"] += 1
            g["bytes"] += e["size"]
        out = sorted(groups.values(), key=lambda g: (-g["bytes"], g["key"]))
        preset = [k for k in ("CL1",) if k in groups]
        return {"groups": out, "preset": preset,
                "total_files": sum(g["count"] for g in out),
                "total_bytes": sum(g["bytes"] for g in out)}
