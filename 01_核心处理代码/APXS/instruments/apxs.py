# -*- coding: utf-8 -*-
"""好奇号 APXS α 粒子 X 射线光谱后端。

数据：MSL-M-APXS-4/5-RDR-V1.0（PDS3）
  文件名：apa_406545731rsp01020050178_______p1.{csv,lbl,xml}
          └仪┘└──SCLK──┘└类型┘                    └版
  rsp = 标定能谱（1024 通道，能量轴 eV 就在 csv 里）；
  rwp = 氧化物成分表（16 项 wt% + 统计误差）——APXS 的核心产出。

csv 第一行是 HEADER（目标名等备注），之后每行「能量 eV, 计数」。
首通道的值是总积分时间（秒）、末通道是溢出计数——都从谱里剔除，
不然首通道一个几千会把 y 轴完全压扁。

观测 = 同 SCLK 前缀的 rsp + rwp 一对；没有 shot 概念。
"""

from __future__ import annotations

import os
import re

import numpy as np

from . import base, lines, pds3

# 前缀 apa = 常规科学测量；apb = 背景测量（结构相同，Target 名为背景点）。
NAME_RX = re.compile(r"^ap[a-z]_(\d{9})(r[a-z]{2})", re.I)

TYPE_LABELS = {
    "RSP": "RSP · 标定能谱（看图用）",
    "RWP": "RWP · 氧化物成分表（wt%）",
    "RSE": "RSE · 工程数据",
}


def _pair_identity(stem):
    """RSP/RWP 按 SCLK 和文件其余标识配对；版本各取最高，不按序列猜配。"""
    key = re.sub(r"r(?:sp|wp)", "", stem.lower(), count=1)
    return re.sub(r"_p\d+$", "", key)


def _version(stem):
    m = re.search(r"_p(\d+)$", stem, re.I)
    return int(m.group(1)) if m else 0


class APXSBackend(base.Backend):
    id = "apxs"
    name = "APXS α 粒子 X 射线光谱"
    rover = "msl"
    rover_label = "好奇号"
    subtitle = "下载 · 看谱 · 元素峰识别 · 氧化物成分"
    xlabel = "能量 (keV)"
    ylabel = "计数"
    accent = "#b3541e"
    steps = [("smooth", "滑动平滑（7 点）", False)]
    pds_roots = {
        "calibrated": {
            "label": "标定数据 RDR（谱 + 成分表）",
            "url": "https://pds-geosciences.wustl.edu/msl/"
                   "msl-m-apxs-4_5-rdr-v1/mslapx_1xxx/data"},
    }

    # -- 本地扫描 -----------------------------------------------------------

    def scan_folder(self, root):
        obs_map = {}
        rwp_by_id = {}
        for path, fname, dirname in self.walk_files(root, r".*\.csv$"):
            m = NAME_RX.match(fname)
            if not m:
                continue
            sclk, ptype = m.group(1), m.group(2).lower()
            stem = os.path.splitext(fname)[0]
            lbl = os.path.join(os.path.dirname(path), stem + ".lbl")
            sol, seq = None, ""
            lab = None
            if os.path.exists(lbl):
                lab = pds3.parse_lbl(lbl)
                sol = lab.get_int("PLANET_DAY_NUMBER")
                seq = (lab.get("SEQUENCE_ID") or "").strip().lower()
            if sol is None:
                sol = self.sol_from_dirname(dirname)
            if sol is None:
                continue
            sol = int(sol)
            if not seq:
                seq = "apa%s" % sclk
            # 跨仪器对照元数据：目标名（rsp 首行 Target:）+ 停车位（label RMC）
            sd = pds3.site_drive_from_lbl(lab)
            header_target = None
            if ptype == "rsp":
                try:
                    with open(path, encoding="utf-8", errors="replace") as f:
                        header_target = pds3.target_from_apxs_header(f.readline())
                except OSError:
                    header_target = None
            if ptype == "rsp":
                obs = obs_map.setdefault(
                    (sol, seq), base.Observation(sol, seq, "RSP"))
                if header_target and not getattr(obs, "target", ""):
                    obs.target = header_target
                if sd and getattr(obs, "site", None) is None:
                    obs.site, obs.drive = sd
                obs.points.append(base.PointFile(
                    stem, {"rsp": path}, sclk=sclk, point_no=1,
                    meta={"seq": seq, "target": header_target or ""}))
            elif ptype == "rwp":
                identity = (sol, seq, _pair_identity(stem))
                previous = rwp_by_id.get(identity)
                if previous is None or _version(stem) > _version(os.path.splitext(previous)[0]):
                    rwp_by_id[identity] = path
        out = []
        for key, obs in sorted(obs_map.items()):
            best = {}
            for pt in obs.points:
                identity = _pair_identity(pt.stem)
                if identity not in best or _version(pt.stem) > _version(best[identity].stem):
                    best[identity] = pt
            obs.points = [best[k] for k in sorted(best)]
            obs.rwp_files = []
            for i, pt in enumerate(obs.points):
                pt.point_no = i + 1
                rwp = rwp_by_id.get((*key, _pair_identity(pt.stem)))
                if rwp:
                    pt.files["rwp"] = rwp
                    obs.rwp_files.append(rwp)
            out.append(obs)
        return out

    # -- 读点 ---------------------------------------------------------------

    # -- 自研峰面积层（Phase 4 产品化；原始信号量，不是浓度） -------------

    def quant_payload(self, obs):
        """观测级自研峰面积 payload（供工作台 quant 面板，APXS 口径）。

        读观测全部 rsp 能谱（与 load_point 同构口径：剔首末通道）→
        逐元素 Kα 高斯拟合峰高/净面积+误差+检出判定。无 rsp 或拟合
        失败显式降级，绝不产出虚构峰面积。
        """
        from . import apxs_fit
        spectra, live_times = [], []
        for pt in obs.points:
            path = pt.files.get("rsp")
            if not path:
                continue
            try:
                data = pds3.read_table(path, skiprows=1)
            except Exception:                            # noqa: BLE001
                continue
            if data.ndim != 2 or data.shape[1] < 2:
                continue
            energy_kev = data[:, 0] / 1000.0
            counts = data[:, 1] if data.shape[1] == 2 else np.median(
                data[:, 1:-1], axis=1)
            live_times.append(float(counts[0]))
            energy_kev, counts = energy_kev[1:-1], counts[1:-1]
            spectra.append(np.asarray(counts, dtype=float))
        if not spectra:
            return {"available": False, "reason": "本观测没有可读的 rsp 能谱"}
        mean = np.mean(spectra, axis=0)
        fits = apxs_fit.spectrum_areas(energy_kev, mean)
        areas = {el: r["area"] for el, r in fits.items()
                 if r["area"] is not None}
        rows = []
        for el, r in fits.items():
            if r["area"] is None:
                continue
            rows.append({"label": el, "value": float(r["amp"]),
                         "err": float(r["error"]),
                         "detected": bool(r["detected"])})
        # FP-lite 自算氧化物 wt%（原始谱独立定量；标定不在册时显式降级）
        official = None
        for pt in obs.points:
            rwp_path = pt.files.get("rwp")
            if rwp_path:
                try:
                    import csv as _csv
                    official = {}
                    for r in _csv.reader(open(rwp_path, encoding="utf-8",
                                              errors="replace")):
                        if len(r) >= 2 and r[0].strip():
                            try:
                                official[r[0].strip()] = float(r[1])
                            except ValueError:
                                continue
                except Exception:                        # noqa: BLE001
                    official = None
                break
        fp = None
        try:
            from . import apxs_fp
            fp = apxs_fp.predict(areas)
            if live_times:
                fp["live_median_s"] = float(np.median(live_times))
            if fp.get("available") and official:
                comp = []
                for ox in apxs_fp.OXIDE_ORDER:
                    pv, ov = fp["oxides"].get(ox), official.get(ox)
                    if pv is None and ov is None:
                        continue
                    rel = ((pv - ov) / ov) if (pv is not None and ov) else None
                    comp.append({"oxide": ox, "fp": pv, "official": ov,
                                 "rel": rel})
                fp["compare"] = comp
        except Exception as exc:                         # noqa: BLE001
            fp = {"available": False, "reason": "FP 层异常：%s" % str(exc)[:60]}
        return {"available": True, "mode": "apxs_peaks", "rows": rows,
                "official": official,
                "fp": fp,
                "n_points_used": len(spectra),
                "scope": "APXS 峰拟合层：元素 Kα 峰高/净面积（原始信号量，"
                         "不是浓度）。自算定量见 FP 层；官方 rwp 氧化物表"
                         "仍是定量引用口径。"}

    def load_point(self, point, shot, steps_on):
        path = point.files["rsp"]
        data = pds3.read_table(path, skiprows=1)      # 首行是备注 HEADER
        log = [("读取", "%s · %d 通道" % (os.path.basename(path), len(data)))]
        note = ""
        try:
            with open(path, encoding="utf-8", errors="replace") as f:
                first = f.readline().strip().strip('"').strip(",")
            note = first if first else ""
        except OSError:
            pass
        x = data[:, 0] / 1000.0                       # eV → keV
        y = data[:, 1]
        log.append(("通道", "已剔除首通道（总积分时间 %.0f s）与末通道（溢出计数）"
                    % y[0]))
        x, y = x[1:-1], y[1:-1]
        if steps_on.get("smooth"):
            y = base.smooth(y, 7)
            log.append(("平滑", "7 点滑动平均"))
        meta = dict(point.meta, sclk=point.sclk)
        if note:
            meta["note"] = note
        return {"x": x, "y": y, "log": log, "n_shots": 1, "meta": meta}

    def annotate(self, x, y):
        return lines.match_lines(x, y, lines.APXS_LINES, tol=0.09, min_lines=1, wavelength_frame='energy_keV')

    # -- 氧化物成分表 --------------------------------------------------------

    def oxides(self, obs, point=None):
        """配对的 rwp csv → 成分行。返回 None 表示这次观测没有成分表。"""
        if point is None:
            # 旧调用只在单点观测时无歧义；多点时不能默认借第一份成分表。
            if len(obs.points) != 1:
                return None
            point = obs.points[0]
        path = point.files.get("rwp")
        if not path:
            return None
        rows, note = [], ""
        try:
            with open(path, encoding="utf-8", errors="replace") as f:
                note = f.readline().strip()
                for line in f:
                    parts = [p.strip().strip('"') for p in line.split(",")]
                    if len(parts) < 3:
                        continue
                    try:
                        wt = float(parts[1])
                        err = float(parts[2])
                    except ValueError:
                        continue
                    rows.append({"oxide": parts[0], "wt": wt, "err": err})
        except OSError:
            return None
        if not rows:
            return None
        return {"file": os.path.basename(path), "rows": rows, "note": note,
                "total": round(sum(r["wt"] for r in rows
                                   if not r["oxide"].strip().isdigit()), 2)}

    # -- 下载面板分组 --------------------------------------------------------

    def parse_listing(self, entries):
        groups = {}
        for e in entries:
            name = e["name"]
            ext = os.path.splitext(name)[1].lower()
            if ext not in (".csv", ".lbl", ".xml"):
                continue
            m = NAME_RX.match(name)
            key = m.group(2).upper() if m else "OTHER"
            label = TYPE_LABELS.get(key, key + " · 其它")
            g = groups.setdefault(key, {"key": key, "label": label,
                                        "files": [], "count": 0, "bytes": 0})
            g["files"].append(e)
            g["count"] += 1
            g["bytes"] += e["size"]
        out = sorted(groups.values(), key=lambda g: (-g["bytes"], g["key"]))
        preset = [k for k in ("RSP", "RWP") if k in groups]
        return {"groups": out, "preset": preset,
                "total_files": sum(g["count"] for g in out),
                "total_bytes": sum(g["bytes"] for g in out)}
