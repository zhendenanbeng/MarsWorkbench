# -*- coding: utf-8 -*-
"""好奇号 ChemCam LIBS 后端。

数据：MSL-M-CHEMCAM-LIBS-4/5-RDR-V1.0（PDS3）
  文件名：cl5_398645626ccs_f0030004ccam02013p3.{csv,lbl}
          └级别┘└─SCLK─┘└类型┘└请求号┘└序列号┘└产品版
  cl5 = 完全标定；ccs / rdr = 主动激光光谱（两种刻度）；
  CL1/CL9 psv = 被动光谱（无激光，配官方波长表看图）；cr0 prc = RMI 远景图。

CSV：若干行 '#' 注释头 + 6144 行 × [WAVELENGTH, SPECTRA_001..NNN,
MEDIAN, AVERAGE]——波长轴就在文件里，每列一发激光。
sol 不在文件名里：优先从 .lbl 的 PLANET_DAY_NUMBER 取，没有 .lbl
就认所在目录（下载时保留 PDS 的 solXXXXX 目录结构）。
"""

from __future__ import annotations

import os
import re
from pathlib import Path

import numpy as np

from . import base, pds3

NAME_RX = re.compile(
    r"^(cl[0-9])_(\d{9})([a-z]{3})_(\w+?)(ccam\d{5})([a-z]\d)\.(csv|tab)$", re.I)
RMI_RX = re.compile(r"^cr[0-9]_\d{9}prc_.*\.tif$", re.I)

TYPE_LABELS = {
    ("CL5", "CCS"): "CL5 CCS · 主动激光光谱（标定，计数刻度）",
    ("CL5", "RDR"): "CL5 RDR · 主动激光光谱（标定，辐射亮度单位）",
    ("CL0", "PSV"): "CL0 PSV · 被动光谱（未标定，本版只下载不看图）",
    ("CL1", "PSV"): "CL1 PSV · 被动光谱（官方波长轴，DN）",
    ("CL9", "PSV"): "CL9 PSV · 被动光谱（官方波长轴，DN）",
    ("CR0", "PRC"): "CR0 PRC · RMI 远景相机图像",
    ("CR1", "PRC"): "CR1 PRC · RMI 远景相机图像",
}
PAIR_EXTS = (".csv", ".tab", ".lbl", ".tif", ".xml")


class ChemCamBackend(base.Backend):
    id = "chemcam"
    name = "ChemCam LIBS 激光击穿光谱"
    rover = "msl"
    rover_label = "好奇号"
    subtitle = "下载 · 看谱 · 元素线识别"
    xlabel = "波长 (nm)"
    ylabel = "强度"
    accent = "#b3541e"
    wavelength_path = None  # 可指定已冻结的离线PDS波长表；默认下载到.cache/。
    steps = [("dust", "丢弃前 5 发（除尘）", True),
             ("despike", "去尖峰（跨发 MAD）", True),
             ("baseline", "连续背景扣除", False)]
    pds_roots = {
        "calibrated": {
            "label": "标定数据 RDR（看图用，推荐）",
            "url": "https://pds-geosciences.wustl.edu/msl/"
                   "msl-m-chemcam-libs-4_5-rdr-v1/mslccm_1xxx/data"},
    }

    # -- 本地扫描 -----------------------------------------------------------

    def scan_folder(self, root):
        obs_map = {}
        for path, fname, dirname in self.walk_files(root, r".*\.(csv|tab)$"):
            m = NAME_RX.match(fname)
            if not m:
                continue
            level, sclk, ptype = m.group(1).upper(), m.group(2), m.group(3).lower()
            passive = ptype=='psv' and level in ('CL1','CL9') and m.group(7).lower()=='tab'
            if not passive and not (ptype in ('ccs','rdr') and m.group(7).lower()=='csv'):continue
            seq = m.group(5).lower()
            stem = os.path.splitext(fname)[0]
            lbl = os.path.join(os.path.dirname(path), stem + ".lbl")
            sol = None
            lab = None
            if os.path.exists(lbl):
                lab = pds3.parse_lbl(lbl)
                sol = lab.get_int("PLANET_DAY_NUMBER")
                seq = (lab.get("SEQUENCE_ID") or "").strip().lower() or seq
            if sol is None:
                sol = self.sol_from_dirname(dirname)
            if sol is None:
                continue
            sol = int(sol)  # LBL与目录兜底必须同类型，避免混合目录排序失败。
            # 跨仪器对照元数据：目标名（csv 注释头）+ 停车位（label RMC）
            target = pds3.target_from_chemcam_csv(path)
            sd = pds3.site_drive_from_lbl(lab)
            key = (sol, seq, level+'-PSV' if passive else 'LIBS')
            obs = obs_map.get(key)
            if obs is None:
                obs = obs_map[key] = base.Observation(sol, seq, "CL5-" + ptype.upper())
                obs.passive = passive
                obs.target = target or ""
                obs.site, obs.drive = (sd if sd else (None, None))
            if target and not obs.target:
                obs.target = target
            obs.points.append(base.PointFile(
                stem, {ptype: path}, sclk=sclk,
                point_no=len(obs.points) + 1,
                meta={"level": level, "seq": seq, "passive": passive,
                      "target": target or ""}))

        out = []
        for (sol, seq, kind), obs in sorted(obs_map.items()):
            # 同一 (sol, seq, sclk) 有 ccs 与 rdr 两个文件 = 同一个点的两种刻度
            by_sclk = {}
            metadata = {}
            for pt in obs.points:
                by_sclk.setdefault(pt.sclk, {}).update(pt.files)
                metadata[pt.sclk]=pt.meta
            pts = []
            for i, (sclk, files) in enumerate(sorted(by_sclk.items())):
                prod = kind if obs.passive else ("RDR" if "rdr" in files else "CCS")
                pts.append(base.PointFile(
                    "%s_%s" % (seq, sclk), files, sclk=sclk, point_no=i + 1,
                    meta=dict(metadata[sclk],product=prod)))
            obs.points = pts
            if pts:
                obs.product = pts[0].meta.get("product", "?")
                out.append(obs)
        from . import rmi
        images=rmi.inventory(self,root)
        for obs in out:obs.rmi_images=[r for r in images if r['sol']==str(obs.sol) and r['seq']==obs.seq]
        return out

    # -- 读点 ---------------------------------------------------------------

    def load_point(self, point, shot, steps_on):
        files = point.files
        if 'psv' in files:
            from . import chemcam_passive
            path=files['psv'];label=os.path.splitext(path)[0]+'.lbl'
            wave=self.wavelength_path
            if wave is None:
                for folder in Path(path).resolve().parents[:2]:
                    candidate=folder/'document/ccam_default_wave.tab'
                    if candidate.is_file():wave=candidate;break
            x,y,meta=chemcam_passive.read_spectrum(path,label,wave)
            log=[('读取',os.path.basename(path)+' · 单列平均被动谱，无激光发数'),
                 ('波长','PDS官方6444通道表；剔除300个无效边缘通道，波长与强度同步重排'),
                 ('模式','被动谱不做除尘、跨发去尖峰或元素发射线判定')]
            if steps_on.get('baseline'):
                y=base.poly_baseline(x,y,log=log)
            return {'x':x,'y':y,'log':log,'n_shots':0,'shot_numbers':[],
                    'passive':True,'xlabel':'波长 (nm)','ylabel':'平均强度 (DN)',
                    'meta':dict(point.meta,sclk=point.sclk,**meta)}
        prod = "rdr" if "rdr" in files else "ccs"
        path = files[prod]
        data = pds3.read_table(path)
        if data.ndim != 2 or data.shape[1] < 3:
            raise ValueError("CSV 结构不对（%d 列，需要波长+发数+统计）"
                             % (data.shape[1] if data.ndim == 2 else -1))
        x = data[:, 0]
        shots = data[:, 1:-2]                 # 丢掉末尾 MEDIAN/AVERAGE 两列
        log = [("读取", "%s · %d 通道 × %d 发（%s 刻度）"
                % (os.path.basename(path), len(x), shots.shape[1],
                   "辐射亮度" if prod == "rdr" else "计数"))]
        return self._process(x, shots, shot, steps_on, log, point)

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

    # -- 自研 PLS 定量（Phase 1 产品化；独立标定，不是官方成分替代） --------

    def quant_payload(self, obs):
        """观测级自研 PLS 成分 payload（供工作台 quant 面板）。

        模型（pls_model.json，由数据准备脚本 03 --fit-model 生成）在位且
        观测含主动 CCS 点时，读全部点 → 与 load_point 同构的预处理链 →
        平均谱 → 逐氧化物预测（精度用留一靶 CV RMSE 标注）。模型缺失、
        通道不齐或文件读取失败时返回 {"available": False, ...}——显式
        降级，绝不产出虚构成分。这是独立复现参考，不是官方成分替代。
        """
        from . import base as _b
        try:
            import project_layout
            model_path = (project_layout.data_root() / "ChemCam" / "缓存"
                          / "pls_model.json")
        except Exception:                                # noqa: BLE001
            return {"available": False, "reason": "数据根不可用"}
        if not model_path.is_file():
            return {"available": False,
                    "reason": "模型未训练：先运行 ChemCam 数据准备脚本 03 --fit-model"}
        from . import libs_quant
        doc = libs_quant.load_model(model_path)
        if doc is None:
            return {"available": False, "reason": "模型文件损坏，请重新训练"}
        spectra = []
        used = 0
        for pt in obs.points:
            path = pt.files.get("ccs")
            if not path:
                continue
            try:
                data = pds3.read_table(path)
            except Exception:                            # noqa: BLE001
                continue
            if data.ndim != 2 or data.shape[1] < 3:
                continue
            wl = data[:, 0]
            shots = data[:, 1:-2]
            shots = _b.drop_dust(shots, 5)
            shots = _b.despike_shots(shots)
            cols = [_b.poly_baseline(wl, shots[:, j]) for j in range(shots.shape[1])]
            spectra.append(np.column_stack(cols).mean(axis=1))
            used += 1
        if not used:
            return {"available": False, "reason": "本观测没有可读的主动 CCS 点"}
        mean = np.mean(spectra, axis=0)
        try:
            values = libs_quant.predict_mean(mean, doc)
        except ValueError as exc:
            return {"available": False, "reason": str(exc)}
        return {"available": True, "values": values, "n_points_used": used,
                "scope": doc["scope"],
                "trained_targets": doc["meta"].get("trained_targets", [])}

    # -- 下载面板分组 --------------------------------------------------------

    def parse_listing(self, entries):
        groups, order = {}, []
        for e in entries:
            name = e["name"]
            ext = os.path.splitext(name)[1].lower()
            if ext not in PAIR_EXTS:
                continue
            m = re.match(r"^((?:cl|cr)\d)_\d{9}([a-z]{3})_", name, re.I)
            if not m:
                key, label = "OTHER", "其它文件"
            else:
                lvl, ptype = m.group(1).upper(), m.group(2).upper()
                key = "%s_%s" % (lvl, ptype)
                label = TYPE_LABELS.get((lvl, ptype), "%s %s · 其它" % (lvl, ptype))
            g = groups.setdefault(key, {"key": key, "label": label,
                                        "files": [], "count": 0, "bytes": 0})
            if not g["files"]:
                order.append(key)
            g["files"].append(e)
            g["count"] += 1
            g["bytes"] += e["size"]
        out = sorted(groups.values(), key=lambda g: (-g["bytes"], g["key"]))
        preset = [k for k in ("CL5_RDR", "CL5_CCS") if k in groups]
        return {"groups": out, "preset": preset,
                "total_files": sum(g["count"] for g in out),
                "total_bytes": sum(g["bytes"] for g in out)}
