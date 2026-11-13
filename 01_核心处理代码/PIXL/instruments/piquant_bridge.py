# -*- coding: utf-8 -*-
"""PIQUANT（PIXL 官方开源定量，github.com/pixlise/piquant）桥接层。

约定：PIQUANT 可执行文件路径由环境变量 MARSWORKBENCH_PIQUANT 或
数据目录 PIXL/缓存/piquant_path.json 指定；未配置/不可执行时所有入口
显式返回 {"available": False, reason}——绝不猜测路径、绝不阻塞工作台。
调用约定：piquant CLI 以 detector config + SFS csv 为输入，输出 JSON
逐元素 wt%；本桥只做 subprocess 编排与结果落盘，不复制其物理算法。

对账口径：PIQUANT 输出 vs 官方 RQC wt%（同 PMC 逐元素偏差）写入
PIXL/缓存/piquant_vs_rqc.json；PIQUANT 是官方开源代码，其输出代表
"官方算法复现"，与本仓库自研峰面积层（sfs_spectra）口径分开。
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path


def piquant_bin():
    """定位 PIQUANT 可执行文件；不可用返回 None。"""
    env = os.environ.get("MARSWORKBENCH_PIQUANT")
    if env and os.path.isfile(env):
        return env
    try:
        import project_layout
        cfg = (project_layout.data_root() / "PIXL" / "缓存" / "piquant_path.json")
        if cfg.is_file():
            doc = json.loads(cfg.read_text(encoding="utf-8"))
            p = doc.get("piquant", "")
            if p and os.path.isfile(p):
                return p
    except Exception:                                    # noqa: BLE001
        pass
    found = shutil.which("piquant")
    return found


def available():
    return piquant_bin() is not None


def status():
    p = piquant_bin()
    if not p:
        return {"available": False,
                "reason": "PIQUANT 未配置。设置环境变量 MARSWORKBENCH_PIQUANT 指向"
                          "可执行文件（github.com/pixlise/piquant，Go 构建），或写入"
                          " PIXL/缓存/piquant_path.json {\"piquant\": \"路径\"}。"
                          "官方开源定量与本仓库自研峰面积层口径分开。"}
    return {"available": True, "bin": p}


def run_on_sfs(sfs_csv, work_dir, config=None, timeout=600):
    """对一份 SFS csv 跑 PIQUANT，返回 (输出JSON路径 or None, reason)。

    config: PIQUANT detector 配置目录（必需参数，来自官方标定数据集）。
    """
    st = status()
    if not st["available"]:
        return None, st["reason"]
    if not config or not Path(config).is_dir():
        return None, "PIQUANT detector 配置目录缺失（官方标定数据集）"
    out_dir = Path(work_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_json = out_dir / (Path(sfs_csv).stem + "_piquant.json")
    try:
        proc = subprocess.run([piquant_bin(), "-c", str(config), "-i", str(sfs_csv),
                               "-o", str(out_json)],
                              capture_output=True, text=True, timeout=timeout)
        if proc.returncode != 0 or not out_json.is_file():
            return None, "PIQUANT 退出码 %d：%s" % (proc.returncode,
                                                    (proc.stderr or "")[:120])
        return str(out_json), None
    except (subprocess.TimeoutExpired, OSError) as exc:
        return None, str(exc)[:120]
