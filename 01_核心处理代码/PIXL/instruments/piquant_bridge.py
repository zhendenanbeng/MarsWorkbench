# -*- coding: utf-8 -*-
"""PIQUANT（PIXL 官方开源定量，github.com/pixlise/piquant）桥接层。

约定：PIQUANT 可执行文件路径由环境变量 MARSWORKBENCH_PIQUANT 或
数据目录 PIXL/缓存/piquant_path.json 指定；未配置/不可执行时所有入口
显式返回 {"available": False, reason}——绝不猜测路径、绝不阻塞工作台。
调用约定（v2，2026-11-19 对齐官方 CLI）：PIQUANT 是 C++/CMake 项目
（MSVC/NMake 本地构建），quant 子命令为位置参数：
    Piquant quant <config.msa> <calibration.csv> <spectrum.msa> "<元素表>" <plot.csv>
wt% 定量表写在 stdout（无 JSON 输出），本桥解析 stdout 定量行。
谱输入用 MSA（EMSA/MAS）格式：sfs_to_msa() 把自研管线的逐 PMC 谱导出。

对账口径：PIQUANT（官方算法复现，FM 配置+FM 标定）vs 官方 RQC wt%
（同 PMC 逐元素偏差）由 数据准备脚本/02_piquant_quantify.py 写入
PIXL/缓存/piquant_vs_rqc.json；与本仓库自研峰面积层（sfs_spectra）
口径分开。FM 标定文件官方 COMMENT 明示用于 ENB（B 探测器）谱。
"""
from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
from pathlib import Path


def _config_doc():
    """读 piquant_path.json（bin/config/calibration/testdata 路径注册）。"""
    try:
        import project_layout
        cfg = (project_layout.data_root() / "PIXL" / "缓存" / "piquant_path.json")
        if cfg.is_file():
            return json.loads(cfg.read_text(encoding="utf-8"))
    except Exception:                                    # noqa: BLE001
        pass
    return None


def piquant_bin():
    """定位 PIQUANT 可执行文件；不可用返回 None。"""
    env = os.environ.get("MARSWORKBENCH_PIQUANT")
    if env and os.path.isfile(env):
        return env
    doc = _config_doc()
    if doc:
        p = doc.get("piquant", "")
        if p and os.path.isfile(p):
            return p
    return shutil.which("piquant")


def available():
    return piquant_bin() is not None


def fm_config():
    doc = _config_doc() or {}
    p = doc.get("config", "")
    return p if p and os.path.isfile(p) else None


def fm_calibration():
    doc = _config_doc() or {}
    p = doc.get("calibration", "")
    return p if p and os.path.isfile(p) else None


def testdata_dir():
    doc = _config_doc() or {}
    p = doc.get("testdata_dir", "")
    return p if p and os.path.isdir(p) else None


def status():
    p = piquant_bin()
    if not p:
        return {"available": False,
                "reason": "PIQUANT 未配置。设置环境变量 MARSWORKBENCH_PIQUANT 指向"
                          "可执行文件（github.com/pixlise/piquant，C++/CMake 构建），"
                          "或写入 PIXL/缓存/piquant_path.json {\"piquant\": \"路径\"}。"
                          "官方开源定量与本仓库自研峰面积层口径分开。"}
    missing = [n for n, f in (("config", fm_config()), ("calibration", fm_calibration()))
               if not f]
    if missing:
        return {"available": False,
                "reason": "PIQUANT 可执行文件在位但缺 %s（PIXL FM 配置/标定，"
                          "见 piquant_path.json 的 config/calibration 字段）" % "/".join(missing)}
    return {"available": True, "bin": p, "config": fm_config(),
            "calibration": fm_calibration()}


# ---- MSA（EMSA/MAS）谱导出 ------------------------------------------------
# 官方 spectrum.md：MSA 头关键字 + #SPECTRUM : + 逐通道计数 + #ENDOFDATA。
# XPERCHAN/OFFSET 写入本扫描自证轴（calibrate_energy v4），PIQUANT 优先
# 读谱文件自身标定；LIVETIME 用该 PMC 实测值。
MSA_HEADER = """#FORMAT : EMSA/MAS spectral data file
#VERSION : TC202v2.0 PIXL
#TITLE : MarsWorkbench SFS export ({title})
#OWNER : MarsWorkbench
#DATE : {date}
#NPOINTS : {npoints}
#NCOLUMNS : 1
#XUNITS : eV
#YUNITS : COUNTS
#DATATYPE : YY
#XPERCHAN : {xperchan:.5f} eV per channel
#OFFSET : {offset:.4f} eV of first channel
#SIGNALTYPE : XRF
#LIVETIME : {livetime:.2f}
#REALTIME : {realtime:.2f}
#SPECTRUM :
"""


def sfs_to_msa(spec_counts, livetime_s, axis, out_path, title="PMC"):
    """把一份 4096 通道计数谱写成 PIQUANT 可读 MSA 文件，返回路径。

    spec_counts: 原始计数（PIQUANT 自己做活时间处理，喂计数不喂 cps）；
    livetime_s: 活时间秒；axis: sfs_spectra.calibrate_energy() 结果 dict。
    """
    spec = list(spec_counts)
    if len(spec) != 4096:
        raise ValueError("SFS 谱长度 %d ≠ 4096" % len(spec))
    livetime_s = float(livetime_s) if livetime_s and livetime_s > 0 else 1.0
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    # MSA 文件为纯 ASCII：非 ASCII 字符（如中文标题）一律替换为下划线
    title = "".join(ch if ord(ch) < 128 else "_" for ch in str(title))
    date = __import__("datetime").datetime.now().strftime("%m-%d-%Y  %H:%M:%S")
    with open(out, "w", encoding="ascii", newline="\n") as fh:
        fh.write(MSA_HEADER.format(
            title=title, date=date, npoints=len(spec),
            xperchan=axis["keV_per_ch"] * 1000.0,
            offset=axis["zero_keV"] * 1000.0,
            livetime=livetime_s, realtime=livetime_s))
        # NCOLUMNS=1：每行一个计数值（read_EMSA_PIXL.cpp 逐行取前
        # numColumns 个值作计数——带通道列会把通道号当计数读入）
        for c in spec:
            fh.write("%.1f\n" % float(c))
        fh.write("#ENDOFDATA : \n")
    return str(out)


# ---- stdout wt% 定量表解析 ------------------------------------------------
# quant 输出形如 "    SiO2     18.87 %      Si_K   int ... coeff ... re_c ..."，
# 氧化物行（xxO/xO-y）与元素行（单元素名）混排；抓 "<名> <数> %" 行，
# 按官方氧化物/元素白名单归类，另抓 "Element sum" 与 O（化学计量）。
WT_ROW = re.compile(
    r"^\s{0,6}(?P<name>[A-Z][a-zA-Z0-9]{0,5}(?:-[A-Z])?)\s+(?P<wt>\d+\.?\d*|\.\d+)\s*%"
    r"(?P<rest>.*)$", re.X)
OXIDE_OF = {"SiO2": "SiO2", "K2O": "K2O", "P2O5": "P2O5", "CaO": "CaO",
            "TiO2": "TiO2", "Cr2O3": "Cr2O3", "MnO": "MnO", "FeO-T": "FeO-T",
            "MgO": "MgO", "Al2O3": "Al2O3", "Na2O": "Na2O", "SO3": "SO3",
            "Cl": "Cl"}
ELEMENT_OF = {"Sr": "Sr", "Ar": "Ar", "Ni": "Ni", "Zn": "Zn", "Y": "Y"}
SUM_RE = re.compile(r"Element sum\s+(\d+\.?\d*)\s*%")
O_RE = re.compile(r"^\s+O\s+(\d+\.?\d*)\s*%\s+\(\s*(\d+\.?\d*)\s*% from stoichiometry",
                  re.M)


def parse_quant_stdout(stdout):
    """从 quant stdout 抽 wt% 表：{"oxides": {...}, "elements": {...},
    "element_sum_pct": x, "o_stoich_pct": y}。"""
    oxides, elements = {}, {}
    for line in stdout.splitlines():
        m = WT_ROW.match(line)
        if not m:
            continue
        name = m.group("name")
        try:
            wt = float(m.group("wt"))
        except ValueError:
            continue
        if name in OXIDE_OF:
            oxides[OXIDE_OF[name]] = wt
        elif name in ELEMENT_OF:
            elements[ELEMENT_OF[name]] = wt
    out = {"oxides": oxides, "elements": elements}
    sm = SUM_RE.search(stdout)
    if sm:
        out["element_sum_pct"] = float(sm.group(1))
    om = O_RE.search(stdout)
    if om:
        out["o_stoich_pct"] = float(om.group(2))
    return out


def run_quant(spectrum_msa, out_dir, elements, config=None, calibration=None,
              timeout=300):
    """对一份 MSA 谱跑 PIQUANT quant，返回 (定量 dict or None, reason)。

    elements: 元素线表（如 "Si_K K_K P_K Ca_K Ti_K Cr_K Mn_K Fe_K Sr_K Ar_I"，
    官方 FM 标定覆盖范围）；config/calibration 缺省取 piquant_path.json 的 FM 文件。
    """
    st = status()
    if not st["available"]:
        return None, st["reason"]
    cfg = config or fm_config()
    cal = calibration or fm_calibration()
    if not cfg or not cal:
        return None, "PIQUANT FM 配置/标定缺失（piquant_path.json）"
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_csv = out_dir / (Path(spectrum_msa).stem + "_plot.csv")
    try:
        proc = subprocess.run(
            [piquant_bin(), "quant", str(cfg), str(cal), str(spectrum_msa),
             elements, str(plot_csv)],
            capture_output=True, text=True, timeout=timeout,
            errors="replace", cwd=str(out_dir))
    except (subprocess.TimeoutExpired, OSError) as exc:
        return None, str(exc)[:160]
    quant = parse_quant_stdout(proc.stdout or "")
    if proc.returncode != 0 or not quant.get("oxides"):
        return None, ("PIQUANT 退出码 %d：%s" % (proc.returncode,
                     (proc.stderr or proc.stdout or "")[-160:]))
    quant["plot_csv"] = str(plot_csv)
    quant["config"] = cfg
    quant["calibration"] = cal
    return quant, None
