# -*- coding: utf-8 -*-
"""生成五台仪器数据处理流程图 HTML（poster 单页，供 IAB 截图为 PNG）。"""
import json
from pathlib import Path

OUT = Path(r"C:/Users/unnat/AppData/Local/Temp/mw_papers/flowcharts")
OUT.mkdir(parents=True, exist_ok=True)

# 语义色（低饱和填充 + 语义描边）：蓝=数据/IO，绿=自研计算，琥珀=官方/文献对账基准，灰=输出/冻结
CSS = """
*{margin:0;padding:0;box-sizing:border-box}::-webkit-scrollbar{width:0;height:0}html{overflow:hidden}
body{margin:0;background:#fff;font-family:'Segoe UI','Microsoft YaHei',sans-serif;color:#1e293b}
.poster{width:1280px;padding:36px 44px 44px}
.head{border-bottom:3px solid #0f172a;padding-bottom:14px;margin-bottom:6px}
.head h1{font-size:26px;font-weight:700;letter-spacing:.5px}
.head h1 .en{color:#64748b;font-weight:600;font-size:15px;margin-left:10px}
.head .sub{font-size:13px;color:#475569;margin-top:6px;line-height:1.5}
.legend{display:flex;gap:14px;margin:14px 0 4px;font-size:12px;color:#475569;align-items:center}
.legend .chip{display:inline-flex;align-items:center;gap:6px}
.legend .sw{width:14px;height:14px;border-radius:3px;display:inline-block;border:1.5px solid}
.phase{margin-top:22px}
.phase .pt{display:flex;align-items:baseline;gap:10px;background:#f1f5f9;border-left:5px solid #0f172a;
  padding:7px 14px;font-size:16px;font-weight:700;border-radius:0 4px 4px 0}
.phase .pt .no{color:#64748b;font-size:13px;font-weight:600}
.phase .body{padding:14px 8px 0 26px}
.node{position:relative;background:#fff;border:1.5px solid #94a3b8;border-left-width:5px;
  border-radius:6px;padding:10px 16px;margin-bottom:0;max-width:900px}
.node .t{font-size:14.5px;font-weight:700;line-height:1.4}
.node .d{font-size:12.5px;color:#475569;line-height:1.55;margin-top:3px}
.node .m{font-size:11.5px;color:#64748b;font-family:Consolas,monospace;margin-top:4px}
.arrow{height:26px;display:flex;align-items:center;justify-content:flex-start;padding-left:60px}
.arrow .stem{width:2px;height:100%;background:#64748b}
.arrow{position:relative}
.arrow::after{content:'';position:absolute;left:55px;bottom:0;border:7px solid transparent;border-top-color:#64748b}
.row{display:flex;gap:14px;margin-top:0}
.row .node{flex:1;max-width:none;margin-bottom:0}
.rowwrap{position:relative;padding:0 8px 0 26px}
.rowwrap .arrow-m{height:26px;display:flex;justify-content:center}
.arrow-m .stem{width:2px;height:100%;background:#64748b}
.arrow-m{position:relative}
.arrow-m::after{content:'';position:absolute;left:50%;margin-left:-7px;bottom:0;border:7px solid transparent;border-top-color:#64748b}
.merge{font-size:12px;color:#475569;text-align:center;margin-top:2px}
.foot{margin-top:26px;border-top:1.5px solid #cbd5e1;padding-top:10px;font-size:12px;color:#475569;line-height:1.6}
.foot b{color:#0f172a}
.blue{background:#eff6ff;border-color:#3b82f6}
.green{background:#f0fdf4;border-color:#10b981}
.amber{background:#fff7ed;border-color:#f59e0b}
.slate{background:#f8fafc;border-color:#64748b}
"""

def node(t, d, m=None, cls="blue"):
    mm = '<div class="m">%s</div>' % m if m else ""
    return '<div class="node %s"><div class="t">%s</div><div class="d">%s</div>%s</div>' % (cls, t, d, mm)

def arrow():
    return '<div class="arrow"><div class="stem"></div></div>'

def arrow_m():
    return '<div class="arrow-m"><div class="stem"></div></div>'

def phase(no, title, body):
    return ('<div class="phase"><div class="pt"><span class="no">%s</span><span>%s</span></div>'
            '<div class="body">%s</div></div>' % (no, title, body))

def row(nodes):
    return '<div class="rowwrap"><div class="row">%s</div>%s</div>' % (
        "".join(nodes), arrow_m())

LEGEND = ('<div class="legend">'
          '<span class="chip"><span class="sw" style="background:#3b82f6"></span>数据获取 / IO</span>'
          '<span class="chip"><span class="sw" style="background:#10b981"></span>自研计算层</span>'
          '<span class="chip"><span class="sw" style="background:#f59e0b"></span>官方 / 文献对账基准</span>'
          '<span class="chip"><span class="sw" style="background:#64748b"></span>输出与冻结</span>'
          '</div>')

def page(fname, h1, en, sub, phases, foot):
    html = ('<!DOCTYPE html><html lang="zh"><head><meta charset="utf-8"><style>%s</style></head>'
            '<body><div class="poster"><div class="head"><h1>%s<span class="en">%s</span></h1>'
            '<div class="sub">%s</div></div>%s%s<div class="foot">%s</div></div></body></html>'
            % (CSS, h1, en, sub, LEGEND, "".join(phases), foot))
    (OUT / fname).write_text(html, encoding="utf-8")
    return fname

pages = []

# ================= SHERLOC =================
p = [
    phase("阶段 1", "数据获取（PDS 归档）",
          node("PDS 归档产品下载与原件校验",
               "SHERLOC RWS/RRS 光谱 + ACI/WATSON 图像；保留 PDS 原名与配对 label/XML，SHA-256 校验",
               "SHERLOC/数据准备脚本/01_read_rrs · 04_fetch_stowed_background · 05_build_stowed_background")),
    arrow(),
    phase("阶段 2", "读数与预处理（自研层）",
          node("产品解析与观测组织",
               "逐产品解析 → 光谱×图像配对成观测（同点多模式聚合）",
               "sherloc/reader.py · organize.py", "green")
          + arrow()
          + node("光谱预处理",
                 "波长轴规整、去噪与区域截取，为峰拟合供谱",
                 "sherloc/preprocess.py", "green")
          + arrow()
          + node("仪器本底扣除（可开关）",
                 "收纳臂暗天空本底谱按形状缩放相减；缩放 M 仅按 800 cm⁻¹ 主带窗估计（Jakubek 2024 口径），覆盖外通道不动",
                 "sherloc/background.py · 数据准备脚本 04/05", "green")),
    arrow(),
    phase("阶段 3", "峰提取与矿物判读（自研层）",
          node("峰提取与峰拟合",
               "逐谱提取拟合峰（峰位/宽度/强度）",
               "sherloc/peaks.py · fitting.py", "green")
          + arrow()
          + node("规则库峰位匹配",
                 "REFERENCE 峰位表匹配打分（库内是峰位表，不含完整代表谱）",
                 "sherloc/library.py", "green")
          + arrow()
          + row([node("矿物判定（规则标签代理）", "峰组 → 化学族/矿物候选；规则标签是代理，不是矿物真值", "sherloc/minerals.py", "green"),
                 node("逐点质控分组", "analyse_spectra 输出 mineral_groups 供报告与判读", "sherloc/qc.py", "green")])),
    arrow(),
    phase("阶段 4", "对账（文献 / 冻结基准）",
          row([node("文献逐靶对账",
                    "Sun 2023 / Corpolongo 2023 逐靶：14 类目 7 复现，3 类经 silicate 粗类映射命中",
                    "verify/paper_crosscheck_check.py", "amber"),
               node("谱库止损诊断（冻结基线）",
                    "Sol 留出 5 折代理 AUC ≤0.648 < 0.70 继续线 → close_no_heldout_path，跨树复现",
                    "verify/phase10_check.py · C_交付诊断.json", "amber")])
          + '<div class="merge">对账结论只进台账与报告注记，不改判读输出</div>'),
    arrow(),
    phase("阶段 5", "输出与复核",
          node("质控报告 / 导出 / 论文对照标注",
               "_qc 批量质控报告、PNG+当前谱 CSV 导出；逐点论文对照（一致/不一致/待查 + 备注）存 _notes，不改自动分类",
               "02_图形界面代码/view.py · SHERLOC 台", "slate")),
]
pages.append(page("f_sherloc.html", "SHERLOC 数据处理流程", "Scanning Habitable Environments with Raman & Luminescence for Organics & Chemicals",
                  "拉曼/荧光光谱 + ACI/WATSON 图像 · 自研判读层与官方/文献对账口径分开", p,
                  "<b>口径边界：</b>分类是判读对照工具，检查通过不代表矿物鉴定准确率；规则库为峰位表（无完整代表谱），判读不替代独立矿物鉴定。"))

# ================= SuperCam =================
p = [
    phase("阶段 1", "数据获取（PDS 归档）",
          node("校准谱 FITS 下载与解析",
               "data_calibrated_spectra 集合（7933 通道/发）；自研 FITS 解析器读谱，SCCT 在轨标定靶训练集组织",
               "SuperCam/数据准备脚本/01_run_pls_cv.py · instruments/fitsio.py")),
    arrow(),
    phase("阶段 2", "预处理（自研层，ChemCam 共享模块）",
          node("除尘 → 去尖峰 → 多项式基线 → 波长裁剪",
               "与画图同构的预处理链；保留发为逐发比例分母，发间独立性未被证明",
               "ChemCam/instruments/base.py（共享）", "green")),
    arrow(),
    phase("阶段 3", "自研定量（PLS，从谱级原始数据自算成分）",
          node("PLS1 回归（纯 numpy NIPALS）",
               "逐氧化物建模，auto 选成分数并随 JSON 持久化；输出 8 主氧化物预测",
               "ChemCam/instruments/libs_quant.py（共享）", "green")
          + arrow()
          + node("留一靶交叉验证（防泄漏）",
                 "按靶分组 LOO → RMSE/bias 逐氧化物冻结 scam_pls_cv.json（17 靶，确定性逐位复现）",
                 "libs_quant.loo_target_cv", "green")),
    arrow(),
    phase("阶段 4", "对账（官方 / 文献基准）",
          row([node("官方 MOC 表并排",
                    "官方成分列并排展示（官方口径），自研列口径分开标注",
                    "instruments/supercam.py · quant_ui", "amber"),
               node("文献区间判定（1/8 落带内）",
                    "TiO2 0.455 落官方 RMSEP 0.3 的 ±60% 带 [0.1,0.5]；锚点=Ytsma 2023 Table 1 转录；其余 7 氧化物带外如实入册",
                    "gold_intervals.json · moc_rmsep", "amber"),
               node("独立数据审计",
                    "留一靶 CV 从训练谱独立复算，与冻结 RMSE 逐位对账",
                    "verify/data_audit.py", "amber")])
          + '<div class="merge">对账结论只进台账/报告注记，不回写官方产品</div>'),
    arrow(),
    phase("阶段 5", "输出",
          node("PLS 成分面板 / CI-WATSON 图像浏览 / 冻结模型",
               "观测级成分（非逐点表）随 RMSE 标注；pls_model.json 冻结；CI/WATSON 图像浏览与质控",
               "02_图形界面代码/instruments/quant_ui.py · SuperCam 台", "slate")),
]
pages.append(page("f_supercam.html", "SuperCam 数据处理流程", "SuperCam LIBS + Passive Spectroscopy",
                  "自研 PLS 成分层（谱级原始数据自算）与官方 MOC/文献对账口径分开", p,
                  "<b>口径边界：</b>17 靶最小训练集精度弱于官方 334 标样为预期；官方 MOC 成分仍是定量引用口径，自研预测不回写官方产品。"))

# ================= ChemCam =================
p = [
    phase("阶段 1", "数据获取（PDS 归档）",
          node("CCS 谱按靶名下载",
               "官方 MOC 表按标定靶筛观测 → 下载 CCS 训练谱（6 在轨标定靶）",
               "ChemCam/数据准备脚本/04_download_ccs_by_target.py")),
    arrow(),
    phase("阶段 2", "预处理（自研层）",
          node("除尘 → 去尖峰 → 多项式基线 → 波长裁剪",
               "与画图同构的预处理链（SuperCam 共享同一模块）",
               "ChemCam/instruments/base.py", "green")),
    arrow(),
    phase("阶段 3", "自研定量（PLS）",
          node("PLS1 回归（6 在轨标定靶训练）",
               "逐氧化物建模（auto 选成分数）；留一靶交叉验证防泄漏 → chemcam_pls_cv.json 冻结",
               "03_run_pls_cv.py · libs_quant.py", "green")),
    arrow(),
    phase("阶段 4", "对账（文献 / 官方基准）",
          row([node("Wiens 2012 Table 8c 区间判定",
                    "5/8 氧化物达区间或更优：TiO2 0.190/MgO/Na2O 达区间、K2O 优于区间、FeOT 经 Fe2O3 当量近似处区间内",
                    "gold_intervals.json · chemcam_pls_cv.json", "amber"),
               node("官方 MOC 对账并排",
                    "官方成分列并排（官方口径）；独立数据审计复算 CV",
                    "verify/data_audit.py", "amber")])
          + '<div class="merge">对账结论只进台账/报告注记，不回写官方产品</div>'),
    arrow(),
    phase("阶段 5", "输出",
          node("PLS 成分面板 / RMI 元素分布图 / 图像浏览",
               "观测级成分随 RMSE 标注；RMI 共帧元素分布投影（CAHV，不做跨图像投影）",
               "quant_ui.py · rmi_overlay.py · ChemCam 台", "slate")),
]
pages.append(page("f_chemcam.html", "ChemCam 数据处理流程", "Chemistry and Camera Complex (LIBS)",
                  "自研 PLS 成分层（谱级原始数据自算）与官方 MOC/文献对账口径分开", p,
                  "<b>口径边界：</b>官方 Clegg 2017 标定含 408 地面标样，官方 MOC 成分为定量引用口径；自研 6 靶模型是独立复现参考，不替代官方。"))

# ================= APXS =================
p = [
    phase("阶段 1", "数据获取（PDS 归档）",
          node("rsp 能谱 + rwp 官方氧化物表下载",
               "APXS 能谱原件与官方衍生表（rwp 只作对账基准）",
               "APXS/instruments/apxs.py · 数据准备脚本/01_multitarget_calib.py")),
    arrow(),
    phase("阶段 2", "自研峰拟合层",
          node("14 元素 Kα 高斯峰拟合",
               "线性背景 + Pu 子体参考线（子体线不计元素）；输出峰面积 + 泊松误差 + 检出判定；峰位自证 argmax vs 理论 Kα ≤0.06 keV",
               "APXS/instruments/apxs_fit.py", "green")),
    arrow(),
    phase("阶段 3", "浓度路线（历史负结果 + 比值层）",
          node("绝对面积→浓度：不可行（负结果入档）",
               "25 靶绝对面积标定斜率≈0——跨靶几何/接触距离因子主导；历史结论冻结，不再走绝对标定路线",
               "数据准备脚本/01_multitask_calib 结论", "green")
          + arrow()
          + node("Si 内标比值层（幂律双参标定）",
                 "元素/SiO2 比值 log-log 标定；Ca 相关 0.976、留出中位相对误差 2.5%；rocknest3_rp 时代外边界验证（负结果入档）",
                 "数据准备脚本/02_ratio_calibration.py", "green")),
    arrow(),
    phase("阶段 4", "对账（官方 / 文献基准）",
          row([node("官方 rwp 氧化物表对账",
                    "元素序与比值级对账（rocknest3_rp 首靶：Fe/Ca 面积比 6.69 vs 官方元素比 4.25，同量级区间）",
                    "verify/apxs_fit_check.py", "amber"),
               node("文献身份与独立审计",
                    "Schmidt 2013 rocknest3_rp 身份对账；峰面积独立复算逐元素对账",
                    "verify/paper_quant_check.py · data_audit.py", "amber")])
          + '<div class="merge">对账结论只进台账/报告注记，不回写官方产品</div>'),
    arrow(),
    phase("阶段 5", "输出",
          node("峰面积表 / 比值标定缓存 / 官方氧化物配对表",
               "apxs_peak_areas.json · apxs_ratio_calib.json · UI 氧化物配对表与峰拟合可视化",
               "APXS 台", "slate")),
]
pages.append(page("f_apxs.html", "APXS 数据处理流程", "Alpha-Particle X-ray Spectrometer",
                  "自研峰拟合 + Si 内标比值层与官方 rwp 对账口径分开", p,
                  "<b>口径边界：</b>比值预测是元素/SiO2 半定量信号，不是绝对浓度；全基本参数法定量不在此层（Campbell 2012 路线仅参数转录备查）。"))

# ================= PIXL =================
p = [
    phase("阶段 1", "数据获取（PDS 归档）",
          row([node("SFS 原始谱下载",
                    "data_raw_spectroscopy 集合：逐 PMC 4096 通道 XRF 计数，ENA/ENB 双探测器 CSV；清单+直链+SHA-256 校验",
                    "PIXL/instruments/pixl_download.py"),
               node("官方产品获取（只作对账基准）",
                    "RQC 逐点成分 / RXL 坐标 / MCC 图像；官方衍生表不作为输出",
                    "随包示例 · pixl_browser.py")])),
    phase("阶段 2", "能量轴自证（v4，自研）",
          node("多点线系能量自证",
               "10 线系支持评分（Si/S/Cl/K/Ca/Ti/Mn/Fe）→ g≈0.00786 keV/ch；发现并修正旧轴 13% 增益偏差；双扫描独立自证一致 ±0.3%",
               "PIXL/instruments/sfs_spectra.py calibrate_energy v4", "green")),
    arrow(),
    phase("阶段 3", "自研谱层（峰面积信号）",
          node("Kα 窗净峰面积 + 检出判定",
               "逐元素窗净面积 + 泊松误差；A/B 探测器择低（Schmidt 2025 Table S4 口径）",
               "sfs_spectra.peak_areas", "green")
          + arrow()
          + row([node("66 对比值矩阵",
                      "全元素两两 log-log 比值相关（Top：Si/Fe 0.994）",
                      "01_sfs_pipeline.py", "green"),
                 node("扫描内半定量标定",
                      "wt = k×area 过原点最小二乘（Ca R²=0.981）；k 仅本扫描有效",
                      "01_sfs_pipeline.py", "green")])),
    phase("阶段 4", "PIQUANT 官方开源定量层",
          node("SFS→MSA 导出 → Piquant quant",
               "MSA 单列计数格式 + 本扫描自证轴；FM 配置 + FM 标定（官方 COMMENT 指定 ENB）→ 逐 PMC wt%；sha256 指纹防静默换版",
               "instruments/piquant_bridge.py · 数据准备脚本/02_piquant_quantify.py", "green")),
    arrow(),
    phase("阶段 5", "对账（三线，官方 / 文献基准）",
          row([node("自研层 × RQC",
                    "逐 PMC 相关 + 盐校正对账（Schmidt 2025：8/10 氧化物落 ±σ）",
                    "sfs_vs_rqc.json", "amber"),
               node("PIQUANT × RQC（双扫描）",
                    "Bellegarde 281 点：FeO-T/SiO2/CaO/Al2O3 中位偏差 -1.6/-2.3/-4.0/-4.3%、相关 ≥0.98；sol0300 334 点复验泛化",
                    "piquant_vs_rqc.json · piquant_vs_rqc_sol0300.json", "amber"),
               node("独立数据审计（三路复算）",
                    "rows↔stats 一致性 / RQC 原表抽 20 点交叉核验 / spot 重 quant 逐位复现",
                    "verify/data_audit.py · piquant_check.py", "amber")])
          + '<div class="merge">三线口径分开：自研峰面积 / PIQUANT 官方算法复现 / 官方 RQC（A/B 择低）</div>'),
    arrow(),
    phase("阶段 6", "输出",
          node("UI 面板 / 对账缓存 / 批量入口",
               "MCC 图像+逐点成分浏览、SFS 面板（能量轴+比值+PIQUANT 状态+批量对账按钮）、独立导出",
               "pixl_view.py · pixl_workflow.py PiquantJob · PIXL 台", "slate")),
]
pages.append(page("f_pixl.html", "PIXL 数据处理流程", "Planetary Instrument for X-ray Lithochemistry",
                  "六层定量链：能量自证 → 自研峰面积 → 比值/标定 → PIQUANT 官方层 → 三线对账 → 输出", p,
                  "<b>口径边界：</b>自研峰面积是原始信号量不是浓度；PIQUANT 输出代表官方算法复现；官方 RQC 为定量引用口径；FM 标定不覆盖 P2O5/MnO/Cr2O3/Sr（无 ECF 条目）如实排除。"))

for f in sorted(OUT.glob("*.html")):
    print("生成", f.name, f.stat().st_size, "bytes")
