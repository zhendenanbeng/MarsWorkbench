# -*- coding: utf-8 -*-
"""矿物类别判定，规则来自 Corpolongo et al. (2023, JGR Planets)。

论文口径：

* 噪声 RMSD 取自 2500-3000 cm-1 的安静区，取逐点值的中位数
* 潜在拉曼信号：FWHM 至少 30 cm-1（约 3 像素），强度至少 3 倍噪声 RMSD
* 判定用平均谱，不是单点谱
* 只到矿物类别，不到种属，除非有明确的次级峰证据

论文里的几条关键经验：

* 橄榄石和辉石的识别峰在 800 cm-1 以上，SHERLOC 能认；长石、层状硅酸盐、
  角闪石的主峰在 800 以下，被激光注入滤光片挡掉，只能在 1050 附近看到一个
  宽带，记为未定硅酸盐。
* 硫酸盐和辉石的主峰都在 970-1050 cm-1，靠次级峰和 O-H 带区分。
* 高氯酸盐和磷酸盐的主峰都在 950-970 cm-1，靠次级峰区分；没有次级峰时
  只能记为高氯酸盐或磷酸盐。
* 暗背景扣除会在 970 和 981 cm-1 留下约两成频率出现的假峰。
"""

from __future__ import annotations

import numpy as np

from . import library
from .peaks import half_maximum
from .fitting import background_context, solve_nonnegative_profiles

MIN_FWHM_CM = 30.0
MAX_FWHM_CM = 200.0
SNR_MIN = 3.0

# 拉曼指纹区。**上界从 1500 改到 1800（2026-09-15）**：两份论文用的都是 800–1800
# ——Corpolongo et al. (2023) 的峰位总图写明是 "the 800–1800 cm-1 region"；
# Hurowitz et al. (2025) 判 Bright Angel 的有机质靠的是 ~1600 cm-1 的 G 带。
# 原来 1500 的上界把 G 带整条通道切掉了（sol 370 的 1576/1604/1632 在类别那一栏
# 从来不出现）。这里是**唯一**的定义处，qc.py 直接引用，免得两处走散。
FINGERPRINT = (800.0, 1800.0)
# 有机质 G 带（碳的 E2g 伸缩）单独一条通道：判矿物类别用不上它，
# 但它是 SHERLOC 判有机质的唯一依据，得单独报。
G_BAND = (1560.0, 1680.0)

# 熔融石英（SHERLOC 光学窗口）的背景特征，两篇论文都点过名：
#   Corpolongo et al. (2023) §2.4：~480 与 ~800 cm-1 的弱峰，来自仪器光学；
#   Hurowitz et al. (2025)：~1600 cm-1 的 Si-O 伸缩泛音，与有机质 G 带重叠，
#   他们是用灰色仪器背景谱去比、才把 G 带认出来的。
# 我们的窗口从 800 起，所以只有 ~800 与 ~1600 这两处需要提醒。
SILICA_BACKGROUND_CM = (800.0, 1600.0)

# 峰的**突出度**判据（prominence）。默认 0 = 关闭，行为和以前完全一致。
#
# 为什么需要它：以前的判据只看"峰高 ≥ 3σ"，但那是在跟**噪声**比。
# SHERLOC 的谱上普遍坐着一个宽荧光包，扣完基线后还剩一个几十到上百
# 波数宽的鼓包；鼓包上的噪声会让 argmax 在相邻通道之间乱跳，于是同一个
# 鼓包在几十个点上被"测"成 976 / 1006 / 1035 / 1103 … 等不同峰位，
# 聚类就把一个宽包劈成一串互不相干的簇（sol 370 Alfalfa 实测：
# 41 个有点的点被拆成 4 簇，其中 42 个点次掉进"未归类"）。
# 突出度把峰高改成跟**自己两侧最近的山谷**比，鼓包上的假极大立刻掉到 1σ 以下。
#
# 取值口径：和论文的显著性判据同源（3×噪声 RMSD），所以默认给 3.0 是
# "有理有据"的起点；但先留一个开关，扫参数看它到底改善了什么、破坏了什么。
PROM_SIGMA = 0.0
# 谱库比对用的峰门槛。**2026-09-15 从 5σ 降到 3σ**：
# 论文判"潜在拉曼信号"用的就是 ≥3σ，而 5σ 会把关键的次级峰砍掉——
# 实测 sol 162（Guillaumes）的高氯酸盐组，论文点名的两个次级峰
# （1091.7、1150.4）在组平均谱上是 3.9σ 和 4.0σ，全被 5σ 挡在门外，
# 于是谱库比对只看得到主峰、任何高氯酸盐标准都对不上（得分只有 0.1）。
# 早期改成 5σ 是因为没有高氯酸盐标准，噪声峰会把得分拉垮（误判成蛋白石）；
# 现在库里有 10 个氧氯酸盐样品，降回 3σ 才是和论文一致的口径。
LIB_SNR_MIN = 3.0
NOISE_WINDOW = (2500.0, 3000.0)
# 谱线伪影。**2026-09-16 补全**：原来只有 970/981 两个，现在换成论文支持材料
# Text S3 给的完整清单（Corpolongo et al. 2023 的 Supporting Information，
# 从 JPL 仓库的记录 10.48577/jpl.IMFOK5 里取到的那份完整稿）。
#
# 论文原话：这些是"单像素在暗/亮谱里背景都偏高（>1000 计数）"造成的，
# 伪影只影响**单个像素**，而 Raman 带 FWHM ≥ 3 像素，所以大多数情况下能区分。
# 唯独 970 和 981 最要命——它们正好落在高氯酸盐/磷酸盐的主峰区，
# 让那一带的低强度峰无法可靠定类（正文 2.4 节专门讲了这件事）。
ARTIFACT_CM = (970.0, 981.0)
ARTIFACT_OTHER_CM = (
    331.0, 694.0, 784.0, 1164.0, 1307.0, 1431.0, 1477.0, 1534.0, 1580.0,
    1636.0, 2563.0, 3002.0, 3360.0, 3504.0, 3528.0, 3763.0, 3802.0, 3847.0,
    3863.0, 3947.0,
)

# 水带（O-H 伸缩）判据，照操作手册 A-4 / A-5：
#   拟合区间 3000-4000，峰位限制 3200-3700，FWHM 40-300，SNR > 3
OH_FIT_WINDOW = (3000.0, 4000.0)
OH_PEAK_WINDOW = (3200.0, 3700.0)
OH_MIN_FWHM = 40.0
# 手册 A-5 给的单峰上限是 300 cm-1。但 SHERLOC 的水带在这个靶点上是
# 3229 和 3415 两个峰叠在一起（Corpolongo et al., 2023 对 Quartier 的报告），
# 两峰包络实测 331 cm-1，按单峰上限会被整条否掉。所以放宽到 400，
# 并同时用"双峰拟合优先"来保证这不是在给荧光鼓包开后门。
OH_MAX_FWHM = 400.0

# 稳健 σ：手册 A-5 说某些光谱受荧光影响、σ 被抬高，
# 这时改用该次扫描所有光谱 σ 的中位数。倍数手册没给，这里取 3。
SIGMA_OUTLIER_FACTOR = 3.0

# 拟合用的物理边界。上限给 250 是留出硅酸盐那种 180 cm-1 的宽带；
# 下限 20 是"比仪器分辨本领（41）还窄一半"的宽容值——真窄到这个程度
# 基本就是噪声尖刺，不该被当峰。
FIT_MIN_FWHM = 20.0
FIT_MAX_FWHM = 250.0
# 拟合出的峰心最多允许偏离"本簇波段位置的中位数"这么远。宽峰上不加这个
# 约束，最小二乘会顺着尾巴把峰心拽出类别区间（Dourbes 的橄榄石就是这么
# 从 839 跑到 879 的）。
FIT_CENTER_TOL = 30.0

# 「一个峰还是两个峰」检验的门槛
SPLIT_MIN_POINTS = 4          # 点太少的组不拆
SPLIT_MIN_SEP_CM = 20.0       # 两个组分至少隔这么远才有意义
SPLIT_RMS_GAIN = 0.80         # 加一个组分至少要把残差降 20% 才算数
SPLIT_RESOLVE = 0.8           # 两个峰心的距离至少要有较窄那个峰宽的 0.8 倍

# ---- 扫描级兜底（2026-09-15 加，被 sol 1201 逼出来的）-----------------------
#
# 主流程要求**单个点**就达到 ≥3σ 且 FWHM ≥ 30 才有资格进聚类。当单点噪声
# 很高时（sol 1201 是 150–220 计数，火山口底那批只有 62–142），100 个点里
# 可能只有 2–7 个过门，聚类就散成 1–3 个点的小簇、落在类别区间的缝里，
# 报告于是给出"没有可用的拉曼信号"——而**全扫描平均谱上其实有 3–4σ 的带**。
#
# 所以补一条：逐点路径数据不足时，改看全扫描平均谱。
SCAN_LEVEL_MAX_POINTS = 10     # 过门的点少于这个数，才启用扫描级检查
SCAN_LEVEL_MIN_SNR = 3.0       # 带峰相对"平均谱自身噪声"的下限
# 0.30 是用全库标定出来的（verify/scan_level_sweep.py）：
# 已知真检出的这个比值在 0.30–4.29；已知"论文说零"的扫描里最高的是
# sol 141 第一次扫描的 0.28。sol 1201 四次扫描是 0.33–0.51。
# 也就是说 0.30 这条线把"真检出 + sol 1201"和"已知的零"分开，余量不大，
# 所以报告里必须把实测数值一起印出来，让人自己判断。
SCAN_LEVEL_MIN_RATIO = 0.30    # 带峰相对"单点噪声中位数"的下限
SCAN_LEVEL_MIN_FWHM = 30.0
SCAN_LEVEL_MAX_FWHM = 200.0

# 扫描级兜底**单独用 800–1500 的窗口 + 普通多项式基线**，不跟着 FINGERPRINT 走。
# 为什么（2026-09-15 实测）：把它扩到 800–1800 之后，平均谱上"最强的带"变成了
# 熔融石英光学窗口的 ~800 背景——它在 sol 141、sol 181 这两台"论文说没有矿物"
# 的对照里同样出现，于是兜底从"找矿物"变成了"找仪器背景"，已知的零全部冒候选。
# 换成掩峰迭代基线更糟：那条基线专门保峰，把背景一起保下来了。
# 而 800–1500 + 普通多项式这套组合是标定过的（真检出 0.30–4.29，已知的零 ≤0.28）。
# 1600 附近另走 g_band 那条独立通道，两者互不干扰。
SCAN_LEVEL_WINDOW = (800.0, 1500.0)

# 有机质 G 带的幅度下界（相对单点噪声中位数）。
# 和扫描级兜底同一个思路：只看显著性的话，上千点的平均谱里噪声也够 3σ。
G_BAND_MIN_RATIO = 0.10

RANGES = {
    "olivine": (815.0, 865.0),
    "perchlorate_chlorate": (920.0, 945.0),
    # 论文正文原话：高氯酸盐和磷酸盐的主峰落在"约 950 到 970 cm-1"。
    # 之前这里写到 995，会把 987 那种辉石峰也吞进来。上限收到 975。
    "perchlorate_or_phosphate": (945.0, 975.0),
    "pyroxene": (985.0, 1012.0),
    "sulfate": (1000.0, 1035.0),
    "silicate": (1030.0, 1075.0),
    "carbonate": (1070.0, 1100.0),
}

LABELS = {
    "olivine": "橄榄石",
    "carbonate": "碳酸盐",
    "sulfate": "硫酸盐",
    "pyroxene": "辉石",
    "silicate": "未定硅酸盐",
    "perchlorate_chlorate": "高氯酸盐/氯酸盐",
    "perchlorate_or_phosphate": "高氯酸盐或磷酸盐",
    "phosphate": "磷酸盐",
    "sodium_perchlorate": "高氯酸钠",
    "sulfate_or_pyroxene": "硫酸盐或辉石（缺次级峰）",
    "unassigned": "未归类",
}

# 逐点图上用的缩写，对齐 Corpolongo et al. (2023) 图里的字母习惯：
# C 碳酸盐、S 硫酸盐、O 橄榄石、Si 硅酸盐、P 高氯酸盐或磷酸盐。
SHORT = {
    "olivine": "O",
    "carbonate": "C",
    "sulfate": "S",
    "pyroxene": "Px",
    "silicate": "Si",
    "perchlorate_chlorate": "Cl",
    "perchlorate_or_phosphate": "P",
    "phosphate": "Ph",
    "sodium_perchlorate": "NaP",
    "sulfate_or_pyroxene": "S/Px",
    "unassigned": "?",
}

# 各矿物类别的配色。质控台的网页和出图脚本都用这一份，
# 免得两边各写一套颜色、过两天就对不上了。
CLASS_COLORS = {
    "olivine": "#1a9e63",
    "carbonate": "#2f7de1",
    "sulfate": "#d4751f",
    "silicate": "#8b5cf6",
    "pyroxene": "#0d8b99",
    "perchlorate_or_phosphate": "#c2185b",
    "perchlorate_chlorate": "#a33ea1",
    "phosphate": "#d32f2f",
    "sodium_perchlorate": "#6a3bbd",
    "sulfate_or_pyroxene": "#b0872f",
    "unassigned": "#6b7683",
}


def baseline_masked(x, y, degree=8, rounds=6, mask_cm=25.0, snr=3.0):
    """掩峰迭代基线 —— 论文 2.3 节第 5 步的做法。

    论文原话：基线拟合在**一组不含拉曼峰的点**上，用最小二乘；有些样品还要做
    第二次。我们原来的做法是把多项式拟合在**所有点**上，于是被 SHERLOC 那种
    宽包背景拽跑：实测组平均谱扣完基线后残差仍是正的、还带一个从 900 衰减到
    1800 的宽包，真正的次级峰全埋在肩上。

    这里按论文的思路迭代：拟合 → 在残差上找显著局部极大 → 把它们周围
    ±mask_cm 的点排除 → 重新拟合，重复若干轮。实测噪声 sigma 从 35 降到 13.5。
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size < degree + 3:
        return np.zeros_like(y)
    keep = np.ones_like(y, dtype=bool)
    coef = np.polyfit(x, y, degree)
    for _ in range(rounds):
        resid = y - np.polyval(coef, x)
        sig = resid[keep].std() if keep.sum() > degree + 2 else resid.std()
        peaks = [x[j] for j in range(1, x.size - 1)
                 if resid[j] > resid[j - 1] and resid[j] >= resid[j + 1]
                 and sig > 0 and resid[j] > snr * sig]
        new = np.ones_like(y, dtype=bool)
        for p in peaks:
            new &= ~(np.abs(x - p) <= mask_cm)
        if new.sum() < degree + 3 or new.sum() == keep.sum():
            break
        keep = new
        coef = np.polyfit(x[keep], y[keep], degree)
    return np.polyval(coef, x)


def baseline(x, y, degree=5, n_iter=10, clip=1.5):
    keep = np.ones_like(y, dtype=bool)
    coef = np.polyfit(x[keep], y[keep], degree)
    for _ in range(n_iter):
        resid = y - np.polyval(coef, x)
        sigma = resid[keep].std()
        keep = resid < clip * sigma
        if keep.sum() < degree + 2:
            break
        coef = np.polyfit(x[keep], y[keep], degree)
    return np.polyval(coef, x)


def prominence(y, j, span):
    """峰 j 的突出度：峰高减去「到左右第一个更高点之间最低谷」中较高的那个。

    span 是左右各看多远（通道数）。用它当半径是必须的——不设半径时，
    一条缓慢下坡的荧光背景会让每个小鼓包都"永不下降"，突出度算成无穷大。
    """
    lo_end = max(0, j - span)
    hi_end = min(len(y) - 1, j + span)
    left, lmin = j, y[j]
    while left > lo_end and y[left] <= y[j]:
        lmin = min(lmin, y[left])
        left -= 1
    right, rmin = j, y[j]
    while right < hi_end and y[right] <= y[j]:
        rmin = min(rmin, y[right])
        right += 1
    return float(y[j] - max(lmin, rmin))


def find_bands(x, y, sigma, min_fwhm=MIN_FWHM_CM, max_fwhm=MAX_FWHM_CM,
               snr_min=SNR_MIN, prom_sigma=PROM_SIGMA):
    """按论文判据找波段。prom_sigma > 0 时额外要求突出度达标（见 PROM_SIGMA 注释）。"""
    out = []
    n = len(y)
    # 突出度的搜索半径按落在这个窗口里的通道数折中：太短会把真峰的山谷
    # 找到自己肩上，太长会把隔壁矿物的峰当成自己的山谷。
    span_ch = max(3, int(round(120.0 / max(1e-9, float(np.median(np.diff(x)))))))
    for j in range(1, n - 1):
        if not (y[j] > y[j - 1] and y[j] >= y[j + 1]):
            continue
        if not np.isfinite(sigma) or sigma <= 0 or y[j] < snr_min * sigma:
            continue
        if prom_sigma > 0 and prominence(y, j, span_ch) < prom_sigma * sigma:
            continue
        # 保持开发集已验证的候选筛选口径。跨格宽度是交点宽度的上界，
        # 不是精确FWHM；插值测量/截断状态单独输出，不能偷偷替换阈值后丢失弱峰。
        half = y[j] / 2.0
        k = m = j
        while k > 0 and y[k] > half:k -= 1
        while m < n - 1 and y[m] > half:m += 1
        fwhm = x[m] - x[k]
        width = half_maximum(x, y, j, allow_censored=True)
        if fwhm < min_fwhm or fwhm > max_fwhm:
            continue
        out.append({"pos": float(x[j]), "height": float(y[j]),
                    "fwhm": float(fwhm), "snr": float(y[j] / sigma),
                    "fwhm_interpolated":width['fwhm'] if width else None,
                    "width_censored":width['censored'] if width else True,
                    "width_method":"sample_bracket_upper_bound"})
    return out


def robust_sigma(sig, factor=SIGMA_OUTLIER_FACTOR):
    """手册 A-5 的 σ 补丁：σ 被荧光抬高的光谱，改用扫描中位 σ。

    返回 (实际使用的 σ 数组, 被替换的点号列表)。
    """
    sig = np.asarray(sig, dtype=float)
    good = sig[np.isfinite(sig) & (sig > 0)]
    if good.size < 3:
        return sig, []
    med = float(np.median(good))
    bad = np.where(~np.isfinite(sig) | (sig > factor * med))[0]
    out = sig.copy()
    for i in bad:
        out[i] = med
    return out, [int(i) for i in bad]


def _pseudo_voigt(x, a, x0, fwhm, eta):
    """伪 Voigt：高斯与洛伦兹按 eta 混合。eta=0 纯高斯，eta=1 纯洛伦兹。"""
    s = max(fwhm, 1e-6) / 2.0
    d = (x - x0) / s
    g = np.exp(-np.log(2.0) * d * d)
    l = 1.0 / (1.0 + d * d)
    return a * ((1.0 - eta) * g + eta * l)


def _voigt_basis(x, cx, fw):
    """高斯分量与洛伦兹分量，分开给出来。"""
    d = (x - np.asarray(cx, dtype=float)) / max(float(fw), 1e-6) * 2.0
    return np.exp(-np.log(2.0) * d * d), 1.0 / (1.0 + d * d)


def _solve_linear(x, y, comps, nonnegative=False, context=None):
    """给定各分量的 (峰心, 峰宽)，线性解出 幅度/基线 并返回残差。

    关键点：伪 Voigt 可以写成 c1*高斯 + c2*洛伦兹，
    所以只要峰心和峰宽定了，剩下的参数（两个分量系数 + 线性底的两个系数）
    都是**线性**的，一次最小二乘就能求出，不需要迭代、也不会发散。
    峰心/峰宽用一个有限网格搜，于是整个拟合是确定性的：
    同样的输入在任何解释器上都给同样的输出。
    """
    cols = []
    for cx, fw in comps:
        g, l = _voigt_basis(x, cx, fw)
        cols.append(g)
        cols.append(l)
    if nonnegative:
        return solve_nonnegative_profiles(x,y,np.column_stack(cols),context=context)
    cols.append(np.asarray(x, dtype=float))
    cols.append(np.ones_like(np.asarray(y, dtype=float)))
    a_mat = np.column_stack(cols)
    coef, *_ = np.linalg.lstsq(a_mat, y, rcond=None)
    model = a_mat @ coef
    resid = np.asarray(y, dtype=float) - model
    return coef, float(np.sum(resid * resid)), model


def _grid_search(x, y, cx0, span, cx_step, widths, exclude=None, fixed=(),
                 within=None, nonnegative=False):
    """在 (峰心 × 峰宽) 网格上搜最优格点。fixed 是已经定下来的其他分量。

    每个分量都必须有正幅度——否则最小二乘会用"减去一个 Voigt 形状"
    去贴残差，解出来的峰是负的，物理上毫无意义（踩过这个坑）。
    """
    best = None
    context=background_context(x,y) if nonnegative else None
    for cx in np.arange(cx0 - span, cx0 + span + 1e-9, cx_step):
        if not (x[0] - 1e-6 <= cx <= x[-1] + 1e-6):
            continue
        if within is not None and not within[0]-1e-9<=cx<=within[1]+1e-9:
            continue
        if exclude and abs(cx - exclude[0]) < exclude[1]:
            continue
        for fw in widths:
            comps = list(fixed) + [(float(cx), float(fw))]
            coef, err, model = _solve_linear(x, y, comps, nonnegative=nonnegative,context=context)
            # 只约束"正在搜的这个分量"必须为正。已固定的分量如果也一并约束，
            # 会在水带那种双峰重叠的场合把所有候选都杀光（踩过）。
            # 固定分量的正负由调用方在最后统一检查。
            if coef[2 * (len(comps) - 1)] + coef[2 * (len(comps) - 1) + 1] <= 0.0:
                continue
            if best is None or err < best[0]:
                best = (err, float(cx), float(fw), coef, model)
    return best


def fit_voigt(x, y, x0, fwhm0, min_fwhm=FIT_MIN_FWHM, max_fwhm=FIT_MAX_FWHM,
              center_span=90.0, within=None, nonnegative=False):
    """用一个伪 Voigt 峰 + 线性底拟合一段光谱，返回峰心/峰宽/幅度。

    先粗网格再细网格，最后检查结果有没有撞到物理边界。撞边界或峰宽跑到
    判据下限以下，都说明这次拟合不可信，调用方应该退回半高穿零法。

    nonnegative=True时，高斯/洛伦兹系数均非负，用于物理约束复核。
    默认保留历史候选定位的有符号拟合，但返回真实eta和profile_valid，
    不能通过截断eta把非物理解伪装成合法伪Voigt；RMS须可由返回参数重建。
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size < 10:
        return None
    if (x.ndim!=1 or x.shape!=y.shape or not np.isfinite(x).all() or not np.isfinite(y).all()
            or not np.all(np.diff(x)>0) or not np.isfinite([x0,fwhm0,min_fwhm,max_fwhm,center_span]).all()
            or min_fwhm<=0 or max_fwhm<=min_fwhm or center_span<0):
        raise ValueError('无效的峰形拟合数据或边界')

    if within is not None:
        if len(within)!=2 or not np.isfinite(within).all() or within[0]>=within[1]:
            raise ValueError('峰心约束必须是有限且递增的区间')
        # 峰心只在这个范围内找——这是"本簇的波段到底在哪"的直接约束
        center_span = max(center_span, 0.0)
        lo_c = max(x[0], within[0])
        hi_c = min(x[-1], within[1])
        if hi_c <= lo_c:
            return None
        x0 = 0.5 * (lo_c + hi_c)
        center_span = 0.5 * (hi_c - lo_c)

    widths = np.array([min_fwhm, 25.0, 30.0, 40.0, 50.0, 65.0, 85.0,
                       110.0, 140.0, 180.0, 220.0, max_fwhm])
    widths = np.unique(widths[(widths >= min_fwhm) & (widths <= max_fwhm)])
    if widths.size == 0:
        return None

    allowed=(max(float(x[0]),x0-center_span),min(float(x[-1]),x0+center_span))
    coarse = _grid_search(x, y, x0, center_span, 5.0, widths,within=allowed,nonnegative=nonnegative)
    if coarse is None:
        return None
    _, cx0, fw0, _, _ = coarse

    fine_w = np.arange(max(min_fwhm, fw0 - 12.0),
                       min(max_fwhm, fw0 + 12.0) + 1e-9, 2.0)
    if fine_w.size == 0:
        fine_w = np.array([fw0])
    fine = _grid_search(x, y, cx0, 6.0, 1.0, fine_w,within=allowed,nonnegative=nonnegative)
    if fine is None:
        fine = coarse
    err, cx, fw, coef, _ = fine

    c1, c2 = float(coef[0]), float(coef[1])
    amp = c1 + c2
    if not np.isfinite(amp) or amp <= 0.0:
        return None
    eta = float(c2 / amp)
    profile_valid=bool(c1>=0 and c2>=0)
    at_bound = (fw <= min_fwhm + 1e-9 or fw >= max_fwhm - 1e-9
                or cx<=allowed[0]+1e-9 or cx>=allowed[1]-1e-9)
    return {"amp": float(amp), "center": float(cx), "fwhm": float(fw),
            "eta": eta, "rms": float(np.sqrt(err / x.size)),
            "slope": float(coef[2]), "intercept": float(coef[3]),
            "at_bound": bool(at_bound), "c1": c1, "c2": c2,
            "profile_valid":profile_valid,
            "model":"pseudo_voigt" if profile_valid else "signed_gaussian_lorentz",
            "coefficient_constraint":"nonnegative" if nonnegative else "signed_candidate"}


def fit_two_voigt(x, y, x0, fwhm0, sigma, min_fwhm=MIN_FWHM_CM,
                  max_fwhm=FIT_MAX_FWHM, snr_min=SNR_MIN, sep=SPLIT_MIN_SEP_CM,
                  within=None, nonnegative=False):
    """检验这一段是不是需要两个峰。结果里 ok=True 才表示第二个峰成立。

    做法：先拟合一个组分，在残差上找第二个组分，再交替精修。
    比四维网格便宜得多，而且全程没有迭代发散的可能。

    和 hydration 一样：总是返回字典，不成立时给出 reject 原因，
    这样"为什么没拆"永远查得到。

    within=(lo, hi) 限制两个分量的峰心必须落在这个范围内——拆出来的峰
    必须是本簇自己的波段支撑得起的，不能借隔壁簇的峰。
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    # 6 个待解参数（两分量各 2 个系数 + 线性底 2 个），给 12 个通道就够；
    # 卡在 20 会把 ±90 cm-1 这种只有 19 个通道的窗口全挡掉（踩过）。
    if x.size < 12 or not (sigma > 0):
        return {"ok": False, "reject": "窗口只有 %d 个通道或 σ 无效" % x.size}

    one = fit_voigt(x, y, x0, fwhm0, min_fwhm=min_fwhm, max_fwhm=max_fwhm,
                    within=within if nonnegative else None,nonnegative=nonnegative)
    if one is None:
        return {"ok": False, "reject": "单组分拟合就失败了"}
    c1, f1 = one["center"], one["fwhm"]

    _, _, model1 = _solve_linear(x, y, [(c1, f1)],nonnegative=nonnegative)
    resid = y - model1
    widths = np.arange(min_fwhm, max_fwhm + 1e-9, 5.0)
    sec = _grid_search(x, resid, x0, 90.0, 5.0, widths,
                       exclude=(c1, sep), fixed=(),nonnegative=nonnegative,
                       within=within if nonnegative else None)
    if sec is None:
        return {"ok": False, "reject": "残差上找不到第二个峰"}
    _, c2, f2, _, _ = sec

    # 交替精修：各自在对方的背景下再搜一遍
    for _ in range(2):
        r = _grid_search(x, y, c1, 15.0, 3.0, widths, fixed=[(c2, f2)],nonnegative=nonnegative,
                         within=within if nonnegative else None)
        if r:
            _, c1, f1 = r[0], r[1], r[2]
        r = _grid_search(x, y, c2, 15.0, 3.0, widths,
                         exclude=(c1, sep), fixed=[(c1, f1)],nonnegative=nonnegative,
                         within=within if nonnegative else None)
        if r:
            _, c2, f2 = r[0], r[1], r[2]

    if abs(c1 - c2) < sep:
        return {"ok": False, "reject": "两峰只差 %.0f cm-1（需要 ≥%.0f）"
                                     % (abs(c1 - c2), sep),
                "one": {"center": one["center"], "fwhm": one["fwhm"]}}
    # 可分辨性：两峰距离要够开。1050 与 1083 只差 33 cm-1、各自峰宽 50 上下，
    # 这种"肩"在自动拟合里分不开，硬拆只会把本来完整的一簇切碎。
    need = SPLIT_RESOLVE * min(f1, f2)
    if abs(c1 - c2) < need:
        return {"ok": False,
                "reject": "两峰差 %.0f cm-1，小于可分辨下限 %.0f（较窄峰的 %.0f 倍）"
                          % (abs(c1 - c2), need, SPLIT_RESOLVE)}
    comps = sorted([(c1, f1), (c2, f2)])
    if within is not None:
        outside = [c for c, _ in comps if not (within[0] <= c <= within[1])]
        if outside:
            return {"ok": False,
                    "reject": "有分量落在本簇波段范围 %s 之外（%.0f）"
                              % ("%.0f-%.0f" % within, outside[0])}
    coef, err, _ = _solve_linear(x, y, comps,nonnegative=nonnegative)
    amp1 = float(coef[0] + coef[1])
    amp2 = float(coef[2] + coef[3])
    fw1, fw2 = comps[0][1], comps[1][1]

    if amp2 <= 0 or amp1 <= 0:
        return {"ok": False, "reject": "有两个组分的幅度是负的"}
    if amp2 < snr_min * sigma:                     # 第二峰太弱，别硬加
        return {"ok": False, "reject": "第二峰只有 %.1fσ（需要 ≥%.0fσ）"
                                     % (amp2 / sigma, snr_min)}
    if not (min_fwhm <= fw1 <= max_fwhm and min_fwhm <= fw2 <= max_fwhm):
        return {"ok": False, "reject": "峰宽 %.0f/%.0f 超出 %.0f-%.0f"
                                     % (fw1, fw2, min_fwhm, max_fwhm)}
    rms2 = float(np.sqrt(err / x.size))
    if rms2 > SPLIT_RMS_GAIN * one["rms"]:         # 加一个峰没能明显改善
        return {"ok": False,
                "reject": "加第二峰只把残差降了 %.0f%%（需要 ≥%.0f%%）"
                          % ((1 - rms2 / one["rms"]) * 100,
                             (1 - SPLIT_RMS_GAIN) * 100)}
    return {"ok": True,
            "comps": [{"center": comps[0][0], "fwhm": fw1, "amp": amp1,
                       'eta':float(coef[1]/amp1),'profile_valid':bool(np.all(coef[:2]>=0))},
                      {"center": comps[1][0], "fwhm": fw2, "amp": amp2,
                       'eta':float(coef[3]/amp2),'profile_valid':bool(np.all(coef[2:4]>=0))}],
            'slope':float(coef[4]),'intercept':float(coef[5]),
            'profile_valid':bool(np.all(coef[:4]>=0)),
            'components_significant':bool(min(amp1,amp2)>=snr_min*sigma),
            "rms": rms2, "rms_one": one["rms"], "reject": None}


def hydration(x, y, sigma, min_fwhm=OH_MIN_FWHM, max_fwhm=OH_MAX_FWHM):
    """O-H 伸缩带。手册口径：拟合 3000-4000，峰位限 3200-3700，FWHM 40-300。

    硫酸盐的关键判据之一：有 O-H 带就支持含水硫酸盐。

    这个函数**总是返回一个字典**，检到时 ok=True，没检到时 ok=False 并给出
    reject 原因。因为"没检到"和"检到但被规则否掉"是完全不同的两件事，
    报告里要能区分。
    """
    out = {"ok": False, "reject": None}
    sel = (x >= OH_FIT_WINDOW[0]) & (x <= OH_FIT_WINDOW[1])
    pk = (x >= OH_PEAK_WINDOW[0]) & (x <= OH_PEAK_WINDOW[1])
    if not sel.any() or not pk.any() or not (sigma > 0):
        out["reject"] = "没有 3000-4000 的数据，或 σ 无效"
        return out
    xs, seg = x[sel], y[sel]
    k = int(np.argmax(np.where(pk[sel], seg, -np.inf)))
    if seg[k] < SNR_MIN * sigma:
        out["reject"] = ("3200-3700 内最强点只有 %.1fσ（需要 ≥3σ）"
                         % (seg[k] / sigma))
        return out

    half = seg[k] / 2.0
    L = R = k
    while L > 0 and seg[L] > half:L -= 1
    while R < len(seg)-1 and seg[R] > half:R += 1
    fw0 = float(xs[R]-xs[L])
    measured = half_maximum(xs, seg, k)

    def usable(comp):
        return (OH_PEAK_WINDOW[0] <= comp["center"] <= OH_PEAK_WINDOW[1]
                and min_fwhm <= comp["fwhm"] <= max_fwhm
                and comp["amp"] >= SNR_MIN * sigma)

    # 先试两个组分：论文里 Quartier 的水带就是 3229 与 3415 两个，单峰会宽到 300+。
    two = fit_two_voigt(xs, seg, float(xs[k]), max(fw0, min_fwhm), sigma,
                        min_fwhm=min_fwhm, max_fwhm=max_fwhm)
    if two and two.get("ok"):
        good = [c for c in two["comps"] if usable(c)]
        if good:
            best = max(good, key=lambda c: c["amp"])
            out.update(pos=best["center"], fwhm=best["fwhm"], height=best["amp"],
                       method="voigt2", comps=[dict(c) for c in good],
                       snr=float(best["amp"] / sigma), rms=two["rms"])
            out['profile_valid']=two.get('profile_valid',True)
            out['requires_review']=not out['profile_valid']
            out["ok"] = True
            return out
        out["reject"] = ("双组分拟合出的两个峰都不合格（"
                         + "、".join("%.0f cm-1/FWHM %.0f" % (c["center"], c["fwhm"])
                                     for c in two["comps"]) + "）")

    fit = fit_voigt(xs, seg, float(xs[k]), max(fw0, min_fwhm),
                    min_fwhm=min_fwhm, max_fwhm=max_fwhm)
    if (fit is not None and not fit["at_bound"] and fit["fwhm"] >= min_fwhm
            and fit["amp"] >= SNR_MIN * sigma):
        center, fw, amp = fit["center"], fit["fwhm"], fit["amp"]
        method = "voigt"
    else:
        center, fw, amp = float(xs[k]), fw0, float(seg[k])
        method = "halfmax"

    out.update(pos=center, fwhm=fw, height=amp, method=method,
               snr=float(amp / sigma), rms=fit["rms"] if fit else float("nan"))
    out['profile_valid']=fit.get('profile_valid',True) if method=='voigt' else None
    out['requires_review']=method=='halfmax' or out['profile_valid'] is False
    out['width_censored'] = method == 'halfmax' and measured is None
    out['fwhm_interpolated'] = measured['fwhm'] if measured else None
    if not (OH_PEAK_WINDOW[0] <= center <= OH_PEAK_WINDOW[1]):
        out["reject"] = "峰位 %.0f 落在 3200-3700 之外" % center
        return out
    if not (min_fwhm <= fw <= max_fwhm):
        out["reject"] = "峰宽 %.0f 超出 %.0f-%.0f" % (fw, min_fwhm, max_fwhm)
        return out
    if amp < SNR_MIN * sigma:
        out["reject"] = "拟合幅度只有 %.1fσ（需要 ≥3σ）" % (amp / sigma)
        return out
    out["ok"] = True
    return out


def g_band(y, x, sigma, min_snr=SNR_MIN):
    """在一条**平均谱**上找有机质 G 带（约 1600 cm⁻¹）。

    为什么要单独开这条通道：判矿物类别用不上 1600 附近的带，但它是 SHERLOC
    判定有机质的唯一依据——Hurowitz et al. (2025) 在 Bright Angel 组（Cheyava
    Falls / Apollo Temple / Walhalla Glades）就是靠这条带认出有机碳的，而
    Corpolongo et al. (2023) 用的分析区间本来也是 800–1800。我们原来把窗口
    卡在 1500，等于把这条通道整个切掉了。

    两个必须说清的边界：
      1. **同一位置有熔融石英 Si–O 伸缩泛音背景**（仪器光学窗口）。论文是用
         灰色仪器背景谱去比才把 G 带分出来的。所以这里只报"候选 + 实测数字"，
         并附上提醒，不下"有有机质"的结论。
      2. G 带比矿物带弱得多，所以判据和扫描级兜底同源：相对**平均谱自身噪声**
         ≥3σ，且幅度 ≥ G_BAND_MIN_RATIO × 单点噪声中位数。

    返回 None 或 {pos, fwhm, height, snr, ratio}。
    """
    if y is None or not np.isfinite(sigma) or sigma <= 0:
        return None
    sel = (x >= G_BAND[0]) & (x <= G_BAND[1])
    if sel.sum() < 10:
        return None
    # 基线必须用**掩峰迭代**那条（和谱库取峰同源）。普通多项式拟合在这个
    # 120 cm-1 的小窗里会把 G 带本身当成基线吸掉——实测：换之前 1576/1604/1632
    # 三个峰一条都报不出来，换之后立刻出来。
    ym = y - baseline_masked(x, y)
    # 找带的下限放到 2σ：3σ 常常刚好卡在噪声边缘上（实测 sol 370 的扫描平均谱
    # 是 2.99σ，一刀切就把整条通道关掉了）。低于 min_snr 的照样报出来，
    # 但 passed=False，报告里会写明"未过 3σ 门槛"，让人自己判断。
    xs, resid = x[sel], ym[sel]
    bands = find_bands(xs, resid, sigma, min_fwhm=MIN_FWHM_CM,
                       max_fwhm=MAX_FWHM_CM, snr_min=2.0)
    if not bands:
        return None
    b = max(bands, key=lambda d: d["snr"])
    return {"pos": float(b["pos"]), "center": float(b["pos"]),
            "fwhm": float(b["fwhm"]), "height": float(b["height"]),
            "snr": float(b["snr"]),
            "passed": bool(b["snr"] >= min_snr and b["fwhm"] >= MIN_FWHM_CM),
            "ratio_to_point_sigma": float("nan"),   # 由调用方按自己的 σ 填
            "width_censored": b.get('width_censored',False),
            "fwhm_interpolated":b.get('fwhm_interpolated'),
            "window": tuple(G_BAND),
            "note": ("该区同时有熔融石英 Si-O 泛音背景（约 1600 cm-1，"
                     "Hurowitz et al. 2025），判读需对照仪器背景谱")}


def robust_scale(y):
    """稳健的噪声尺度：1.4826 × MAD。局部一条带只占窗口一小部分，不影响它。"""
    y = np.asarray(y, dtype=float)
    if y.size < 8:
        return float("nan")
    return float(1.4826 * np.median(np.abs(y - np.median(y))))


def scan_level(x, spectra, sig_pt, noise_window=NOISE_WINDOW):
    """扫描级弱信号检查：在**全扫描平均谱**上找带。

    什么时候会走到这里：逐点路径过门的点太少（见 SCAN_LEVEL_MAX_POINTS）。
    那时逐点聚类已经没有统计意义，但整片扫描合起来可能还有带。

    两个条件必须同时满足，缺一不可：

    1. **统计显著**：带峰 ≥ 3 × 平均谱自身噪声（稳健尺度）。
    2. **幅度下界**：带峰 ≥ SCAN_LEVEL_MIN_RATIO（0.30）× 单点噪声中位数。

    第 2 条不是装饰：平均谱的噪声按 1/√N 缩小，1296 点的 Survey 扫描里
    区区 7 个计数就能凑出 4σ，只看第 1 条会把"论文说是零"的扫描报成候选。
    阈值用全库标定过（`verify/scan_level_sweep.py`）：已知真检出的这个比值
    在 0.30–4.29，已知"论文说零"的扫描最高只有 0.28。阈值必须卡在两者
    之间，取 0.30——**注意不能照旧文档改成 0.25**：0.25 < 0.28，会把
    已知的阴性（0.28 那次）放进候选。

    返回的候选带**不算矿物类别**，只作为"这里可能有弱带，去看一眼"的提示，
    报告里单独一块列出，不参与与论文的逐点对照。
    """
    spectra = np.asarray(spectra)
    n = spectra.shape[0]
    band = (x >= SCAN_LEVEL_WINDOW[0]) & (x <= SCAN_LEVEL_WINDOW[1])
    xs = x[band]
    if xs.size < 20 or n < 2:
        return None

    sig_med = float(np.nanmedian(sig_pt)) if np.size(sig_pt) else float("nan")
    if not np.isfinite(sig_med) or sig_med <= 0:
        return None

    # 扣基线要按**整条谱**做（和逐点路径、成组路径保持一致），
    # 之后再把 800-1500 切出来找带。水带（3200-3700）必须用整条谱算，
    # 所以这里同时留一份全谱平均值——踩过坑：只给切片会 IndexError。
    sub = np.empty((n, x.size))
    for i in range(n):
        sub[i] = spectra[i] - baseline(x, spectra[i])
    mean_full = sub.mean(axis=0)
    mean = mean_full[band]
    resid = mean - baseline(xs, mean)
    sig_mean = robust_scale(resid)
    if not np.isfinite(sig_mean) or sig_mean <= 0:
        return None

    cands = []
    for b in find_bands(xs, resid, sig_mean, min_fwhm=SCAN_LEVEL_MIN_FWHM,
                        max_fwhm=SCAN_LEVEL_MAX_FWHM, snr_min=SCAN_LEVEL_MIN_SNR):
        ratio = b["height"] / sig_med
        if ratio < SCAN_LEVEL_MIN_RATIO:
            continue
        refined = {"pos": b["pos"], "center": b["pos"], "height": b["height"],
                   "fwhm": b["fwhm"], "snr": b["snr"], "method": "scan_mean",
                   "ratio_to_point_sigma": ratio}
        # 峰心用伪 Voigt 拟合精化（2026-09-20）：argmax 落在离散通道上，
        # 在鼓包肩上会系统性偏 10-20 cm-1；而次级峰的窗口判定（硫酸盐
        # 1110-1160、高氯酸盐 1075-1160）恰恰吃这个精度。采信条件与
        # 逐点主峰同一套：拟合成功、未撞边界、峰宽达标、幅度 ≥3σ_mean，
        # 且峰心离 argmax ≤30 cm-1；任一不满足退回 argmax（与旧版一致）。
        win = ((xs >= b["pos"] - b["fwhm"] - 40.0)
               & (xs <= b["pos"] + b["fwhm"] + 40.0))
        if win.sum() >= 10:
            fit = fit_voigt(xs[win], resid[win], b["pos"], b["fwhm"],
                            min_fwhm=SCAN_LEVEL_MIN_FWHM,
                            max_fwhm=SCAN_LEVEL_MAX_FWHM,
                            within=(b["pos"] - 30.0, b["pos"] + 30.0))
            fit_ok = (fit is not None and not fit["at_bound"]
                      and fit["fwhm"] >= SCAN_LEVEL_MIN_FWHM
                      and fit["amp"] >= 3.0 * sig_mean)
            if fit_ok:
                refined["center"] = fit["center"]
                refined["fwhm"] = fit["fwhm"]
                refined["height"] = fit["amp"]
                refined["method"] = "voigt_scan_mean"
                refined['profile_valid']=fit['profile_valid']
                refined['requires_review']=not fit['profile_valid']
        cands.append(refined)
    if not cands:
        return None

    # 一条宽带上会冒出好几个极大值（1016/1045/1074/1093 其实是同一个鼓包），
    # 按"半高区间是否重叠"合并，只留每簇里最强的那个，免得把一条带报成四条。
    cands.sort(key=lambda c: c["pos"])
    merged = []
    for c in cands:
        lo = c["pos"] - c["fwhm"] / 2.0
        if merged and lo <= merged[-1]["_hi"]:
            merged[-1]["_hi"] = max(merged[-1]["_hi"], c["pos"] + c["fwhm"] / 2.0)
            merged[-1]["_n"] += 1
            if c["snr"] > merged[-1]["snr"]:
                keep = merged[-1]["_n"]
                merged[-1] = dict(c)
                merged[-1]["_hi"] = c["pos"] + c["fwhm"] / 2.0
                merged[-1]["_n"] = keep
        else:
            d = dict(c)
            d["_hi"] = c["pos"] + c["fwhm"] / 2.0
            d["_n"] = 1
            merged.append(d)
    for d in merged:
        d.pop("_hi", None)
        d["merged_count"] = d.pop("_n", 1)
    cands = merged

    # 半谱一致性：把 100 个点对半拆开，两边各自的平均谱上都要看得见这条带。
    # 真带是相干的，平均两半都会留下；几个点的偶然突起只在一边有。
    half = n // 2
    halves = []
    if half >= 4:
        for sl in (slice(0, half), slice(half, n)):
            m = sub[sl].mean(axis=0)[band]
            r2 = m - baseline(xs, m)
            halves.append((r2, robust_scale(r2)))

    sel = (x >= noise_window[0]) & (x <= noise_window[1])
    sig_noise = (float(np.sqrt(np.mean(mean_full[sel] ** 2)))
                 if sel.any() else sig_mean)
    oh = hydration(x, mean_full, sig_noise)
    for c in cands:
        c["half"] = [float("nan"), float("nan")]
        if halves:
            for k, (r2, s2) in enumerate(halves):
                if not np.isfinite(s2) or s2 <= 0:
                    continue
                i0 = int(np.argmin(np.abs(xs - c["center"])))
                lo = max(0, i0 - 3)
                hi = min(xs.size, i0 + 4)
                c["half"][k] = float(r2[lo:hi].max() / s2)
        c["half_min"] = (min(c["half"]) if np.all(np.isfinite(c["half"]))
                         else float("nan"))
        c["classes"] = classify_bands([c], oh)
    cands.sort(key=lambda c: -c["snr"])
    return {"bands": cands, "sigma_mean": sig_mean, "sigma_point": sig_med,
            "n_points": n}


def classify_bands(bands, oh=None, context=None):
    """把波段映射到矿物类别。

    `bands` 是要**定类别**的峰（主峰），`context` 是**只用来查次级峰**的其他峰。
    分开是因为两者角色不同：类别由主峰的位置定，而"能不能定死"（硫酸盐 vs 辉石、
    高氯酸盐 vs 磷酸盐）要看同一张平均谱上别处还有没有配套的次级峰。

    为什么要这个参数（2026-09-15 发现）：以前每次都只传一个主峰进来，
    于是下面 `sulf_secondary` / `naclo4_secondary` 两个判断**永远不可能成立**
    ——"高氯酸钠"这个类别等于死代码，硫酸盐也只能靠水带才判得出来。
    而论文的判据恰恰是"靠次级峰区分"（Corpolongo 2023 §3）。
    """
    out = []
    if not bands:
        return out
    if oh is not None and not oh.get("ok", True):
        oh = None                      # 被规则否掉的水带，不能当证据
    r = lambda k: RANGES[k]
    pool = list(bands) + list(context or [])

    def has(lo, hi):
        return any(lo <= b.get("center", b["pos"]) <= hi for b in pool)

    sulf_secondary = has(1110.0, 1160.0) or has(1200.0, 1240.0)

    # 高氯酸盐的次级峰：论文支持材料 4.8 节原话——高氯酸盐要"有 1000-1150 cm⁻¹
    # 的**低强度**次级峰"才能与磷酸盐区分。窗口收紧到 1075-1160：Clark et al.
    # (2026) 表 1 实测的 7 种高氯酸盐，**有分量的**反对称伸缩都落在这里
    # （Mn 1077/1115、Fe 1080/1121、K 1082/1119、Ca 1083/1117、Na 1091/1146、
    # Al 1107、Mg 1109；只有 Al 的 1052 是 0.038 的弱峰，排在窗外）
    #
    # 为什么不是论文正文说的 1000-1150：那样太松。实测 sol 186（Bellegarde）
    # 的 967 组会被 1035/1064 两个峰带成"钠高氯酸盐"，而论文对 Bellegarde
    # 的结论恰恰是"磷酸盐与高氯酸盐无法区分"。收到 1075-1160 之后，
    # 三个靶点分别落到论文的三种结论上：Guillaumes → 钠高氯酸盐；
    # Bellegarde / Montpezat → 高氯酸盐或磷酸盐（不定）；Dourbes/Quartier → 氯酸盐。
    # "低强度"也要卡：Na perchlorate 实测是主峰的 0.19 和 0.10 倍，
    # 所以要求次级峰的 SNR 落在主峰的 2%-60% 之间——不卡的话，
    # 组平均谱上蹭进来的碳酸盐峰（比主峰还强）会把类别带偏。
    naclo4_secondary = False
    naclo4_used = []
    _ms = (bands[0].get("obs_snr") or 0) if bands else 0
    if _ms > 0:
        for b in pool:
            q = b.get("center", b["pos"])
            if 1075.0 <= q <= 1160.0 and 0.02 * _ms <= (b.get("snr", 0) or 0) <= 0.6 * _ms:
                naclo4_secondary = True
                naclo4_used.append(q)

    for b in bands:
        # 论文用高斯拟合定峰心；带宽内的强度质心最接近那个值
        p = b.get("center", b["pos"])
        ev = []
        if r("olivine")[0] <= p <= r("olivine")[1]:
            key = "olivine"
            ev.append("峰位落在橄榄石区间 815-865 cm-1")
        elif r("perchlorate_chlorate")[0] <= p <= r("perchlorate_chlorate")[1]:
            key = "perchlorate_chlorate"
            ev.append("峰位 920-945，与多种氧氯酸盐主峰一致")
        elif r("perchlorate_or_phosphate")[0] <= p <= r("perchlorate_or_phosphate")[1]:
            if naclo4_secondary:
                key = "sodium_perchlorate"
                ev.append("主峰 945-975，且有 1050-1160 的低强度次级峰（实测 %s，"
                          "为主峰的 2%%%%-60%%%%）——论文支持材料 4.8 节的钠高氯酸盐判据"
                          % "、".join("%.0f" % v for v in sorted(set(naclo4_used))[:4]))
            else:
                key = "perchlorate_or_phosphate"
                ev.append("主峰 945-975，无可用的低强度次级峰，"
                          "高氯酸盐与磷酸盐无法区分")
        elif r("pyroxene")[0] <= p <= r("sulfate")[1]:
            if sulf_secondary or oh is not None:
                key = "sulfate"
                if sulf_secondary:
                    used = [b["pos"] for b in pool
                            if 1110.0 <= b.get("center", b["pos"]) <= 1160.0
                            or 1200.0 <= b.get("center", b["pos"]) <= 1240.0]
                    ev.append("主峰 1000-1035，且有 1110-1160 或 1200-1240 次级峰"
                              "（组平均谱上实测 %s）"
                              % "、".join("%.0f" % v for v in sorted(set(used))[:4]))
                if oh is not None:
                    ev.append("观察到 O-H 带 %.0f cm-1，支持含水硫酸盐" % oh["pos"])
            else:
                key = "sulfate_or_pyroxene"
                ev.append("主峰 985-1035，既无硫酸盐次级峰也无含水特征，与辉石无法区分")
        elif r("silicate")[0] <= p <= r("silicate")[1]:
            if b["fwhm"] >= 120.0:
                key = "silicate"
                ev.append("1050 附近宽带，FWHM %.0f cm-1，与未定硅酸盐一致" % b["fwhm"])
            else:
                key = "carbonate" if p >= 1070.0 else "silicate"
                ev.append("峰位在硅酸盐与碳酸盐的过渡带")
        elif r("carbonate")[0] <= p <= r("carbonate")[1]:
            key = "carbonate"
            ev.append("峰位 1070-1100，与方解石标准 1084 cm-1 一致")
        else:
            key = "unassigned"
            ev.append("峰位不在论文报告的任何类别区间内")

        for a in ARTIFACT_CM:
            if abs(p - a) <= 4.0:
                ev.append("注意：峰位在 %.0f cm-1 附近，暗背景扣除常在此留假峰" % a)
        for a in ARTIFACT_OTHER_CM:
            if abs(p - a) <= 4.0:
                ev.append("注意：%.0f cm-1 是论文列出的单像素伪影位置"
                          "（Supporting Information Text S3），不要当矿物带" % a)
        # 两处熔融石英背景（仪器光学窗口），判读前必须对照仪器背景谱
        if abs(p - 800.0) <= 30.0:
            ev.append("注意：该区有熔融石英光学窗口的背景带（约 800 cm-1，"
                      "Corpolongo 2023 §2.4），与硅酸盐/橄榄石的低波数侧重叠；"
                      "可在工作台预处理勾选「仪器本底扣除」后对照复看")
        if abs(p - 1600.0) <= 60.0:
            ev.append("注意：该区有熔融石英 Si-O 伸缩泛音背景（约 1600 cm-1，"
                      "Hurowitz 2025），与有机质 G 带重叠，判读需对照仪器背景谱")
        out.append({"key": key, "class": LABELS.get(key, key),
                    "peak": b, "evidence": ev})
    return sorted(out, key=lambda c: -c["peak"]["snr"])


def analyse_scan(x, spectra, noise_window=NOISE_WINDOW, cluster_cm=15.0, min_points=2,
                 min_fwhm=MIN_FWHM_CM, robust=True, sigma_factor=SIGMA_OUTLIER_FACTOR,
                 snr_min=SNR_MIN, prom_sigma=PROM_SIGMA):
    """按论文流程分析一次扫描。

    第一步逐点找合格波段；第二步把峰位相近的点分组，只对有信号的点求平均谱，
    用平均谱定类别。第二步很关键：把所有点一起平均会把信号稀释掉。
    """
    spectra = np.asarray(spectra)
    n = spectra.shape[0]
    sel = (x >= noise_window[0]) & (x <= noise_window[1])

    # **两套残差，各司其职**（2026-09-16）：
    #   subs      —— 在**整条谱**上扣基线。**逐点找峰必须用它**：拉曼区的宽带
    #                （未定硅酸盐 ~1050）会被"只在拉曼区拟合的基线"吸掉，
    #                单点就一个峰都过不了门（实测 sol 162 直接变成 0 个组）。
    #   subs_band —— 只在**拉曼区**（FINGERPRINT）内扣基线。**组平均谱上的
    #                观测峰表用它**：论文支持材料的流程就是"background
    #                subtraction in the Raman region"，而且只有这样弱次级峰
    #                才露得出来——实测 sol 162（Guillaumes）的高氯酸盐组，
    #                论文点名的 1091.7 与 1150.4 在拉曼区残差上有 3.9σ/4.0σ，
    #                在整条谱残差上连 3σ 都够不到。
    bnd = (x >= FINGERPRINT[0]) & (x <= FINGERPRINT[1])
    xb = x[bnd]
    subs = np.empty_like(spectra)
    subs_band = np.empty((n, xb.size))
    full = np.empty_like(spectra)
    sig = np.empty(n)
    for i in range(n):
        ys = spectra[i] - baseline(x, spectra[i])
        full[i] = ys
        subs[i] = ys
        subs_band[i] = spectra[i, bnd] - baseline(xb, spectra[i, bnd])
        sig[i] = float(np.sqrt(np.mean(ys[sel] ** 2))) if sel.any() else np.nan

    # 手册 A-5：σ 被荧光抬高的光谱改用扫描中位 σ，否则那些点会被误判成"没信号"
    if robust:
        sig_used, sigma_replaced = robust_sigma(sig, sigma_factor)
    else:
        sig_used, sigma_replaced = sig, []
    sigma_scan = float(np.median(sig_used)) if n else float("nan")

    mean_sub = spectra.mean(axis=0) - baseline(x, spectra.mean(axis=0))
    med_sub = np.median(subs, axis=0)

    # 第一步：逐点找波段
    rw = (x >= FINGERPRINT[0]) & (x <= FINGERPRINT[1])
    hits = []
    for i in range(n):
        for b in find_bands(x[rw], subs[i][rw], sig_used[i],
                            min_fwhm=min_fwhm, snr_min=snr_min,
                            prom_sigma=prom_sigma):
            b["pt"] = i
            hits.append(b)

    # 第二步：一维间隙聚类。相邻峰位差不超过 cluster_cm 归一组，
    # 同时限制整组跨度，否则会一路串连把相隔几百波数的峰混在一起。
    span_cm = 4.0 * cluster_cm
    hits.sort(key=lambda b: b["pos"])
    clusters = []
    for b in hits:
        if clusters:
            cur = clusters[-1]
            if (b["pos"] - cur[-1]["pos"]) <= cluster_cm and (b["pos"] - cur[0]["pos"]) <= span_cm:
                cur.append(b)
                continue
        clusters.append([b])

    def build_group(cl, allow_small=False):
        """把一簇波段折叠成一个「组」：求平均谱、拟合、定类别。"""
        pts = sorted({b["pt"] for b in cl})
        if len(pts) < (1 if allow_small else min_points):
            return None
        pos0 = float(np.median([b["pos"] for b in cl]))
        gm = subs[pts].mean(axis=0)              # 整条谱残差的组平均（拟合用）
        gm_full = full[pts].mean(axis=0)         # 同上，噪声/水带用
        gm_band = subs_band[pts].mean(axis=0)    # 拉曼区残差的组平均（观测峰表用）
        g_sigma = (float(np.sqrt(np.mean(gm_full[sel] ** 2)))
                   if sel.any() else np.nan)

        # 拟合窗口按本簇自己的波段范围定，而不是固定的 ±90。
        # 固定窗口会"看见"隔壁簇的峰，于是用邻居的峰来解释自己——
        # 实测簇心 1141 的组被拟到 1127，就是因为窗口伸进了 1083 的碳酸盐峰。
        band_lo = float(np.min([b["pos"] for b in cl]))
        band_hi = float(np.max([b["pos"] for b in cl]))
        win_lo = min(band_lo - 40.0, pos0 - 60.0)
        win_hi = max(band_hi + 40.0, pos0 + 60.0)
        near = (x >= win_lo) & (x <= win_hi)

        # Voigt 拟合定峰心/峰宽（手册 A-4）。拟合失败、撞边界、
        # 或峰宽比判据下限还窄，都算失败，退回半高穿零法。
        fw0 = float(np.median([b["fwhm"] for b in cl]))
        fit = fit_voigt(x[near], gm[near], pos0, fw0, min_fwhm=min_fwhm,
                        max_fwhm=FIT_MAX_FWHM,
                        within=(pos0 - FIT_CENTER_TOL, pos0 + FIT_CENTER_TOL))
        # 拟合要同时满足四件事才采信：不是 None、没撞边界、峰宽达标、
        # 而且**幅度显著**（≥3σ）。最后一条最关键——平均谱上一个 SNR<3 的
        # 包络不是真峰，采信它就会把峰心带到尾巴上去。
        fit_ok = (fit is not None and not fit["at_bound"] and fit["fwhm"] >= min_fwhm
                  and fit["amp"] >= snr_min * g_sigma)
        if fit_ok:
            refined = {"pos": pos0, "center": fit["center"], "height": fit["amp"],
                       "fwhm": fit["fwhm"], "eta": fit["eta"], "rms": fit["rms"],
                       "method": "voigt",
                       "snr": float(fit["amp"] / g_sigma) if g_sigma else float("nan")}
            refined.update(profile_valid=fit['profile_valid'],profile_model=fit['model'],
                           gaussian_coefficient=fit['c1'],lorentz_coefficient=fit['c2'])
        else:
            if fit is None:
                why = "拟合不出结果"
            elif fit["at_bound"]:
                why = "撞边界"
            elif fit["fwhm"] < min_fwhm:
                why = "峰宽 %.0f 低于判据下限" % fit["fwhm"]
            else:
                why = "拟合幅度只有 %.1fσ，不显著" % (fit["amp"] / g_sigma
                                                if g_sigma else float("nan"))
            # 退路：不采信拟合时，峰心直接取本簇波段位置的中位数。
            # 这比"平均谱的 argmax"稳得多——波段分布是直接证据，
            # 平均谱上的包络还可能被基线/相邻峰带偏。
            seg, xs = gm[near], x[near]
            k = int(np.argmin(np.abs(xs - pos0)))
            while k < len(seg) - 1 and seg[k + 1] > seg[k]:
                k += 1
            half = seg[k] / 2.0
            L = R = k
            while L > 0 and seg[L] > half:L -= 1
            while R < len(seg)-1 and seg[R] > half:R += 1
            measured = half_maximum(xs,seg,k)
            refined = {"pos": pos0, "center": float(pos0), "height": float(seg[k]),
                       "fwhm": float(xs[R]-xs[L]), "method": "halfmax",
                       "fwhm_interpolated": measured['fwhm'] if measured else None,
                       "width_censored":measured is None,
                       "width_method":"sample_bracket_upper_bound",
                       "fallback": why,
                       "snr": float(seg[k] / g_sigma) if g_sigma else float("nan")}
        oh = hydration(x, gm_full, g_sigma)

        # 给谱库比对用的观测峰集合：在**组平均谱**上找 800-1800 内的全部合格波段。
        # 为什么不能只用聚类簇里的那些峰：簇内峰都聚在主峰附近（±15 cm-1 内），
        # 而判类需要的次级峰（碳酸盐 1432/1745、硫酸盐 1134）离得很远，不在簇里。
        # 取谱库比对用的峰之前，先在组平均谱上做一次**掩峰迭代基线**（论文 2.3 节
        # 第 5 步）。原来那版基线被宽包带偏，残差里只有主峰一个鼓包，
        # 次级峰（碳酸盐 ν3 ~1432、硫酸盐 ν3 ~1134）根本露不出来。
        # **只扣一次基线**。gm 本身已经是逐点扣完基线之后的平均，再减一次
        # 掩峰基线等于扣两次——实测（2026-09-15）后果很严重：sol 349 那个
        # 22 点的高氯酸盐组，二次基线之后 800–1800 内**一个合格峰都不剩**，
        # sol 162 的两个组也从 4-5 个峰掉到 1-2 个，谱库比对直接变成空白。
        # 二次基线当初是为了让碳酸盐 1432、硫酸盐 1134 这些次级峰从宽包肩上
        # 露出来；现在次级峰能直接看到（1091.7、1150.4、1111.3 都在 3σ 以上），
        # 那一步就不再需要了。
        gm_lib = gm_band
        lib_sel = np.ones(xb.size, dtype=bool)      # gm_band 已经是拉曼区切片
        # 噪声必须用**标准口径**那个 sigma（g_sigma，逐点扣基线后在安静区取 RMS），
        # 不能用掩峰基线之后的残差算——掩峰基线会把安静区也压平，
        # sigma 掉到极小，噪声就冒充成十几个 sigma 的假峰（实测踩过：
        # sol 370 的组平均谱真实噪声 885 计数、1600 cm-1 那个"峰"其实只有 0.3σ，
        # 但用掩峰后的 sigma 一算就成了 18σ）。
        lib_bands = find_bands(xb[lib_sel], gm_lib[lib_sel], g_sigma,
                               min_fwhm=min_fwhm, snr_min=LIB_SNR_MIN)
        # 峰心用伪 Voigt 精化（2026-09-20）：find_bands 给的是 argmax，
        # 在鼓包肩上会系统性偏 10-20 cm-1——而次级峰窗口判定（硫酸盐
        # 1110-1160、高氯酸盐 1075-1160）和谱库比对都吃这个精度。
        # 采信条件沿用逐点主峰那套：拟合成功、未撞边界、峰宽达标、
        # 幅度 ≥3σ；失败退回 argmax（ctx 里的 center 同步受益）。
        for _b in lib_bands:
            _b.setdefault("center", _b["pos"])
            _w = ((xb[lib_sel] >= _b["pos"] - _b["fwhm"] - 40.0)
                  & (xb[lib_sel] <= _b["pos"] + _b["fwhm"] + 40.0))
            if _w.sum() < 10:
                continue
            _f = fit_voigt(xb[lib_sel][_w], gm_lib[lib_sel][_w],
                           _b["pos"], _b["fwhm"], min_fwhm=min_fwhm,
                           within=(_b["pos"] - 30.0, _b["pos"] + 30.0))
            if (_f is not None and not _f["at_bound"]
                    and _f["fwhm"] >= min_fwhm
                    and _f["amp"] >= 3.0 * g_sigma):
                _b["center"] = _f["center"]
                _b["fwhm"] = _f["fwhm"]
                _b["height"] = _f["amp"]
                _b['profile_valid']=_f['profile_valid']
        # 最多留最强的 4 个峰：组平均谱上弱峰多半是噪声，带进比对只会拉低得分。
        # 带上 SNR 当权重交给比对——主峰和杂峰的份量不一样。
        lib_bands.sort(key=lambda b: -b["snr"])
        seen_cm = set()
        observed = []
        for b in lib_bands[:4]:
            key_cm = round(b["center"])
            if key_cm in seen_cm:
                continue
            seen_cm.add(key_cm)
            observed.append((key_cm, round(float(b["snr"]), 1)))
        observed.sort()

        # 次级峰上下文：组平均谱上除主峰以外的那些峰。
        # 只交给 classify_bands 用来查"配套峰"，不参与定类别本身。
        # （主峰由 classify_bands 自己的 b 提供，两者的 SNR 口径要通过
        #  refined["obs_snr"] 对齐——见下面）
        ctx = [{"pos": float(p), "center": float(p), "fwhm": 40.0,
                "height": 1.0, "snr": float(s)}
               for p, s in observed
               if abs(float(p) - refined["center"]) > 15.0]
        # 主峰在**同一份观测峰表**里的 SNR：高氯酸盐那条判据要拿次级峰跟主峰比
        # 相对强度，用拟合幅度去比原始峰高是不自洽的。
        # **已知伪影位置不算主峰**（论文支持材料 Text S3 的那份清单）：
        # 实测 sol 186 的 Survey 组就是被 977（981 伪影 ±4 内）顶成了"主峰"，
        # 于是拿 1141 当次级峰、判成钠高氯酸盐；而论文对那个扫描的结论是硫酸盐。
        _art = ARTIFACT_CM + ARTIFACT_OTHER_CM
        refined["obs_snr"] = max(
            (s for p, s in observed
             if 940.0 <= p <= 980.0
             and not any(abs(p - a) <= 4.0 for a in _art)), default=None)
        classes = classify_bands([refined], oh, context=ctx)
        review_reasons=[]
        if refined['method']=='halfmax':
            review_reasons.append('峰形拟合未通过；当前类别依赖近似半高宽和聚类中心')
        if refined.get('profile_valid') is False:
            review_reasons.append('高斯/洛伦兹系数含负值：仅为经验候选定位，不是物理伪Voigt，需做非负约束复核')
        if oh.get('ok') and oh.get('requires_review'):
            review_reasons.append('用于含水判据的O-H带依赖近似或非物理拟合，不能单独作为确定物相证据')
        if any(b.get('profile_valid') is False for b in lib_bands[:4]):
            review_reasons.append('观测峰表中含负分量拟合，次级峰归属需要复核')
        if refined.get('width_censored'):
            review_reasons.append('测量窗口未包含两侧完整半高交点，峰宽不能当精确测量值')
        if refined.get('snr',0)<snr_min:
            review_reasons.append('组平均峰未达到当前SNR门槛')
        if review_reasons:
            for c in classes:
                c['requires_review']=True
                c['evidence'].extend(review_reasons)

        return {
            "pos": pos0, "n_points": len(pts), "pts": pts,
            "selection_basis": "按逐点候选峰筛选成员后计算组平均谱；组平均SNR不是独立验证，也不是矿物置信度",
            # 真正进了这个簇的波段位置。诊断聚类是否乱合并就看这一串，
            # 不能看"每个点的最佳波段"——那是另一个量。
            "bands": cl,
            "band_pos": [float(b["pos"]) for b in cl],
            "band_lo": band_lo, "band_hi": band_hi,
            "sigma_mean": g_sigma, "band": refined, "oh": oh,
            "classes": classes,
            "requires_review":bool(review_reasons),"review_reasons":review_reasons,
            "class_context": [float(b["pos"]) for b in ctx],
            # 有机质 G 带：单开一条通道，不参与矿物判类
            "g_band": g_band(gm_full, x, g_sigma),
            "observed_peaks": observed,
            "library": library.best_per_class(observed)[:3] if observed else [],
            "split_into": None,
        }

    groups = []
    for cl in clusters:
        g = build_group(cl)
        if g is not None:
            groups.append(g)

    # 「一个峰还是两个峰」检验（手册 A-4：用所需的最小峰数）。
    # 这一步专门解决 1050 的硅酸盐宽带和 1083 的碳酸盐峰被并成一组的问题。
    # 判据是形状而不是宽度——硫酸盐那组也宽 58 cm-1，但它是货真价实的单峰。
    split_count = 0
    final_groups = []
    for g in groups:
        if g["n_points"] < SPLIT_MIN_POINTS:
            g["split_test"] = {"ok": False, "reject": "只有 %d 个点，不拆" % g["n_points"]}
        elif g["band"]["method"] != "voigt":
            g["split_test"] = {"ok": False, "reject": "单组分拟合已退回半高法，不拆"}
        else:
            near = ((x >= min(g["band_lo"] - 40.0, g["pos"] - 60.0))
                    & (x <= max(g["band_hi"] + 40.0, g["pos"] + 60.0)))
            gm = subs[g["pts"]].mean(axis=0)
            two = fit_two_voigt(x[near], gm[near], g["pos"],
                                g["band"]["fwhm"], g["sigma_mean"],
                                min_fwhm=min_fwhm,
                                within=(g["band_lo"] - 15.0, g["band_hi"] + 15.0))
            g["split_test"] = two
            if two is not None and two.get("ok"):
                c_a, c_b = two["comps"][0]["center"], two["comps"][1]["center"]
                cl_a, cl_b = [], []
                for b in g["bands"]:
                    (cl_a if abs(b["pos"] - c_a) <= abs(b["pos"] - c_b) else cl_b).append(b)
                ga = build_group(cl_a) if cl_a else None
                gb = build_group(cl_b) if cl_b else None
                if ga is not None and gb is not None:
                    if not two.get('profile_valid',True) or not two.get('components_significant',True):
                        warning='拆分模型含负峰形分量或弱组分，拆分仅作候选，需物理约束复核'
                        for child in (ga,gb):
                            child['requires_review']=True;child['review_reasons'].append(warning)
                            for c in child['classes']:
                                c['requires_review']=True;c['evidence'].append(warning)
                    ga["split_into"] = [c_a, c_b]
                    gb["split_into"] = [c_a, c_b]
                    final_groups.extend([ga, gb])
                    split_count += 1
                    continue
        final_groups.append(g)
    groups = final_groups
    groups.sort(key=lambda g: -g["n_points"])

    # 逐点明细：每个点自己的最优合格波段，加上它落在哪个聚类组、被判成什么类。
    # 这张表是和论文逐点对照的唯一依据，所以每个字段都要能追到源头。
    hit_count = {}
    best_hit = {}
    for b in hits:
        i = int(b["pt"])
        hit_count[i] = hit_count.get(i, 0) + 1
        cur = best_hit.get(i)
        if cur is None or (b["fwhm"], b["height"]) > (cur["fwhm"], cur["height"]):
            best_hit[i] = b

    point_classes = []
    point_detail = []
    for i in range(n):
        gl = [g for g in groups if i in g["pts"]]
        uniq = []
        seen = set()
        for g in gl:
            for c in g["classes"]:
                if c["key"] in seen:
                    continue
                seen.add(c["key"])
                uniq.append({"key": c["key"], "class": c["class"],
                             "group_pos": float(g["pos"]),
                             "group_center": float(g["band"]["center"]),
                             "requires_review":g.get('requires_review',False)})
        b = best_hit.get(i)
        point_detail.append({
            "point": i,
            "sigma": float(sig[i]),
            "n_bands": int(hit_count.get(i, 0)),
            "pos": float(b["pos"]) if b else None,
            "fwhm": float(b["fwhm"]) if b else None,
            "snr": float(b["snr"]) if b else None,
            "height": float(b["height"]) if b else None,
            "groups": [float(g["pos"]) for g in gl],
            "bands": [{"pos": float(b["pos"]), "fwhm": float(b["fwhm"]),
                       "snr": float(b["snr"])} for b in hits if int(b["pt"]) == i],
            "classes": uniq,
            "labeled": bool(gl),
        })
        point_classes.append(sorted({c["class"] for c in uniq}))

    allpts = sorted({p for g in groups for p in g["pts"]})
    if allpts:
        scan_mean = subs[allpts].mean(axis=0)
        scan_sigma = float(np.sqrt(np.mean(scan_mean[sel] ** 2))) if sel.any() else np.nan
        scan_bands = find_bands(x[rw], scan_mean[rw], scan_sigma,
                                min_fwhm=min_fwhm, snr_min=snr_min)
        scan_oh = hydration(x, scan_mean, scan_sigma)
        scan_classes = classify_bands(scan_bands, scan_oh)
    else:
        scan_mean, scan_sigma, scan_bands, scan_oh, scan_classes = None, np.nan, [], None, []

    # 扫描级兜底：进入峰组的点太少时，再看一眼**全部点**的平均谱。
    # 注意和上面的 scan_* 不是一回事：那三个量只平均"成了组的点"，
    # 这里平均全部点——正因为过门的点太少，才需要退到全部点上来。
    n_with_band = int(sum(1 for p in point_classes if p))
    n_candidate_points = len(hit_count)
    n_named_points = sum(any(c['key'] != 'unassigned' for c in p['classes']) for p in point_detail)
    scan_lv = None
    if n_with_band < SCAN_LEVEL_MAX_POINTS:
        try:
            scan_lv = scan_level(x, spectra, sig_used, noise_window=noise_window)
        except Exception:                            # noqa: BLE001
            scan_lv = None
    if scan_lv is not None:
        scan_lv["why"] = ("只有 %d 个点进入峰组（少于 %d），"
                          "所以补看全扫描平均谱" % (n_with_band, SCAN_LEVEL_MAX_POINTS))

    # 有机质 G 带通道：在**全扫描平均谱**上找 ~1600 的带，每次扫描都做。
    # 与扫描级兜底分开：兜底是"矿物判不出来时的替补"，这条是有机质的独立通道。
    scan_mean_all = full.mean(axis=0)            # 全谱长度，G 带和水带都要用
    sig_mean_all = (float(np.sqrt(np.mean(scan_mean_all[sel] ** 2)))
                    if sel.any() else float("nan"))
    gb_scan = g_band(scan_mean_all, x, sig_mean_all)
    if gb_scan is not None:
        sig_pt_med = float(np.nanmedian(sig_used)) if n else float("nan")
        gb_scan["ratio_to_point_sigma"] = (gb_scan["height"] / sig_pt_med
                                           if sig_pt_med else float("nan"))
        gb_scan["passed"] = bool(gb_scan["passed"]
                                 and gb_scan["ratio_to_point_sigma"] >= G_BAND_MIN_RATIO)
    for g in groups:
        gb = g.get("g_band")
        if gb is not None:
            sig_pt_med = float(np.nanmedian(sig_used)) if n else float("nan")
            gb["ratio_to_point_sigma"] = (gb["height"] / sig_pt_med
                                          if sig_pt_med else float("nan"))
            gb["passed"] = bool(gb["passed"]
                                and gb["ratio_to_point_sigma"] >= G_BAND_MIN_RATIO)

    return {
        "x": x, "mean_sub": mean_sub, "median_sub": med_sub,
        "sigma_scan": sigma_scan, "sigma_per_point": sig,
        "sigma_used": sig_used, "sigma_replaced": sigma_replaced,
        "noise_window": tuple(noise_window), "min_fwhm": float(min_fwhm),
        "robust_sigma": bool(robust), "sigma_factor": float(sigma_factor),
        "snr_min": float(snr_min),
        "groups": groups,
        "n_split": int(split_count),
        "scan_mean_sub": scan_mean, "scan_sigma": scan_sigma,
        "scan_bands": scan_bands, "scan_oh": scan_oh, "scan_classes": scan_classes,
        "scan_level": scan_lv,
        "g_band_scan": gb_scan,
        "point_classes": point_classes,
        "point_detail": point_detail,
        "n_points_with_band": n_with_band,
        # 旧字段兼容历史调用/兜底门槛，其实际含义为进入峰组的点数。
        "n_points_in_groups": n_with_band,
        "n_points_with_candidate_bands": n_candidate_points,
        "n_points_with_named_classes": int(n_named_points),
    }
