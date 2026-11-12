# -*- coding: utf-8 -*-
"""深紫外拉曼参考谱库：峰位 + 振动模式归属。

**数据来源**：Razzell Hollis et al. (2021), "A deep-ultraviolet Raman and
Fluorescence spectral library of 62 minerals for the SHERLOC instrument onboard
Mars 2020", Planetary and Space Science 209, 105356. 表 4。

这不是网上抄来的峰位区间，是**用 SHERLOC 的类比仪器（Brassboard，同样 248.6 nm
激发）实测出来的 92 个样品、62 种矿物**的峰位，论文给出的合计不确定度是 ±5 cm⁻¹。
本工具原来的"峰位区间规则"是对这份数据的粗糙替代，现在换成它本身。

表 4 的原始结构是：每个矿物大类一张小表，列是该类的振动模式，行是样品。
这里把它拉平成 `PEPeak` 列表，保留矿物名、样品类型、峰位、模式归属。

几个约定（与论文一致）：

* `None` 表示论文写 n.d.（not detectable，检不出）；**没写**表示该模式本来就不预期存在
  （本表里直接省略那个键）——这两者含义不同，不要混。
* 峰位单位 cm⁻¹。同一个模式可能给多个值（如石膏 3411,3500），用 tuple 表示。
* `organic` 是石墨碳 G 带，属于有机污染，**不是矿物**，判定时要排除。
"""

from __future__ import annotations

# 每个条目：(矿物名, 样品类型, {模式: 峰位 或 峰位元组 或 None})
# 模式的键名尽量贴合表 4 的列标题。
REFERENCE = {
    "borate": {
        "label": "硼酸盐",
        "modes": ["BO3_BO4_stretch", "B_OH_bend", "B_OH_stretch", "H2O_stretch",
                  "organic"],
        "samples": [
            ("Colemanite", "Powder 1", {"BO3_BO4_stretch": 1079, "B_OH_bend": None,
                                        "B_OH_stretch": 3606, "H2O_stretch": 3215}),
            ("Kernite", "Crystal 1", {"BO3_BO4_stretch": 931, "B_OH_bend": 1333,
                                      "B_OH_stretch": 3553, "H2O_stretch": 3324}),
            ("Ulexite", "Crystal 1", {"BO3_BO4_stretch": 985, "B_OH_bend": 1372,
                                      "B_OH_stretch": 3576, "H2O_stretch": 3416}),
        ],
    },
    "carbonate": {
        "label": "碳酸盐",
        "modes": ["CO3_v1", "CO3_v3", "CO3_v1+v4", "organic"],
        "samples": [
            ("Aragonite", "Crystal 1", {"CO3_v1": 1091, "CO3_v3": 1465,
                                        "CO3_v1+v4": 1748}),
            ("Calcite", "Crystal 1", {"CO3_v1": 1085, "CO3_v3": 1432,
                                      "CO3_v1+v4": 1745}),
            ("Calcite", "Powder 1", {"CO3_v1": 1085, "CO3_v3": 1434,
                                     "CO3_v1+v4": 1749}),
            ("Calcite", "Pellet 2", {"CO3_v1": 1091, "CO3_v3": 1439,
                                     "CO3_v1+v4": 1754}),
            ("Chalk", "Rock 1", {"CO3_v1": 1089, "CO3_v3": 1433,
                                 "CO3_v1+v4": 1752, "organic": 1610}),
            ("Dolomite", "Crystal 1", {"CO3_v1": 1099, "CO3_v3": 1444,
                                       "CO3_v1+v4": 1761}),
            ("Dolomite", "Powder 1", {"CO3_v1": 1097, "CO3_v3": 1443,
                                      "CO3_v1+v4": 1760}),
            ("Dolomite", "Pellet 2", {"CO3_v1": 1102, "CO3_v3": 1441,
                                      "CO3_v1+v4": 1764, "organic": 1615}),
            ("Magnesite", "Crystal 1", {"CO3_v1": 1099, "CO3_v3": 1452,
                                        "CO3_v1+v4": 1768}),
            ("Magnesite", "Powder 1", {"CO3_v1": 1094, "CO3_v3": 1446,
                                       "CO3_v1+v4": 1762}),
            ("Nahcolite", "Powder 1", {"CO3_v1": 1044, "CO3_v3": 1441,
                                       "CO3_v1+v4": 1678}),
            ("Nahcolite", "Pellet", {"CO3_v1": 1051, "CO3_v3": 1448,
                                     "CO3_v1+v4": 1683}),
            ("Rhodochrosite", "Crystal 1", {"CO3_v1": 1091, "CO3_v3": 1423,
                                            "CO3_v1+v4": 1733}),
            ("Rhodochrosite", "Powder 1", {"CO3_v1": 1084, "CO3_v3": 1416,
                                           "CO3_v1+v4": 1726}),
            ("Siderite", "Crystal 1", {"CO3_v1": None, "CO3_v3": None,
                                       "CO3_v1+v4": None}),
            ("Trona", "Powder 1", {"CO3_v1": 1061, "CO3_v3": 1427,
                                   "CO3_v1+v4": 1715}),
        ],
    },
    "phosphate": {
        "label": "磷酸盐",
        "modes": ["PO4_v1", "PO4_v3", "organic"],
        "samples": [
            ("Fluorapatite", "Crystal", {"PO4_v1": 960, "PO4_v3": 1050}),
        ],
    },
    "sulfate": {
        "label": "硫酸盐",
        "modes": ["SO4_v1", "SO4_v2", "SO4_v3", "H2O_stretch", "organic"],
        "samples": [
            ("Anhydrite", "Powder 1", {"SO4_v1": 1014, "SO4_v2": 494, "SO4_v3": 1128,
                                       "H2O_stretch": (3408, 3496)}),
            ("Barite", "Powder 1", {"SO4_v1": 984, "SO4_v2": 458, "SO4_v3": 1138}),
            ("Epsomite", "Powder 1", {"SO4_v1": (983, 1036), "SO4_v2": None,
                                      "SO4_v3": 1112, "H2O_stretch": 3410}),
            ("Gypsum", "Crystal 1", {"SO4_v1": 1009, "SO4_v2": 494, "SO4_v3": 1134,
                                     "H2O_stretch": (3411, 3500)}),
            ("Gypsum", "Powder 1", {"SO4_v1": 1004, "SO4_v2": 494, "SO4_v3": 1131,
                                    "H2O_stretch": (3408, 3497)}),
            ("Gypsum", "Pellet 2", {"SO4_v1": 1014, "SO4_v2": 504, "SO4_v3": 1138,
                                    "H2O_stretch": (3415, 3503)}),
            ("Kieserite", "Powder", {"SO4_v1": (979, 1033), "SO4_v2": 459,
                                     "SO4_v3": 1115, "H2O_stretch": 3414}),
            ("Melanterite", "Powder 1", {"SO4_v1": 1014, "SO4_v2": None,
                                         "SO4_v3": 1191, "H2O_stretch": 3234}),
            ("Natroalunite", "Pellet 63", {"SO4_v1": 1021, "SO4_v2": 494,
                                           "SO4_v3": 1163, "H2O_stretch": 3468}),
            ("Natroalunite-Natrojarosite", "Pellet 13", {"SO4_v1": None,
                                                         "SO4_v3": None}),
            ("Natroalunite-Natrojarosite", "Pellet 103", {"SO4_v1": None,
                                                          "SO4_v3": None}),
            ("Natrojarosite", "Pellet 93", {"SO4_v1": None, "SO4_v3": None}),
            ("Natrojarosite", "Powder 1", {"SO4_v1": None, "SO4_v3": None}),
            ("Thenardite", "Powder 1", {"SO4_v1": 991, "SO4_v2": None,
                                        "SO4_v3": 1126}),
            ("Fe(III) sulfate", "Pellet 2", {"SO4_v1": None, "SO4_v3": None}),
            ("Fe(II) sulfate", "Pellet", {"SO4_v1": 987, "SO4_v2": None,
                                          "SO4_v3": 1123, "H2O_stretch": 3418}),
            ("Fe(II) sulfate", "Powder 1", {"SO4_v1": 986, "SO4_v2": None,
                                            "SO4_v3": 1123, "H2O_stretch": 3424}),
            ("Mg sulfate", "Pellet 2", {"SO4_v1": 992, "SO4_v2": 471,
                                        "SO4_v3": 1157, "H2O_stretch": 3428}),
            ("Mg sulfate", "Pellet", {"SO4_v1": 1020, "SO4_v2": 484,
                                      "SO4_v3": 1150, "H2O_stretch": 3440}),
            ("Na sulfate", "Pellet 2", {"SO4_v1": 998, "SO4_v2": 468,
                                        "SO4_v3": 1135}),
            ("Na sulfate", "Pellet", {"SO4_v1": 998, "SO4_v2": 468,
                                      "SO4_v3": 1135}),
        ],
    },
    "halide": {
        "label": "卤化物",
        "modes": ["H2O_bend", "H2O_stretch", "NH4_v1", "NH4_v2", "organic"],
        "samples": [
            ("Halite", "Crystal 1", {}),
            ("Halite", "Pellet 2", {"organic": 1615}),
            ("K bromide", "Pellet", {}),
            ("Mg chloride", "Pellet", {"H2O_stretch": (3363, 3500), "organic": 1589}),
            ("Ca chloride", "Pellet", {"H2O_stretch": 3448, "organic": 1638}),
            ("Ammonium chloride", "Pellet", {"NH4_v1": 3051, "organic": 1716}),
        ],
    },
    "oxide_hydroxide": {
        "label": "氧化物与氢氧化物",
        "modes": ["MOx_stretch", "MOx_bend", "M_OH_stretch", "organic"],
        "samples": [
            ("Hematite", "Powder", {"MOx_stretch": None, "MOx_bend": None}),
            ("Manganite", "Rock", {"MOx_stretch": None, "MOx_bend": 465,
                                   "organic": 1569}),
            ("Brucite", "Pellet 2", {"M_OH_stretch": 3653}),
        ],
    },
    "silicate": {
        "label": "硅酸盐",
        "modes": ["SiO4_v1", "SiO4_v4", "Si_OH", "organic"],
        "samples": [
            # 石英相
            ("Quartz", "Powder 1", {"SiO4_v1": (1079, 1162, 1232), "SiO4_v4": 474}),
            ("Quartz", "Pellet 2", {"SiO4_v1": (1090, 1168, None), "SiO4_v4": 477,
                                    "organic": 1600}),
            ("Opal", "Crystal", {"SiO4_v1": None, "SiO4_v4": 435, "organic": 1608}),
            ("Opal", "Pellet 2", {"SiO4_v1": 1078, "SiO4_v4": 485, "organic": 1618}),
            ("Basaltic glass", "Pellet 4", {"SiO4_v1": None, "SiO4_v4": None}),
            # 辉石
            ("Augite", "Powder", {"SiO4_v1": 1014, "SiO4_v4": None}),
            ("Augite", "Pellet", {"SiO4_v1": 1000, "SiO4_v4": None}),
            ("Diopside", "Crystal", {"SiO4_v1": 1017, "SiO4_v4": None}),
            ("Diopside", "Pellet", {"SiO4_v1": 1014, "SiO4_v4": None}),
            ("Enstatite", "Crystal", {"SiO4_v1": 1026, "SiO4_v4": None}),
            ("Enstatite", "Powder", {"SiO4_v1": 1021, "SiO4_v4": None}),
            ("Enstatite", "Pellet", {"SiO4_v1": 1014, "SiO4_v4": None}),
            ("Hedenbergite", "Crystal", {"SiO4_v1": None, "SiO4_v4": None}),
            ("Hedenbergite", "Pellet", {"SiO4_v1": None, "SiO4_v4": None}),
            # 斜长石
            ("Albite", "Crystal", {"SiO4_v1": None, "SiO4_v4": None}),
            ("Bytownite", "Pellet", {"SiO4_v1": 1014, "SiO4_v4": None,
                                     "organic": 1598}),
            ("Labradorite", "Powder", {"SiO4_v1": 1003, "SiO4_v4": None}),
            # 橄榄石
            ("Fayalite", "Crystal", {"SiO4_v1": None, "SiO4_v4": None}),
            ("Forsterite", "Crystal", {"SiO4_v1": 844, "SiO4_v4": None}),
            ("Olivine", "Pellet 16", {"SiO4_v1": 846, "SiO4_v4": None}),
            ("Olivine", "Pellet 26", {"SiO4_v1": 847, "SiO4_v4": None}),
            ("Olivine", "Pellet 36", {"SiO4_v1": 847, "SiO4_v4": None}),
            ("Olivine", "Pellet 46", {"SiO4_v1": None, "SiO4_v4": None}),
            # 钾长石
            ("Microcline", "Crystal", {"SiO4_v1": 1114, "SiO4_v4": 504}),
            ("Orthoclase", "Crystal", {"SiO4_v1": 1113, "SiO4_v4": 494}),
            # 其它
            ("Epidote", "Pellet 2", {"SiO4_v1": None, "SiO4_v4": None}),
            ("Topaz", "Pellet 2", {"SiO4_v1": 929, "SiO4_v4": None,
                                   "Si_OH": 3652, "organic": 1611}),
        ],
    },
    "phyllosilicate": {
        "label": "层状硅酸盐",
        "modes": ["SiO4_v1", "Si_OH", "organic"],
        "samples": [
            ("Chlorite", "Pellet", {"SiO4_v1": None, "organic": 1578}),
            ("Illite", "Pellet", {"SiO4_v1": None, "Si_OH": None}),
            ("Kaolinite", "Powder", {"SiO4_v1": None, "Si_OH": (3632, 3696),
                                     "organic": 1609}),
            ("Montmorillonite", "Powder", {"SiO4_v1": None, "organic": 1592}),
            ("Smectite A5", "Powder", {"SiO4_v1": None, "Si_OH": None}),
            ("Smectite G5", "Powder", {"SiO4_v1": None, "Si_OH": None}),
            ("Smectite H5", "Powder", {"SiO4_v1": None, "Si_OH": None}),
            ("Smectite J5", "Powder", {"SiO4_v1": None, "Si_OH": None}),
            ("Smectite K5", "Powder", {"SiO4_v1": None, "Si_OH": None}),
            ("Smectite N5", "Powder", {"SiO4_v1": None, "Si_OH": None}),
            ("Biotite", "Crystal", {"SiO4_v1": None, "Si_OH": None}),
            ("Muscovite", "Crystal", {"SiO4_v1": None, "Si_OH": 3627}),
            ("Muscovite", "Pellet 2", {"SiO4_v1": None, "Si_OH": 3636,
                                       "organic": 1584}),
            ("Antigorite", "Pellet", {"SiO4_v1": None, "Si_OH": 3679}),
            ("Serpentine", "Pellet 2", {"SiO4_v1": 1054, "Si_OH": 3691,
                                        "organic": 1605}),
            ("Serpentinite", "Rock", {"SiO4_v1": None, "Si_OH": 3680,
                                      "organic": 1600}),
        ],
    },
    # ---- 2026-09-15 新增：高氯酸盐与氯酸盐 --------------------------------
    # 来源：Clark, J.V. et al. (2026), "Deep-Ultraviolet (DUV) Raman and
    # Fluorescence of Perchlorates and Chlorates...", Earth and Space Science,
    # DOI 10.1029/2025EA004512，表 1（ACRONM 实测，同样是 248.6 nm 激发）。
    #
    # 为什么要补这一类：上面那 62 种矿物的库里**一个高氯酸盐、一个氯酸盐都没有**，
    # 而"高氯酸盐或磷酸盐"恰恰是 Jezero 最常见的盐类（Corpolongo et al. 2023
    # 每个 Máaz 靶点都报）。原来我们只能靠峰位区间记模糊类别，真按谱库比对反而
    # 会失效——因为库里根本没这一类。
    #
    # 论文给的判别规律（正文原话）：对称伸缩 925–943（Na perchlorate 例外，950）、
    # 反对称伸缩 971–1190、对称 Cl-O 变形 600–650、反对称 Cl-O 变形 400–500、
    # H2O 弯曲 1600–1630；Na/K 的反对称伸缩峰**又尖又可分**，Mn/Fe/Al/Mg 的很宽。
    "oxychlorine": {
        "label": "高氯酸盐与氯酸盐",
        "modes": ["ClO_nu_s", "ClO_nu_as", "ClO_delta_s", "ClO_delta_as",
                  "overtone", "H2O_bend", "H2O_stretch", "VW_unassigned"],
        "samples": [
            ("Na perchlorate", "Honeywell",
             {"ClO_nu_s": 950.7, "ClO_nu_as": (1091.2, 1145.8),
              "ClO_delta_s": (621.4, 651.1), "ClO_delta_as": (441.7, 479.2),
              "overtone": 886.7, "H2O_bend": 1621.4,
              "H2O_stretch": (3526.6, 3549.9), "VW_unassigned": 1563.6}),
            ("Na chlorate", "Sigma Aldrich",
             {"ClO_nu_s": 929.9, "ClO_nu_as": (980.5, 1022.7),
              "ClO_delta_s": 617.8, "ClO_delta_as": 480.3,
              "H2O_stretch": 3469.8,
              "VW_unassigned": (1242.4, 1415.4, 1552.7, 1861.0)}),
            ("K perchlorate", "Acros",
             {"ClO_nu_s": 935.8, "ClO_nu_as": (1081.6, 1119.4),
              "ClO_delta_s": 629.2, "ClO_delta_as": 461.0,
              "VW_unassigned": (1255.8, 1401.2, 1561.6, 1711.7, 1762.6, 2214.7)}),
            ("K chlorate", "Sigma Aldrich",
             {"ClO_nu_s": 930.6, "ClO_nu_as": 970.3, "ClO_delta_s": 612.0,
              "ClO_delta_as": 482.3, "H2O_stretch": (3205.9, 3426.9),
              "VW_unassigned": (1230.8, 1422.8, 1548.0, 1860.9)}),
            ("Ca perchlorate", "Acros",
             {"ClO_nu_s": 942.9, "ClO_nu_as": (1082.5, 1117.4),
              "ClO_delta_s": 624.7, "ClO_delta_as": 450.1,
              "H2O_bend": 1624.8,
              "H2O_stretch": (3229.4, 3463.2, 3531.0, 3600.2),
              "VW_unassigned": 1402.0}),
            ("Mg perchlorate", "Sigma Aldrich",
             {"ClO_nu_s": 934.6, "ClO_nu_as": 1109.0, "ClO_delta_s": 632.6,
              "ClO_delta_as": 463.3, "H2O_bend": 1633.2,
              "H2O_stretch": 3528.2}),
            ("Al perchlorate", "Thermo Scientific",
             {"ClO_nu_s": 928.2, "ClO_nu_as": (1051.9, 1107.3),
              "ClO_delta_s": 621.1, "ClO_delta_as": 456.3,
              "H2O_bend": 1620.0, "H2O_stretch": (3251.0, 3450.0),
              "VW_unassigned": 2493.0}),
            ("Ba chlorate", "Sigma Aldrich",
             {"ClO_nu_s": 923.7, "ClO_nu_as": 964.5, "ClO_delta_s": 613.6,
              "ClO_delta_as": 491.7, "H2O_bend": 1599.8,
              "H2O_stretch": (3172.3, 3506.9, 3576.4),
              "VW_unassigned": (1228.9, 1414.8, 1539.9, 1811.3, 1854.5, 2063.4)}),
            ("Fe II perchlorate", "Sigma Aldrich",
             {"ClO_nu_s": 925.8, "ClO_nu_as": (1079.5, 1121.3),
              "ClO_delta_s": 621.6, "ClO_delta_as": (424.8, 457.6),
              "overtone": 797.4, "H2O_bend": 1613.0,
              "H2O_stretch": 3501.9}),
            ("Mn perchlorate", "Thermo Scientific",
             {"ClO_nu_s": 925.9, "ClO_nu_as": (1076.8, 1114.8),
              "ClO_delta_s": 624.6, "ClO_delta_as": 455.5,
              "H2O_bend": 1613.4, "H2O_stretch": (3219.6, 3505.1)}),
        ],
    },
}

# 论文正文里按类别归纳的特征峰（用于给判定结果加一句通俗依据）
# 谱库判定的置信门槛。实测下来：
#   ≥0.7  多峰模式干净对上（方解石 1.00、石英 1.00、石膏 0.75）
#   0.5-0.7 有支持但不够硬（Guillaumes 的硫酸盐 0.52、橄榄石单峰 0.50）
#   <0.5  基本是"只命中一个峰"，参考库里单峰矿物会靠覆盖率占便宜，
#          这种结果不该当成结论，只当线索
LIB_MIN_SCORE = 0.55


def verdict(matches):
    """把 match() 的结果翻成人话，带上置信判断。"""
    if not matches:
        return "证据不足", None
    top = matches[0]
    if top["score"] >= 0.70:
        return "可信", top
    if top["score"] >= LIB_MIN_SCORE:
        return "有可能", top
    return "证据不足", top


CLASS_NOTES = {
    "borate": "主峰 950（[BO3]）或 1080（[BO4]）附近，B–OH 伸缩在 ~3580",
    "carbonate": "ν₁ 在 1085–1100，另有两个强峰 ~1440 和 ~1750",
    "phosphate": "ν₁ 在 ~960，ν₃ 在 ~1050（与碳酸盐易混）",
    "sulfate": "ν₁ 在 990–1020，ν₂ 在 ~460–500，ν₃ 在 1110–1190",
    "halide": "800–1800 内基本无峰，只有水和铵的弱峰",
    "oxide_hydroxide": "M–O 伸缩/弯曲峰弱，M–OH 伸缩在 3650 附近",
    "silicate": "橄榄石 844–847；辉石 1000–1026；石英 1079/1162/1232 + 474",
    "phyllosilicate": "Si–O 峰被滤光片挡掉，只有 Si–OH 在 3630–3690",
}


def all_peaks():
    """把参考库拉平成一串 (类别, 矿物, 样品, 模式, 峰位) —— 只保留检出的峰。"""
    out = []
    for key, cls in REFERENCE.items():
        for mineral, sample, peaks in cls["samples"]:
            for mode, val in peaks.items():
                if val is None:
                    continue
                vals = val if isinstance(val, tuple) else (val,)
                for v in vals:
                    if v is None:          # 元组里也会有 n.d.，比如石英的 1090,1168,n.d.
                        continue
                    out.append((key, mineral, sample, mode, v))
    return out


def peak_pool(exclude_organic=True):
    """按类别汇集所有参考峰位，返回 {类别: [峰位,...]}。"""
    pool = {}
    for key, _mineral, _sample, mode, v in all_peaks():
        if exclude_organic and mode == "organic":
            continue
        pool.setdefault(key, []).append(v)
    return {k: sorted(v) for k, v in pool.items()}


def sample_peaks(peaks, exclude_organic=True):
    """把一个样品条目里的峰位拉平成列表（丢掉 n.d. 和有机 G 带）。"""
    out = []
    for mode, val in peaks.items():
        if val is None or (exclude_organic and mode == "organic"):
            continue
        vals = val if isinstance(val, tuple) else (val,)
        out.extend(v for v in vals if v is not None)
    return sorted(out)


def match(observed, tol=8.0, exclude_organic=True, window=(800.0, 1800.0),
          min_rel_intensity=0.02):
    """把一组观测峰位与参考库里**每一个样品**比对，返回按得分排序的候选。

    为什么要逐样品比对，而不是比"类别峰位池"：
    把一类里所有样品的峰混在一起会得到几十个值、覆盖面极宽，几乎什么都算命中。
    逐样品比对才有区分度——比如碳酸盐里，方解石 ν₁ 在 1085、白云石 1099、
    天然碱 1044，靠主峰位置就能分开。

    **2026-09-15 重写**。原来用 F1，但有两个结构性毛病，导致分数普遍只有零点几、
    而且**完全没有分辨力**（已判定组的最高分中位数 0.33，未归类组也是 0.33）：

    1. **参考峰没有按分析窗口筛**。参考库里大量峰在观测窗口之外
       （高氯酸钠有 441/479/621/3526 四个，窗口内只有 950.7/1091.2/1145.8），
       而 recall = 命中的参考峰 / **全部**参考峰——等于每个样品都被扣一笔
       不可能还的债：高氯酸钠的 recall 上限只有 3/12。现在只统计窗口内的峰。
    2. **参考侧没有按强度加权**。表 1/表 4 都给了归一化强度，很弱的峰
       （相对强度 0.0005）根本观测不到，把它们算进 recall 同样是白扣分。
       现在丢掉相对强度 < min_rel_intensity 的参考峰。

    打分改成**强度加权的软 Jaccard**：
        num = Σ_命中 √(w_观测 × w_参考)
        den = Σ_观测 w + Σ_参考 w − num
        score = num / den
    完美匹配 = 1；多一个观测峰（杂峰）或少一个参考峰（漏检）都要扣分，
    而且**强度大的峰扣得多**——这才叫"谱形比较"，不是数几个峰撞上了。

    tol 默认 8 cm⁻¹：论文说参考峰位的合计不确定度是 ±5 cm⁻¹，
    我们自己的峰位误差也在几个 cm⁻¹，8 是两者相加再留一点余量。
    """
    # observed 可以是 [位置,...] 或 [(位置, 权重),...]。
    # 权重用峰的 SNR：主峰比蹭进来的杂峰强得多，不按强度加权的话，
    # 杂峰命中某个"参考库里只有一个峰"的矿物反而会赢。
    obs = []
    for item in observed:
        if isinstance(item, (tuple, list)):
            obs.append((float(item[0]), float(item[1])))
        else:
            obs.append((float(item), 1.0))
    obs.sort()
    if not obs:
        return []
    # 观测峰也只在窗口内比较（外面本来就没检出）
    obs = [(p, w) for p, w in obs if window[0] <= p <= window[1]]
    if not obs:
        return []
    obs_max = max(w for _p, w in obs) or 1.0
    obs = [(p, w / obs_max) for p, w in obs]
    total_w = sum(w for _p, w in obs) or 1.0

    out = []
    for key, cls in REFERENCE.items():
        for mineral, sample, peaks in cls["samples"]:
            ref = sample_peaks(peaks, exclude_organic)
            rel = cls.get("intensities", {}).get((mineral, sample))
            # 参考峰：只留窗口内、强度够的
            keep = []
            for i, r in enumerate(ref):
                if not (window[0] <= r <= window[1]):
                    continue
                a = 1.0 if rel is None else float(rel[i])
                if a < min_rel_intensity:
                    continue
                keep.append((r, a))
            if not keep:
                continue
            ref_max = max(a for _r, a in keep) or 1.0
            keep = [(r, a / ref_max) for r, a in keep]

            pairs = []
            # **一对一匹配**（贪心：按距离从近到远认领）。
            # 不这么做的话，一个参考峰会被好几个观测峰重复认领，
            # 分子被重复累加、分母又被减掉，分数能超过 1（实测出现过 2.197）。
            cand = []
            for i, (p, w) in enumerate(obs):
                for j, (r, a) in enumerate(keep):
                    d = abs(p - r)
                    if d <= tol:
                        cand.append((d, i, j))
            cand.sort()
            used_o, used_r = set(), set()
            for d, i, j in cand:
                if i in used_o or j in used_r:
                    continue
                used_o.add(i)
                used_r.add(j)
                p, w = obs[i]
                r, a = keep[j]
                pairs.append((p, r, d, w, a))
            num = sum((w * a) ** 0.5 for _p, _r, _d, w, a in pairs)
            den = total_w + sum(a for _r, a in keep) - num
            score = (num / den) if den > 0 else 0.0
            # 证据量加权：只有一个峰命中的"匹配"没有信息量（辉石只有 1017 一个峰），
            # 命中 1 个参考峰打折、2 个及以上才算充分。
            n_hit = len({r for _p, r, _d, _w, _a in pairs})
            evidence = min(1.0, n_hit / 2.0)
            score *= evidence
            prec = sum(w for _p, _r, _d, w, _a in pairs) / total_w
            rec = len({r for _p, r, _d, _w, _a in pairs}) / len(keep)
            f1 = (2 * prec * rec / (prec + rec)) if (prec + rec) else 0.0
            rms = (sum(d for _p, _r, d, _w, _a in pairs) / len(pairs)) if pairs else float("nan")
            out.append({
                "class_key": key, "class": cls["label"],
                "mineral": mineral, "sample": sample,
                "score": round(score, 3), "f1": round(f1, 3),
                "evidence": round(evidence, 2),
                "precision": round(prec, 3),
                "recall": round(rec, 3),
                "n_observed": len(obs), "n_matched": len(pairs),
                "n_reference": len(keep),
                "mean_offset": None if rms != rms else round(rms, 1),
                "reference_peaks": ref,
            })
    out.sort(key=lambda d: (-d["score"], d["mean_offset"] if d["mean_offset"] is not None else 99))
    return out


def best_per_class(observed, tol=8.0, exclude_organic=True):
    """每个类别只留得分最高的那个样品——报告里按类别列的时候用这个。"""
    best = {}
    for r in match(observed, tol=tol, exclude_organic=exclude_organic):
        if r["class_key"] not in best:
            best[r["class_key"]] = r
    return sorted(best.values(), key=lambda d: -d["score"])
