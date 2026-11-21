# -*- coding: utf-8 -*-
"""把散落的文件整理成扫描。

分组规则全部由 PDS4 label 判定，不靠猜文件名：

* 活动(activity) = 同一个 sol + 同一个 sequence_id
* 活动内部，同类型、同处理、expID 依次递增且点数相同的多个光谱产品，
  属于同一次扫描被拆成的若干段（例如 expID 1 与 expID 2 各 50 条，
  合起来才是完整的 100 点扫描）
* 其余同类型产品各自成为独立扫描

坐标表 RLS、激光能量 RLI 等辅助产品按 SCLK 或 expID 归入对应的段。
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from .reader import Product

SPECTRAL_TYPES = ("RRS", "RRA", "RRB", "RCS", "RCA", "RCB",
                  "ERA", "ERB", "ERP", "ECA", "ECB")
AUX_TYPES = ("RLS", "RLP", "RLI", "RMO", "RCC")


@dataclass
class Scan:
    """一次扫描：由一到多个光谱产品段组成，外加它的辅助文件。"""

    sol: str
    sequence_id: str
    sclk: str
    parts: list = field(default_factory=list)
    products: list = field(default_factory=list)
    # 同一次观测的原始版本备选组（RAW 版与 processed 版合并时保留，
    # view 层据此列出「ERA 原始」条目——用户可从暗背景扣除开始自行预处理）
    alts: list = field(default_factory=list)

    @property
    def key(self) -> str:
        return "%s|%s|%s" % (self.sol, self.sequence_id, self.sclk)

    @property
    def expids(self):
        return [str(p.pid.flags.get("expid", "")) for p in self.parts]

    @property
    def all_sclks(self):
        return [p.pid.sclk for p in self.parts]

    @property
    def label_time(self) -> str:
        for p in list(self.parts) + list(self.products):
            if p.label.start_time:
                return p.label.start_time
        return ""

    def by_type(self, *types) -> list:
        return [p for p in self.products if p.product_type in types]

    def pick(self, *types) -> Optional[Product]:
        """取某类型里处理最完整的一个。"""
        got = self.by_type(*types)
        if not got:
            return None

        def rank(p):
            f = p.pid.flags if not p.pid.is_image else {}
            score = sum(1 for k in ("wavelength", "cosmic", "gain", "normalize")
                        if f.get(k) in (True, "C", "G", "N"))
            return (-score, -p.size_mb)

        return sorted(got, key=rank)[0]

    @property
    def spectra(self) -> Optional[Product]:
        return self.parts[0] if self.parts else None

    @property
    def is_split(self) -> bool:
        return len(self.parts) > 1

    @property
    def n_points(self) -> Optional[int]:
        total = sum(p.n_spectra or 0 for p in self.parts)
        return total or None

    def read_regions(self):
        """读 Region 1/2/3；被拆分的扫描按段顺序拼接。"""
        import numpy as np
        wl = None
        chunks = [[], [], []]
        for p in self.parts:
            if wl is None:
                wl = p.wavelength_nm()
            rs = p.regions()
            for i in range(3):
                if rs[i] is not None:
                    chunks[i].append(rs[i])
        out = []
        for c in chunks:
            out.append(np.vstack(c) if c else None)
        return wl, out[0], out[1], out[2]

    def read_shot_positions(self):
        """拼接各段的 RLS 坐标；第二段的点号整体后移。

        RLS 的行数 = 点数 × 该次扫描用到的 ACI 图像张数，所以**按点号
        去重**取坐标，不能按行数切（`pts[:n]` 隐含"行序图像优先"的
        假设，一旦顺序不是这样就会丢点）。同一点在多张图上各有一行，
        保留第一次出现的那张图的坐标。
        """
        allpts = []
        offset = 0
        for p in self.parts:
            rls = None
            for q in self.products:
                if q.product_type in ("RLS", "RLP") and q.pid.sclk == p.pid.sclk:
                    rls = q
                    break
            if rls is None:
                rls = self.pick("RLS", "RLP")
            pts = rls.shot_positions() if rls else None
            n = p.n_spectra or (len(pts) if pts else 0)
            if pts:
                seen = set()
                for x in pts:
                    num = int(x["number"])
                    if num < 0 or num >= n or num in seen:
                        continue
                    seen.add(num)
                    y = dict(x)
                    y["number"] = num + offset
                    allpts.append(y)
            offset += n
        return allpts or None


def _point_range(p: Product, items):
    """找同 SCLK 的 RLS，返回 (起始点号, 结束点号, 点数)；没有则 None。

    注意：RLS 的行数 = 点数 x 该次扫描用到的 ACI 图像张数，
    所以点数要按**不同的点号**去数，不能数行数。
    """
    rls = None
    for q in items:
        if q.product_type in ("RLS", "RLP") and q.pid.sclk == p.pid.sclk:
            rls = q
            break
    if rls is None:
        return None
    try:
        pts = rls.shot_positions()
    except Exception:                                # noqa: BLE001
        return None
    if not pts:
        return None
    nums = sorted({x["number"] for x in pts})
    return (nums[0], nums[-1], len(nums))


def _is_half(a, b):
    """两段是否同属一次被拆分的扫描。

    判据：两段都从第 0 点开始、点数相同，且点数**未达完整扫描规模**。
    HDR / Detail 的完整扫描是 100 点；只有不足 100 点、从 0 开始的坐标表
    才是「半段」（对应指南里的 108 / 208 现象）。
    """
    if not a or not b:
        return False
    if a[0] != 0 or b[0] != 0 or a[2] != b[2]:
        return False
    if a[1] != a[2] - 1 or b[1] != b[2] - 1:
        return False
    return a[2] < FULL_SCAN_MIN_POINTS


FULL_SCAN_MIN_POINTS = 100


def _norm_sol(sol: str) -> str:
    return str(int(sol)) if sol and sol.isdigit() else (sol or "")


def _has_wavelength(p: Product) -> bool:
    try:
        return p.find_table(r"WAVELENGTH") is not None
    except Exception:                                # noqa: BLE001
        return False


def _merge_same_scan(scans: list) -> list:
    """把"同一次观测的不同版本"合并成一次扫描。

    两条判据，满足其一就合并：

    1. **SCLK 相同** —— 同一条命令产出的两份产品。
    2. **RLS 坐标逐点相同，且处理级别不同** —— 同一次扫描的 RAW 版与
       processed 版。PDS 会把同一次观测按不同处理级别各发一份，文件名
       后缀不同（`w108cgn` = 已做波长校正/去宇宙射线/增益/归一化；
       `b108zpz` = 全都没做），但它们的 SCLK 不一样，只靠第 1 条合不到一起。

       实测踩过：sol 186 的两次 HDR_250 扫描各发了 3 份产品，被当成 6 次
       独立扫描；论文里 Bellegarde 只有 2 次 HDR 扫描（共 200 条谱）。
       加这条判据后归一成 2 次。

       为什么要求"处理级别不同"：sol 207（Garde）里 expID 1 出现了两次
       （SCLK 0685346165 和 0685347508），处理级别**完全相同**，那是两次
       真实的不同观测，坐标也不同，不能被并掉。

    合并后优先保留带波长表的那一组作为主光谱，其余只作为文件列出。
    """
    used = [False] * len(scans)
    sig_cache = {}
    out = []
    for i, a in enumerate(scans):
        if used[i]:
            continue
        used[i] = True
        group = [a]
        sids = set(a.all_sclks)
        sig = _rls_signature(a, sig_cache)
        for j in range(i + 1, len(scans)):
            b = scans[j]
            if used[j] or b.sol != a.sol or b.sequence_id != a.sequence_id:
                continue
            same_sclk = bool(sids & set(b.all_sclks))
            same_points = False
            if not same_sclk and sig is not None:
                sig_b = _rls_signature(b, sig_cache)
                same_points = (sig_b == sig and _proc_flags(b) != _proc_flags(a))
            if same_sclk or same_points:
                group.append(b)
                used[j] = True
                sids |= set(b.all_sclks)
        if len(group) == 1:
            out.append(a)
            continue
        best = next((g for g in group if all(_has_wavelength(p) for p in g.parts)), group[0])
        merged = Scan(sol=a.sol, sequence_id=a.sequence_id, sclk=best.sclk,
                      parts=list(best.parts))
        seen = set()
        for g in group:
            for p in g.products:
                if p.data_path not in seen:
                    seen.add(p.data_path)
                    merged.products.append(p)
        # 原始版本组挂为 alts（不并入返回列表——phase10/data_audit 的计数
        # 契约不受影响；view 层据此追加「ERA 原始」条目）
        for g in group:
            if g is not best:
                merged.alts.append(g)
        out.append(merged)
    return out


def _rls_signature(scan, cache):
    """这次扫描的坐标指纹：RLS 里每个点的 (samp, line)。

    坐标逐点相同 = 打在石块上的同一批位置 = 同一次扫描。
    这是对**数据本身**的核对，不是靠文件名猜。
    """
    if scan.key in cache:
        return cache[scan.key]
    sig = None
    rls = scan.pick("RLS", "RLP")
    if rls is not None:
        try:
            pts = rls.shot_positions()
        except Exception:                            # noqa: BLE001
            pts = None
        if pts:
            sig = tuple((round(float(p["samp"]), 3), round(float(p["line"]), 3))
                        for p in pts)
    cache[scan.key] = sig
    return sig


def _proc_flags(scan):
    """这次扫描用到的处理级别（波长校正 / 去宇宙射线 / 增益 / 归一化）。"""
    out = set()
    for p in scan.parts:
        f = p.pid.flags
        out.add((f.get("wavelength"), f.get("cosmic"), f.get("gain"),
                 f.get("normalize")))
    return frozenset(out)


def build_scans(products) -> list:
    """把产品列表整理成扫描列表。"""
    activities = {}
    for p in products:
        sol = _norm_sol(p.label.sol or p.pid.sol)
        seq = (p.label.sequence_id or p.pid.sequence_id).lower()
        if not sol or not seq:
            continue
        activities.setdefault((sol, seq), []).append(p)

    scans = []
    for (sol, seq), items in sorted(activities.items()):
        spectral = [p for p in items if p.product_type in SPECTRAL_TYPES]
        others = [p for p in items if p.product_type not in SPECTRAL_TYPES]

        ranges = {p.name: _point_range(p, items) for p in spectral}

        by_type = {}
        for p in spectral:
            by_type.setdefault(p.product_type, []).append(p)

        made = []
        for _, group in by_type.items():
            group.sort(key=lambda p: p.pid.sclk)
            i = 0
            while i < len(group):
                run = [group[i]]
                j = i + 1
                while j < len(group):
                    prev, nxt = group[j - 1], group[j]
                    if (_is_half(ranges.get(prev.name), ranges.get(nxt.name))
                            and (prev.n_spectra or -1) == (nxt.n_spectra or -2)):
                        run.append(nxt)
                        j += 1
                    else:
                        break
                made.append(run)
                i = j

        for parts in made:
            sc = Scan(sol=sol, sequence_id=seq, sclk=parts[0].pid.sclk, parts=list(parts))
            sids = {p.pid.sclk for p in parts}
            eids = {str(p.pid.flags.get("expid", "")) for p in parts}
            for q in others:
                if q.pid.sclk in sids or str(q.pid.flags.get("expid", "")) in eids:
                    sc.products.append(q)
            sc.products.extend(parts)
            # 原始主动帧（ERA/ECA）：把配对的暗帧（ERB/ECB，同 SCLK）折进
            # products——preprocess._match 在 products 里找 ERB 触发暗背景
            # 扣除步；再把同序列 RRS 也折进去，波长轴借用才有着落。
            if parts[0].product_type in ("ERA", "ECA"):
                for q in spectral:
                    if q.product_type in ("ERB", "ECB") and q.pid.sclk in sids:
                        sc.products.append(q)
                    elif q.product_type in ("RRS", "RRA"):
                        sc.products.append(q)
            scans.append(sc)
    return _merge_same_scan(scans)


def scan_folder(folder: str):
    """扫描文件夹，返回 (产品列表, 扫描列表)。"""
    from .reader import collect_products
    products = collect_products(folder)
    return products, build_scans(products)


def summary_rows(scans) -> list:
    rows = []
    for ob in scans:
        s = ob.spectra
        rows.append({
            "sol": ob.sol,
            "sequence_id": ob.sequence_id,
            "sclk": ob.sclk,
            "time": ob.label_time,
            "n_files": len(ob.products),
            "n_points": ob.n_points,
            "spectra": s.product_type if s else "-",
            "spectra_file": s.name if s else "-",
            "has_rls": bool(ob.by_type("RLS", "RLP")),
            "split": ("%d 段" % len(ob.parts)) if ob.is_split else "单段",
            "expids": "".join(ob.expids),
            "types": " ".join(sorted({p.product_type for p in ob.products})),
        })
    return rows
