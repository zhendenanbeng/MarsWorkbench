# -*- coding: utf-8 -*-
"""仪器后端的公共部分：数据模型、LIBS 通用预处理、下载器。

每个「操作台」就是一个 Backend：
  - 描述自己（名字、坐标轴、预处理开关、PDS 数据源、产品类型表）
  - scan_folder(root)：把本地文件夹整理成「观测 → 点 → 文件」
  - load_point(...)：读一个点的原始数据（波长轴 + 各发激光的强度矩阵）
  - annotate(...)：对显示中的光谱做元素线判定

观测/点的口径与 SHERLOC 工作台一致：观测 = 同 sol 同序列号的一次活动，
点 = 观测里的一个打点位置（ChemCam/SuperCam 一个点一个文件，
文件里再分激光发数 shot；APXS 一次观测就是一个单点谱）。
"""

from __future__ import annotations

import os
import re

import numpy as np

from . import pds_listing

# --------------------------------------------------------------------------
# 数据模型
# --------------------------------------------------------------------------

class PointFile:
    """一个打点位置对应的本地文件（或一组同义产品文件）。"""

    def __init__(self, stem, files, sclk="", point_no=0, meta=None):
        self.stem = stem              # 不含扩展名的主名
        self.files = files            # {产品类型: 绝对路径}
        self.sclk = sclk
        self.point_no = point_no
        self.meta = meta or {}


class Observation:
    def __init__(self, sol, seq, product):
        self.sol = str(sol)
        self.seq = str(seq)
        self.product = product        # 显示用的主产品类型（如 CL5-RDR）
        self.points = []

    @property
    def key(self):
        return "%s|%s" % (self.sol, self.seq)

    @property
    def n_points(self):
        return len(self.points)


# --------------------------------------------------------------------------
# LIBS 通用预处理
# --------------------------------------------------------------------------

def drop_dust(shots, n=5, log=None):
    """丢弃前 n 发（除尘）：打在表面尘土上的头几发谱会完全不同。"""
    if shots.shape[1] <= n + 3:       # 剩得太少就不动了
        return shots
    if log is not None:
        log.append(("除尘", "丢弃前 %d 发（表面尘土），剩 %d 发"
                    % (n, shots.shape[1] - n)))
    return shots[:, n:]


def despike_shots(shots, log=None):
    """跨发去尖峰：每个通道对各发取中位数/MAD，超 8σ 的单发值替换为中位数。

    宇宙射线 / 热像素表现为「某一个发数、某一个通道」的孤立尖峰，
    沿发数方向做稳健统计正好抓住它们。
    """
    if shots.shape[1] < 6:
        return shots
    med = np.median(shots, axis=1, keepdims=True)
    mad = np.median(np.abs(shots - med), axis=1, keepdims=True)
    sig = np.maximum(1.4826 * mad, 1e-9)
    bad = np.abs(shots - med) > 8.0 * sig
    n_bad = int(bad.sum())
    if n_bad:
        shots = shots.copy()
        shots[bad] = np.broadcast_to(med, shots.shape)[bad]
    if log is not None:
        log.append(("去尖峰", "跨发 MAD 判定，替换 %d 个孤立尖峰（共 %d 值）"
                    % (n_bad, shots.size)))
    return shots


def poly_baseline(x, y, degree=3, rounds=5, log=None):
    """迭代多项式连续背景扣除（对总谱做，不管发数）。

    用 SHERLOC 侧同款的**对称稳健裁剪**（sherloc/background.symmetric_baseline）：
    以残差中位数为中心双侧裁剪、每轮用全残差重新估计尺度。原先「只压上包络」
    的单侧裁剪会让基线整体偏高，扣出来的谱带系统性负偏置——发射谱线是正峰，
    对称裁剪同样能保护峰不受拟合侵蚀，还不引入偏置。

    输入不满足严格前提（波长非递增、含 NaN 等）时跳过扣除并在日志说明，
    绝不输出一条被 NaN 洗掉的假谱。
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    # 对称算法要求严格递增：排序 + 去掉重复波长（分光计接缝处可能出现），
    # 算完再把基线插值回原始通道位置
    order = np.argsort(x, kind="stable")
    xs, ys = x[order], y[order]
    uniq = np.concatenate(([True], np.diff(xs) > 0))
    dropped = int((~uniq).sum())
    if dropped:
        xs, ys = xs[uniq], ys[uniq]
    try:
        from sherloc.background import symmetric_baseline
        base = symmetric_baseline(xs, ys, degree=degree, n_iter=rounds, clip=1.5)
    except Exception:                                   # noqa: BLE001
        if log is not None:
            log.append(("基线", "输入不满足拟合前提（非递增/含无效值），跳过扣除"))
        return y
    if log is not None:
        log.append(("基线", "对称稳健 %d 次多项式（迭代 %d 轮%s）"
                    % (degree, rounds,
                       "，去重 %d 个重复波长" % dropped if dropped else "")))
    return y - np.interp(x, xs, base)


def smooth(y, k=7):
    """滑动平均（APXS 可选平滑）。"""
    if len(y) < k:
        return y
    kernel = np.ones(k) / k
    pad = np.concatenate([np.full(k // 2, y[0]), y, np.full(k // 2, y[-1])])
    return np.convolve(pad, kernel, mode="valid")


# --------------------------------------------------------------------------
# Backend 基类
# --------------------------------------------------------------------------

class Backend:
    id = ""
    name = ""
    rover = ""
    rover_label = ""
    subtitle = ""
    xlabel = "波长 (nm)"
    ylabel = "强度"
    accent = "#1a63c8"               # 界面主色
    steps = []                       # [(key, 中文名, 默认开)]
    pds_roots = {}                   # {kind: {"label": str, "url": str}}
    sol_pattern = re.compile(r"sol_?(\d+)$", re.I)
    walk_skip = {"_qc", "_png", "_maps", "_aci", ".git", ".cache", "_notes", "__pycache__", "_libraries"}

    # -- 子类要实现的接口 ---------------------------------------------------

    def scan_folder(self, root):
        raise NotImplementedError

    def load_point(self, point, shot, steps_on):
        """返回 {"x": …, "shots": 矩阵(n通道×n发), "log": [(步骤, 文字)],
        "n_shots": int, "meta": {...}}"""
        raise NotImplementedError

    def annotate(self, x, y):
        from . import lines
        return lines.match_lines(x, y, lines.LIBS_LINES, tol=0.35, min_lines=2, wavelength_frame='vacuum_nm')

    def annotate_data(self,data):
        if data.get('passive'):return [],[]
        marks,elements=self.annotate(np.asarray(data['x']),np.asarray(data['y']))
        if self.id in ('chemcam','supercam') and 'shot_matrix' in data:
            from .evidence import enrich
            marks=enrich(data['x'],data['shot_matrix'],marks)
        return marks,elements

    def parse_listing(self, entries):
        """PDS sol 目录文件列表 → 下载面板的分组结构。"""
        raise NotImplementedError

    def point_label(self, obs, point, shot, n_shots):
        return "第 %d 点 · 第 %s 发" % (point.point_no, shot)

    # -- 公共实现 -----------------------------------------------------------

    def walk_files(self, root, pattern):
        """递归找匹配文件，返回 [(路径, 文件名, 上层sol)]。跳过输出目录。"""
        out = []
        rx = re.compile(pattern, re.I)
        for dirpath, dirnames, filenames in os.walk(root):
            dirnames[:] = [d for d in dirnames
                           if d not in self.walk_skip and not d.startswith("下载_")]
            for f in filenames:
                if rx.match(f):
                    out.append((os.path.join(dirpath, f), f, os.path.basename(dirpath)))
        return out

    def sol_from_dirname(self, dirname):
        m = self.sol_pattern.match(dirname or "")
        return str(int(m.group(1))) if m else None


# --------------------------------------------------------------------------
# 下载（PDS → 本地），规则与 SHERLOC 工作台一致
# --------------------------------------------------------------------------

def download_files(entries, dest, job):
    """逐个下载 entries [{name,size,url}] 到 dest。

    job 是共享状态 dict：cancel 标志、done/total/bytes_*、log 列表。
    已存在且大小一致跳过；.part 先写后改名；每个文件失败重试 3 次。
    """
    import time
    import urllib.request

    os.makedirs(dest, exist_ok=True)
    job.update(total=len(entries), done=0,
               bytes_total=sum(e["size"] for e in entries), bytes_done=0,
               out=dest, log=[])
    for e in entries:
        if job.get("cancel"):
            job["log"].append("已取消")
            break
        from .download_queue import target_path,check_length
        path = target_path(dest,e['name'])
        if os.path.exists(path) and os.path.getsize(path) == e["size"]:
            job["done"] += 1
            job["bytes_done"] += e["size"]
            continue
        tmp = path + ".part"
        err = None
        for attempt in range(3):
            if job.get("cancel"):
                break
            try:
                req = urllib.request.Request(
                    e["url"], headers={"User-Agent": pds_listing.USER_AGENT})
                got = 0
                with urllib.request.urlopen(req, timeout=180) as r, open(tmp, "wb") as fh:
                    while True:
                        if job.get("cancel"):
                            break
                        chunk = r.read(262144)
                        if not chunk:
                            break
                        fh.write(chunk)
                        got += len(chunk)
                if job.get("cancel"):
                    break
                check_length(r,got)
                job["bytes_done"] += got
                os.replace(tmp, path)
                job["done"] += 1
                err = None
                break
            except Exception as exc:                # noqa: BLE001
                err = exc
                job["log"].append("第 %d 次失败，重试 %s：%s"
                                  % (attempt + 1, e["name"], str(exc)[:50]))
                time.sleep(2.0 * (attempt + 1))
        if job.get("cancel"):
            if os.path.exists(tmp):
                try:
                    os.remove(tmp)
                except OSError:
                    pass
            break
        if err is not None:
            job["log"].append("放弃 %s：%s" % (e["name"], str(err)[:70]))
    if not job.get('cancel') and job.get('done',0)<job.get('total',0):job['error']='未完成%d个文件'%(job['total']-job['done'])
    if not job.get("cancel") and not job.get("error"):
        job["log"].append("完成：%d 个文件，%.1f MB"
                          % (job["done"], job.get("bytes_done", 0) / 1048576))
