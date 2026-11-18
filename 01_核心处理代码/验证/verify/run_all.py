# -*- coding: utf-8 -*-
"""一键跑全部护栏。任何一项挂了就报红、退出码非零。

**为什么要这个东西**：护栏原来有两条，改完代码得分别手动敲，很容易漏一条。
漏掉的代价是实打实的——"碳酸盐 48 打到 3"（2026-09-13）和
"processed 数据被二次归一"（2026-09-16）这两件事，都是靠人记得跑护栏才发现的。
这个脚本的作用就是把"我记得跑"变成"机制会跑"。

两条护栏各管什么：

* `regression.py`（下限 72；processed 练习树恢复后实跑 328 项）——**判类结果不许跑偏**：
  对照组必须零类别、金标准点数与峰位、单元级结论（Séítah 橄榄石 vs Máaz 盐类）、
  扫描分组结构、已知假阳性上限、扫描级兜底、有机质 G 带通道。
* `duplicate_steps.py`（下限 30；processed 练习树恢复后实跑 235 项）——**processed 数据不许被二次处理**：
  已归一的产品不许再归一、不许再清尖峰；没归一的产品必须照做（反向断言，
  防止"别重复"修成"该做也不做"）。

用法：

    python verify\\run_all.py          # 跑全部，每个脚本只打印末尾几行
    python verify\\run_all.py -v       # 打印完整输出
    python verify\\run_all.py --auc    # 额外跑谱库比对的分辨力量（慢；当前
                                      # AUC 0.60 不达标，只报告、不当失败）
"""

from __future__ import annotations

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout

import os
import subprocess
import sys
import json
import hashlib
import tempfile
from pathlib import Path
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = str(_layout.project_layout.APP)
sys.path.insert(0, ROOT)
import appenv
LAST_RUN = str(appenv.location('reports_selfcheck','last_run.json'))


def source_hashes():
    return _layout.project_layout.source_hashes()


def write_snapshot(results, started, auc_result=None, tested_hashes=None, scope='developer'):
    for row in results:
        row['ok'] = row['code'] == 0 and row['n'] is not None and row['n'] >= row['expect']
    data = {'version': 1, 'scope':scope, 'started_utc': started,
            'finished_utc': datetime.now(timezone.utc).isoformat(), 'results': results,
            'suite_count': len(results), 'assertions': sum(r['n'] or 0 for r in results),
            'failures': sum(not r['ok'] for r in results), 'source_sha256': tested_hashes if tested_hashes is not None else source_hashes(),
            'auc_informational': auc_result}
    fd, temp = tempfile.mkstemp(dir=os.path.dirname(LAST_RUN), suffix='.part')
    try:
        with os.fdopen(fd, 'w', encoding='utf-8') as f:
            json.dump(data, f, ensure_ascii=False, indent=2)
        os.replace(temp, LAST_RUN)
    finally:
        if os.path.exists(temp):os.unlink(temp)
    return data

# (脚本, 一句话说明, 预期断言条数)
CHECKS = [
    ("regression.py", "判类结果护栏（processed 练习树已从 PDS 16.0 恢复；328 项含 7 项数据版本 WARN）", 328),
    ("duplicate_steps.py", "不重复预处理护栏（processed 不许再处理；下限绑定当前练习树）", 30),
    ("integration_guards.py", "集成路径护栏（地图冒烟需上色缓存；下限绑定当前树）", 40),
    ("ui_workflows_check.py", "启动/下载/发号/摘要护栏", 35),
    ("exploration_check.py", "连续特征/光谱分型/任务接口", 52),
    ("instruments_check.py", "ChemCam/SuperCam/APXS 后端与官方成分表", 6),
    ("multi_instrument_check.py", "混合产品报告/导出/缓存/会话", 43),
    ("data_integrity_check.py", "FITS物理值缩放/APXS逐谱成分配对", 13),
    ("scientific_evidence_check.py", "真空线表/独立谱线/重叠/分析失败语义", 45),
    ("fit_profile_check.py", "拟合参数/RMS一致性、非负约束、峰心边界与复核产物", 32),
    ("repeatability_check.py", "候选点统计口径/全窗噪声重放/共同结构/复核产物", 54),
    ("chemcam_passive_check.py", "ChemCam PSV原件/通道对应/失败语义与RMI闸门", 40),
    ("phase4_check.py", "第四期逐发证据与工作流契约", 90),
    ("gold_ledger.py --check", "金标准差异台账与GOLD对账", 2),
    ("phase5_r3_gaps.py", "KNOWN_GAP 四条缺口复核断言", 4),
    ("phase5_cross_check.py", "跨仪器对照（ChemCam×APXS）", 6),
    ("phase6_rmi_check.py", "第六/七期 RMI 共帧投影与元素分布图（投影/产品/高亮/坐标导出）", 25),
    ("phase7_pixl_check.py", "第七/八期 SuperCam×PIXL 对照（RMC/映射/负载/降级/报告页）", 8),
    ("phase7_r2bg_experiment.py --check", "R2 续作：背景竞争判别失败事实锁定", 3),
    ("phase9_spatial_experiment.py --check", "第九期 R2 空间维度判别失败事实锁定", 3),
    ("phase10_check.py", "第十期图像浏览、报告闭环与谱库止损事实（library 指纹段已绑 2026-11-19 重算基线，止损结论跨树复现）", 16),
    ("phase11_check.py", "第十一期产品状态与发布契约", 23),
    ("pixl_browser_check.py", "PIXL独立台、官方坐标与误差、缺值、导出与隔离", 26),
    ("pixl_download_check.py", "PIXL清单、官方直链、下载校验、复用/取消与失败隔离", 20),
    ("workflow_alignment_check.py", "工作台操作统一：真实PMC导航、独立导出、全目录报告、PIQUANT对账入口与共享偏好", 17),
    ("instbkg_check.py", "仪器本底扣除：合成谱M恢复、峰位保持与缺文件/坏JSON安全侧、部分覆盖越界契约", 20),
    ("published_targets_check.py", "已发表靶点台账：结构、来源DOI、定位与只增不缩", 12),
    ("gold_intervals_check.py", "区间数值基准与标定台账：四仪器自主定量对账契约", 15),
    ("libs_quant_check.py", "LIBS自研PLS定量：合成谱真值恢复、CV结构、auto选择、模型往返与安全降级", 21),
    ("pixl_sfs_check.py", "PIXL SFS谱层：读取、能量自证、峰面积恢复与安全降级", 10),
    ("piquant_check.py", "PIQUANT官方开源定量桥：FM配置/标定齐备、stdout解析、MSA导出与官方谱基线复现（PIQUANT未配置则跳过）", 13),
    ("apxs_fit_check.py", "APXS峰拟合：合成能谱峰高恢复、检出判定与越界不猜值", 9),
    ("paper_crosscheck_check.py", "SHERLOC逐靶对账摘要：复现数只增不缩、判定结构与类目映射注记", 7),
    ("theme_check.py", "整套主题、导航、配色偏好与原控件/图表/报告", 18),
    ("layout_check.py", "目录分离、数据原件、手册与打包路径", 12),
    ("data_audit.py", "独立数据审计：训练谱SHA/CV真留一复算/SFS与APXS峰面积复算/台账指纹/第二扫描一致性/PIQUANT复算", 30),
    ("paper_quant_check.py", "文献定量对账：Bellegarde盐校正+Rocknest3身份+比值矩阵+扫描内标定+APXS峰位自证+APXS比值层+PIQUANT×RQC双扫描对账", 32),
]

DIST_SCRIPTS={'multi_instrument_check.py','data_integrity_check.py','scientific_evidence_check.py',
              'fit_profile_check.py','repeatability_check.py','chemcam_passive_check.py','phase11_check.py','pixl_browser_check.py','pixl_download_check.py','workflow_alignment_check.py','instbkg_check.py','published_targets_check.py','gold_intervals_check.py','libs_quant_check.py','pixl_sfs_check.py','piquant_check.py','apxs_fit_check.py','paper_crosscheck_check.py','theme_check.py','layout_check.py'}


def run(script, args=()):
    """跑一个护栏脚本，返回 (退出码, 输出行列表)。"""
    env = dict(os.environ, PYTHONIOENCODING="utf-8")
    script, *extra = script.split()
    path = os.path.join(HERE, script)
    if not os.path.exists(path):
        return 2, ["找不到脚本：%s" % path]
    p = subprocess.run([sys.executable, path, *extra, *args],
                       cwd=ROOT, env=env,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    text = p.stdout.decode("utf-8", "replace")
    return p.returncode, text.splitlines()


def main(argv):
    global LAST_RUN
    started = datetime.now(timezone.utc).isoformat()
    tested_hashes = source_hashes()
    verbose = "-v" in argv or "--verbose" in argv
    with_auc = "--auc" in argv
    distribution = '--dist' in argv
    if distribution:
        if with_auc:raise ValueError('--dist只跑离线自包含检查，不与--auc混用')
        LAST_RUN=str(appenv.location('reports_selfcheck','dist_last_run.json'))
    checks=[row for row in CHECKS if not distribution or row[0] in DIST_SCRIPTS]

    print("=" * 78)
    print("MarsWorkbench 发行自包含检查（不依赖02/03开发数据）" if distribution else "SHERLOC 护栏总检")
    print("=" * 78)

    results = []
    for script, what, expect in checks:
        print("\n" + "-" * 78)
        print("跑 %s —— %s" % (script, what))
        print("-" * 78)
        code, lines = run(script,['--dist'] if distribution and script in ('phase11_check.py','layout_check.py') else [])
        if not lines:
            lines = ["（这个脚本没有输出）"]
        if verbose:
            for ln in lines:
                print(ln)
        else:
            # 失败的行一定在末尾（脚本最后才汇总），所以打印尾巴就够
            tail = 12 if code else 4
            if len(lines) > tail:
                print("  …（前面 %d 行略，加 -v 看全部）" % (len(lines) - tail))
            for ln in lines[-tail:]:
                print(ln)
        # 从输出里捞"共 N 项断言"，用来跟预期条数对账——
        # 断言被谁删掉了几条，这里会直接暴露。
        n_assert = None
        for ln in lines:
            if "项断言" in ln and "共" in ln:
                # "共 146 项断言，0 项失败" —— 只看"项断言"前面那一段，
                # 不然会把后面的失败数也当成数字拼进来（146 + 0 = 1460）。
                head = ln.split("项断言", 1)[0].split("共", 1)[1]
                digits = "".join(ch for ch in head if ch.isdigit())
                if digits:
                    n_assert = int(digits)
                    break
        skipped = (any(ln.startswith("SKIP") for ln in lines)
                      and not n_assert)   # 只有 0 断言才算整守卫跳过
        results.append({"script": script, "what": what, "code": code,
                        "n": n_assert, "expect": expect, "skipped": skipped})

    auc_result = None
    if with_auc:
        print("\n" + "-" * 78)
        print("附带：规则组与谱库分数的内部一致性（代理AUC，非独立矿物准确率）")
        print("-" * 78)
        code, lines = run("library_separation.py", ["8", "12", "16", "20"])
        auc_result = {'code': code, 'lines': lines, 'scope': '规则标签代理AUC，只报告，不参与工程通过判据'}
        for ln in lines:
            print(ln)

    print("\n" + "=" * 78)
    print("总检结果")
    print("=" * 78)
    bad = 0
    for r in results:
        n = "?" if r["n"] is None else r["n"]
        # 退出码为0但没有汇总，或实际检查数缩水，也不能显示为通过。
        if r.get("skipped"):
            ok = True            # 数据不在位的显式跳过：不计失败、免条数对账
        else:
            ok = (r["code"] == 0 and r["n"] is not None and r["n"] >= r["expect"])
        if not ok:
            bad += 1
        # 条数只设下限：这两个护栏都是**遍历数据文件夹**的，往 processed/ 里
        # 多放几个 sol，条数自然就涨了（2026-09-17：加了 sol 865 / 1853，
        # 从 146+152 涨到 165+182）。真正要报警的是**变少**——那说明有断言
        # 被删掉了，是很危险的信号。
        drift = ""
        if r["n"] is not None:
            if r["n"] < r["expect"]:
                drift = "  <== 比预期少（预期至少 %d）：有断言被删掉了？" % r["expect"]
            elif r["n"] > r["expect"]:
                drift = "  （比预期多 %d 条，多半是数据文件夹里多了 sol）" % (r["n"] - r["expect"])
        print("  %-6s %-22s %s 项%s" % (
            "跳过" if r.get("skipped") else ("通过" if ok else "失败"),
            r["script"], n, "" if r.get("skipped") else drift))

    print()
    try:
        write_snapshot(results, started, auc_result, tested_hashes,'distribution' if distribution else 'developer')
    except OSError as exc:
        print('无法保存本次总检记录：%s' % exc)
        return 1
    if source_hashes() != tested_hashes:
        print('总检期间源码发生变化，请重跑；本次记录不对应当前源码。')
        return 1
    if bad:
        print("有 %d 个护栏没通过——**不要提交，先查清楚是哪一步坏了。**" % bad)
        return 1
    skipped = sum(1 for r in results if r.get("skipped"))
    if skipped:
        print("全部通过。%d 条护栏、共 %s 项断言（%d 条因数据不在位显式跳过）。"              % (len(results) - skipped, sum(r["n"] or 0 for r in results), skipped))
    else:
        print("全部通过。%d 条护栏、共 %s 项断言。"
              % (len(results), sum(r["n"] or 0 for r in results)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
