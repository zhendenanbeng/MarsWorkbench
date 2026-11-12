"""Generate an auditable trust/boundary report from completed local checks. No analysis rerun."""

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout
import argparse
import ast
from collections import Counter
from datetime import datetime, timezone, timedelta
import hashlib
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
APP=_layout.project_layout.APP
ROOT=_layout.project_layout.ROOT
DOCS=_layout.project_layout.MANUALS
sys.path.insert(0,str(HERE))
import run_all
import appenv


def constants(path, names):
    result={}
    for node in ast.parse(path.read_text(encoding='utf-8')).body:
        if isinstance(node,ast.Assign):
            for target in node.targets:
                if isinstance(target,ast.Name) and target.id in names:
                    result[target.id]=ast.literal_eval(node.value)
    return result


def load_run(path, now):
    """Never present missing, stale, failed, or source-mismatched runs as current green."""
    try:
        data=json.loads(Path(path).read_text(encoding='utf-8'))
        stamp=datetime.fromisoformat(data['finished_utc'])
        if stamp.tzinfo is None:raise ValueError('时间戳缺少时区')
        age=now-stamp
        rows=data['results']
        expected={script:n for script,_what,n in run_all.CHECKS}
        valid=(data.get('version')==1 and isinstance(rows,list) and bool(rows)
               and data['suite_count']==len(rows)
               and len(rows)==len(expected) and {r['script'] for r in rows}==set(expected)
               and all(r['expect']>=expected[r['script']] for r in rows)
               and data['assertions']==sum(r['n'] or 0 for r in rows)
               and all(r.get('ok')==(r['code']==0 and r['n'] is not None and r['n']>=r['expect']) for r in rows)
               and data['failures']==sum(not r['ok'] for r in rows))
        if not valid:raise ValueError('总检结构/条数不一致')
        if age>timedelta(days=7) or age<timedelta(minutes=-5):
            return data,'数据过期，请重跑 run_all（时间戳超期或来自未来）'
        if data['failures']:return data,'总检存在失败，请修复后重跑 run_all'
        if data.get('source_sha256')!=run_all.source_hashes():
            return data,'源码已变化或记录未绑定源码，请重跑 run_all'
        return data,''
    except (OSError,ValueError,KeyError,TypeError,OverflowError):
        return None,'数据过期，请重跑 run_all（缺少或无法核验 _last_run.json）'


def generate(last_run=None, now=None):
    now=now or datetime.now(timezone.utc)
    if now.tzinfo is None:raise ValueError('now须带时区')
    data,warning=load_run(last_run or Path(run_all.LAST_RUN),now)
    config=constants(HERE/'regression.py',{'GOLD','NEW_GOLD','KNOWN_GAP','KNOWN_FALSE_POSITIVE'})
    ledger=json.loads((HERE/'gold_ledger.json').read_text(encoding='utf-8'))
    gold=sum(len(v) for v in config['GOLD'].values())
    new=len(config['NEW_GOLD'])
    rows=ledger['rows']
    text=['# 软件可信度与已知边界报告',
          '\n生成时间：'+now.astimezone(timezone(timedelta(hours=8))).isoformat(timespec='seconds')+'（Asia/Shanghai）。此报告是运行记录与已知边界的自动汇编，不是科学准确率证书。']
    if warning:text.append('\n<span style="color:#c62828;font-weight:bold">'+warning+'</span>')
    text+=['\n## 护栏总览']
    if data:
        text+=['\n最近记录：'+str(data.get('finished_utc','未知'))+'；是否当前可用：'+('否' if warning else '是')+'。',
               '\n记录内总计：%s组 / %s项；失败组数：%s。数据无效或过期时，这些数字只作历史记录。'%(data.get('suite_count','?'),data.get('assertions','?'),data.get('failures','?')),
               '\n| 脚本 | 实际/预期下限 | 结果 |\n|---|---|---|']
        for r in data.get('results',[]):
            text.append('| `%s` | %s / %s | %s |'%(r.get('script','?'),r.get('n','?'),r.get('expect','?'),'通过' if r.get('ok') else '失败'))
    else:text.append('\n没有可核验的总检记录；未把缺失当作零失败。')
    text+=['\n独立数值审计不由此脚本自动重跑，也不沿用旧条数冒充本次结果。可查看各期审计原件。',
           '\n## 金标准与差异台账',
           '\n当前GOLD：%d个类别区间；NEW_GOLD：%d条扫描/类别区间；合计%d条。金标准防结果漂移，区间不是矿物真值。'%(gold,new,gold+new),
           '\n`gold_ledger.json`共%d行；它是既有对照报告的历史快照，工具点数不是当前运行结果。'%(len(rows)),
           '\n| 类别 | 台账行数 | 差异性质摘要 |\n|---|---|---|']
    for category in sorted({r['category'] for r in rows}):
        counts=Counter(r.get('kind','未记录') for r in rows if r['category']==category)
        text.append('| %s | %d | %s |'%(category,sum(counts.values()),'；'.join('%s %d'%(k,n) for k,n in sorted(counts.items()))))
    text+=['\nKNOWN_FALSE_POSITIVE：'+ '；'.join('Sol%s / %s，上限%d点（工程容忍上限，不是科学放行）'%(sol,seq,cap) for sol,seq,cap,_ in config['KNOWN_FALSE_POSITIVE']),
           '\n## 已知边界',
           '\n### R3：四条缺口及不自动补的原因']
    reasons={('293','未定硅酸盐'):'弱带落入碳酸盐主导组或未归类，强拆组会改变已通过区间。',
             ('293','高氯酸盐或磷酸盐'):'同族已用更细族名检出，不能为贴论文词汇而改名。',
             ('257','辉石'):'与硫酸盐保留为缺次级峰的模糊类，证据不足以自动拆分。',
             ('257','硫酸盐'):'同一模糊组不能按论文点数分摊，需要独立判读证据。'}
    for sol,names in config['KNOWN_GAP'].items():
        for name in names:text.append('\n- Sol%s %s：%s'%(sol,name,reasons.get((sol,name),'新增缺口：尚未补写自动处理依据，须复核。')))
    text+=['\n依据：[R3缺口复核](第五期验证_2026-10-24/R3_缺口复核结论.md)。',
           '\n### Sol349：三类已尝试的判别量均未给出可靠分离',
           '\n峰形、背景竞争与空间一致性实验未获得兼顾保留样例和排除该样例的可靠阈值。这是所测开发样例和方法下的失败，不是证明所有未来方法都不可能。',
           '\n依据：[峰形实验](第五期验证_2026-10-24/R2_假阳性实验结论.md)、[背景竞争](第七期验证_2026-10-24/T3_背景竞争实验.md)、[空间实验](第九期验证_2026-10-24/第九期结论_主线止损与空间实验.md)。',
           '\n### SuperCam投影：两期止损与再战条件',
           '\n第八期相机几何研究完成部分闭合核验，灼蚀斑闸门未能判定；第九期跨仪器帧/共靶关系未闭合，未接产品。需官方帧关系文档或两仪器共视的标定目标，并重新核验对应表面及误差；图像能显示不等于可以投影。',
           '\n依据：[第八期闸门](第八期验证_2026-10-24/T1T2_CI投影与闸门结论.md)、[第九期结论](第九期验证_2026-10-24/第九期结论_主线止损与空间实验.md)。',
           '\n明确不做：ChemCam定量PLS、为凑论文类名修改分类器、未验证的SuperCam叠图。',
           '\n## 各期验证文档索引',
           '\n以下为当前文档目录树内的Markdown记录；历史结论有适用时点，不能自动当作当前功能状态。',
           '\n| 目录 | 记录 |\n|---|---|']
    for directory in sorted(p for p in DOCS.iterdir() if p.is_dir()):
        for path in sorted(directory.rglob('*.md')):
            rel=path.relative_to(DOCS).as_posix()
            text.append('| %s | [%s](<%s>) |'%(directory.name,path.stem,rel))
    text+=['\n## 口径申明',
           '\n快速浏览不等于定量。检查通过不等于矿物鉴定准确率；当前缺少独立盲标真值，不能由工程护栏数计算假阳性/假阴性率。',
           '\n逐发证据是当前预处理后保留发的经验检出/支持比例，发次独立性未获验证；SNR和FWHM是局部统计量，截断†是宽度下界。经验比例不是置信概率或浓度。',
           '\n用户论文标注只记录人工复核，不反向训练默认判定。谱库代理AUC比较软件规则组，标签不是独立矿物真值；无论AUC高低，均不能单独证明矿物学正确性。']
    return '\n'.join(text)+'\n'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path)
    parser.add_argument('--last-run',type=Path,default=Path(run_all.LAST_RUN))
    args=parser.parse_args()
    date=datetime.now(timezone(timedelta(hours=8))).date().isoformat()
    out=args.out or (_layout.project_layout.INSTALL/'构建缓存/诊断')/('可信度报告_'+date+'.md')
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(generate(args.last_run),encoding='utf-8')
    print(out)
    return 0

if __name__=='__main__':raise SystemExit(main())
