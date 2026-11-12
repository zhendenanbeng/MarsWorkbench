"""Time-boxed score-only experiments; rule labels are proxies, never mineral truth.

Three predeclared probes: remove relative SNR weights, mode-wise reference consensus,
and observed local-pair strength. Five Sol-blocked outer folds; choose the method
using training AUC only, then evaluate held-out Sols. No production implementation.
"""

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout
import argparse
from collections import Counter,defaultdict
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import sys
from unittest.mock import patch

import numpy as np

HERE=Path(__file__).resolve().parent
APP=_layout.project_layout.APP
sys.path.insert(0,str(APP))
from sherloc import library,qc
from sherloc.organize import scan_folder
from library_separation import auc

METHODS=('baseline','unweighted_snr','mineral_consensus','local_pair_ratio')


def score_auc(rows,method):
    pos=[r['scores'][method] for r in rows if r['named']]
    neg=[r['scores'][method] for r in rows if not r['named']]
    return auc(pos,neg) if pos and neg else None


def consensus_table():
    table=deepcopy(library.REFERENCE)
    for cls in table.values():
        by_mineral=defaultdict(list)
        for mineral,sample,modes in cls['samples']:
            by_mineral[mineral].append(modes)
        samples=[]
        for mineral,variants in by_mineral.items():
            modes={}
            for mode in {m for variant in variants for m in variant}:
                entries=[]
                for variant in variants:
                    value=variant.get(mode)
                    values=value if isinstance(value,tuple) else (value,)
                    values=sorted(float(v) for v in values if v is not None)
                    if values:entries.append(values)
                if not entries:continue
                centers=[float(np.median([v[i] for v in entries if len(v)>i])) for i in range(max(map(len,entries)))]
                modes[mode]=centers[0] if len(centers)==1 else tuple(centers)
            samples.append((mineral,'mode-wise median (diagnostic)',modes))
        cls['samples']=samples
        cls.pop('intensities',None)
    return table


def local_pair_score(observed, matches):
    """Observed-only pair-strength probe; no invented reference intensity ratios."""
    values=[]
    for candidate in matches:
        used=set();strengths=[]
        refs=[v for v in candidate['reference_peaks'] if 800<=v<=1800]
        for center,snr in sorted(observed,key=lambda r:-r[1]):
            possible=[(abs(center-v),j) for j,v in enumerate(refs) if j not in used and abs(center-v)<=8]
            if possible:
                _,j=min(possible);used.add(j);strengths.append(float(snr))
        pair=(max(0,strengths[1]/strengths[0])**.5) if len(strengths)>=2 and strengths[0]>0 else 0
        values.append(candidate['score']*(.5+.5*pair))
    return max(values,default=0.)


def collect():
    rows=[];reference=consensus_table();invariant=True;inputs={}
    _,observations=scan_folder(str(APP.parent/'03_下载数据/SHERLOC/练习数据'))
    for ob in observations:
        if ob.spectra is None:continue
        wl,r1,_r2,_r3=ob.read_regions()
        if wl is None or r1 is None:continue
        source=[]
        for part in ob.parts:
            for name in (part.data_path,part.label.path):
                if name:
                    path=Path(name).resolve();relative=path.relative_to(APP.parent).as_posix()
                    inputs[relative]=hashlib.sha256(path.read_bytes()).hexdigest();source.append(relative)
        result=qc.analyse_spectra(wl,r1) # exactly library_separation.collect's R1-only semantics
        for index,group in enumerate(result.get('mineral_groups') or []):
            peaks=group.get('observed_peaks') or []
            if not peaks:continue
            key=group['classes'][0]['key'] if group['classes'] else 'unassigned'
            baseline=library.match(peaks,tol=8)
            scaled=library.match([(p,w*100) for p,w in peaks],tol=8)
            base_scores={(r['class_key'],r['mineral'],r['sample']):r['score'] for r in baseline}
            scale_scores={(r['class_key'],r['mineral'],r['sample']):r['score'] for r in scaled}
            invariant=invariant and base_scores==scale_scores
            with patch.object(library,'REFERENCE',reference):
                combined=library.match(peaks,tol=8)
            top=baseline[0] if baseline else {}
            rows.append({'sol':str(ob.sol),'seq':ob.sequence_id,'clock':ob.spectra.sclk,'sources':source,
                         'group':index,'named':key!='unassigned','rule_key':key,'peaks':peaks,
                         'library_class':top.get('class_key','none'),'library_mineral':top.get('mineral','none'),
                         'scores':{'baseline':max((r['score'] for r in baseline),default=0),
                                   'unweighted_snr':max((r['score'] for r in library.match([p for p,w in peaks],tol=8)),default=0),
                                   'mineral_consensus':max((r['score'] for r in combined),default=0),
                                   'local_pair_ratio':local_pair_score(peaks,baseline)}})
        print('collected Sol',ob.sol,ob.sequence_id,len(rows),flush=True)
    return rows,invariant,inputs


def diagnose():
    rows,invariant,inputs=collect()
    sols=sorted({r['sol'] for r in rows},key=int)
    folds=[]
    for index in range(5):
        held=set(sols[index::5]);train=[r for r in rows if r['sol'] not in held];test=[r for r in rows if r['sol'] in held]
        train_auc={m:score_auc(train,m) for m in METHODS}
        # Baseline is eligible; tie order favors the unchanged method.
        eligible=[m for m in METHODS if train_auc[m] is not None]
        selected=max(eligible,key=lambda m:train_auc[m]) if eligible else 'baseline'
        folds.append({'test_sols':sorted(held,key=int),'train_sols':sorted(set(sols)-held,key=int),
                      'train_n':len(train),'test_n':len(test),'selected_from_training_only':selected,
                      'train_auc':train_auc,'test_auc':{m:score_auc(test,m) for m in METHODS}})
    macro={m:float(np.mean([f['test_auc'][m] for f in folds if f['test_auc'][m] is not None]))
           if any(f['test_auc'][m] is not None for f in folds) else None for m in METHODS}
    chosen=[f['test_auc'][f['selected_from_training_only']] for f in folds]
    nested=float(np.mean([v for v in chosen if v is not None])) if any(v is not None for v in chosen) else None
    cross=Counter((r['rule_key'],r['library_class']) for r in rows)
    inversions=Counter()
    for pos in (r for r in rows if r['named']):
        for neg in (r for r in rows if not r['named']):
            p,n=pos['scores']['baseline'],neg['scores']['baseline']
            if p<=n:inversions[(pos['rule_key'],neg['library_class'])]+=1 if p<n else .5
    visible=[len([v for v in library.sample_peaks(modes) if 800<=v<=1800]) for cls in library.REFERENCE.values() for _,_,modes in cls['samples']]
    improved=[m for m in METHODS[1:] if macro[m] is not None and macro[m]>.70]
    decision='close_no_heldout_path' if not improved or nested is None or nested<=.70 else 'candidate_requires_separate_validation'
    return {'scope':'规则标签代理AUC，不是矿物真值混淆矩阵；全部为开发样例的Sol分组留出复用，不是独立外部盲测。',
            'row_count':len(rows),'named':sum(r['named'] for r in rows),'unnamed':sum(not r['named'] for r in rows),
            'global_snr_scale_invariant':invariant,'reference_samples':len(visible),
            'reference_window_empty':sum(n==0 for n in visible),'reference_window_single_peak':sum(n==1 for n in visible),
            'full_data_auc':{m:score_auc(rows,m) for m in METHODS},'sol_blocked_macro_auc':macro,
            'train_selected_heldout_macro_auc':nested,'folds':folds,
            'rule_vs_library_counts':[{'rule':a,'library':b,'n':n} for (a,b),n in cross.most_common()],
            'baseline_ranking_inversion_burden':[{'positive_rule':a,'negative_library':b,'weighted_pairs':n} for (a,b),n in inversions.most_common()],
            'decision':decision,'default_changed':False,
            'limitations':['库内是峰位表，不含完整代表谱；mode-wise median不是实测光谱均值。',
                           '输入只保留最多四个拟合峰；统一噪声尺度可改变上游选峰，但在当前打分的最大值归一化后不改变相对权重。',
                           '局部峰比试验只用观测侧，库中缺成对强度标定，不能当物理峰比模型。',
                           '不同化学族或混合/模糊规则类不能一一对应；交叉表和排序倒置不是矿物误分类率。'],
            'production_library_sha256':hashlib.sha256((APP/'sherloc/library.py').read_bytes()).hexdigest(),
            'production_source_sha256':{p.relative_to(APP).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((APP/'sherloc').glob('*.py'))},
            'input_sha256':inputs,'rows':rows}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    if args.out.exists():parser.error('请用新路径保留既有诊断证据')
    data=diagnose();args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps({k:data[k] for k in ['row_count','named','unnamed','full_data_auc','sol_blocked_macro_auc','train_selected_heldout_macro_auc','decision']},ensure_ascii=False))

if __name__=='__main__':main()
