# -*- coding: utf-8 -*-
"""论文靶点的前后对照；记录偏差，不把开发集包装成盲测准确率。"""

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout
from collections import Counter
from datetime import datetime, timezone
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT=_layout.project_layout.APP
sys.path.insert(0,str(ROOT))
from sherloc import qc
from sherloc.organize import scan_folder
from instruments import REGISTRY


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',required=True,type=Path)
    args=parser.parse_args()
    if args.out.exists():parser.error('请用新的输出文件名保留前后证据')
    _,scans=scan_folder(str(ROOT.parent/'03_下载数据/SHERLOC/练习数据'))
    rows=[]
    wanted={('141','srlc10203'),('349','srlc11420'),('293','srlc11360'),
            ('257','srlc11360'),('370','srlc11360'),('162','srlc10500'),('186','srlc10500')}
    for scan in scans:
        if (str(scan.sol),scan.sequence_id) not in wanted or scan.spectra is None:continue
        wl,r1,r2,r3=scan.read_regions()
        if wl is None or r1 is None:continue
        res=qc.analyse_spectra(wl,r1,r2,r3)
        classes=Counter(c['class'] for p in res['point_detail'] for c in p['classes'] if c['key']!='unassigned')
        entry={'sol':str(scan.sol),'seq':scan.sequence_id,'points':len(r1),'classes':dict(classes),
               'named_point_assignments':sum(classes.values()),
               'raw_band_points':res['n_points_with_band'],  # 兼容旧审计；旧名实际指峰组点数
               'points_in_groups':res['n_points_in_groups'],
               'points_with_candidate_bands':res['n_points_with_candidate_bands'],
               'points_with_named_classes':res['n_points_with_named_classes'],
               'groups':[{'n':g['n_points'],'center':g['band']['center'],'fwhm':g['band']['fwhm'],
                          'snr':g['band']['snr'],'classes':[c['key'] for c in g['classes']],
                          'method':g['band']['method']} for g in res['mineral_groups']]}
        rows.append(entry)
        print(entry['sol'],entry['seq'],entry['points'],entry['classes'],flush=True)
    other=[]
    for name,b in REGISTRY.items():
        for obs in b.scan_folder(str(_layout.project_layout.VERIFY/'fixtures'/name)):
            for pt in obs.points:
                d=b.load_point(pt,None,{k:v for k,_label,v in b.steps})
                if d.get('passive'):continue
                marks,elements=b.annotate(np.asarray(d['x']),np.asarray(d['y']))
                other.append({'instrument':name,'sol':obs.sol,'product':obs.product,'sclk':pt.sclk,
                              'elements':elements,'marks':marks})
    result={'generated_utc':datetime.now(timezone.utc).isoformat(),'sherloc':rows,'other_instruments':other,
            'code_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                           for folder in ('sherloc','instruments') for p in (ROOT/folder).glob('*.py')},
            'scope':'开发集上的论文靶点对照与数值变更审计；不计算独立科学准确率。'}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    with args.out.open('x',encoding='utf-8') as f:json.dump(result,f,ensure_ascii=False,indent=2,allow_nan=False)
    return 0


if __name__=='__main__':raise SystemExit(main())
