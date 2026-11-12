# -*- coding: utf-8 -*-
"""对称MAD背景的受控噪声与真实开发样例试验，不修改默认算法。"""

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout
import argparse
import hashlib
import json
from pathlib import Path
import sys
from unittest.mock import patch
import numpy as np

ROOT=_layout.project_layout.APP;sys.path.insert(0,str(ROOT));sys.path.insert(0,str(_layout.project_layout.VERIFY))
from sherloc import minerals,qc
from sherloc.background import symmetric_baseline
import science_benchmark


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',required=True,type=Path)
    parser.add_argument('--mode',choices=['symmetric','quiet_centered'],default='symmetric')
    args=parser.parse_args()
    noise=args.out.with_name(args.out.stem+'_noise.json')
    if args.out.exists() or noise.exists():parser.error('请用新的输出名，保留已有试验及噪声证据')
    args.out.parent.mkdir(parents=True,exist_ok=True)
    original=minerals.baseline
    def quiet_centered(x,y,**kw):
        base=original(x,y,**kw);quiet=(x>=2500)&(x<=3000)
        return base+np.median((y-base)[quiet]) if quiet.sum()>=3 else base
    selected=symmetric_baseline if args.mode=='symmetric' else quiet_centered
    rng=np.random.default_rng(5817);x=np.linspace(800,4000,340)
    data=rng.normal(0,1,(500,len(x)));records={}
    for name,fn in [('original',original),(args.mode,selected)]:
        residual=np.array([row-fn(x,row) for row in data])
        records[name]={'mean_residual':float(residual.mean()),'rms':float(np.sqrt(np.mean(residual**2))),
                       'fraction_above_3':float(np.mean(residual>3)),
                       'mean_baseline':float((data-residual).mean())}
    print('Controlled noise:',records,flush=True)
    noise.write_text(json.dumps({'seed':5817,'rows':len(data),'channels':len(x),'truth':'zero background, independent unit Gaussian noise','results':records},indent=2),encoding='utf-8')
    with patch.object(minerals,'baseline',side_effect=selected),patch.object(qc,'poly_baseline',side_effect=selected),patch.object(sys,'argv',['science_benchmark.py','--out',str(args.out)]):
        science_benchmark.main()
    result=json.loads(args.out.read_text(encoding='utf-8'))
    result['experimental_override']={'functions':['minerals.baseline','qc.poly_baseline'],'method':args.mode,
                                     'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')


if __name__=='__main__':main()
