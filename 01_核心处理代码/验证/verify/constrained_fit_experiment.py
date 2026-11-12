# -*- coding: utf-8 -*-
"""先在受控上下文替换线性峰形求解器，评估物理约束是否破坏已知信号。"""

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

ROOT=_layout.project_layout.APP
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(_layout.project_layout.VERIFY))
from sherloc import minerals
from sherloc.fitting import solve_nonnegative_profiles
import science_benchmark


def constrained(x,y,comps,nonnegative=False,context=None):
    profiles=np.column_stack([part for center,width in comps for part in minerals._voigt_basis(x,center,width)])
    return solve_nonnegative_profiles(x,y,profiles,context=context)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',required=True,type=Path)
    parser.add_argument('--mode',choices=['constrained','reject_invalid'],default='constrained')
    args=parser.parse_args()
    original=minerals.fit_voigt
    count=[0]
    def reject_invalid(*values,**kw):
        result=original(*values,**kw)
        if result and (result['c1']<0 or result['c2']<0):
            count[0]+=1;return None
        return result
    hook=patch.object(minerals,'_solve_linear',side_effect=constrained) if args.mode=='constrained' else patch.object(minerals,'fit_voigt',side_effect=reject_invalid)
    with hook,patch.object(sys,'argv',['science_benchmark.py','--out',str(args.out)]):
        science_benchmark.main()
    result=json.loads(args.out.read_text(encoding='utf-8'))
    result['experimental_override']={'mode':args.mode,'invalid_profiles_rejected':count[0],
                                     'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')


if __name__=='__main__':main()
