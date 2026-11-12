"""Bind the successful frozen run to working and Git-index source bytes."""

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[2]/'02_图形界面代码'))
import source_paths as _layout
import hashlib
import json
from pathlib import Path
import subprocess
import sys

PROJECT=_layout.project_layout.ROOT
APP=_layout.project_layout.APP
sys.path.insert(0,str(APP))
import appenv

def main():
    snapshot=json.loads((appenv.resources()/'_last_run.json').read_text(encoding='utf-8'))
    bad_work=[];bad_index=[]
    for relative,sha in snapshot['source_sha256'].items():
        path=PROJECT/relative
        if hashlib.sha256(path.read_bytes()).hexdigest()!=sha:bad_work.append(relative)
        staged=subprocess.check_output(['git','show',':'+path.relative_to(PROJECT).as_posix()],cwd=PROJECT)
        if hashlib.sha256(staged).hexdigest()!=sha:bad_index.append(relative)
    result={'frozen_assertions':snapshot['assertions'],'source_files':len(snapshot['source_sha256']),
            'working_mismatch':bad_work,'git_index_mismatch':bad_index,'all_match':not bad_work and not bad_index}
    evidence=_layout.project_layout.INSTALL/'机器验收'/('v'+appenv.VERSION);evidence.mkdir(exist_ok=True)
    (evidence/'源码与Git索引字节核对.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    assert result['all_match'] and snapshot['failures']==0,result
    print(json.dumps(result,ensure_ascii=False))

if __name__=='__main__':main()
