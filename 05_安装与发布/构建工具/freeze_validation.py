"""Freeze a successful current complete run into installation metadata."""
from pathlib import Path
import json,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'02_图形界面代码'))
import source_paths,appenv
sys.path.insert(0,str(source_paths.project_layout.VERIFY))
import run_all

def main():
    source=Path(run_all.LAST_RUN);data=json.loads(source.read_text(encoding='utf-8'))
    if data['scope']!='developer' or data['failures'] or data['assertions']<992 or data['source_sha256']!=run_all.source_hashes():raise RuntimeError('需要与当前源码一致的完整通过记录')
    target=source_paths.project_layout.INSTALL/'资源模板/_last_run.json'
    target.write_bytes(source.read_bytes())
    appenv.resources()
    print('Frozen %d assertions from %s'%(data['assertions'],data['finished_utc']))
if __name__=='__main__':main()
