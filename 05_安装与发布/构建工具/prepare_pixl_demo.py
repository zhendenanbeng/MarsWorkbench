"""Validate original PIXL example CSV/XML/MCC bytes."""
from pathlib import Path
import sys,json,hashlib
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'02_图形界面代码'))
import source_paths,appenv

def main():
    root=appenv.resources()/'pixl-demo'
    rows=json.loads((root/'manifest.json').read_text(encoding='utf-8'))['files']
    for row in rows:
        if hashlib.sha256((root/row['file']).read_bytes()).hexdigest()!=row['sha256']:raise ValueError('PIXL示例校验失败：'+row['file'])
    print('Verified',len(rows),'PIXL original files')
if __name__=='__main__':main()
