"""Validate canonical public examples and assemble installation assets; never alter originals."""
from pathlib import Path
import sys,json,hashlib
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'02_图形界面代码'))
import source_paths,appenv

def main():
    resources=appenv.resources()
    rows=json.loads((resources/'sample-data/manifest.json').read_text(encoding='utf-8'))['files']
    for row in rows:
        path=resources/row['path']
        if hashlib.sha256(path.read_bytes()).hexdigest()!=row['sha256']:raise ValueError('示例校验失败：'+row['path'])
    print('Verified',len(rows),'sample/cache files')
if __name__=='__main__':main()
