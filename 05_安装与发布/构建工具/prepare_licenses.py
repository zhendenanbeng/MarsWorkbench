"""Verify retained notices and their checksum ledger before installation builds."""
from pathlib import Path
import sys,json,hashlib
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'02_图形界面代码'))
import source_paths,appenv

def main():
    root=appenv.resources()/'licenses';manifest=json.loads((root/'manifest.json').read_text(encoding='utf-8'))
    for row in manifest['files']:
        if hashlib.sha256((root/row['file']).read_bytes()).hexdigest()!=row['sha256']:raise ValueError('许可校验失败：'+row['file'])
    for name in ('LICENSE','THIRD_PARTY_NOTICES.md','DATA_SOURCES.md','Saira-OFL.txt'):
        if not (root/name).is_file():raise FileNotFoundError(name)
    print('Verified license ledger and notices')
if __name__=='__main__':main()
