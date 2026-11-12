# -*- coding: utf-8 -*-
"""只读从 PDS 重新下载冻结样例，核对字节哈希；不覆盖本地产品。"""

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import urllib.request

ROOT=_layout.project_layout.APP
WORKSPACE=_layout.project_layout.ROOT
HOST='https://pds-geosciences.wustl.edu'
PREFIX={
    'apxs':'/msl/msl-m-apxs-4_5-rdr-v1/mslapx_1xxx/data/',
    'chemcam':'/msl/msl-m-chemcam-libs-4_5-rdr-v1/mslccm_1xxx/data/',
    'supercam':'/m2020/urn-nasa-pds-mars2020_supercam/data_calibrated_spectra/',
    'sherloc':'/m2020/urn-nasa-pds-mars2020_sherloc/data_processed/',
}


def source(rec):
    p=Path(rec['path'])
    if '_moc' in p.parts:return None  # 本地三行摘录不是官方原文件；全表另行比对。
    if 'fixtures' in p.parts:
        index=p.parts.index('fixtures');instrument=p.parts[index+1]
    else:instrument='sherloc'
    return HOST+PREFIX[instrument]+p.parent.name+'/'+p.name.lower()


def fetch(rec):
    url=source(rec)
    result={'path':rec['path'],'url':url,'expected_sha256':rec['sha256']}
    try:
        with urllib.request.urlopen(url,timeout=60) as r:
            raw=r.read();modified=r.headers.get('Last-Modified')
        digest=hashlib.sha256(raw).hexdigest()
        local=WORKSPACE/rec['path']
        result.update(remote_sha256=digest,bytes=len(raw),last_modified=modified,
                      passed=digest==rec['sha256']==hashlib.sha256(local.read_bytes()).hexdigest())
    except Exception as exc:result.update(passed=False,error=str(exc))
    print(('PASS ' if result['passed'] else 'FAIL ')+rec['path'],flush=True)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',required=True,type=Path)
    args=parser.parse_args()
    if args.out.exists():parser.error('输出已存在；换用新的证据文件名。')
    frozen=_layout.project_layout.VERIFY/'independent_manifest.json'
    manifest=json.loads(frozen.read_text(encoding='utf-8'))
    records=[r for r in manifest['files'] if source(r)]
    with ThreadPoolExecutor(max_workers=4) as pool:results=list(pool.map(fetch,records))
    failures=sum(not r['passed'] for r in results)
    output={'generated_utc':datetime.now(timezone.utc).isoformat(),
            'method':'PDS重新下载SHA256 = 冻结SHA256 = 本地当前SHA256；不写回源文件',
            'manifest_sha256':hashlib.sha256(frozen.read_bytes()).hexdigest(),
            'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'results':results,'checks':len(results),'failures':failures,
            'limits':'字节一致性证明来源与版本一致，不验证仪器团队标定或矿物结论本身。'}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    with args.out.open('x',encoding='utf-8') as f:json.dump(output,f,ensure_ascii=False,indent=2)
    print('共 %d 项断言，%d 项失败'%(len(results),failures))
    return int(failures>0)


if __name__=='__main__':raise SystemExit(main())
