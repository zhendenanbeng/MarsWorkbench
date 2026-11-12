"""Download and audit the proposed Sol1000 raster without assuming shared RMI coordinates."""

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import re
import sys
import time
import urllib.parse
import urllib.request
import numpy as np
from PIL import Image

ROOT=_layout.project_layout.APP;sys.path.insert(0,str(ROOT))
from instruments import pds3
ORIGIN='https://pds-geosciences.wustl.edu'
VOLUME=ORIGIN+'/msl/msl-m-chemcam-libs-4_5-rdr-v1/mslccm_1xxx/'

def fetch(url,dest):
    dest=Path(dest)
    if not dest.exists():
        for attempt in range(3):
            try:
                data=urllib.request.urlopen(url,timeout=90).read();break
            except OSError:
                if attempt==2:raise
                time.sleep(1+attempt)
        dest.parent.mkdir(parents=True,exist_ok=True)
        with dest.open('xb') as stream:stream.write(data)
    raw=dest.read_bytes()
    return {'url':url,'path':str(dest.resolve()),'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache',type=Path,default=_layout.project_layout.legacy_cache_root()/'phase3_probe')
    parser.add_argument('--out',type=Path,required=True)
    a=parser.parse_args()
    if a.out.exists():parser.error('Use a new audit filename')
    listing=urllib.request.urlopen(VOLUME+'data/sol01000/',timeout=45).read().decode()
    links=re.findall(r'href=["\x27]([^"\x27]+)',listing,re.I)
    names=[Path(u).name for u in links if re.match(r'(cl5_\d+ccs_.*ccam02000p3\.(csv|lbl)|cr0_\d+prc_.*ccam02000l1\.(tif|lbl))$',Path(u).name,re.I)]
    tasks=[(VOLUME+'data/sol01000/'+n,a.cache/'sol01000'/n) for n in names]
    for n in ['ccam_default_wave.tab','ccam_default_wave.lbl']:
        tasks.append((VOLUME+'document/'+n,a.cache/'document'/n))
    for directory,stem in [('sol00010','cl1_398380611psv_f0030000ccam04010p3'),('sol01000','cl9_486263266psv_f0481570ccam02000p3')]:
        for ext in ['tab','lbl']:tasks.append((VOLUME+'data/'+directory+'/'+stem+'.'+ext,a.cache/directory/(stem+'.'+ext)))
    with ThreadPoolExecutor(max_workers=2) as pool:files=list(pool.map(lambda t:fetch(*t),tasks))
    points=[];images=[]
    for n in names:
        p=a.cache/'sol01000'/n
        if n.startswith('cl5_') and n.endswith('.lbl'):
            text=p.read_text();state=re.search(r'GROUP\s*=\s*INSTRUMENT_STATE_PARMS(.*?)END_GROUP\s*=\s*INSTRUMENT_STATE_PARMS',text,re.S).group(1)
            lab=pds3.parse_lbl(p)
            coord=lambda key:float(re.search(r'MSL:'+key+r'\s*=\s*([-+0-9.]+)',state).group(1))
            points.append({'label':n,'sclk':n[4:13],'line':coord('SPECIAL_LINE'),'sample':coord('SPECIAL_SAMPLE'),
                           'azimuth':lab.get('INSTRUMENT_AZIMUTH'),'elevation':lab.get('INSTRUMENT_ELEVATION')})
        elif n.endswith('.tif'):
            with Image.open(p) as im:images.append({'file':n,'width':im.width,'height':im.height,'mode':im.mode})
    unique={(p['line'],p['sample']) for p in points}
    gate={'nine_points':len(points)==9,'two_images':len(images)==2,'unique_coordinates':len(unique)==len(points),
          'inside_images':all(0<=p['sample']<im['width'] and 0<=p['line']<im['height'] for p in points for im in images)}
    wave=np.loadtxt(a.cache/'document/ccam_default_wave.tab',ndmin=2)
    tables={p.name:list(pds3.read_table(p,delimiter=None,skiprows=29).shape) for p in a.cache.glob('sol*/cl*psv*.tab')}
    result={'generated_utc':datetime.now(timezone.utc).isoformat(),'files':files,'points':sorted(points,key=lambda p:p['sclk']),
            'images':images,'gate':gate,'shared_overlay_authorized_by_data':all(gate.values()),
            'wave_shape':list(wave.shape),'passive_shapes':tables}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    with a.out.open('x',encoding='utf-8') as f:json.dump(result,f,ensure_ascii=False,indent=2)
    print(json.dumps({k:v for k,v in result.items() if k!='files'},ensure_ascii=False,indent=2))

if __name__=='__main__':main()
