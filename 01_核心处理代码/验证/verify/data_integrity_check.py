# -*- coding: utf-8 -*-
"""不依赖 Astropy 的数据完整性护栏：缩放整数和 APXS 成分配对。"""

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout
import csv
from pathlib import Path
import struct
import sys
import tempfile
from unittest.mock import patch

import numpy as np

ROOT=_layout.project_layout.APP
sys.path.insert(0,str(ROOT))
from instruments import fitsio
from instruments.apxs import APXSBackend
from instruments.report import build_report
import console

TOTAL=FAIL=0
def check(ok,name):
    global TOTAL,FAIL
    TOTAL+=1
    if not ok:FAIL+=1;print(' FAIL '+name)


def header(cards):
    data=''.join((k.ljust(8)+'= '+v).ljust(80) for k,v in cards)+ 'END'.ljust(80)
    return data.ljust((len(data)+2879)//2880*2880).encode('ascii')


def scaling(temp):
    # 用 FITS 标准的字节布局直接写入已知存储值，避免待测解析器生成自己的真值。
    primary=header([('SIMPLE','T'),('BITPIX','8'),('NAXIS','0'),('EXTEND','T')])
    table=header([('XTENSION',"'BINTABLE'"),('BITPIX','8'),('NAXIS','2'),('NAXIS1','4'),
                  ('NAXIS2','3'),('PCOUNT','0'),('GCOUNT','1'),('TFIELDS','2'),
                  ('EXTNAME',"'TEST'"),('TTYPE1',"'Unsigned'"),('TFORM1',"'I'"),('TZERO1','32768'),
                  ('TTYPE2',"'Scaled'"),('TFORM2',"'I'"),('TSCAL2','2.5D0'),('TZERO2','-3D0')])
    raw=b''.join(struct.pack('>hh',a,b) for a,b in [(-32768,1),(0,-2),(32767,10)])
    path=temp/'scaled.fits';path.write_bytes(primary+table+raw.ljust(2880,b'\0'))
    with patch.dict(sys.modules,{'astropy':None,'astropy.io':None,'astropy.io.fits':None}):
        result=fitsio.read_bintable(str(path),'TEST')
    check(np.array_equal(result['Unsigned'],[0,32768,65535]),'TZERO无符号整数须正确还原')
    check(result['Unsigned'].dtype.kind=='u','无符号整数保留整数精度')
    check(np.array_equal(result['Scaled'],[-.5,-8.,22.]),'TSCAL/TZERO和D指数须正确还原')


def paired_apxs(temp):
    source=_layout.project_layout.VERIFY/'fixtures/apxs/sol00102'
    folder=temp/'sol00102';folder.mkdir()
    # 三个同序列观测：不同FeO、缺成分表、高版本覆盖旧成分，防止按目录/顺序借值。
    for sclk,wt,has_rwp in [('406545731',18.39,True),('406545732',31.25,True),('406545733',0,False)]:
        for p in source.glob('*'):
            if p.suffix not in ('.csv','.lbl') or ('rwp' in p.name and not has_rwp):continue
            text=p.read_text(encoding='utf-8')
            if 'rwp' in p.name and p.suffix=='.csv':text=text.replace('18.39','%.2f'%wt)
            (folder/p.name.replace('406545731',sclk)).write_text(text,encoding='utf-8')
    old=next(folder.glob('*406545732rwp*.csv'))
    (folder/old.name.replace('_p1','_p2')).write_text(old.read_text().replace('31.25','32.50'))
    old_label=old.with_suffix('.lbl')
    (folder/old_label.name.replace('_p1','_p2')).write_bytes(old_label.read_bytes())
    rsp=next(folder.glob('*406545732rsp*.csv'))
    (folder/rsp.name.replace('_p1','_p2')).write_bytes(rsp.read_bytes())
    old_label=rsp.with_suffix('.lbl')
    (folder/old_label.name.replace('_p1','_p2')).write_bytes(old_label.read_bytes())
    b=APXSBackend();obs=b.scan_folder(str(folder))[0]
    check(len(obs.points)==3,'同点多版本RSP只保留一个谱点')
    pts={pt.sclk:pt for pt in obs.points}
    def fe(pt):
        ox=b.oxides(obs,pt)
        return None if not ox else next(r['wt'] for r in ox['rows'] if r['oxide']=='FeO')
    check(fe(pts['406545731'])==18.39,'首点成分正确')
    check(fe(pts['406545732'])==32.5,'第二点取自身最高版本成分')
    check(fe(pts['406545733']) is None,'缺成分表不可借用同序列其他点的值')
    check(b.oxides(obs) is None,'未指定点的多点观测不可给出歧义成分')
    payload=console.plot_for({'observations':[obs]},b,0,1,None,{})
    check(next(r['wt'] for r in payload['oxides']['rows'] if r['oxide']=='FeO')==32.5,'API返回当前点成分')
    missing=console.plot_for({'observations':[obs]},b,0,2,None,{})
    check('oxides' not in missing,'API缺表时不残留其他点成分')
    out=temp/'report';result=build_report(b,[obs],str(out),steps_on={})
    with (out/result['csv'][0]).open(encoding='utf-8-sig',newline='') as f:
        rows=list(csv.DictReader(f))
    values={r['SCLK']:r['APXS FeO wt%'] for r in rows}
    check(values=={'406545731':'18.39','406545732':'32.5','406545733':''},'报告逐点数值与缺值必须保持配对')
    html=(out/'index.html').read_text(encoding='utf-8')
    check('406545731rwp' in html and '406545732rwp' in html,'HTML同时展示两份来源成分表')
    check('_p2.csv' in html,'HTML成分来源标记最新版本')


def main():
    with tempfile.TemporaryDirectory(prefix='sherlock-integrity-') as name:
        temp=Path(name);scaling(temp);paired_apxs(temp)
    print('共 %d 项断言，%d 项失败'%(TOTAL,FAIL))
    return int(FAIL>0)


if __name__=='__main__':raise SystemExit(main())
