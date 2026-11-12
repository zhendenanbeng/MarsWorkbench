"""Offline real-product contracts for PSV wavelength pairing and the failed RMI gate."""

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout
import csv
import hashlib
import json
from pathlib import Path
import re
import sys
import tempfile
from unittest.mock import patch
import numpy as np
from PIL import Image

ROOT=_layout.project_layout.APP;sys.path.insert(0,str(ROOT))
from instruments.chemcam import ChemCamBackend
from instruments import chemcam_passive as passive,report
import console

FIX=_layout.project_layout.VERIFY/'fixtures/chemcam'
TOTAL=FAIL=0
def check(ok,name):
    global TOTAL,FAIL
    TOTAL+=1
    if not ok:FAIL+=1;print(' FAIL '+name)
def rejected(fn,name):
    try:fn()
    except ValueError:check(True,name)
    else:check(False,name)

def main():
    b=ChemCamBackend();obs=b.scan_folder(str(FIX));active=[o for o in obs if not getattr(o,'passive',False)]
    check(len(active)==1 and set(active[0].points[0].files)=={'ccs','rdr'},'主动CCS/RDR仍合并为同一点')
    spectra=[o for o in obs if getattr(o,'passive',False)]
    check({o.product for o in spectra}=={'CL1-PSV','CL9-PSV'},'扫描必须纳入CL1/CL9 TAB，产品独立')
    grid=(FIX/'document/ccam_default_wave.tab').read_bytes()
    check(hashlib.sha256(grid).hexdigest()==passive.WAVE_SHA256,'波长表必须保持冻结原件')
    wave=np.array([float(line) for line in grid.splitlines()])
    valid=wave>0;order=np.argsort(wave[valid],kind='stable')
    check(len(wave)==6444 and valid.sum()==6144,'区分6444原始通道和6144有效通道')
    check(np.all(wave[~valid]==-999.9999999),'300个无效边缘波长哨兵值须可解释')
    for o in spectra:
        pt=o.points[0];path=Path(pt.files['psv']);label=path.with_suffix('.lbl')
        raw=np.array([float(line) for line in path.read_bytes().splitlines()[29:]])
        d=b.load_point(pt,9,{'dust':True,'despike':True})
        check(np.array_equal(d['x'],wave[valid][order]),o.product+'波长轴须按原始通道掩码成对排序')
        check(np.array_equal(d['y'],raw[valid][order]),o.product+'强度不得截错300行或调换谱段')
        check(d['n_shots']==0 and not d['shot_numbers'] and d['passive'],o.product+'没有激光发数')
        check(np.all(np.diff(d['x'])>0) and 240<d['x'][0]<250 and 900<d['x'][-1]<910,o.product+'排序后覆盖官方UV/VIS/VNIR范围')
        check(d['ylabel']=='平均强度 (DN)',o.product+'不得把DN冒充反射率或辐射亮度')
        plain=b.load_point(pt,None,{})
        check(np.array_equal(plain['y'],d['y']),o.product+'除尘和跨发去尖峰选项不得更改被动谱')
        with patch.object(b,'annotate',side_effect=AssertionError('passive annotation forbidden')):
            payload=console.plot_for({'observations':[o]},b,0,0,3,{'dust':True})
        check(not payload['annotate']['marks'] and not payload['annotate']['elements'],o.product+'不得执行元素发射线判定')
        check(payload['shots']==[{'i':-1,'label':'全谱'}],o.product+'界面只提供全谱且标记被动模式')
    with tempfile.TemporaryDirectory(prefix='chemcam-passive-') as tmp:
        dest=Path(tmp);o=spectra[0];pt=o.points[0];path=Path(pt.files['psv']);label=path.with_suffix('.lbl')
        bad=dest/'wrong.tab';bad.write_bytes(grid+b'0\n')
        rejected(lambda:passive.ensure_wavelengths(bad),'错误缓存须拒绝，不能回退到猜测波长')
        rejected(lambda:passive.read_spectrum(path,dest/'absent.lbl',FIX/'document/ccam_default_wave.tab'),'缺少LBL须明确拒绝')
        for old,new,name in [('ROWS                            = 6444','ROWS                            = 6144','错误行数'),
                             ('COLUMNS                         = 1','COLUMNS                         = 2','多列PSV'),
                             ('DIGITAL_NUMBER','REFLECTANCE','未支持强度单位'),
                             (',30)',',29)','错误表偏移')]:
            text=label.read_text();check(old in text,'负例修改必须确实命中 '+name)
            badlabel=dest/'bad.lbl';badlabel.write_text(text.replace(old,new),encoding='ascii')
            rejected(lambda:passive.read_spectrum(path,badlabel,FIX/'document/ccam_default_wave.tab'),'拒绝 '+name)
        target=dest/path.name;target.write_bytes(path.read_bytes()[:-44]);labelcopy=target.with_suffix('.lbl');labelcopy.write_bytes(label.read_bytes())
        rejected(lambda:passive.read_spectrum(target,labelcopy,FIX/'document/ccam_default_wave.tab'),'截断TAB须拒绝')
        with patch.object(b,'annotate',side_effect=AssertionError('passive annotation forbidden')):
            report.build_report(b,spectra,str(dest/'qc'))
            console.render_point_png(b,o,0,None,{},str(dest/'passive.png'))
        with Image.open(dest/'passive.png') as im:check(im.width>100 and im.height>100,'真实PSV可以导出谱图PNG')
        page=(dest/'qc/index.html').read_text(encoding='utf-8')
        check('被动谱不执行元素发射线识别' in page and '平均强度 (DN)' in page,'报告须保留被动模式与正确单位')
        check(not list((dest/'qc').glob('*_lines.csv')),'被动报告不能生成元素发射线证据表')
        rows=[row for p in (dest/'qc').glob('*_points.csv') for row in csv.DictReader(p.open(encoding='utf-8-sig'))]
        check(all(not row['发数'] and not row['检出元素'] for row in rows),'被动CSV不报告发数和元素检出')
        # 同一sol/序列的主动、被动仍是两个观测，不能按SCLK合并跨模式数据。
        mixed_dir=dest/'sol00013';mixed_dir.mkdir()
        # 主动文件故意缺LBL而使用目录sol；被动文件从LBL取sol，须统一键类型。
        for p in (FIX/'sol00013').glob('*.csv'): (mixed_dir/p.name).write_bytes(p.read_bytes())
        for p in (FIX/'sol01000').glob('*psv*'): (mixed_dir/p.name).write_bytes(p.read_bytes())
        for p in mixed_dir.glob('*psv*.lbl'):
            text=re.sub(r'(PLANET_DAY_NUMBER\s*=\s*)1000',r'\g<1>13',p.read_text())
            p.write_text(text.replace('ccam02000','ccam02013').replace('CCAM02000','CCAM02013'))
        # Explicit scan grouping test uses filename and label sequence consistently.
        mixed=b.scan_folder(str(mixed_dir))
        matching=[o for o in mixed if o.sol=='13' and o.seq=='ccam02013']
        check({o.product for o in matching}=={'RDR','CL9-PSV'},'同一sol/序列须独立保留主动和被动观测')
    entries=[{'name':'cl9_486263266psv_f0481570ccam02000p3.'+ext,'size':1} for ext in ['tab','lbl']]
    entries += [{'name':'cr0_486263086prc_f0481570ccam02000l1.tif','size':1}]
    listing=b.parse_listing(entries)
    check(sum(g['count'] for g in listing['groups'] if g['key']=='CL9_PSV')==2,'下载面板必须成对列出TAB和LBL')
    check(not listing['preset'],'PSV和RMI仍为可选下载，不能自动加入默认预选')
    labels=sorted((FIX/'sol01000').glob('cl5_*ccs*.lbl'))
    coordinates=[tuple(int(re.search(r'MSL:'+key+r'\s*=\s*(\d+)',p.read_text()).group(1)) for key in ['SPECIAL_LINE','SPECIAL_SAMPLE']) for p in labels]
    check(len(coordinates)==9 and set(coordinates)=={(497,532)},'真实9点均为同一相机内光斑位置，必须阻止免配准叠图假设')
    print('共 %d 项断言，%d 项失败'%(TOTAL,FAIL))
    return int(FAIL>0)

if __name__=='__main__':raise SystemExit(main())
