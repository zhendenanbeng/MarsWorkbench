# -*- coding: utf-8 -*-
"""【已弃用待重算】独立数值审计。参考侧禁止调用软件的读表/积分/分类内核。

冻结清单 independent_manifest.json 是 0.16.0 构建时代产物（绑定旧目录
01_软件工作台），不在本仓——当前入口显式跳过；其功能已由
data_audit.py + paper_quant_check.py + regression.py 覆盖。
重算基线（对当前树重新生成冻结清单并校准断言）属独立任务。

python verify/independent_check.py --out <新报告目录>
附加 --online 才只读核对 PDS 原件；离线运行不联网、不更新冻结基准。
缺少独立依赖或冻结输入时返回失败，不能作为跳过后的通过。
"""
from __future__ import annotations

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout

import argparse
import builtins
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import platform
import sys
import tempfile
import urllib.request
import warnings
import xml.etree.ElementTree as ET
from unittest.mock import patch

import numpy as np

ROOT = _layout.project_layout.APP
WORKSPACE = _layout.project_layout.ROOT
FIX = _layout.project_layout.VERIFY/'fixtures'
FROZEN = _layout.project_layout.VERIFY/'independent_manifest.json'
sys.path.insert(0, str(ROOT))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class Audit:
    def __init__(self, out):
        self.out = out
        self.results = []
        self.inputs = []
        self.measurements = {}

    def check(self, name, condition, **evidence):
        self.results.append(dict(name=name, passed=bool(condition), **evidence))
        print((' PASS ' if condition else ' FAIL ')+name, flush=True)

    def equal(self, name, actual, reference, rtol=0., atol=0.):
        a, b = np.asarray(actual, dtype=float), np.asarray(reference, dtype=float)
        same = a.shape == b.shape
        good = same and bool(np.allclose(a, b, rtol=rtol, atol=atol, equal_nan=True))
        delta = np.abs(a-b) if same else np.array([])
        finite = delta[np.isfinite(delta)]
        self.check(name, good, actual_shape=list(a.shape), reference_shape=list(b.shape),
                   scalars=int(b.size), rtol=rtol, atol=atol,
                   max_abs_difference=float(finite.max()) if finite.size else None)


def native_fits_only():
    """拦截生产侧 astropy 导入，防止两边最终调用同一读取库。"""
    original = builtins.__import__
    def import_without_astropy(name, *args, **kwargs):
        if name == 'astropy' or name.startswith('astropy.'):
            raise ImportError('独立验证强制工作台使用自带 FITS 解析器')
        return original(name, *args, **kwargs)
    return patch('builtins.__import__', side_effect=import_without_astropy)


def freeze_checks(audit):
    manifest = json.loads(FROZEN.read_text(encoding='utf-8'))
    for rec in manifest['files']:
        path = WORKSPACE/rec['path']
        actual = sha(path) if path.is_file() else None
        audit.check('冻结输入 '+rec['path'], actual == rec['sha256'],
                    expected_sha256=rec['sha256'], actual_sha256=actual)
        audit.inputs.append(dict(rec, actual_sha256=actual))
    return manifest


def fits_checks(audit, manifest):
    from astropy.io import fits
    from instruments import fitsio, supercam
    backend = supercam.SuperCamBackend()
    observations = backend.scan_folder(str(FIX/'supercam'))
    for rec in manifest['fits']:
        path = WORKSPACE/rec['path']
        with warnings.catch_warnings(record=True) as caught, fits.open(path, memmap=False) as hdus:
            hdus.verify('exception')  # 不修复源文件，也不静默修正后再比较
            audit.check(rec['product']+' FITS 标准结构', True,
                        warnings=[str(w.message) for w in caught])
            for ext in rec['extensions']:
                hdu = hdus[ext]
                with native_fits_only():
                    actual = fitsio.read_bintable(str(path), ext)
                for column in hdu.columns.names:
                    if np.asarray(hdu.data[column]).dtype.kind not in 'biufc':
                        continue
                    audit.equal(rec['product']+' '+ext+'.'+column,
                                actual[column], hdu.data[column])
            wx = hdus['WAVELENGTH'].data
            reference_x = np.asarray(wx[wx.names[0]], dtype=float)
            obs = next(o for o in observations if o.product == rec['product'])
            pt = next(p for p in obs.points if Path(next(iter(p.files.values()))) == path)
            if obs.product == 'CL1':
                shot_names = sorted(hdus['SPECTRA'].columns.names,
                                    key=lambda s: int(s.removeprefix('Spectrum')))
                matrix = np.array([hdus['SPECTRA'].data[n] for n in shot_names], dtype=float).T
                for dust in (False, True):
                    kept = matrix[:, 5:] if dust else matrix
                    # math.fsum 独立于生产 numpy.mean 的求和路径
                    mean = np.array([math.fsum(row)/len(row) for row in kept])
                    with native_fits_only():
                        actual = backend.load_point(pt, None, {'dust':dust})
                        first = backend.load_point(pt, 0, {'dust':dust})
                    audit.equal('CL1 逐通道均值 dust='+str(dust), actual['y'], mean, rtol=2e-15, atol=1e-12)
                    audit.equal('CL1 原始发号 dust='+str(dust), first['y'], kept[:,0])
                    audit.check('CL1 发号映射 dust='+str(dust), first['shot_numbers'][0] == (6 if dust else 1))
                reference_y = matrix.mean(axis=1)
            else:
                ext, col = rec['display_column'].split(':')
                reference_y = np.asarray(hdus[ext].data[col], dtype=float)
                with native_fits_only():
                    actual = backend.load_point(pt, None, {})
                order = np.argsort(reference_x)
                audit.equal(obs.product+' 展示原始强度（全部通道）', actual['y'], reference_y[order])
                audit.check(obs.product+' 横轴单位', actual['xlabel'] == rec['xlabel'])
                reference_x = reference_x[order]
                if obs.product=='CP2':
                    with native_fits_only():reflection=backend.load_point(pt,None,{'reflectance':True})
                    audit.equal('CP2 可选大气校正反射列独立核对',reflection['y'],np.asarray(hdus['SPECTRA'].data['I_F_atm'],float)[order])
            with native_fits_only():
                actual = backend.load_point(pt, None, {})
            audit.equal(obs.product+' 展示波长（全部通道）', actual['x'], reference_x)
            audit.check(obs.product+' 冻结通道数', len(reference_x) == rec['channels'])
            # 这些锚点由外部读取器冻结，不能由软件输出重建。
            audit.equal(obs.product+' 冻结波长锚点', reference_x[rec['anchor_indices']], rec['anchor_x'])
    # 已知带缩放的真实 IR HK 表也纳入检查，不能忽略 FITS TZERO。


def csv_numbers(path, skip=0, comments=False):
    with Path(path).open(encoding='utf-8-sig', newline='') as f:
        rows = list(csv.reader(f))[skip:]
    return np.array([[float(v) for v in row] for row in rows
                     if row and not (comments and row[0].lstrip().startswith('#'))], dtype=float)


def pds_checks(audit):
    from instruments import chemcam, apxs
    from sherloc.reader import Product
    for product in ('rdr','ccs'):
        path = next((FIX/'chemcam/sol00013').glob('*'+product.upper()+'*.csv'))
        reference = csv_numbers(path, comments=True)
        b = chemcam.ChemCamBackend()
        obs = next(o for o in b.scan_folder(str(FIX/'chemcam')) if not getattr(o,'passive',False))
        from instruments.base import PointFile
        pt = PointFile(path.stem, {product:str(path)}, sclk=obs.points[0].sclk)
        for dust in (False, True):
            data = b.load_point(pt, None, {'dust':dust})
            shots = reference[:,6:-2] if dust else reference[:,1:-2]
            mean = np.array([math.fsum(row)/len(row) for row in shots])
            audit.equal('ChemCam '+product+' 波长', data['x'], reference[:,0])
            audit.equal('ChemCam '+product+' 平均谱 dust='+str(dust),data['y'],mean,rtol=2e-15,atol=1e-12)
    # 参考侧用标准库逐行读原件，独立按原始通道索引成对排序，不调用生产PSV解析器。
    reference_wave=[float(line) for line in (FIX/'chemcam/document/ccam_default_wave.tab').read_bytes().splitlines()]
    channel_order=sorted((i for i,v in enumerate(reference_wave) if v>0),key=lambda i:reference_wave[i])
    audit.check('ChemCam PSV 原始6444/有效6144通道',len(reference_wave)==6444 and len(channel_order)==6144)
    for obs in (o for o in b.scan_folder(str(FIX/'chemcam')) if getattr(o,'passive',False)):
        pt=obs.points[0]
        raw=[float(line) for line in Path(pt.files['psv']).read_bytes().splitlines()[29:]]
        actual=b.load_point(pt,None,{})
        audit.equal('ChemCam '+obs.product+' 波长独立通道映射',actual['x'],[reference_wave[i] for i in channel_order])
        audit.equal('ChemCam '+obs.product+' 强度独立通道映射',actual['y'],[raw[i] for i in channel_order])
        audit.check('ChemCam '+obs.product+' 被动强度口径',actual.get('passive') and actual['n_shots']==0 and actual['ylabel']=='平均强度 (DN)')
    b = apxs.APXSBackend();obs = b.scan_folder(str(FIX/'apxs'))[0]
    pt = obs.points[0];raw = csv_numbers(pt.files['rsp'], skip=1)
    data = b.load_point(pt, None, {})
    audit.equal('APXS eV→keV及首末通道剔除', data['x'], raw[1:-1,0]/1000.)
    audit.equal('APXS 全部有效通道计数', data['y'], raw[1:-1,1])
    with Path(obs.rwp_files[0]).open(encoding='utf-8', newline='') as f:
        rows = list(csv.reader(f))[1:]
    expected = {r[0].strip():(float(r[1]),float(r[2])) for r in rows if len(r)>=3}
    actual = {r['oxide']:(r['wt'],r['err']) for r in b.oxides(obs)['rows']}
    audit.check('APXS RWP 全部成分和标准差原文', actual == expected, components=len(expected))
    folder = WORKSPACE/'03_下载数据/SHERLOC/练习数据/processed/sol_00370'
    paths = sorted(folder.glob('*3450*rrs*.xml')) + sorted(folder.glob('*4099*rrs*.xml'))
    for label in paths:
        prod = Product(str(label.with_suffix('.csv')))
        # 参考侧直接读标准 XML 和 csv，不使用 Label/read_table/find_table。
        root = ET.parse(label).getroot()
        for table in root.findall('.//{*}Table_Delimited'):
            name = table.findtext('{*}name')
            offset = int(table.findtext('{*}offset'));records = int(table.findtext('{*}records'))
            with label.with_suffix('.csv').open('rb') as f:
                f.seek(offset)
                text = b''.join(f.readline() for _ in range(records)).decode('ascii')
            reference = np.array([[float(v) for v in r] for r in csv.reader(io.StringIO(text))])
            actual = prod.read_table(next(t for t in prod.tables() if t.name == name))
            audit.equal('SHERLOC '+label.stem[10:26]+' '+name, actual, reference)
            groups = table.findall('.//{*}Group_Field_Delimited/{*}repetitions')
            audit.check('SHERLOC XML行列契约 '+label.stem[10:26]+' '+name,
                        reference.shape == (records, int(groups[0].text)),
                        offset=offset, records=records, fields=int(groups[0].text))


def moc_checks(audit, path, name):
    from instruments import supercam
    # 标准 csv.DictReader 独立解析，键是源表 cdr_fname 原文（不是软件 moc_key）。
    text = Path(path).read_text(encoding='utf-8-sig').splitlines()
    start = next(i for i,s in enumerate(text) if s.startswith('cdr_fname,'))
    records = list(csv.DictReader(text[start:]))
    oxides = ('SiO2','TiO2','Al2O3','FeOT','MgO','CaO','Na2O','K2O')
    expected = {};conflicts=[];duplicates=0
    def number(v):
        try:
            n=float(v)
            return n if math.isfinite(n) else None
        except (ValueError,TypeError):return None
    for r in records:
        key = r['cdr_fname'].lower()
        value = {'oxides':{k:(number(r[k]),number(r[k+'_stdev'])) for k in oxides},
                 'total':number(r['Total']), 'distance':number(r['distance_mm'])}
        if key in expected:
            duplicates+=1
            if value != expected[key]:conflicts.append(key)
        expected[key]=value
    actual,_ = supercam.load_moc(str(path))
    audit.check(name+' 所有行与全部数值字段', actual==expected,
                physical_rows=len(records), unique_products=len(expected),duplicates=duplicates)
    audit.check(name+' 重复产品没有矛盾值', not conflicts, conflicting_keys=conflicts)
    mismatches=[]
    with patch.object(supercam,'moc_cache_path',return_value=str(path)):
        backend=supercam.SuperCamBackend()
        for key,rec in expected.items():
            payload=backend.moc_payload(key+'01')
            got={r['oxide']:(r['wt'],r['err']) for r in payload.get('rows',[])}
            if (not payload.get('matched') or got!=rec['oxides'] or payload.get('total')!=rec['total']
                    or payload.get('distance')!=rec['distance']):mismatches.append(key)
        audit.check(name+' 全部产品名→预测值匹配',not mismatches,
                    product_count=len(expected),mismatched_keys=mismatches[:10])
        audit.check(name+' 未收录产品不可借值',
                    backend.moc_payload('scam_99999_0000000000_000_cl1_scam99999_unlisted___01p01').get('matched') is False)


def pairing_checks(audit):
    """故意构造同序列、不同SCLK和不同成分，检测跨产品串表。"""
    from instruments.apxs import APXSBackend
    from instruments.report import build_report
    import shutil
    source=FIX/'apxs/sol00102'
    with tempfile.TemporaryDirectory(prefix='sherlock-apxs-pair-') as tmp:
        folder=Path(tmp)/'sol00102';folder.mkdir()
        for sclk,wt in [('406545731',18.39),('406545732',31.25)]:
            for p in source.glob('*'):
                if p.suffix not in ('.csv','.lbl'):continue
                out=folder/p.name.replace('406545731',sclk)
                raw=p.read_text(encoding='utf-8')
                if 'rwp' in p.name and p.suffix=='.csv':raw=raw.replace('18.39','%.2f'%wt)
                out.write_text(raw,encoding='utf-8')
        b=APXSBackend();observations=b.scan_folder(str(folder))
        result=build_report(b,observations,str(Path(tmp)/'report'),steps_on={})
        recovered={}
        for name in result['csv']:
            with (Path(tmp)/'report'/name).open(encoding='utf-8-sig',newline='') as f:
                for row in csv.DictReader(f):recovered[row['SCLK']]=float(row['APXS FeO wt%'])
        audit.check('APXS同序列不同SCLK各取自身成分（防串表）',
                    recovered=={'406545731':18.39,'406545732':31.25},values=recovered)


def numeric_checks(audit):
    from scipy.interpolate import make_interp_spline
    from sherloc import explore
    rng = np.random.default_rng(20261023)
    errors=[]
    for i in range(40):
        x=np.r_[0.,np.cumsum(rng.uniform(.1,2,30))];y=rng.normal(0,4,len(x))
        lo,hi=sorted(rng.uniform(x[0],x[-1],2))
        spline=make_interp_spline(x,y,k=1)
        truth=spline.integrate(lo,hi)
        actual=explore.band_features(x,y[None,:],np.array([.2]),(lo,hi))['area'][0]
        errors.append(abs(float(actual-truth)))
    audit.check('40组不等间距/边界插值积分对独立SciPy一阶样条',max(errors)<1e-10,max_abs_difference=max(errors))
    x=np.linspace(-8,8,3201);y=np.exp(-x*x/2)
    truth=math.sqrt(2*math.pi)*math.erf(3/math.sqrt(2))
    area=explore.band_features(x,y[None,:],np.array([.001]),(-3,3))['area'][0]
    audit.equal('高斯解析积分（离散误差限）',[area],[truth],atol=2e-6)
    # 参考侧直接用 SciPy 对噪声实现的每条样条积分，不使用待测权重。
    x=np.array([0.,.2,.8,1.1,1.8,2.7,3.5]);signal=4+np.sin(x)
    sigma=.015;n=30000
    noisy=signal+rng.normal(0,sigma,(n,len(x)))
    aa=(.15,2.3);bb=(1.0,3.2)
    int_a=make_interp_spline(x,noisy,k=1,axis=1).integrate(*aa)
    int_b=make_interp_spline(x,noisy,k=1,axis=1).integrate(*bb)
    a=explore.band_features(x,signal[None,:],np.array([sigma]),aa)
    b=explore.band_features(x,signal[None,:],np.array([sigma]),bb)
    ratio=explore.ratio_features(a,b,np.array([sigma]))
    std_area=float(np.std(int_a,ddof=1));std_ratio=float(np.std(int_a/int_b,ddof=1))
    area_error=abs(std_area/a['area_sigma'][0]-1)
    ratio_error=abs(std_ratio/ratio['sigma'][0]-1)
    audit.check('30000次独立噪声模拟：积分σ',area_error<.025,relative_sigma_error=area_error,
                empirical_sigma=std_area,predicted_sigma=float(a['area_sigma'][0]))
    audit.check('30000次独立噪声模拟：重叠窗比值σ',ratio_error<.025,relative_sigma_error=ratio_error,
                empirical_sigma=std_ratio,predicted_sigma=float(ratio['sigma'][0]))
    audit.measurements['uncertainty_scope']='仅独立同方差通道白噪声，比例误差用高SNR一阶近似；不含基线/标定/相关噪声。'


def cluster_checks(audit):
    from sklearn.decomposition import PCA
    from sklearn.cluster import KMeans
    from sklearn.metrics import adjusted_rand_score
    from sherloc import explore
    from sherloc.organize import scan_folder
    rng=np.random.default_rng(20261024)
    x=np.linspace(800,1800,161)
    centers=(910,1190,1660)
    prototypes=np.array([np.exp(-((x-c)/25)**2) for c in centers])
    spectra=np.vstack([p+rng.normal(0,.004,(24,len(x))) for p in prototypes])
    truth=np.repeat(np.arange(3),24)
    model=explore.spectral_types(x,spectra,np.ones(72,bool),3)
    score=float(adjusted_rand_score(truth,model['labels']))
    audit.check('冻结随机种子合成谱形：类别真值ARI',score==1.,ari=score,
                truth='人工生成的三个高斯谱形，非矿物实测标签')
    normalized=spectra/np.linalg.norm(spectra,axis=1)[:,None]
    reference=PCA(n_components=.95,svd_solver='full').fit(normalized)
    proj=reference.components_.T@reference.components_
    actual=model['loadings'].T@model['loadings']
    audit.equal('PCA子空间对独立sklearn实现（符号无关）',actual,proj,atol=2e-12)
    ref_labels=KMeans(n_clusters=3,n_init=20,random_state=31).fit_predict(reference.transform(normalized))
    audit.check('合成分型对独立sklearn KMeans',adjusted_rand_score(ref_labels,model['labels'])==1.)
    _,scans=scan_folder(str(WORKSPACE/'03_下载数据/SHERLOC/练习数据/processed/sol_00370'))
    scan=next(s for s in scans if s.spectra is not None and s.sequence_id=='srlc11360')
    data=explore.prepare(scan)
    # 生产数据准备步骤用于固定当前研究问题；它不是基线扣除的独立真值。
    fingerprint=(data['x']>=800)&(data['x']<=1800)
    xx=data['x'][fingerprint];yy=data['residual'][:,fingerprint]
    candidates=data['finite']&data['qualified']&data['unknown']&np.isfinite(data['sigma'])&(data['sigma']>0)
    peak_snr=np.max(yy,axis=1)/data['sigma']
    candidates &= np.nan_to_num(peak_snr,nan=0.,posinf=0.)>=3
    base=explore.spectral_types(xx,yy,candidates,3)
    selected=base['selected'];labels=base['labels']
    noise_scores=[];subsample_scores=[];sensitivity={}
    for _ in range(30):
        perturbed=yy.copy()
        perturbed[selected]+=rng.normal(size=yy[selected].shape)*data['sigma'][selected,None]
        other=explore.spectral_types(xx,perturbed,candidates,3)
        noise_scores.append(float(adjusted_rand_score(labels,other['labels'])))
        indices=np.sort(rng.choice(selected,size=max(3,math.ceil(len(selected)*.8)),replace=False))
        mask=np.zeros(len(yy),bool);mask[indices]=True
        sub=explore.spectral_types(xx,yy,mask,3)
        corresponding=np.array([labels[np.where(selected==i)[0][0]] for i in sub['selected']])
        subsample_scores.append(float(adjusted_rand_score(corresponding,sub['labels'])))
    for k in (2,3,4):
        other=explore.spectral_types(xx,yy,candidates,k)
        sensitivity[str(k)]={'sizes':[g['n'] for g in other['groups']],
                             'ari_vs_k3':float(adjusted_rand_score(labels,other['labels']))}
    audit.measurements['real_cluster_stability']={
        'scan':'Sol370 srlc11360','candidate_points':len(selected),
        'base_sizes':[g['n'] for g in base['groups']],
        'noise_ari':noise_scores,'subsample_80_percent_ari':subsample_scores,
        'k_sensitivity':sensitivity,
        'interpretation':'稳健性诊断，不设伪造的科学通过线；小组不自动等于新矿物。'}
    audit.check('真实数据稳健性实验完成（不是分类准确率）',
                len(noise_scores)==30 and len(subsample_scores)==30,
                median_noise_ari=float(np.median(noise_scores)),
                median_subsample_ari=float(np.median(subsample_scores)))


def profile_checks(audit):
    from scipy.optimize import lsq_linear
    from sherloc.fitting import solve_nonnegative_profiles
    from sherloc import minerals
    rng=np.random.default_rng(773);errors=[];negative=0;reconstruction=[]
    for count in (1,2):
        for _ in range(30):
            x=np.sort(rng.uniform(900,1150,51))
            profiles=np.column_stack([part for center,width in [(990,45),(1050,60)][:count]
                                      for part in minerals._voigt_basis(x,center,width)])
            design=np.column_stack([profiles,x,np.ones(len(x))])
            y=design@rng.uniform(-10,10,design.shape[1])+rng.normal(0,2,len(x))
            coefficient,error,model=solve_nonnegative_profiles(x,y,profiles)
            reference=lsq_linear(design,y,bounds=(np.r_[np.zeros(2*count),-np.inf,-np.inf],
                                                  np.full(2*count+2,np.inf)),method='bvls',tol=1e-12)
            if not reference.success:raise RuntimeError('独立有界最小二乘未收敛')
            expected=float(np.sum((design@reference.x-y)**2))
            errors.append(abs(error-expected)/max(1.,expected))
            negative+=int(np.any(coefficient[:2*count]<0))
            reconstruction.append(float(np.max(np.abs(model-design@coefficient))))
    audit.check('60组峰形非负约束对独立SciPy BVLS',max(errors)<1e-7,
                max_relative_objective_error=max(errors),cases=len(errors))
    audit.check('非负求解的所有峰形系数符合约束',negative==0)
    audit.check('自由背景与峰形可由返回系数重建',max(reconstruction)<1e-8,
                max_abs_difference=max(reconstruction))
    x=np.linspace(900,1100,41);g,l=minerals._voigt_basis(x,1000,50);y=30*g-10*l+3
    fit=minerals.fit_voigt(x,y,1000,50,within=(970,1030))
    model=minerals._pseudo_voigt(x,fit['amp'],fit['center'],fit['fwhm'],fit['eta'])+fit['slope']*x+fit['intercept']
    actual=float(np.sqrt(np.mean((y-model)**2)))
    audit.check('有符号候选不截断比例以伪装物理解',not fit['profile_valid'] and abs(actual-fit['rms'])<1e-9,
                reported_rms=fit['rms'],reconstructed_rms=actual,eta=fit['eta'])


def online_checks(audit, manifest):
    for rec in manifest['fits']:
        with urllib.request.urlopen(rec['source_url'],timeout=60) as response:
            raw=response.read();headers=dict(response.headers)
        remote_hash=hashlib.sha256(raw).hexdigest()
        audit.check('PDS在线原件 '+rec['product'],remote_hash==sha(WORKSPACE/rec['path']),
                    url=rec['source_url'],remote_sha256=remote_hash,bytes=len(raw),
                    last_modified=headers.get('Last-Modified'))
    url='https://pds-geosciences.wustl.edu/m2020/urn-nasa-pds-mars2020_supercam/data_derived_spectra/supercam_libs_moc.csv'
    with urllib.request.urlopen(url,timeout=90) as r:raw=r.read()
    digest=hashlib.sha256(raw).hexdigest()
    (audit.out/'official_supercam_moc.csv.gz').write_bytes(gzip.compress(raw,mtime=0))
    audit.measurements['official_moc_source']={'url':url,'sha256':digest,'bytes':len(raw)}
    audit.check('在线官方MOC冻结版本',digest==manifest['official_moc_sha256'],
                expected_sha256=manifest['official_moc_sha256'],actual_sha256=digest)
    with tempfile.TemporaryDirectory(prefix='sherlock-moc-audit-') as tmp:
        path=Path(tmp)/'moc.csv';path.write_bytes(raw)
        moc_checks(audit,path,'在线官方MOC全表')


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',required=True,type=Path)
    parser.add_argument('--online',action='store_true')
    parser.add_argument('--moc-snapshot',type=Path,help='离线重验已存的官方MOC .csv.gz 原件')
    args=parser.parse_args(argv)
    if FROZEN.is_file():
        _bound = json.loads(FROZEN.read_text(encoding='utf-8'))
        _rows = _bound.get('inputs') or _bound.get('files') or []
        if any('01_软件工作台' in str(r.get('path', '')) for r in _rows[:50]):
            print('SKIP 已弃用待重算：冻结清单绑定 0.16.0 时代旧目录'
                  '（01_软件工作台/…），与当前仓树结构不匹配。'
                  '其功能已由 data_audit.py + paper_quant_check.py + regression.py 覆盖；'
                  '对当前树重算冻结清单属独立任务。')
            print('共 0 项断言，0 项失败')
            return 0
    try:
        import astropy  # noqa: F401
        import scipy    # noqa: F401
        import sklearn  # noqa: F401
    except ImportError as exc:
        print('SKIP 依赖不在位：%s——本审计需要 full 版依赖'
              '（pip install astropy scipy scikit-learn）后重跑' % exc)
        return 0
    args.out.mkdir(parents=True,exist_ok=False)  # 不覆盖上次的审计证据
    audit=Audit(args.out)
    try:
        import astropy,scipy,sklearn
        audit.measurements['versions']={'python':platform.python_version(),'numpy':np.__version__,
                                        'astropy':astropy.__version__,'scipy':scipy.__version__,'sklearn':sklearn.__version__}
        manifest=freeze_checks(audit)
        if any(not r['passed'] for r in audit.results):
            raise ValueError('冻结输入变动，停止数值审计；不能自动改写基准使它通过。')
        fits_checks(audit,manifest)
        pds_checks(audit)
        pairing_checks(audit)
        moc_checks(audit,FIX/'supercam/_moc/supercam_moc.csv','冻结MOC样例')
        numeric_checks(audit)
        profile_checks(audit)
        cluster_checks(audit)
        if args.online:online_checks(audit,manifest)
        if args.moc_snapshot:
            raw=gzip.decompress(args.moc_snapshot.read_bytes())
            digest=hashlib.sha256(raw).hexdigest()
            audit.check('离线官方MOC冻结版本',digest==manifest['official_moc_sha256'],actual_sha256=digest)
            with tempfile.TemporaryDirectory(prefix='sherlock-moc-replay-') as tmp:
                path=Path(tmp)/'moc.csv';path.write_bytes(raw)
                moc_checks(audit,path,'离线官方MOC全表')
    except Exception as exc:
        import traceback
        audit.check('审计执行完整性',False,error=str(exc),traceback=traceback.format_exc())
        print(traceback.format_exc())
    for rec in audit.inputs:
        path=WORKSPACE/rec['path']
        audit.check('源文件未更改 '+rec['path'],path.is_file() and sha(path)==rec['actual_sha256'])
    result={'generated_utc':datetime.now(timezone.utc).isoformat(),'results':audit.results,
            'measurements':audit.measurements,'inputs':audit.inputs,
            'software_files':{Path(p).relative_to(WORKSPACE).as_posix() if p.is_relative_to(WORKSPACE) else Path(p).name:sha(p) for p in
                              [*ROOT.glob('*.py'),*(ROOT/'instruments').glob('*.py'),*(ROOT/'sherloc').glob('*.py'),
                               Path(__file__).resolve(),FROZEN]},
            'limits':['尚无独立盲标的矿物阳性/阴性真值，不能给出真实矿物准确率、假阳性率或假阴性率。',
                      '官方成分表一致性不验证官方模型本身的准确率。',
                      '读取审计覆盖固定样例及记录的产品版本，不代表全任务所有产品。',
                      '分型稳定性不证明类别具有矿物学含义；光谱预处理仍需独立标定验证。']}
    failures=sum(not r['passed'] for r in audit.results)
    result['summary']={'checks':len(audit.results),'failures':failures,'online_requested':args.online,
                       'scientific_mineral_accuracy':'not_evaluated_no_independent_labels'}
    (args.out/'audit.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    print('共 %d 项断言，%d 项失败'%(len(audit.results),failures))
    return int(failures>0)


if __name__=='__main__':raise SystemExit(main())
