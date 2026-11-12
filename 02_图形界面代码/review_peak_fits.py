# -*- coding: utf-8 -*-
"""对组平均峰作物理约束复核：保留默认候选，输出独立的非负拟合证据。

python review_peak_fits.py <数据目录> --sol 293 --seq srlc11360 -o <新输出目录>
不更改原始数据、默认类别或历史报告；复核通过不等于矿物学证实。
"""
from __future__ import annotations

import source_paths  # source roots; no-op for bundled imports

import argparse
import csv
from datetime import datetime,timezone
import hashlib
import html
import json
import re
from pathlib import Path

import numpy as np
from sherloc import minerals,qc
from sherloc.organize import scan_folder


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def review(scan,out):
    wl,r1,r2,r3=scan.read_regions()
    result=qc.analyse_spectra(wl,r1,r2,r3)
    # 精确复用候选分析的R1/裁剪/宽区基线口径；物理参数求解是另一条路径。
    shift=qc.nm_to_raman(wl);order=np.argsort(shift);sl=slice(qc.PIXEL_LO,qc.PIXEL_HI)
    x=shift[order][sl];y=np.asarray(r1,float)[:,order][:,sl]
    wide=(x>=minerals.FINGERPRINT[0])&(x<=4000);x=x[wide];y=y[:,wide]
    residual=np.array([row-minerals.baseline(x,row) for row in y])
    entries=[]
    for index,group in enumerate(result['mineral_groups'],1):
        mean=residual[group['pts']].mean(axis=0);center=group['pos']
        lo=min(group['band_lo']-40,center-60);hi=max(group['band_hi']+40,center+60)
        selected=(x>=lo)&(x<=hi);xx,yy=x[selected],mean[selected]
        candidate=minerals.fit_voigt(xx,yy,center,group['band']['fwhm'],min_fwhm=minerals.MIN_FWHM_CM,
                                    within=(center-minerals.FIT_CENTER_TOL,center+minerals.FIT_CENTER_TOL))
        physical=minerals.fit_voigt(xx,yy,center,group['band']['fwhm'],min_fwhm=minerals.MIN_FWHM_CM,
                                   within=(center-minerals.FIT_CENTER_TOL,center+minerals.FIT_CENTER_TOL),nonnegative=True)
        reasons=[]
        if candidate is None:reasons.append('候选拟合失败')
        elif not candidate['profile_valid']:reasons.append('候选需要负峰形系数')
        if physical is None:reasons.append('无物理约束拟合结果')
        else:
            if physical['at_bound']:reasons.append('物理拟合触及峰心/峰宽边界')
            sigma=group['sigma_mean']
            if not np.isfinite(sigma) or sigma<=0 or physical['amp']<minerals.SNR_MIN*sigma:
                reasons.append('物理拟合幅度未通过当前噪声门槛')
        fits={}
        for name,fit in [('candidate',candidate),('physical',physical)]:
            if fit is None:fits[name]=np.full(len(xx),np.nan)
            else:
                fits[name]=minerals._pseudo_voigt(xx,fit['amp'],fit['center'],fit['fwhm'],fit['eta'])+fit['slope']*xx+fit['intercept']
                fit['reconstructed_rms']=float(np.sqrt(np.mean((yy-fits[name])**2)))
        filename='group_%03d.csv'%index
        with (out/filename).open('w',encoding='utf-8-sig',newline='') as f:
            writer=csv.writer(f);writer.writerow(['shift_cm-1','mean_residual','candidate_model','physical_model'])
            writer.writerows(zip(xx,yy,fits['candidate'],fits['physical']))
        entries.append({'group':index,'points':group['pts'],'default_classes':[c['class'] for c in group['classes']],
                        'candidate':candidate,'physical':physical,'review_reasons':reasons,
                        'status':'需复核' if reasons else '本次物理拟合检查通过（非矿物证实）','csv':filename})
    inputs={Path(p.data_path) for p in scan.products if Path(p.data_path).is_file()}
    inputs.update(Path(p.label.path) for p in scan.products if p.label.path and Path(p.label.path).is_file())
    return {'sol':str(scan.sol),'seq':scan.sequence_id,'sclk':scan.sclk,'groups':entries,
            'input_files':[{ 'path':str(p),'sha256':digest(p)} for p in sorted(inputs)]}


def run(folder,out,sol=None,seq=None):
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    _,scans=scan_folder(str(folder));records=[]
    for scan in scans:
        if scan.spectra is None or (sol is not None and str(scan.sol)!=str(sol)) or (seq and scan.sequence_id!=seq):continue
        name=re.sub(r'[^A-Za-z0-9_-]','_','sol%s_%s_%s'%(scan.sol,scan.sequence_id,scan.sclk))
        target=out/(name+'_'+hashlib.sha256(scan.key.encode()).hexdigest()[:8]);target.mkdir()
        record=review(scan,target);record['directory']=target.name;records.append(record)
        print(record['sol'],record['seq'],len(record['groups']),'groups',flush=True)
    if not records:raise ValueError('没有匹配的可分析扫描')
    manifest={'generated_utc':datetime.now(timezone.utc).isoformat(),'method':'有符号经验候选与非负峰形系数+自由线性背景对照',
              'limits':'固定组平均谱复核；非独立矿物真值验证，不自动改类别。两条路径共用候选点和基线。',
              'scans':records,'code_sha256':{p.name:digest(p) for p in [Path(__file__),Path(minerals.__file__),Path(qc.__file__),Path(minerals.__file__).with_name('fitting.py')]}}
    (out/'review.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    rows=[]
    for record in records:
        for g in record['groups']:
            candidate,physical=g['candidate'] or {},g['physical'] or {}
            rows.append('<tr><td>Sol %s / %s</td><td>%d</td><td>%s</td><td>%s</td><td>%s</td><td><a href="%s/%s">曲线CSV</a></td></tr>'%
                        (record['sol'],html.escape(record['seq']),g['group'],html.escape('、'.join(g['default_classes'])),
                         html.escape(str(candidate.get('rms','无'))),html.escape(g['status']+'；'+'；'.join(g['review_reasons'])),
                         record['directory'],g['csv']))
    body='<!doctype html><meta charset="utf-8"><title>峰形物理约束复核</title><style>body{font:15px/1.6 sans-serif;max-width:1200px;margin:30px auto}table{border-collapse:collapse}td,th{border:1px solid #ccc;padding:6px}</style><h1>峰形物理约束复核</h1><p>保留默认候选；物理约束检查不等于矿物学证实。参数、RMS和代码/输入哈希见 <a href="review.json">review.json</a>。每组CSV包含可重建的完整模型。</p><table><tr><th>扫描</th><th>组</th><th>默认候选</th><th>候选RMS</th><th>复核状态</th><th>数值</th></tr>'+''.join(rows)+'</table>'
    (out/'index.html').write_text(body,encoding='utf-8')
    return manifest


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('folder',type=Path)
    parser.add_argument('-o','--out',required=True,type=Path);parser.add_argument('--sol');parser.add_argument('--seq')
    args=parser.parse_args();run(args.folder,args.out,args.sol,args.seq)


if __name__=='__main__':main()
