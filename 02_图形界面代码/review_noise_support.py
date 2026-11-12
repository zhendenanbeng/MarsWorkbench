# -*- coding: utf-8 -*-
"""检查候选峰组重复性是否超出扫描噪声模型，不修改默认候选或原始产品。"""

import source_paths  # source roots; no-op for bundled imports
import argparse
import csv
from datetime import datetime,timezone
import hashlib
import html
import json
from pathlib import Path
import numpy as np
from sherloc import minerals,qc,repeatability
from sherloc.organize import scan_folder


def examine(scan,simulations=199,seed=20261023,model='gaussian_ar1'):
    wl,r1,r2,r3=scan.read_regions();result=qc.analyse_spectra(wl,r1,r2,r3)
    shift=qc.nm_to_raman(wl);order=np.argsort(shift);sl=slice(qc.PIXEL_LO,qc.PIXEL_HI)
    x=shift[order][sl];y=np.asarray(r1,float)[:,order][:,sl]
    wide=(x>=minerals.FINGERPRINT[0])&(x<=4000);x=x[wide];y=y[:,wide]
    residual=np.array([row-minerals.baseline(x,row) for row in y])
    quiet=(x>=2500)&(x<=3000);finger=(x>=minerals.FINGERPRINT[0])&(x<=minerals.FINGERPRINT[1])
    raw=np.sqrt(np.mean(residual[:,quiet]**2,axis=1));used,_=minerals.robust_sigma(raw)
    if np.any(~np.isfinite(used)) or np.any(used<=0):raise ValueError('噪声RMS无效')
    calibration=repeatability.calibrate(x[finger],residual[:,quiet]/used[:,None],result['mineral_groups'],simulations,seed,model=model)
    for row,g in zip(calibration['groups'],result['mineral_groups']):
        row.update(center=float(g['band']['center']),classes=[c['class'] for c in g['classes']],
                   default_mean_snr=float(g['band']['snr']),requires_fit_review=bool(g.get('requires_review')))
    files={Path(p.data_path) for p in scan.products if Path(p.data_path).is_file()}
    files.update(Path(p.label.path) for p in scan.products if p.label.path and Path(p.label.path).is_file())
    return {'sol':str(scan.sol),'seq':scan.sequence_id,'sclk':scan.sclk,'points':len(y),'calibration':calibration,
            'point_counts':{k:result[k] for k in ('n_points_with_candidate_bands','n_points_in_groups','n_points_with_named_classes')},
            'input_files':[{'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(files)]}


def run(folder,out,sol=None,seq=None,simulations=199,model='both'):
    out=Path(out);out.mkdir(parents=True,exist_ok=False);_,scans=scan_folder(str(folder));records=[]
    for scan in scans:
        if scan.spectra is None or (sol is not None and str(scan.sol)!=str(sol)) or (seq and scan.sequence_id!=seq):continue
        for method in (['gaussian_ar1','joint_blocks'] if model=='both' else [model]):
            record=examine(scan,simulations,model=method);records.append(record)
            print(record['sol'],record['seq'],method,'above reference',sum(g['exceeds_simulated_reference'] for g in record['calibration']['groups']),
                  'diagnostics',record['calibration'].get('diagnostic_warnings'),flush=True)
    if not records:raise ValueError('没有匹配扫描')
    manifest={'generated_utc':datetime.now(timezone.utc).isoformat(),'scans':records,
              'scope':'后处理的模型条件重复性检验，不修改默认矿物候选，也不提供物相正确率',
              'code_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                             [Path(__file__),Path(repeatability.__file__),Path(minerals.__file__),Path(qc.__file__)]}}
    (out/'noise_support.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    rows=[]
    with (out/'groups.csv').open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.writer(f);writer.writerow(['sol','seq','sclk','noise_model','group','points','center_cm-1','default_classes','mean_snr','simulation_exceedance_fraction','model_status','repeatability_status'])
        for s in records:
            for g in s['calibration']['groups']:
                rate=g['simulation_exceedance_fraction'];rate_text='—' if rate is None else '%.3g'%rate
                method=s['calibration'].get('model','unsupported')
                writer.writerow([s['sol'],s['seq'],s['sclk'],method,g['group_index'],g['observed_points'],g['center'],'；'.join(g['classes']),g['default_mean_snr'],rate,s['calibration']['status'],g['status']])
                rows.append('<tr><td>%s/%s</td><td>%s</td><td>%s</td><td>%d</td><td>%.1f</td><td>%s</td><td>%s</td></tr>'%
                            (s['sol'],html.escape(s['seq']),html.escape(method),html.escape('、'.join(g['classes'])),g['observed_points'],g['default_mean_snr'],rate_text,html.escape(g['status'])))
    body='<!doctype html><meta charset="utf-8"><title>峰组噪声复现对照</title><style>body{font:15px/1.6 sans-serif;max-width:1200px;margin:30px auto}table{border-collapse:collapse}td,th{border:1px solid #ccc;padding:6px}</style><h1>峰组噪声复现对照</h1><p>模拟全指纹区挑峰与分组后的最大组点数，避免把先筛选再平均的SNR当独立证明。模拟超越率不是已校准p值或物相置信度；超出参考范围仍可能来自仪器背景，未超出也不能否定弱信号。两种模型结果可用于评估模型敏感性。</p><p><a href="groups.csv">逐组CSV</a> · <a href="noise_support.json">模型参数、模拟分布与输入/代码哈希</a></p><table><tr><th>扫描</th><th>噪声模型</th><th>默认候选</th><th>点数</th><th>组平均SNR</th><th>模拟超越率</th><th>参考状态</th></tr>'+''.join(rows)+'</table>'
    (out/'index.html').write_text(body,encoding='utf-8')
    return manifest


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('folder',type=Path);p.add_argument('-o','--out',required=True,type=Path)
    p.add_argument('--sol');p.add_argument('--seq');p.add_argument('--simulations',type=int,default=199)
    p.add_argument('--model',choices=['both','gaussian_ar1','joint_blocks'],default='both')
    a=p.parse_args();run(a.folder,a.out,a.sol,a.seq,a.simulations,a.model)


if __name__=='__main__':main()
