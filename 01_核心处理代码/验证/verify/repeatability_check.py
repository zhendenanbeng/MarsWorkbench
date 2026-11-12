# -*- coding: utf-8 -*-
"""候选统计口径、全窗噪声重放与报告的行为检查；不声称验证物相准确率。"""

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout
import csv
import json
from pathlib import Path
import sys
import tempfile
from unittest.mock import patch
import numpy as np

ROOT=_layout.project_layout.APP;sys.path.insert(0,str(ROOT))
from sherloc import minerals as m,qc,repeatability as r
import review_noise_support as report
import view

TOTAL=FAIL=0
def check(ok,name):
    global TOTAL,FAIL
    TOTAL+=1
    if not ok:FAIL+=1;print(' FAIL '+name)

def rejected(call,name,error=ValueError):
    try:call()
    except error:check(True,name)
    else:check(False,name)

def detector():
    rng=np.random.default_rng(6401);x=np.arange(800.,1801.,5.)
    y=rng.normal(size=(12,len(x)))
    y+=7*np.exp(-4*np.log(2)*((x-1080)/45)**2)
    y[0,80]=30;y[1,0]=40;y[2,50:53]=8
    expected=sorted((b['pos'],i) for i,row in enumerate(y) for b in m.find_bands(x,row,1.,prom_sigma=0))
    check(sorted(r.peak_hits(x,y))==expected,'重放候选须逐项等于默认逐点筛选（含尖刺/边界/平台）')
    check(r.max_cluster_count([])==0,'空候选最大组为0')
    check(r.max_cluster_count([(900,0),(910,0),(915,1)])==2,'同点多个候选不能重复计点')
    check(r.max_cluster_count([(900,0),(916,1)])==1,'超过相邻间隙须断组')
    chain=[(900+15*i,i) for i in range(6)]
    check(r.max_cluster_count(chain)==5,'链式聚类须遵守总跨度60')
    groups=[]
    for pos,point in sorted(expected):
        if not groups or pos-groups[-1][-1][0]>15 or pos-groups[-1][0][0]>60:groups.append([])
        groups[-1].append((pos,point))
    check(r.max_cluster_count(expected)==max(len({i for _,i in g}) for g in groups),'最大组数与直接分组独立写法一致')

def replay():
    rng=np.random.default_rng(8157);x=np.arange(800.,1801.,5.);q=rng.normal(size=(16,100))
    original=q.copy();groups=[{'pts':list(range(16))},{'pts':[0]}]
    for model in ['gaussian_ar1','joint_blocks']:
        result=r.calibrate(x,q,groups,simulations=39,model=model)
        same=r.calibrate(x,q,groups,simulations=39,model=model)
        check(result==same,model+'固定种子须可重现全部输出')
        check(result['status']=='reference_only',model+'正常噪声也只能得到模型参考结论')
        check(len(result['null_max_cluster_points'])==39,model+'实际模拟次数不得缩水')
        check(result['groups'][0]['exceeds_simulated_reference'],model+'16点共同峰组应超出本次独立噪声参考')
        expected_single=(1+sum(v>=1 for v in result['null_max_cluster_points']))/40
        check(result['groups'][1]['simulation_exceedance_fraction']==expected_single,model+'单点也须对照全窗最大组，不能预设结论')
        expected=(1+sum(v>=16 for v in result['null_max_cluster_points']))/40
        check(result['groups'][0]['simulation_exceedance_fraction']==expected,model+'超越率含有限模拟修正')
        check(result['monte_carlo_resolution']==1/40,model+'分辨率必须受模拟次数限制')
        check('不是已校准p值' in result['limits'] and '未独立验证' in result['limits'],model+'限制必须明示模型和迁移假设')
        check('p_value' not in result['groups'][0],model+'不能把模拟频率命名为已校准p值')
        json.dumps(result,allow_nan=False)
    check(np.array_equal(q,original),'估计和模拟不得更改原始安静区')
    shared=.3*rng.normal(size=100)+7*np.exp(-4*np.log(2)*((np.arange(100)-50)/10)**2)
    common=np.tile(shared,(16,1))
    joint=r.calibrate(x,common,groups,simulations=39,model='joint_blocks')
    independent=r.calibrate(x,common,groups,simulations=39,model='gaussian_ar1')
    check(np.median(joint['null_max_cluster_points'])>np.median(independent['null_max_cluster_points']),'联合块重放必须保留独立点模型丢失的共同结构')
    check(not joint['groups'][0]['exceeds_simulated_reference'],'可由共同安静区结构重放的全点峰组不能当独立异常')
    heavy=q.copy();heavy[:,::10]=100
    failed=r.calibrate(x,heavy,groups,simulations=19)
    check(failed['status']=='diagnostics_failed' and failed['diagnostic_warnings'],'厚尾不符合Gaussian时必须提示')
    check(not any(g['exceeds_simulated_reference'] for g in failed['groups']),'模型诊断失败时不能输出超出参考的肯定状态')
    zero=r.calibrate(x,np.zeros_like(q),groups,simulations=19,model='joint_blocks')
    check(zero['status']=='unsupported' and zero['model']=='joint_blocks','零噪声应保留所请求模型及失败语义')
    check(len(zero['groups'])==2 and all(g['simulation_exceedance_fraction'] is None for g in zero['groups']),'无法建模不能静默丢组或输出0超越率')
    for args in [dict(simulations=18),dict(simulations=2000),dict(simulations=True),dict(alpha=0),dict(model='unknown')]:
        rejected(lambda:r.calibrate(x,q,groups,**args),'无效模拟参数须拒绝 '+str(args))
    rejected(lambda:r.calibrate(x[::-1],q,groups),'倒序轴须拒绝')
    rejected(lambda:r.calibrate(x,q[:1],groups),'仅1点无法校准扫描噪声')
    rejected(lambda:r.calibrate(x,q,[{'pts':[len(q)]}]),'越界成员须拒绝')

def counts_and_products():
    x=np.arange(800.,4001.,5.);y=np.tile(np.where(np.arange(len(x))%2,1.,-1.),(4,1))
    for i,center in enumerate([950,1080,1080,1300]):y[i]+=20*np.exp(-4*np.log(2)*((x-center)/50)**2)
    original=y.copy()
    with patch.object(m,'baseline',side_effect=lambda x,y:np.zeros_like(y)):
        result=m.analyse_scan(x,y)
    check(result['n_points_with_candidate_bands']==4,'统计须包括未成组的孤立候选点')
    check(result['n_points_in_groups']==2,'孤立候选不能计入峰组点数')
    check(result['n_points_with_band']==result['n_points_in_groups'],'旧字段保留历史峰组含义及兜底行为')
    check(result['n_points_with_named_classes']==2,'命名点数按唯一原始点计数')
    check(all('不是独立验证' in g['selection_basis'] for g in result['groups']),'峰组须携带先筛选后平均的证据限制')
    check(np.array_equal(y,original),'统计更正不能更改光谱')
    rng=np.random.default_rng(17);shift=np.linspace(0,4500,2148);wl=1/(1/qc.LASER_NM-shift*1e-7)
    r1=np.array([30+100*np.exp(-.5*((shift-1084)/25)**2)+rng.normal(0,1,len(shift)) for _ in range(6)])
    class Scan:
        sol='9999';sequence_id='synthetic';sclk='123';products=[];spectra=True;key='synthetic/123'
        def read_regions(self):return wl,r1,None,None
    res=qc.analyse_spectra(wl,r1)
    check(res['n_points_with_candidate_bands']==sum(p['n_bands']>0 for p in res['point_detail']),'QC输出须保持逐点候选统计一致')
    check(res['n_points_in_groups']==len({i for g in res['mineral_groups'] for i in g['pts']}),'QC峰组点数须去重')
    summary=view.mineral_summary(view.scan_minerals({'mineral_cache':{}},Scan()))
    check(summary['n_candidate_points']==res['n_points_with_candidate_bands'] and summary['n_in_groups']==res['n_points_in_groups'],'工作台摘要须传递不同点数口径')
    with tempfile.TemporaryDirectory(prefix='sherlock-replay-') as folder:
        out=Path(folder)/'report'
        with patch.object(report,'scan_folder',return_value=(None,[Scan()])):
            result=report.run(folder,out,simulations=19)
        check(len(result['scans'])==2,'默认报告必须包含两种模型')
        rows=list(csv.DictReader((out/'groups.csv').open(encoding='utf-8-sig')))
        check({row['noise_model'] for row in rows}=={'gaussian_ar1','joint_blocks'},'CSV必须保留模型差异')
        check(all('simulation_exceedance_fraction' in row for row in rows),'CSV必须用描述性超越率名称')
        body=(out/'index.html').read_text(encoding='utf-8')
        check('不是已校准p值' in body and '物相置信度' in body,'HTML须直接解释诊断边界')
        check(bool(json.loads((out/'noise_support.json').read_text(encoding='utf-8'))['code_sha256']),'报告须记录执行代码哈希')
        rejected(lambda:report.run(folder,out,simulations=19),'不得覆盖已有复核证据',FileExistsError)

def main():
    detector();replay();counts_and_products()
    print('共 %d 项断言，%d 项失败'%(TOTAL,FAIL))
    return int(FAIL>0)

if __name__=='__main__':raise SystemExit(main())
