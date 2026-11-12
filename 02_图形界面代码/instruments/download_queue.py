"""Shared serial Sol queue: isolated progress, continue failures, stop after cancellation."""
from pathlib import Path
import re

def normalize(values):
    if isinstance(values,str):values=values.split(',')
    result=[];seen=set()
    for value in values:
        value=str(value).strip().lower()
        if not re.fullmatch(r'sol_?\d+',value):raise ValueError('火星日需为PDS目录名，例如sol_00700')
        number=int(re.search(r'\d+',value).group())
        if number not in seen:result.append(value);seen.add(number)
    if not result:raise ValueError('请至少选择一个火星日')
    return result

from .io_integrity import check_length

def target_path(folder,name):
    root=Path(folder).resolve();path=(root/name).resolve()
    if Path(name).name!=name or not path.is_relative_to(root):raise ValueError('PDS文件名不能越出下载目录')
    return str(path)

def run(job,sols,dest,download_one):
    try:
        sols=normalize(sols);root=Path(dest).resolve()
        job.update(running=True,sol_total=len(sols),sol_done=0,sol_index=0,results=[],history=[],log=[],error=None,queue_status='running',queue_out=str(root))
        for index,sol in enumerate(sols):
            if job.get('cancel'):break
            out=root/sol if len(sols)>1 else root
            job.update(sol=sol,sol_index=index+1,done=0,total=0,bytes_done=0,bytes_total=0,error=None,log=[])
            job['history'].append('[%d/%d] %s 开始'%(index+1,len(sols),sol))
            try:
                if not out.resolve().is_relative_to(root):raise ValueError('火星日目录越出下载根目录')
                download_one(sol,str(out))
            except Exception as exc:job['error']=str(exc)
            incomplete=(job.get('done',0)<job.get('total',0))
            error=job.get('error') or ('部分文件未完成' if incomplete else None)
            status='cancelled' if job.get('cancel') else ('failed' if error else 'completed')
            job['results'].append({'sol':sol,'status':status,'done':job.get('done',0),'total':job.get('total',0),'out':str(out),'error':error})
            job['history'].extend(job.get('log',[]));job['log']=[]
            job['history'].append('%s：%s%s'%(sol,{'completed':'完成','failed':'失败','cancelled':'取消'}[status],('；'+error) if error else ''))
            if status!='cancelled':job['sol_done']+=1
            if job.get('cancel'):break
        job['failed_sols']=[r['sol'] for r in job['results'] if r['status']=='failed']
        job['completed_sols']=[r['sol'] for r in job['results'] if r['status']=='completed']
        job['files_done']=sum(r['done'] for r in job['results'])
        job['queue_status']='cancelled' if job.get('cancel') else ('completed_with_errors' if job['failed_sols'] else 'completed')
        job['error']=None;job['out']=str(root)
    except Exception as exc:job['error']=str(exc);job['queue_status']='failed';job.setdefault('history',[]).append('队列错误：'+str(exc))
    finally:job['running']=False

def snapshot(job):
    return dict(job,log=list(job.get('history',[]))+list(job.get('log',[])))

JS='''
function selectedSols(){return Array.from($('dlSol').selectedOptions).map(function(o){return o.value;});}
function queueStatus(d){
 var prefix='Sol '+d.sol_index+'/'+d.sol_total+' '+(d.sol||'')+' · ';
 if(d.running){var fraction=d.bytes_total?d.bytes_done/d.bytes_total:(d.total?d.done/d.total:0);$('dlBar').style.width=(100*((d.sol_index-1)+Math.max(0,Math.min(1,fraction)))/d.sol_total)+'%';$('dlMsg').className='msg';$('dlMsg').textContent=prefix+(d.cancel?'正在取消，等待当前网络读取结束；':'')+'本段 '+(d.done||0)+'/'+(d.total||0)+' 个文件';}
 else {clearInterval(DL.timer);DL.timer=null;$('dlStart').disabled=false;$('dlCancel').hidden=true;$('dlBar').style.width=(100*(d.sol_done||0)/d.sol_total)+'%';
 $('dlMsg').className='msg'+((d.failed_sols||[]).length||d.error?' err':' ok');$('dlMsg').textContent=(d.error?('队列失败：'+d.error):(d.queue_status==='cancelled'?'已取消，保留已完成文件':'队列结束'))+'：'+(d.completed_sols||[]).length+' 天成功，'+(d.failed_sols||[]).length+' 天失败；已保留'+(d.files_done||0)+'个文件 → '+(d.queue_out||d.out);
 var button=$('dlOpenFolder');button.hidden=false;button.onclick=function(){fetch('/api/openfolder?path='+encodeURIComponent(d.queue_out||d.out)+(typeof APP==='undefined'?'':('&app='+APP)));};}
 if(d.log&&d.log.length){$('dlLog').hidden=false;$('dlLog').textContent=d.log.join('\\n');}
}
'''
HTML='<button id="dlOpenFolder" hidden>打开下载目录</button>'
