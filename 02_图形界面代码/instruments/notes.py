"""Dataset-local paper comparison annotations; never used as classifier truth."""
from datetime import datetime,timezone
import json
import os
from pathlib import Path
import tempfile
import threading
import time
from contextlib import contextmanager
LABELS={'consistent':'与论文一致','inconsistent':'与论文不一致','pending':'待查'}
_LOCK=threading.RLock()

@contextmanager
def directory_lock(folder,name='notes.lock'):
    """Cross-process lock: SHERLOC and the common console may save the same dataset."""
    folder.mkdir(exist_ok=True)
    with (folder/name).open('a+b') as stream:
        if stream.tell()==0:stream.write(b'0');stream.flush()
        deadline=time.monotonic()+5
        while True:
            try:
                stream.seek(0)
                if os.name=='nt':
                    import msvcrt
                    msvcrt.locking(stream.fileno(),msvcrt.LK_NBLCK,1)
                else:
                    import fcntl
                    fcntl.flock(stream.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
                break
            except OSError:
                if time.monotonic()>deadline:raise ValueError('其他工作台正在写入，请稍后重试')
                time.sleep(.02)
        try:yield
        finally:
            stream.seek(0)
            if os.name=='nt':msvcrt.locking(stream.fileno(),msvcrt.LK_UNLCK,1)
            else:fcntl.flock(stream.fileno(),fcntl.LOCK_UN)

def obs_key(instrument,obs):
    return json.dumps([instrument,str(obs.sol),getattr(obs,'seq',getattr(obs,'sequence_id','')),
                       getattr(obs,'product','SHERLOC'),getattr(obs,'sclk','')],ensure_ascii=False,separators=(',',':'))

def point_key(instrument,obs,point):
    ident=str(point) if isinstance(point,int) else json.dumps([point.sclk,point.stem],ensure_ascii=False,separators=(',',':'))
    return obs_key(instrument,obs)+'|point:'+ident

def replace_file(source,dest):
    # Windows readers opened without FILE_SHARE_DELETE can briefly block replace.
    for attempt in range(50):
        try:os.replace(source,dest);return
        except PermissionError:
            if attempt==49:raise
            time.sleep(.02)

def _load(root):
    file=Path(root)/'_notes/notes.json'
    if not file.exists():return {'version':1,'annotations':{}}
    try:data=json.loads(file.read_text(encoding='utf-8'))
    except (ValueError,OSError) as exc:raise ValueError('论文标注文件损坏，保留原件；请修复后再保存') from exc
    if not isinstance(data,dict) or data.get('version')!=1 or not isinstance(data.get('annotations'),dict):raise ValueError('论文标注文件版本或结构不支持')
    for r in data['annotations'].values():
        if not isinstance(r,dict) or r.get('status') not in LABELS or not isinstance(r.get('note'),str):raise ValueError('论文标注条目无效，不能静默覆盖')
    return data

def load(root):
    folder=Path(root)/'_notes'
    if not (folder/'notes.json').exists():return {'version':1,'annotations':{}}
    with _LOCK,directory_lock(folder):return _load(root)

def read(root,instrument,obs,point):
    data=load(root)['annotations']
    return {'observation':data.get(obs_key(instrument,obs)), 'point':data.get(point_key(instrument,obs,point)),'source':'用户论文对照标注，非自动真值'}

def marked(root,instrument,obs,data=None):
    entries=(data if data is not None else load(root))['annotations'];key=obs_key(instrument,obs)
    return any(k==key or k.startswith(key+'|point:') for k in entries)

def save(root,instrument,obs,point,scope,status,note):
    if scope not in ('point','observation') or status not in LABELS:raise ValueError('无效标注范围或状态')
    if not isinstance(note,str) or len(note)>2000:raise ValueError('标注备注需为不超过2000字的文本')
    root=Path(root).resolve()
    if not root.is_dir():raise ValueError('数据目录不存在')
    key=obs_key(instrument,obs) if scope=='observation' else point_key(instrument,obs,point)
    with _LOCK,directory_lock(root/'_notes'):
        data=_load(root);data['annotations'][key]={'status':status,'note':note,'updated_utc':datetime.now(timezone.utc).isoformat()}
        folder=root/'_notes';folder.mkdir(exist_ok=True)
        fd,tmp=tempfile.mkstemp(prefix='notes-',suffix='.part',dir=folder)
        try:
            with os.fdopen(fd,'w',encoding='utf-8',newline='\n') as stream:json.dump(data,stream,ensure_ascii=False,indent=2)
            replace_file(tmp,folder/'notes.json')
        finally:
            if os.path.exists(tmp):os.remove(tmp)
    return read(root,instrument,obs,point)

def label(record):return LABELS.get((record or {}).get('status'),'')

def attach_sherloc(result,root,scan):
    entries=load(root)['annotations'];result['paper_observation']=entries.get(obs_key('sherloc',scan))
    result['paper_points']={d['point']:entries.get(point_key('sherloc',scan,int(d['point']))) for d in result.get('point_detail',[])}

HTML='''<div class="paperNotes" style="margin:8px 0;font-size:12px"><b>论文对照标注</b>
<select id="paperScope"><option value="point">当前点</option><option value="observation">当前观测</option></select>
<button type="button" data-paper="consistent">与论文一致</button><button type="button" data-paper="inconsistent">不一致</button><button type="button" data-paper="pending">待查</button>
<input id="paperNote" type="text" maxlength="2000" placeholder="一句备注（随数据目录保存）" style="width:100%;margin-top:5px">
<button type="button" id="paperSave">保存备注</button><span id="paperMessage"></span></div>'''
JS='''
var PAPER=null, paperToken=0, PAPER_ID=null;
function paperApi(){return '/api/paper?app='+(typeof APP==='undefined'?'sherloc':APP)+'&scan='+scanSel.value+'&point='+pt.value;}
function renderPaper(){var r=PAPER&&PAPER[$('paperScope').value];$('paperNote').value=r?r.note:'';$('paperMessage').textContent=r?({'consistent':'与论文一致','inconsistent':'与论文不一致','pending':'待查'}[r.status]):'尚未标注';}
function loadPaperNotes(){var id=($('folderPath')?$('folderPath').textContent:'')+'|'+scanSel.value+'|'+pt.value;if(PAPER_ID===id&&PAPER)return;PAPER_ID=id;PAPER=null;$('paperNote').value='';$('paperMessage').textContent='读取标注中…';document.querySelectorAll('.paperNotes button').forEach(function(b){b.disabled=true;});var token=++paperToken;fetch(paperApi()).then(function(r){return r.json();}).then(function(d){if(token!==paperToken)return;if(d.error){$('paperMessage').textContent=d.error;return;}PAPER=d;renderPaper();document.querySelectorAll('.paperNotes button').forEach(function(b){b.disabled=false;});}).catch(function(e){$('paperMessage').textContent='读取标注失败：'+e;});}
function savePaper(status){
  var token=++paperToken, scope=$('paperScope').value;
  document.querySelectorAll('.paperNotes button').forEach(function(b){b.disabled=true;});
  fetch(paperApi(),{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({scope:scope,status:status,note:$('paperNote').value})})
  .then(function(r){return r.json();}).then(function(d){
    if(token!==paperToken)return;
    document.querySelectorAll('.paperNotes button').forEach(function(b){b.disabled=false;});
    if(d.error){$('paperMessage').textContent=d.error;return;}
    PAPER=d;renderPaper();var opt=scanSel.options[scanSel.selectedIndex];if(opt&&!opt.textContent.startsWith('● '))opt.textContent='● '+opt.textContent;
  }).catch(function(e){if(token!==paperToken)return;document.querySelectorAll('.paperNotes button').forEach(function(b){b.disabled=false;});$('paperMessage').textContent='保存失败：'+e;});
}
document.querySelectorAll('[data-paper]').forEach(function(button){button.addEventListener('click',function(){savePaper(this.dataset.paper);});});
$('paperScope').addEventListener('change',renderPaper);
$('paperSave').addEventListener('click',function(){var r=PAPER&&PAPER[$('paperScope').value];savePaper(r?r.status:'pending');});
'''
