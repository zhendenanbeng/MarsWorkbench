"""Recent-folder settings and shared, editing-safe navigation shortcuts."""
from copy import deepcopy
import json
import os
from pathlib import Path
import tempfile
import threading
from .notes import directory_lock,replace_file
_LOCK=threading.RLock()

def remember(config,path,instrument=None):
    config=deepcopy(config);path=os.path.abspath(os.path.expanduser(path))
    stored=config.get('recent_folders',{} if instrument else [])
    if instrument:
        if not isinstance(stored,dict):stored={}
        values=stored.get(instrument,[])
    else:values=stored
    if not isinstance(values,list):values=[]
    values=[path]+[p for p in values if isinstance(p,str) and os.path.normcase(os.path.abspath(p))!=os.path.normcase(path)]
    if instrument:stored[instrument]=values[:5];config['recent_folders']=stored
    else:config['recent_folders']=values[:5]
    return config

def rows(config,instrument=None):
    values=config.get('recent_folders',{} if instrument else [])
    if instrument:values=values.get(instrument,[]) if isinstance(values,dict) else []
    if not isinstance(values,list):return []
    return [{'path':p,'exists':os.path.isdir(p)} for p in values[:5] if isinstance(p,str)]

def _read(path):
    try:data=json.loads(path.read_text(encoding='utf-8')) if path.exists() else {}
    except (ValueError,OSError) as exc:raise ValueError('设置文件损坏，未覆盖原件') from exc
    if not isinstance(data,dict):raise ValueError('设置文件结构无效')
    import project_layout
    return project_layout.remap_preferences(data)

def read(path):
    path=Path(path).resolve()
    if not path.exists():return {}
    with _LOCK,directory_lock(path.parent,path.name+'.lock'):return _read(path)

def update(path,mutator):
    path=Path(path).resolve()
    with _LOCK,directory_lock(path.parent,path.name+'.lock'):
        data=_read(path)
        result=mutator(data)
        if result is None:return False
        fd,tmp=tempfile.mkstemp(prefix='.settings-',suffix='.part',dir=path.parent)
        try:
            with os.fdopen(fd,'w',encoding='utf-8',newline='\n') as stream:json.dump(result,stream,ensure_ascii=False,indent=1)
            replace_file(tmp,path)
        finally:
            if os.path.exists(tmp):os.remove(tmp)
    return True

def store_folder(path,folder,instrument=None):
    folder=os.path.abspath(folder)
    def changed(data):
        result=remember(data,folder,instrument)
        if instrument:
            if not isinstance(result.get('folders'),dict):result['folders']={}
            result['folders'][instrument]=folder
        else:result['folder']=folder
        return result
    return update(path,changed)

NAV_JS=r'''
function retainedPointId(previous,ids,preferred){return ids.includes(previous)?previous:ids.includes(preferred)?preferred:ids.length?ids[0]:null;}
function navigationCommand(key,buffer,elapsed,selected,count,last,base,digitsLimit){
 if(count<1)return null;
 if(key==='PageUp'||key==='PageDown')return {type:'scan',index:Math.max(0,Math.min(count-1,selected+(key==='PageUp'?-1:1)))};
 if(key==='Home'||key==='End')return {type:'point',index:key==='Home'?0:last};
 if(/^\d$/.test(key)){buffer=(elapsed>700||buffer.length>=(digitsLimit||3)?'':buffer)+key;var index=Number(buffer)-base;return {type:'number',buffer:buffer,index:index>=0&&index<=last?index:null};}
 return null;
}
'''
JS=NAV_JS+r'''
var NUMBER_BUFFER='',NUMBER_TIME=-Infinity,RECENT_TOKEN=0;
function workflowKey(e){
 var tag=e.target&&e.target.tagName;if(e.ctrlKey||e.altKey||e.metaKey||e.isComposing||tag==='TEXTAREA'||tag==='SELECT'||(tag==='INPUT'&&e.target!==pt))return false;
 var command=navigationCommand(e.key,NUMBER_BUFFER,e.timeStamp-NUMBER_TIME,scanSel.selectedIndex,scanSel.options.length,Number(pt.max)||0,__POINT_BASE__);
 if(!command){NUMBER_BUFFER='';return false;}e.preventDefault();
 if(command.type==='scan'){NUMBER_BUFFER='';if(command.index!==scanSel.selectedIndex){scanSel.selectedIndex=command.index;onScanChange();}}
 else {if(command.type==='number'){NUMBER_BUFFER=command.buffer;NUMBER_TIME=e.timeStamp;}else NUMBER_BUFFER='';if(command.index!==null&&command.index!==Number(pt.value)){pt.value=command.index;load();}}
 return true;
}
function loadRecentFolders(host){var token=++RECENT_TOKEN;fetch('/api/recent'+(typeof APP==='undefined'?'':('?app='+APP))).then(function(r){return r.json();}).then(function(data){
 if(token!==RECENT_TOKEN)return;var old=$('recentFolders');if(old)old.remove();var block=document.createElement('div');block.id='recentFolders';block.style.marginBottom='8px';var label=document.createElement('b');label.textContent='最近打开';block.appendChild(label);
 (data.folders||[]).forEach(function(item){var button=document.createElement('button');button.textContent=item.path;button.disabled=!item.exists;button.style.cssText='display:block;max-width:100%;overflow-wrap:anywhere;margin:4px 0';button.addEventListener('click',function(){loadFolder(item.path);});block.appendChild(button);});
 if(!(data.folders||[]).length){var empty=document.createElement('span');empty.textContent=' 暂无记录';block.appendChild(empty);}host.insertBefore(block,host.firstChild);
 });}
'''
