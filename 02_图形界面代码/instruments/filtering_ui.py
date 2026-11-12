"""filtering 的工作台显示块；处理算法位于核心目录。"""
HTML='<input id="scanFilter" type="search" placeholder="搜索目标 / Sol号 / 元素如Fe" style="width:100%;margin-bottom:5px"><p id="filterInfo" class="sub"></p>'

JS='''
var FILTER_TOKEN=0,FILTER_TIMER=null,FILTER_ABORT=null;
function observationControls(on){['pt','ptPrev','ptNext','ptBack10','ptFwd10','shot','csvBtn','lockBtn'].forEach(function(id){if($(id))$(id).disabled=!on;});}
function resetScanFilter(){FILTER_TOKEN++;if(FILTER_TIMER)clearTimeout(FILTER_TIMER);if(FILTER_ABORT)FILTER_ABORT.abort();$('scanFilter').value='';$('filterInfo').textContent='';observationControls(true);$('paperNote').disabled=false;if($('rmiPanel'))$('rmiPanel').hidden=false;if($('scamImagePanel'))$('scamImagePanel').hidden=false;}
function noFilteredObservation(){if(typeof resetScamImages==='function')resetScamImages(null);loadToken++;paperToken++;RMI_TOKEN++;PAPER=null;PAPER_ID=null;CHART=null;curData=null;plot.textContent='没有匹配的观测';hint.textContent='';logEl.replaceChildren();$('ptlabel').textContent='没有当前观测';$('paperMessage').textContent='没有当前观测';$('paperNote').value='';$('csvBtn').disabled=true;$('paperNote').disabled=true;
observationControls(false);document.querySelectorAll('.paperNotes button').forEach(function(b){b.disabled=true;});['elements','oxides','mocbox','absorptionBox','rmiPanel','scamImagePanel'].forEach(function(id){if($(id))$(id).hidden=true;});}
function applyFiltered(list){var previous=scanSel.value;SCANS=list;$('nscan').textContent=list.length;scanSel.replaceChildren();list.forEach(function(s){var o=document.createElement('option');o.value=s.index;o.textContent=s.label;scanSel.appendChild(o);});
 if(!list.length){noFilteredObservation();return;}
 observationControls(true);$('paperNote').disabled=false;if($('rmiPanel'))$('rmiPanel').hidden=false;if($('scamImagePanel'))$('scamImagePanel').hidden=false;
 if(list.some(function(s){return String(s.index)===previous;})){scanSel.value=previous;return;}
 scanSel.value=list[0].index;onScanChange();}
function searchScans(){var token=++FILTER_TOKEN;if(FILTER_ABORT)FILTER_ABORT.abort();FILTER_ABORT=new AbortController();$('filterInfo').textContent='按当前处理配置检索中…';
 fetch(api('scans?filter=1&q='+encodeURIComponent($('scanFilter').value)+'&configured=1&steps='+encodeURIComponent(stepParam())),{signal:FILTER_ABORT.signal}).then(function(r){return r.json();}).then(function(d){if(token!==FILTER_TOKEN||d.cancelled)return;if(d.error){$('filterInfo').textContent=d.error;return;}applyFiltered(d.scans);$('filterInfo').textContent='匹配 '+d.scans.length+' 个观测；'+(d.unresolved&&d.unresolved.length?d.unresolved.length+' 个读取失败未判定；':'')+'过滤只影响列表，批量报告仍处理当前目录全部观测。';}).catch(function(e){if(token!==FILTER_TOKEN)return;if(e.name!=='AbortError')$('filterInfo').textContent='过滤失败：'+e;});}
$('scanFilter').addEventListener('input',function(){if(FILTER_TIMER)clearTimeout(FILTER_TIMER);FILTER_TIMER=setTimeout(searchScans,350);});
'''

