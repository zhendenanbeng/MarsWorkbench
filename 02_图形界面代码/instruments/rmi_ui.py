"""rmi 的工作台显示块；处理算法位于核心目录。"""
HTML='''<details id="rmiPanel" class="elebox"><summary>RMI图像（浏览 + 元素分布图）</summary>
<div style="margin:4px 0"><button id="rmiMapBtn">生成元素分布图</button>
<label style="margin-left:10px"><input type="checkbox" id="rmiOverlay" disabled> 叠加光斑投影</label>
<label style="margin-left:10px">高亮元素 <select id="rmiElem" disabled><option value="">（无）</option></select></label></div>
<p id="rmiMapMsg" class="ev"></p>
<select id="rmiSelect"></select><button id="rmiReset">还原图像</button><p id="rmiInfo" class="ev"></p>
<div id="rmiViewport" style="height:400px;overflow:auto;position:relative;background:#222"><img id="rmiImage" alt="RMI上下文图像" style="display:none;max-width:none"></div></details>'''

JS=r'''
var RMI_ROWS=[],RMI_ID=null,RMI_TOKEN=0,RMI_SCALE=1,RMI_IMAGE_GENERATION=0;
var RMI_MAPS=null,RMI_MAP_SCAN=null,RMI_MAP_TOKEN=0;
function rmiRawURL(index){return '/api/rmi/image?app=chemcam&scan='+scanSel.value+'&image='+index+'&generation='+RMI_IMAGE_GENERATION;}
function rmiMapURL(index){var u='/api/rmi/overlay?app=chemcam&scan='+scanSel.value+'&image='+index;
  var e=$('rmiElem').value;if(e)u+='&elem='+encodeURIComponent(e);return u;}
function showRMI(){if(!RMI_ROWS.length)return;var row=RMI_ROWS[Number($('rmiSelect').value)],image=$('rmiImage');
  var useMap=RMI_MAPS&&RMI_MAP_SCAN===scanSel.value&&$('rmiOverlay').checked;
  var mapRow=useMap?(RMI_MAPS.maps.filter(function(m){return String(m.sclk)===String(row.sclk);})[0]||null):null;
  var url=mapRow?rmiMapURL(row.index):rmiRawURL(row.index);
  $('rmiInfo').textContent='SCLK '+row.sclk+' · '+row.relationship+'；'+row.marker_reason+
    (mapRow?'。已叠加光斑投影（'+mapRow.basis+'，'+mapRow.points+' 点'+($('rmiElem').value?'，高亮 '+$('rmiElem').value:'')+'；未做地面核验，无比例尺）':'。')+
    '滚轮缩放、双击还原；1–99%灰度显示拉伸，不改变原图。';
  if(image.getAttribute('src')===url&&image.naturalWidth){image.style.display='block';resizeRMI();return;}
  image.onload=function(){RMI_SCALE=1;resizeRMI();image.style.display='block';};
  image.onerror=function(){if(mapRow){$('rmiInfo').textContent='叠加图读取失败，先点「生成元素分布图」；当前显示原图。';image.onerror=function(){image.style.display='none';$('rmiInfo').textContent='RMI图像读取失败；未显示旧图。';};image.src=rmiRawURL(row.index);return;}
    image.style.display='none';$('rmiInfo').textContent='RMI图像读取失败；未显示旧图。';};
  image.style.display='none';image.src=url;
}
function resizeRMI(){var image=$('rmiImage');image.style.width=(Math.min(image.naturalWidth,$('rmiViewport').clientWidth||image.naturalWidth)*RMI_SCALE)+'px';image.style.height='auto';}
function resetRmiMap(){RMI_MAPS=null;RMI_MAP_SCAN=null;$('rmiOverlay').checked=false;$('rmiOverlay').disabled=true;
 $('rmiElem').value='';$('rmiElem').disabled=true;$('rmiElem').replaceChildren(new Option('（无）',''));$('rmiMapMsg').textContent='';}
function loadRMIs(){if(!$('rmiPanel'))return;var id=$('folderPath').textContent+'|'+scanSel.value+'|'+pt.value;if(id===RMI_ID)return;RMI_ID=id;var token=++RMI_TOKEN;
 if(RMI_MAP_SCAN!==scanSel.value)resetRmiMap();
 fetch('/api/rmi?app=chemcam&scan='+scanSel.value+'&point='+pt.value).then(function(r){return r.json();}).then(function(d){if(token!==RMI_TOKEN)return;
 if(d.error){$('rmiInfo').textContent=d.error;$('rmiImage').style.display='none';return;}
 var previous=$('rmiSelect').value;RMI_ROWS=d.images;$('rmiSelect').replaceChildren();
 RMI_ROWS.forEach(function(r){var option=document.createElement('option');option.value=r.index;option.textContent=r.name;$('rmiSelect').appendChild(option);});
 if(!RMI_ROWS.length){$('rmiInfo').textContent='没有同Sol/序列的本地RMI图像；可在下载面板选择CR0 PRC。';$('rmiImage').style.display='none';return;}
 if(RMI_ROWS.some(function(r){return String(r.index)===previous;}))$('rmiSelect').value=previous;showRMI();
 }).catch(function(e){$('rmiImage').style.display='none';$('rmiInfo').textContent='RMI列表读取失败：'+e;});}
if($('rmiPanel')){
 $('rmiPanel').addEventListener('toggle',function(){if(this.open)resizeRMI();});
 $('rmiSelect').addEventListener('change',showRMI);$('rmiReset').addEventListener('click',function(){RMI_SCALE=1;resizeRMI();});
 $('rmiOverlay').addEventListener('change',showRMI);
 $('rmiElem').addEventListener('change',function(){if($('rmiOverlay').checked)showRMI();});
 $('rmiMapBtn').addEventListener('click',function(){
   var token=++RMI_MAP_TOKEN;$('rmiMapMsg').textContent='正在生成元素分布图…（含按需下载 EDR label/原片，首次较慢）';
   fetch('/api/rmi/map?app=chemcam&scan='+scanSel.value).then(function(r){return r.json();}).then(function(d){
     if(token!==RMI_MAP_TOKEN)return;
     if(d.error){$('rmiMapMsg').textContent='元素分布图失败：'+d.error;RMI_MAPS=null;$('rmiOverlay').disabled=true;$('rmiOverlay').checked=false;return;}
     RMI_MAPS=d;RMI_MAP_SCAN=scanSel.value;
     var msg='已生成 '+d.maps.length+' 帧元素分布图 → _maps/：'+d.maps.map(function(m){return m.file+'（'+m.basis+'，'+m.points+'点，坐标表 '+m.csv+'）';}).join('；');
     if(d.warnings&&d.warnings.length)msg+=' ⚠ '+d.warnings.join('；');
     msg+='。'+d.note;
     $('rmiMapMsg').textContent=msg;$('rmiOverlay').disabled=false;$('rmiOverlay').checked=true;
     var sel=$('rmiElem');sel.replaceChildren(new Option('（无）',''));
     (d.elements||[]).forEach(function(e){sel.add(new Option(e,e));});
     sel.disabled=!(d.elements||[]).length;sel.value='';
     showRMI();
   }).catch(function(e){if(token===RMI_MAP_TOKEN)$('rmiMapMsg').textContent='元素分布图请求失败：'+e;});
 });
 $('rmiViewport').addEventListener('dblclick',function(){RMI_SCALE=1;resizeRMI();});
 $('rmiViewport').addEventListener('wheel',function(e){if(!$('rmiImage').naturalWidth)return;e.preventDefault();var v=this,r=v.getBoundingClientRect(),x=e.clientX-r.left,y=e.clientY-r.top,old=RMI_SCALE;
 RMI_SCALE=Math.max(.25,Math.min(8,old*(e.deltaY>0?1/1.2:1.2)));resizeRMI();v.scrollLeft=(v.scrollLeft+x)*RMI_SCALE/old-x;v.scrollTop=(v.scrollTop+y)*RMI_SCALE/old-y;
 },{passive:false});
}
'''

