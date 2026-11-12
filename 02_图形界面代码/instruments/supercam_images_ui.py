"""supercam_images 的工作台显示块；处理算法位于核心目录。"""
HTML = '''<details id="scamImagePanel" class="elebox"><summary>SuperCam 图像浏览（CI / WATSON）</summary>
<p class="ev">仅浏览：不叠加、不投影、不做光斑判定。同Sol不表示共靶，WATSON是另一台仪器。</p>
<label>图像Sol <input id="scamImageSol" type="number" min="0" max="9999" style="width:85px"></label>
<button id="scamImageLoad">加载图像清单</button><button id="scamImageRefresh">刷新网络CI清单</button>
<p id="scamImageStatus" class="ev"></p>
<label>CI帧 <select id="scamCI"></select></label><button id="scamCIShow">查看CI</button>
<p id="scamCIStatus" class="ev"></p>
<label>WATSON序列/定位 <select id="scamWatsonGroup"></select></label><br>
<label>WATSON帧 <select id="scamWatson"></select></label><button id="scamWatsonShow">查看WATSON</button>
<p id="scamWatsonStatus" class="ev"></p><button id="scamImageReset">还原图像</button>
<p id="scamImageInfo" class="ev">按需下载，首次CI约25MB；取波段0作1–99%灰度显示拉伸，不改变科学数据，不是彩色合成。</p>
<div id="scamImageViewport" style="height:400px;overflow:auto;background:#222"><img id="scamImage" alt="CI或WATSON上下文图像" style="display:none;max-width:none"></div></details>'''

JS = r'''
var SCAM_IMAGES=null,SCAM_SOL=null,SCAM_LIST_TOKEN=0,SCAM_IMG_TOKEN=0,SCAM_SCALE=1,SCAM_OBJECT_URL=null;
function resetScamImages(sol){if(!$('scamImagePanel')||String(sol)===SCAM_SOL)return;
 SCAM_SOL=String(sol);SCAM_IMAGES=null;SCAM_LIST_TOKEN++;SCAM_IMG_TOKEN++;$('scamImageSol').value=sol==null?'':sol;
 ['scamCI','scamWatsonGroup','scamWatson'].forEach(function(id){$(id).replaceChildren();});
 ['scamCIStatus','scamWatsonStatus'].forEach(function(id){$(id).textContent='';});
 $('scamImageStatus').textContent='当前Sol '+SCAM_SOL+'，点击加载清单（不自动下载图像）';
 $('scamImage').style.display='none';if(SCAM_OBJECT_URL){URL.revokeObjectURL(SCAM_OBJECT_URL);SCAM_OBJECT_URL=null;}
}
function scamWatsonOptions(){var group=SCAM_IMAGES&&SCAM_IMAGES.watson[Number($('scamWatsonGroup').value)];$('scamWatson').replaceChildren();
 if(group)group.frames.forEach(function(r,i){$('scamWatson').add(new Option(r.sclk+' · '+r.kind+' / '+r.encoding,i));});}
function scamImageList(refresh){var sol=$('scamImageSol').value;if(!sol)return;var token=++SCAM_LIST_TOKEN;
 SCAM_IMG_TOKEN++;SCAM_IMAGES=null;$('scamImage').style.display='none';
 ['scamCI','scamWatsonGroup','scamWatson'].forEach(function(id){$(id).replaceChildren();});
 $('scamImageInfo').textContent='等待选择Sol '+sol+'的图像；未沿用上一张图。';
 $('scamImageStatus').textContent='正在加载Sol '+sol+'清单…';
 fetch('/api/scamimg?app=supercam&sol='+sol+(refresh?'&refresh=1':'')).then(function(r){return r.json();}).then(function(d){if(token!==SCAM_LIST_TOKEN)return;
 if(d.error)throw Error(d.error);SCAM_IMAGES=d;$('scamCI').replaceChildren();$('scamWatsonGroup').replaceChildren();
 d.ci.forEach(function(r,i){$('scamCI').add(new Option(r.seq+' · '+r.target+' · 点'+r.point+' · '+r.sclk,i));});
 d.watson.forEach(function(r,i){$('scamWatsonGroup').add(new Option(r.id+' · '+r.count+'帧',i));});scamWatsonOptions();
 $('scamCIStatus').textContent=d.ci_error||('CI '+d.ci.length+'帧');$('scamWatsonStatus').textContent=d.watson_error||(d.watson.length?'WATSON '+d.watson.length+'组':'本地清单无本Sol WATSON；可手动修改图像Sol，不作邻日配对');
 $('scamImageStatus').textContent='Sol '+sol+'清单已处理；请选择帧再按需取图。';
 }).catch(function(e){if(token===SCAM_LIST_TOKEN)$('scamImageStatus').textContent='图像清单失败：'+e.message;});}
function resizeScamImage(){var image=$('scamImage');if(image.naturalWidth)image.style.width=(Math.min(image.naturalWidth,$('scamImageViewport').clientWidth||image.naturalWidth)*SCAM_SCALE)+'px';image.style.height='auto';}
function showScamImage(kind){if(!SCAM_IMAGES)return;var row,query;
 if(kind==='ci'){row=SCAM_IMAGES.ci[Number($('scamCI').value)];if(!row)return;query='name='+encodeURIComponent(row.name)+'&sol='+SCAM_IMAGES.sol;}
 else{var group=SCAM_IMAGES.watson[Number($('scamWatsonGroup').value)];row=group&&group.frames[Number($('scamWatson').value)];if(!row)return;query='stem='+encodeURIComponent(row.stem);}
 var token=++SCAM_IMG_TOKEN,info=$('scamImageInfo'),image=$('scamImage');image.style.display='none';
 info.textContent='正在按需读取 '+kind.toUpperCase()+' '+row.sclk+'…首次下载可能较慢。';
 fetch('/api/scamimg/image?app=supercam&kind='+kind+'&'+query).then(function(r){if(!r.ok)return r.json().then(function(d){throw Error(d.error||'图像请求失败');});return r.blob();}).then(function(blob){if(token!==SCAM_IMG_TOKEN)return;
 if(SCAM_OBJECT_URL)URL.revokeObjectURL(SCAM_OBJECT_URL);SCAM_OBJECT_URL=URL.createObjectURL(blob);
 image.onload=function(){if(token!==SCAM_IMG_TOKEN)return;SCAM_SCALE=1;image.style.display='block';resizeScamImage();};
 image.onerror=function(){if(token===SCAM_IMG_TOKEN){image.style.display='none';info.textContent='PNG解码失败，未显示旧图';}};image.src=SCAM_OBJECT_URL;
 info.textContent=kind.toUpperCase()+' · 图像Sol '+SCAM_IMAGES.sol+' · SCLK '+row.sclk+(kind==='ci'?' · '+row.seq+' · '+row.target+' · 点'+row.point:' · '+row.kind+' / '+row.encoding)+'；仅浏览，无共靶/光斑判定。波段0，1–99%拉伸；滚轮缩放、双击还原。';
 }).catch(function(e){if(token===SCAM_IMG_TOKEN){image.style.display='none';info.textContent='图像读取失败：'+e.message;}});}
if($('scamImagePanel')){
 $('scamImageLoad').onclick=function(){scamImageList(false);};$('scamImageRefresh').onclick=function(){scamImageList(true);};
 $('scamWatsonGroup').onchange=scamWatsonOptions;$('scamCIShow').onclick=function(){showScamImage('ci');};$('scamWatsonShow').onclick=function(){showScamImage('watson');};
 $('scamImageReset').onclick=function(){SCAM_SCALE=1;resizeScamImage();};$('scamImagePanel').addEventListener('toggle',resizeScamImage);
 $('scamImageViewport').ondblclick=function(){SCAM_SCALE=1;resizeScamImage();};
 $('scamImageViewport').addEventListener('wheel',function(e){if($('scamImage').style.display==='none')return;e.preventDefault();var box=this,r=box.getBoundingClientRect(),x=e.clientX-r.left,y=e.clientY-r.top,old=SCAM_SCALE;
 SCAM_SCALE=Math.max(.25,Math.min(8,old*(e.deltaY>0?1/1.2:1.2)));resizeScamImage();box.scrollLeft=(box.scrollLeft+x)*SCAM_SCALE/old-x;box.scrollTop=(box.scrollTop+y)*SCAM_SCALE/old-y;},{passive:false});
}
'''

