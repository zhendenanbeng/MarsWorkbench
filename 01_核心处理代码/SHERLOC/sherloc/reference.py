"""Small pure browser helpers for immutable, coordinate-aware SHERLOC references."""
JS=r'''
function freezeReference(data,label){return JSON.parse(JSON.stringify({x:data.x,series:data.series,xlabel:data.xlabel,product:data.product,label:label}));}
function referenceSelection(reference,data){
  if(!reference)return {series:[],message:''};
  if(reference.xlabel!==data.xlabel||reference.product!==data.product)return {series:[],message:'参考产品或横轴不同，暂不叠加'};
  var keys=data.series.map(function(s){return s.key;});
  var found=reference.series.filter(function(s){return keys.indexOf(s.key)>=0;});
  var missing=keys.filter(function(k){return !reference.series.some(function(s){return s.key===k;});});
  return {series:found,message:missing.length?('参考不含 '+missing.join('、')+'，这些序列暂不叠加'):''};
}
function referenceRead(reference,series,value){
  var x=reference.x;if(!x.length||value<x[0]||value>x[x.length-1])return null;
  var lo=0,hi=x.length-1;while(hi-lo>1){var m=(lo+hi)>>1;if(x[m]<=value)lo=m;else hi=m;}
  var i=Math.abs(x[lo]-value)<=Math.abs(x[hi]-value)?lo:hi;
  var local=Math.min(i>0?x[i]-x[i-1]:Infinity,i+1<x.length?x[i+1]-x[i]:Infinity);
  var y=series.y[i];if(y===null||!isFinite(y)||Math.abs(x[i]-value)>1.5*local)return null;
  return {x:x[i],y:y};
}
function referenceEscape(text){return String(text).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;');}
'''
