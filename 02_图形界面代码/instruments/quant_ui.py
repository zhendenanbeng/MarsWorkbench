"""自研 PLS 成分的工作台显示块；处理算法位于核心目录（libs_quant）。"""
HTML = '<div id="quantBox" class="elebox" hidden></div>'

JS = '''
function renderQuant(d){var box=$('quantBox'),q=d.quant;if(!q){box.hidden=true;return;}
if(!q.available){box.innerHTML='<b>PLS 成分（独立标定）</b><p class="ev">'+q.reason+'</p>';box.hidden=false;return;}
var h='<b>PLS 成分（独立标定，在轨靶训练）</b><p class="ev">'+q.scope+'</p>';
var official=q.official;
h+='<table class="t"><tr><th>氧化物</th><th>预测 wt%</th><th>留一靶 RMSE</th>'+(official?'<th>官方 MOC</th>':'')+'</tr>';
if(q.rows){
  var oc=q.official||null;
  h+='<table class="t"><tr><th>元素</th><th>峰高 (计数)</th><th>±误差</th><th>检出</th>'+(oc?'<th>官方 wt (氧化物)</th>':'')+'</tr>';
  q.rows.forEach(function(r){var ov=(oc&&oc[r.label]!==undefined)?oc[r.label]:null;
    h+='<tr><td>'+r.label+'</td><td>'+r.value.toFixed(0)+'</td><td>±'+r.err.toFixed(0)+'</td><td>'+(r.detected?'检出':'—')+'</td>'+(oc?('<td>'+(ov===null?'—':ov)+'</td>'):'')+'</tr>';});
  h+='</table><p class="ev">本观测使用 '+q.n_points_used+' 个点。'+(oc?'官方列=配对 rwp 的元素 wt（对照参考）。':'')+'</p>';
  if(q.fp){
    if(q.fp.available){
      h+='<p style="margin-top:10px"><b>基本参数法自算 wt%（原始谱独立定量）</b></p><p class="ev">'+(q.fp.scope||'')+' 标定时代：'+(q.fp.calib_era||'—')+'；Sherman 迭代 '+q.fp.iterations+' 轮'+(q.fp.converged?'收敛':'未收敛')+'。'+(q.fp.live_median_s?('本观测积分中位 '+Math.round(q.fp.live_median_s)+' s'+(q.fp.live_median_s<2400?'（低统计——轻元素 MgO/Al2O3/Na2O 偏差可超 50%，对账以官方列为准）':'')+''):'')+'</p>';
      h+='<table class="t"><tr><th>氧化物</th><th>FP wt%</th><th>官方 rwp</th><th>相对偏差</th></tr>';
      (q.fp.compare||[]).forEach(function(r){
        var fv=(r.fp===null||r.fp===undefined)?'—':r.fp.toFixed(2);
        var ov=(r.official===null||r.official===undefined)?'—':r.official.toFixed(2);
        var rv=(r.rel===null||r.rel===undefined)?'—':(100*r.rel).toFixed(1)+'%';
        h+='<tr><td>'+r.oxide+'</td><td>'+fv+'</td><td>'+ov+'</td><td>'+rv+'</td></tr>';});
      h+='</table>';
    }else{
      h+='<p class="ev"><b>基本参数法自算：</b>'+q.fp.reason+'</p>';
    }
  }
  box.innerHTML=h;box.hidden=false;return;}
var order=['SiO2','TiO2','Al2O3','FeOT','MgO','CaO','Na2O','K2O'];
var neg=[];
order.forEach(function(ox){var v=q.values[ox];if(!v)return;
var vv=v.value.toFixed(2);if(v.value<0){neg.push(ox);vv+='↑';}
var oc='—';
if(official){if(official.oxides&&official.oxides[ox]!==undefined)oc=official.oxides[ox][0].toFixed(2);else if(official[ox]&&Array.isArray(official[ox]))oc=official[ox][0].toFixed(2);}
h+='<tr><td>'+ox+'</td><td>'+vv+'</td><td>'+(v.rmse==null?'—':('±'+v.rmse.toFixed(2)))+'</td>'+(official?('<td>'+oc+'</td>'):'')+'</tr>';});
h+='</table><p class="ev">训练靶：'+q.trained_targets.join('、')+'；本观测使用 '+q.n_points_used+' 个点。独立复现参考——官方成分（MOC）仍是定量引用口径。'+(official?('官方列为 MOC 对本观测的预测（靶：'+official.target+'）。'):'')+(neg.length?('负值预测（'+neg.join('、')+'）↑ =低于检出参考零（|值|在 ±RMSE 内属正常波动），不截断以保真。'):'')+'</p>';
box.innerHTML=h;box.hidden=false;}
'''
