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
  box.innerHTML=h;box.hidden=false;return;}
var order=['SiO2','TiO2','Al2O3','FeOT','MgO','CaO','Na2O','K2O'];
order.forEach(function(ox){var v=q.values[ox];if(!v)return;
var oc=(official&&official.oxides[ox])?official.oxides[ox][0].toFixed(2):'—';
h+='<tr><td>'+ox+'</td><td>'+v.value.toFixed(2)+'</td><td>'+(v.rmse==null?'—':('±'+v.rmse.toFixed(2)))+'</td>'+(official?('<td>'+oc+'</td>'):'')+'</tr>';});
h+='</table><p class="ev">训练靶：'+q.trained_targets.join('、')+'；本观测使用 '+q.n_points_used+' 个点。独立复现参考——官方成分（MOC）仍是定量引用口径。'+(official?('官方列为 MOC 对本观测的预测（靶：'+official.target+'）。'):'')+'</p>';
box.innerHTML=h;box.hidden=false;}
'''
