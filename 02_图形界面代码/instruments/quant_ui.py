"""自研 PLS 成分的工作台显示块；处理算法位于核心目录（libs_quant）。"""
HTML = '<div id="quantBox" class="elebox" hidden></div>'

JS = '''
function renderQuant(d){var box=$('quantBox'),q=d.quant;if(!q){box.hidden=true;return;}
if(!q.available){box.innerHTML='<b>PLS 成分（独立标定）</b><p class="ev">'+q.reason+'</p>';box.hidden=false;return;}
var h='<b>PLS 成分（独立标定，在轨靶训练）</b><p class="ev">'+q.scope+'</p>';
h+='<table class="t"><tr><th>氧化物</th><th>预测 wt%</th><th>留一靶 RMSE</th></tr>';
var order=['SiO2','TiO2','Al2O3','FeOT','MgO','CaO','Na2O','K2O'];
order.forEach(function(ox){var v=q.values[ox];if(!v)return;
h+='<tr><td>'+ox+'</td><td>'+v.value.toFixed(2)+'</td><td>'+(v.rmse==null?'—':('±'+v.rmse.toFixed(2)))+'</td></tr>';});
h+='</table><p class="ev">训练靶：'+q.trained_targets.join('、')+'；本观测使用 '+q.n_points_used+' 个点。独立复现参考——官方成分（MOC）仍是定量引用口径。</p>';
box.innerHTML=h;box.hidden=false;}
'''
