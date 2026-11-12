"""absorption 的工作台显示块；处理算法位于核心目录。"""
HTML='<div id="absorptionBox" class="elebox" hidden></div>'

JS='''
function renderAbsorption(d){var box=$('absorptionBox'),a=d.absorption;if(!a){box.hidden=true;return;}
var h='<b>反射吸收特征参考（非矿物鉴定）</b><p class="ev">'+a.scope+'</p><p class="ev">'+(a.enabled?('谱列：'+a.column+'；浏览门槛：深度≥1%、局部分数≥5，非已校准检出限。'):'当前仅显示参考范围；启用“CP2反射吸收参考”后使用原产品I/F数据评分。')+'</p>';
a.references.forEach(function(r){h+='<div class="row2"><span>'+r.label+'</span><span class="ev">'+(r.depth==null?'':('深度 '+(100*r.depth).toFixed(2)+'% · 分数 '+r.snr.toFixed(1)+' · '))+r.reason+'</span><a target="_blank" href="'+r.source+'">依据</a></div>';});box.innerHTML=h;box.hidden=false;}
'''

