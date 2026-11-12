"""observation_summary 的工作台显示块；处理算法位于核心目录。"""
SORT_JS='''<script>
document.querySelectorAll('#observationSummary th[data-sort]').forEach(function(th){function sortRows(){var key=th.dataset.sort,body=document.querySelector('#observationSummary tbody'),rows=Array.from(body.rows),direction=th.dataset.direction==='asc'?-1:1;rows.sort(function(a,b){return direction*(Number(a.dataset[key])-Number(b.dataset[key]));});rows.forEach(function(r){body.appendChild(r);});th.dataset.direction=direction===1?'asc':'desc';}th.addEventListener('click',sortRows);th.addEventListener('keydown',function(e){if(e.key==='Enter'||e.key===' '){e.preventDefault();sortRows();}});});
</script>'''

