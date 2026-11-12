"""The two retained handbooks are the single source of in-app help."""
import html
import json
import re
import appenv

NAMES=('操作手册.md','代码作用手册.md')

def read(name):
    if name not in NAMES:raise ValueError('未知手册')
    path=appenv.resources()/'manuals'/name
    return path.read_text(encoding='utf-8') if path.exists() else '# 手册尚未生成\n请完成当前版本的安装准备。'

def section(text,title):
    match=re.search(r'^## '+re.escape(title)+r'\s*$',text,re.M)
    if not match:return text[:5000]
    end=re.search(r'^## ',text[match.end():],re.M)
    return text[match.start():match.end()+end.start() if end else len(text)]

def inline(text):
    text=html.escape(text)
    text=re.sub(r'`([^`]+)`',r'<code>\1</code>',text)
    text=re.sub(r'\*\*([^*]+)\*\*',r'<strong>\1</strong>',text)
    def link(match):
        label,target=match.groups()
        if target in NAMES:target='/manuals/'+target
        if not target.startswith(('https://','http://','/manuals/','#')):return label
        return '<a href="'+html.escape(html.unescape(target),quote=True)+'">'+label+'</a>'
    return re.sub(r'\[([^\]]+)\]\(([^)]+)\)',link,text)

def render(text):
    output=[];code=[];in_code=False;in_list=False;in_table=False
    for line in text.splitlines():
        if line.strip().startswith('|') and line.strip().endswith('|') and not in_code:
            cells=[cell.strip() for cell in line.strip().strip('|').split('|')]
            if all(re.fullmatch(r'[:\- ]+',cell or '-') for cell in cells):continue
            if not in_table:
                output.append('<table><thead><tr>'+''.join('<th>'+inline(cell)+'</th>' for cell in cells)+'</tr></thead><tbody>');in_table=True
            else:output.append('<tr>'+''.join('<td>'+inline(cell)+'</td>' for cell in cells)+'</tr>')
            continue
        if in_table:output.append('</tbody></table>');in_table=False
        if line.startswith('```'):
            if in_code:output.append('<pre><code>'+html.escape('\n'.join(code))+'</code></pre>');code=[]
            in_code=not in_code;continue
        if in_code:code.append(line);continue
        image=re.match(r'^!\[([^\]]*)\]\((插图/[^)]+)\)$',line)
        if image:
            output.append('<img alt="'+html.escape(image[1],quote=True)+'" src="/manuals/'+html.escape(image[2],quote=True)+'">');continue
        if not line.strip():
            if in_list:output.append('</ul>');in_list=False
            continue
        heading=re.match(r'^(#{1,4})\s+(.+)$',line)
        bullet=re.match(r'^(?:[-*]|\d+\.)\s+(.+)$',line)
        if heading:
            if in_list:output.append('</ul>');in_list=False
            level=len(heading[1]);output.append('<h%d>%s</h%d>'%(level,inline(heading[2]),level))
        elif bullet:
            if not in_list:output.append('<ul>');in_list=True
            output.append('<li>'+inline(bullet[1])+'</li>')
        else:output.append('<p>'+inline(line)+'</p>')
    if in_list:output.append('</ul>')
    if in_table:output.append('</tbody></table>')
    return ''.join(output)

def help_sections():
    operation=read(NAMES[0])
    return {
        'quickstart':render(section(operation,'第一次使用'))+render(section(operation,'各工作台统一的操作方式')),
        'faq':render(section(operation,'故障排查')),
        'scope':render(section(operation,'科学口径')),
        'about':'<h2>MarsWorkbench '+appenv.VERSION+'</h2><p>本地五仪器工作台；不上传分析数据。</p><p><a href="/manuals/操作手册.md">完整操作手册</a> · <a href="/manuals/代码作用手册.md">代码作用手册</a></p>',
    }

def trust_text():
    path=appenv.resources()/'_last_run.json'
    try:
        record=json.loads(path.read_text(encoding='utf-8'))
        lines=['MarsWorkbench '+appenv.VERSION+' · 发行前冻结记录',
               '完成时间（UTC）：'+record.get('finished_utc','未知'),
               '工程检查：%s组 / %s项 / %s失败。'%(record.get('suite_count','?'),record.get('assertions','?'),record.get('failures','?')),
               '此记录来自发行前环境，不表示当前电脑重新跑过全部检查；检查通过不等于矿物鉴定准确率。']
    except (OSError,ValueError):lines=['此开发版本暂无可核验冻结记录。请运行完整检查。']
    lines.append(section(read(NAMES[0]),'科学口径'))
    return '\n\n'.join(lines)
