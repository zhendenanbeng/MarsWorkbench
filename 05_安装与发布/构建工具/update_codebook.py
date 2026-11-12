from pathlib import Path
import ast,re
R=Path(__file__).resolve().parents[2]
C=R/'01_核心处理代码';G=R/'02_图形界面代码';M=R/'04_使用说明';I=R/'05_安装与发布'
roles={
'launch.py':'启动器与仪器导航；启动本机服务、选择窗口及实际端口。',
'console.py':'SuperCam/ChemCam/APXS共用工作台；调度核心后端并呈现光谱、读数、下载与导出。',
'view.py':'SHERLOC工作台；处理界面事件、参数、后台任务和显示，不定义分类阈值。',
'pixl_view.py':'PIXL页面/API；打开目录、选择实际PMC、图像与成分查询、导出及报告调度。',
'product_ui.py':'五台共用导航、帮助、示例、主题与配色交互。',
'product_theme.py':'黄黑视觉样式、离线字体、导航图标与报告阅读样式。',
'product_api.py':'共用帮助、健康检查、示例复制、路径与偏好接口。',
'runtime_http.py':'本机HTTP封装、请求编号、日志和来源检查。',
'source_paths.py':'源码开发入口的路径启动器；安装包使用随包模块。',
'appenv.py':'版本、资源、数据/缓存/设置/日志位置的统一入口。',
'project_layout.py':'读取目录清单，连接各仪器命名空间，定位数据及构建所需资源。',
'manuals.py':'从两本手册生成程序帮助，不再维护重复的历史帮助文档。',
'samples.py':'检查和复制公开示例；原件与工作副本分开。',
'health.py':'组件、权限、缓存、网络及当前版本的可复制诊断。',
'pixl_output.py':'PIXL PNG、CSV与HTML报告输出，数值来自核心读取结果。',
}
text=['# MarsWorkbench 代码作用手册','', '适用版本：0.17.0 / 本地源码与Windows发行。','',
'## 目录与运行关系','',
'核心处理代码按五个仪器分目录；共享基础设施和验证脚本单独列出。图形界面负责交互、HTTP服务、绘图界面及导出调度。数据目录不存程序源码。','',
'为保持原有模块接口，保留了Python导入名“sherloc”和“instruments”。它们是跨源目录的命名空间：例如“sherloc.pixl”实际归PIXL；“sherloc.workbench”归图形界面。真正文件位置以下表为准。','',
'开发启动链：启动入口 → 图形界面/source_paths.py → 共享/project_layout.py → 各仪器模块。安装时将相同文件按导入名组合到发行包；不复制开发数据和旧研究报告。','',
'## 五个仪器','',
'- SHERLOC：PDS标签/光谱读取、死区/基线等预处理、选峰、拟合、分类候选、重复性及ACI/WATSON数据处理。','- SuperCam：LIBS/被动FITS读取、官方MOC查表、吸收特征参考、CI/WATSON图像读取与清单。','- ChemCam：PDS3 LIBS/被动谱、波长表、RMI和已实现的相机几何/元素地图数据处理。','- APXS：能谱读取、配对官方氧化物表及误差。','- PIXL：RQC/RXL/MCC读取、PMC配对、官方成分/误差及描述统计；原始XRF能谱尚未接入。','',
'## 代码索引','',
'实验和数据准备脚本用于研究/维护，不是软件日常启动入口；其中公开论文原件按用户要求移除，重做实验时需另备输入。固定验证样本和机器基准继续保留。','']
for group,root in [('核心处理',C),('图形界面',G),('构建与安装',I/'构建工具')]:
 text+=['### '+group,'','| 文件 | 作用 | 主要入口/函数 |','|---|---|---|']
 for path in sorted(root.rglob('*.py')):
  if any(part in ('__pycache__','.cache','build','lib','work-standard','work-full') for part in path.relative_to(root).parts):continue
  try:tree=ast.parse(path.read_text(encoding='utf-8-sig'))
  except (OSError,SyntaxError):continue
  doc=ast.get_docstring(tree) or '模块接口与辅助实现；见函数名及代码内说明。'
  purpose=roles.get(path.name,re.sub(r'\s+',' ',doc).split('。')[0].split('\n')[0][:150])
  names=[node.name for node in tree.body if isinstance(node,(ast.FunctionDef,ast.AsyncFunctionDef,ast.ClassDef))]
  text.append('| '+path.relative_to(R).as_posix()+' | '+purpose.replace('|','/')+' | '+', '.join(names[:7])+' |')
text+=['','## 开发与安装','',
'源码使用Python 3.12。依赖锁文件在“05_安装与发布/依赖”；普通使用者完整解压安装包后运行启动.bat，无需Python。','',
'开发启动：python 02_图形界面代码/launch.py','',
'完整检查：python 01_核心处理代码/验证/verify/run_all.py','',
'发行独立检查：python 05_安装与发布/构建工具/dist_probe.py','',
'冻结记录：python 05_安装与发布/构建工具/freeze_validation.py','',
'双ZIP：python 05_安装与发布/构建工具/make_zip.py --variant both','',
'开发用wheel：python -m build --wheel --no-isolation --outdir 05_安装与发布/安装包/wheels','',
'构建工具保留许可证、示例来源与SHA，检查打包源码与冻结记录相符。源代码和wheel并未做算法保密加固，不要把目录分离当成防逆向保护。','',
'## 验证与数据','',
'原165项结果回归及182项不重复预处理检查继续使用原件；科学核心10个文件迁移前后字节相同。历史报告用于测试的数值转为“验证/verify/fixtures/历史数值基准”的JSON/CSV，与当前判定独立。','',
'数据迁移的逐文件校验、当前运行记录和安装包校验属于本地构建的机器数据，不随源码仓库分发。源代码历史保留在本仓库的Git中。','']
(M/'代码作用手册.md').write_text('\n'.join(text),encoding='utf-8',newline='\n')

print(M/'代码作用手册.md')
