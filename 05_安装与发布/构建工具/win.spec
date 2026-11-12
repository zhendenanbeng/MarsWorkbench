from pathlib import Path
import os,sys
PROJECT=Path(SPECPATH).parents[1]
sys.path.insert(0,str(PROJECT/'02_图形界面代码'))
import source_paths,appenv
from PyInstaller.utils.hooks import collect_submodules,copy_metadata
layout=source_paths.project_layout
variant=os.environ.get('MW_VARIANT','standard')
assert variant in ('standard','full')
roots=[PROJECT/p for p in layout.CONFIG['namespace_roots']]
hidden=[p.stem for _,p in layout.runtime_sources() if p.parent in roots and p.stem!='__init__']
hidden+=collect_submodules('sherloc')+collect_submodules('instruments')
hidden+=[logical[:-3].replace('/','.') for logical,_ in layout.runtime_sources() if logical.endswith('.py') and Path(logical).stem!='__init__']
hidden=sorted(set(hidden))
datas=[(str(appenv.resources()),'resources')]
for distribution in ('numpy','pillow')+(('opencv-python','astropy','pyerfa','astropy-iers-data','packaging','pyyaml') if variant=='full' else ()):datas+=copy_metadata(distribution)
excluded=['scipy','pandas','matplotlib','pytest','IPython','tkinter','test','pip','setuptools']
if variant=='standard':excluded+=['cv2','astropy','erfa']
else:hidden+=['cv2','astropy.io.fits','erfa']
a=Analysis([str(layout.APP/'launch.py')],pathex=[str(p) for p in roots],binaries=[],datas=datas,hiddenimports=hidden,hookspath=[str(PROJECT/'05_安装与发布/构建工具/pyinstaller_hooks')],hooksconfig={},runtime_hooks=[],excludes=excluded,noarchive=False)
pyz=PYZ(a.pure)
exe=EXE(pyz,a.scripts,[],exclude_binaries=True,name='MarsWorkbench',debug=False,bootloader_ignore_signals=False,strip=False,upx=False,console=True,icon=str(appenv.resources()/'brand/mark.ico'))
coll=COLLECT(exe,a.binaries,a.datas,strip=False,upx=False,name='MarsWorkbench')
