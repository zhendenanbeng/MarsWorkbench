"""Merge documented source roots only when producing an installation wheel."""
from pathlib import Path
import shutil
import sys
from setuptools import setup
from setuptools.command.build_py import build_py

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'02_图形界面代码'))
import source_paths
import appenv
layout=source_paths.project_layout

class MergedBuildPy(build_py):
    def find_package_modules(self,package,package_dir):
        suffix='' if package=='mars_workbench' else package.rsplit('.',1)[-1]
        return [(package,path.stem,path.relative_to(ROOT).as_posix()) for logical,path in layout.runtime_sources()
                if (Path(logical).parent.as_posix() if '/' in logical else '')==suffix]

    def build_package_data(self):
        super().build_package_data()
        target=Path(self.build_lib)/'mars_workbench/resources'
        shutil.copytree(appenv.resources(),target,dirs_exist_ok=True)

cache=layout.INSTALL/'构建缓存'
(cache/'egg-info').mkdir(parents=True,exist_ok=True)
setup(version=appenv.VERSION,cmdclass={'build_py':MergedBuildPy},options={
    'build':{'build_base':str(cache/'wheel')},
    'egg_info':{'egg_base':str(cache/'egg-info')},
})
