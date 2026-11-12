"""First-run bootstrap: create the local .venv and install dependencies.

Called by 启动.bat / 安装依赖.bat so that Chinese paths and pinned versions
never pass through cmd.exe batch parsing. Usage:

    py -3.12 bootstrap.py            # standard deps (numpy + pillow)
    py -3.12 bootstrap.py full       # add opencv-python + astropy
"""
from __future__ import annotations
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
VENV = ROOT / '.venv'
PY = VENV / 'Scripts' / 'python.exe'
LOCK = ROOT / '05_安装与发布' / '依赖' / 'requirements-lock.txt'

STANDARD = ['numpy>=2.2,<2.5', 'pillow>=11,<14', 'setuptools>=77']
FULL_EXTRA = ['opencv-python>=4.10,<6', 'astropy>=6.1,<9']


def pip(*args: str) -> int:
    return subprocess.run([str(PY), '-m', 'pip', '--disable-pip-version-check', *args]).returncode


def main(argv: list[str]) -> int:
    if not PY.exists():
        print('正在创建虚拟环境 .venv ...')
        if subprocess.run([sys.executable, '-m', 'venv', str(VENV)]).returncode:
            print('虚拟环境创建失败，请确认已安装 Python 3.10 - 3.13。')
            return 1
    if 'full' in argv:
        print('正在安装完整版依赖（约 200 MB）...')
        specs = STANDARD + FULL_EXTRA
    else:
        print('正在安装标准版依赖（numpy + pillow）...')
        specs = STANDARD
    # Python 3.12 与锁定文件在同一组合上验证过：优先使用完全一致的版本钉。
    if 'full' not in argv and LOCK.is_file() and sys.version_info[:2] == (3, 12):
        code = pip('install', '-r', str(LOCK), 'setuptools>=77', '--quiet')
    else:
        code = pip('install', *specs, '--quiet')
    if code:
        print('依赖安装失败，请检查网络后重试。')
        return 1
    print('环境就绪。')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
