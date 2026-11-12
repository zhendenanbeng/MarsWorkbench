@echo off
setlocal
chcp 65001 >nul
cd /d "%~dp0"
rem 依赖锁定与打包资源位于 05_安装与发布\依赖；首次运行由 bootstrap.py 准备环境。
if exist ".venv\Scripts\python.exe" (
  start "" /min ".venv\Scripts\python.exe" "02_图形界面代码\launch.py"
  exit /b 0
)
echo [MarsWorkbench] 首次运行：正在准备 Python 环境（约 1-2 分钟，仅一次）。
echo.
where py >nul 2>nul
if errorlevel 1 (
  echo   未找到 Python。请先安装 Python 3.12：
  echo   https://www.python.org/downloads/ （安装时勾选 "Add python.exe to PATH"）
  echo   然后重新双击本文件。
  pause
  exit /b 1
)
set "PYSPEC=-3.12"
py -3.12 --version >nul 2>nul
if errorlevel 1 set "PYSPEC=-3"
py %PYSPEC% "bootstrap.py"
if errorlevel 1 (
  echo   环境准备失败：请检查网络后重新双击本文件。
  pause
  exit /b 1
)
start "" /min ".venv\Scripts\python.exe" "02_图形界面代码\launch.py"
exit /b 0
