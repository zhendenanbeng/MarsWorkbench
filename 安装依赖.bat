@echo off
setlocal
chcp 65001 >nul
cd /d "%~dp0"
rem 标准版依赖见 pyproject.toml；完全一致的锁定版本在 05_安装与发布\依赖。
where py >nul 2>nul
if errorlevel 1 (
  echo 未找到 Python。请先安装 Python 3.12：
  echo https://www.python.org/downloads/ （安装时勾选 "Add python.exe to PATH"）
  pause
  exit /b 1
)
set "PYSPEC=-3.12"
py -3.12 --version >nul 2>nul
if errorlevel 1 set "PYSPEC=-3"
set "MODE=standard"
if "%~1"=="full" set "MODE=full"
py %PYSPEC% "bootstrap.py" %MODE%
if errorlevel 1 (
  echo 依赖安装失败，请检查网络后重试。
  pause
  exit /b 1
)
echo.
echo 依赖安装完成。双击 启动.bat 即可运行 MarsWorkbench。
pause
exit /b 0
