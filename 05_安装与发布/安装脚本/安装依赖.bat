@echo off
chcp 65001 >nul
cd /d "%~dp0..\.."
py -3.12 -m venv "05_安装与发布\运行环境\用户源码环境"
if errorlevel 1 exit /b 1
"05_安装与发布\运行环境\用户源码环境\Scripts\python.exe" -m pip install -r "05_安装与发布\依赖\requirements-full-lock.txt"
pause
