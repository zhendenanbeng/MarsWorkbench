@echo off
chcp 65001 >nul
cd /d "%~dp0..\.."
call "%cd%\启动.bat"
