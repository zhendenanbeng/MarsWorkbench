# -*- coding: utf-8 -*-
"""打开界面窗口的小工具。

用 Edge 的「应用模式」可以开一个没有地址栏、没有标签页的独立窗口，
任务栏里也是单独条目 —— 看起来就是个普通桌面应用，不需要装任何东西。
找不到 Edge 就退回系统默认浏览器。
"""

from __future__ import annotations

import os
import json
import subprocess
import tempfile
import webbrowser


def remember_instrument(instrument, settings_path=None):
    """原子更新最近使用的仪器；保留已有目录与其他设置。"""
    if instrument not in ("sherloc", "supercam", "chemcam", "apxs"):
        return False
    import appenv
    path = settings_path or appenv.settings('console')
    try:
        from instruments.history import update
        return update(path,lambda data:None if data.get('last_app')==instrument else dict(data,last_app=instrument))
    except (OSError,ValueError):
        return False


def find_edge():
    cands = [
        os.path.join(os.environ.get("ProgramFiles(x86)", r"C:\Program Files (x86)"),
                     r"Microsoft\Edge\Application\msedge.exe"),
        os.path.join(os.environ.get("ProgramFiles", r"C:\Program Files"),
                     r"Microsoft\Edge\Application\msedge.exe"),
        os.path.join(os.environ.get("LOCALAPPDATA", ""),
                     r"Microsoft\Edge\Application\msedge.exe"),
        "/usr/bin/microsoft-edge",
        "/Applications/Microsoft Edge.app/Contents/MacOS/Microsoft Edge",
    ]
    for p in cands:
        if p and os.path.exists(p):
            return p
    return None


def open_app_window(url, prefer_window=True, size=(1300, 900)):
    """打开界面。返回 'window' 或 'browser'。"""
    if os.environ.get('MARSWORKBENCH_HEADLESS')=='1':return 'headless'
    if prefer_window:
        edge = find_edge()
        if edge:
            profile = os.path.join(tempfile.gettempdir(), "sherloc_viewer_profile")
            try:
                subprocess.Popen([edge, "--app=" + url,
                                  "--window-size=%d,%d" % size,
                                  "--user-data-dir=" + profile])
                return "window"
            except Exception:                         # noqa: BLE001
                pass
    webbrowser.open(url)
    return "browser"
