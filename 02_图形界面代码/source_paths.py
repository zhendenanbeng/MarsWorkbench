"""Activate the documented source roots; installed binaries use their bundled modules."""
from pathlib import Path
import sys
_gui=Path(__file__).resolve().parent
_shared=_gui.parent/'01_核心处理代码/共享'
if (_shared/'project_layout.py').is_file() and not getattr(sys,'frozen',False):
    sys.path.insert(0,str(_shared))
import project_layout
project_layout.activate()
