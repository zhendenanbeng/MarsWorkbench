"""MarsWorkbench package; compatibility imports preserve the validated science modules."""
from pathlib import Path
import sys
_root=str(Path(__file__).resolve().parent)
if _root not in sys.path:sys.path.insert(0,_root)
import source_paths
if source_paths.project_layout.ROOT:
    __path__.append(str(source_paths.project_layout.CORE/'共享'))
from .appenv import VERSION as __version__
