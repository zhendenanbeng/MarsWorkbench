"""Expose all canonical instrument package paths to PyInstaller's static graph."""
from pathlib import Path
import project_layout

def pre_safe_import_module(api):
    for relative in project_layout.CONFIG['namespace_roots']:
        package=project_layout.ROOT/relative/'instruments'
        if package.is_dir():api.append_package_path(str(package))
