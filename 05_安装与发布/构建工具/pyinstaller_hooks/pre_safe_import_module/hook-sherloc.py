"""Expose the shared SHERLOC/PIXL/UI namespace paths during static analysis."""
import project_layout

def pre_safe_import_module(api):
    for relative in project_layout.CONFIG['namespace_roots']:
        package=project_layout.ROOT/relative/'sherloc'
        if package.is_dir():api.append_package_path(str(package))
