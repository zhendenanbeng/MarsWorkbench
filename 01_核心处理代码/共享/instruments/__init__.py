# -*- coding: utf-8 -*-
"""仪器后端注册表：火星车 → 仪器。

SHERLOC 不在这里——它有自己的完整工作台（view.py），
启动器把 SHERLOC 卡片直接指过去。这里只登记走通用工作台（console.py）
的新仪器。
"""

from __future__ import annotations

from pkgutil import extend_path
__path__ = extend_path(__path__, __name__)

from .apxs import APXSBackend
from .chemcam import ChemCamBackend
from .supercam import SuperCamBackend

REGISTRY = {
    "chemcam": ChemCamBackend(),
    "supercam": SuperCamBackend(),
    "apxs": APXSBackend(),
}


def by_rover():
    """{火星车 id: [(backend), …]}，供启动器分组显示。"""
    out = {}
    for b in REGISTRY.values():
        out.setdefault(b.rover, []).append(b)
    return out
