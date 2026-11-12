# -*- coding: utf-8 -*-
r"""把论文 PDF 转成可以 grep 的纯文本。

用法::

    python verify\extract_paper_text.py "C:\...\Corpolongo2023.pdf"
    python verify\extract_paper_text.py "C:\...\Corpolongo2023.pdf" 输出目录

输出写到 PDF 同目录、同名加 `_text.txt`。为什么要这么一步：
质控台里那些"论文值"必须能逐条回查，而 PDF 里的数字没法直接搜。
转成文本之后，`Select-String "Montpezat"` 就能把 3.2.3 节整段拉出来。

只用 pypdf（内置运行环境里已经有）。扫描版 PDF 取不出文字，会明确报错。

给了第二个参数就写到那个目录——论文放在工作区外面时（比如桌面上的
文献阅读 目录），沙箱不允许往那儿写，只能指到工作区里来。
"""

from __future__ import annotations

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout

import os
import sys


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    src = sys.argv[1]
    if not os.path.isfile(src):
        print("找不到文件：%s" % src)
        return 2
    try:
        from pypdf import PdfReader
    except ImportError:
        print("缺少 pypdf。装法：python -m pip install pypdf")
        return 2

    reader = PdfReader(src)
    chunks = []
    for i, page in enumerate(reader.pages):
        chunks.append("\n\n===== PAGE %d =====\n" % (i + 1))
        chunks.append(page.extract_text() or "")
    text = "".join(chunks)
    if len(text.strip()) < 500:
        print("只取到 %d 个字符——这份 PDF 大概是扫描件，取不出文字层。" % len(text.strip()))
        return 1

    out_dir = sys.argv[2] if len(sys.argv) > 2 else os.path.dirname(os.path.abspath(src))
    os.makedirs(out_dir, exist_ok=True)
    dst = os.path.join(out_dir, os.path.splitext(os.path.basename(src))[0] + "_text.txt")
    try:
        with open(dst, "w", encoding="utf-8") as f:
            f.write(text)
    except OSError as exc:
        print("写不进去：%s\n提示：论文在工作区外时，把输出目录指到工作区里，例如：\n"
              r"  python verify\extract_paper_text.py <论文.pdf> ..\_libraries"
              % exc)
        return 1
    print("共 %d 页，写出 %d 个字符 →\n  %s" % (len(reader.pages), len(text), dst))
    return 0


if __name__ == "__main__":
    sys.exit(main())
