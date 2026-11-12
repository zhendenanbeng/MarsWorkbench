# 第三方许可与致谢

MarsWorkbench 自有代码采用 MIT；这不改变第三方代码、数据、论文和图片的许可。发行包 `_internal/resources/licenses/` 保留本次安装版本的许可原文及 SHA 清单。源码中为 `01_软件工作台/resources/licenses/`。

| 组件 | 本次版本 | 许可概要 | 范围 |
|---|---|---|---|
| CPython | 3.12.10 | PSF 及其内含第三方许可 | 两版运行时 |
| NumPy | 2.4.3 | BSD-3-Clause，内含组件另见原文 | 两版 |
| Pillow | 13.0.0 | MIT-CMU/HPND，图像库另见原文 | 两版 |
| OpenCV / opencv-python | 5.0.0.93 | Apache-2.0，内含 FFmpeg 等另见 LICENSE-3RD-PARTY | 全功能版 |
| Astropy | 8.0.2 | BSD-3-Clause，内含组件另见 licenses | 全功能版 |
| PyERFA | 2.0.1.5 | BSD-3-Clause，ERFA/SOFA 通知见原文 | 全功能版 |
| astropy-iers-data | 0.2026.10.19.0.52.1 | BSD-3-Clause及数据说明 | 全功能版 |
| packaging | 26.3 | Apache-2.0 或 BSD-2-Clause | 全功能版 |
| PyYAML | 6.0.3 | MIT | 全功能版 |
| PyInstaller | 6.22.3 | GPL-2.0 带打包例外，部分文件 Apache-2.0 | 构建工具/引导器 |

完整通知也包含 OpenBLAS、LAPACK、libjpeg、libpng、zlib、FreeType 等实际 wheel 附带的许可；以随包原文及各组件条款为准。标准版携带可选组件的许可文本便于完整审计，不表示安装了这些组件。PyInstaller 打包例外允许按应用自身与依赖许可分发可执行包。

核对来源：[NumPy](https://numpy.org/doc/stable/license.html)、[Python](https://docs.python.org/3.12/license.html)、[PyInstaller 打包例外](https://pyinstaller.org/en/stable/license.html)。本软件未获得 NASA/JPL、CNES 或其他机构背书。

界面内嵌Saira可变字体的Latin子集（Copyright 2020 The Saira Project Authors），SIL Open Font License 1.1，原文随resources/licenses/Saira-OFL.txt提供。工业风UI参考endfield-ui设计规则，未使用游戏官方素材。
