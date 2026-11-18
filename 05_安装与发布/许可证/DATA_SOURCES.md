# 示例数据来源与引用

随包示例来自 NASA Planetary Data System（PDS）公开科学归档。它们不适用本项目代码的 MIT 授权；请遵守具体归档产品的使用、引用和致谢说明，不将公开可访问等同于所有第三方材料均为公共领域。NASA 标识和机构名称不代表对本工具的认可。

逐文件产品名、原始下载 URL、bundle/collection URL、字节数和 SHA-256 见随包 `_internal/resources/sample-data/SOURCES.md` 与 `manifest.json`。缓存种子也逐项列出；WATSON 清单是供演示筛选的历史快照，不代表任务最新完整目录。

- SHERLOC：Mars 2020 SHERLOC Raman/fluorescence 产品及 ACI/WATSON 图像；感谢 NASA/JPL-Caltech 与 SHERLOC 团队。
- SuperCam：Mars 2020 SuperCam LIBS 校准谱及官方 MOC 成分表；感谢 LANL、CNES/IRAP、NASA/JPL-Caltech 与合作机构。
- ChemCam：MSL ChemCam RDR/EDR、sol01000 ccam02000 九点 LIBS 与两帧 RMI；感谢 LANL、CNES/IRAP 和 MSL 团队。
- APXS：MSL sol999 Wallace_DRT 能谱与氧化物表；感谢 APXS 仪器与 PDS Geosciences 团队。
- WATSON 缓存清单：[PDS imgops collection](https://pds.nasa.gov/ds-view/pds/viewCollection.jsp?identifier=urn%3Anasa%3Apds%3Amars2020_imgops%3Adata_watson_imgops)，父 bundle DOI [10.17189/1522846](https://doi.org/10.17189/1522846)。

第三方软件组件（可选层，不入安装包）：
- PIQUANT（PIXL 官方开源定量软件）：github.com/pixlise/piquant，W. T. Elam 等，Caltech/University of Washington 许可（版权与免责声明随 PIXL/缓存/piquant/LICENSE.md 再分发）；本工具的 PIQUANT 桥接层只是 subprocess 编排与结果解析，不复制其物理算法。PIXL FM 配置/标定文件来自该仓库 test data（Chris Heirwegh 标准文件）。

论文引用请按 [PDS 引用指南](https://pds.nasa.gov/datastandards/citing/) 引用实际使用的数据集 DOI、产品 LIDVID/产品名与版本，另记 MarsWorkbench 版本、预处理设置及辅助判读边界。SHERLOC 的规则分类、LIBS 元素线提示、投影图和代理 AUC 均不能替代独立矿物鉴定或定量标定。
