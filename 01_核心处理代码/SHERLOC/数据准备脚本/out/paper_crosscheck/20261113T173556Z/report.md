# 软件分类输出 vs 已发表结论（逐靶对比检查）

- 台账：published_targets.json
- 判定口径：类别候选复现性（非准确率）

## Garde（03_下载数据\SHERLOC\练习数据\raw\sol_00207）

| 论文类目 | 判定 | 软件键 | 说明 |
|---|---|---|---|
| carbonate | reproduced | carbonate | 检出点数 76 |
| feldspar | reproduced | silicate | 检出点数 3 |
| olivine | reproduced | olivine | 检出点数 65 |
| perchlorate | not_reproduced |  | 本 sol 观测中无该类候选（记录差异，不判对错） |
| perchlorate or phosphate | unmapped |  | 映射表未覆盖（请人工核对后补充 PAPER_TO_KEYS） |
| phosphate | not_reproduced |  | 本 sol 观测中无该类候选（记录差异，不判对错） |
| pyroxene | not_reproduced |  | 本 sol 观测中无该类候选（记录差异，不判对错） |
| silicate | unmapped |  | 映射表未覆盖（请人工核对后补充 PAPER_TO_KEYS） |

额外候选：unassigned

## Alfalfa（03_下载数据\SHERLOC\练习数据\raw\sol_00370）

| 论文类目 | 判定 | 软件键 | 说明 |
|---|---|---|---|
| carbonate | not_reproduced |  | 本 sol 观测中无该类候选（记录差异，不判对错） |
| chlorite | reproduced | silicate | 检出点数 30 |
| fe-oxides | excluded |  | SHERLOC 深紫外拉曼类别体系无对应（拉曼惰性/无单矿物键） |
| feldspar | reproduced | silicate | 检出点数 30 |
| olivine | reproduced | olivine | 检出点数 6 |
| organics | excluded |  | SHERLOC 深紫外拉曼类别体系无对应（拉曼惰性/无单矿物键） |
| perchlorate | weak | perchlorate_or_phosphate | 仅模糊键候选（缺次级峰区分），计部分复现 |
| perchlorate or phosphate | unmapped |  | 映射表未覆盖（请人工核对后补充 PAPER_TO_KEYS） |
| phosphate | weak | perchlorate_or_phosphate | 仅模糊键候选（缺次级峰区分），计部分复现 |
| phyllosilicates | reproduced | silicate | 检出点数 30 |
| pyroxene | weak | sulfate_or_pyroxene | 仅模糊键候选（缺次级峰区分），计部分复现 |
| quartz | excluded |  | SHERLOC 深紫外拉曼类别体系无对应（拉曼惰性/无单矿物键） |
| silicate | unmapped |  | 映射表未覆盖（请人工核对后补充 PAPER_TO_KEYS） |

额外候选：unassigned
