# PE 单细胞数据集 Inventory（v1，2026-09-28）

状态：✅=GSE 已确认可下 | 🔶=研究已确认、GSE 号待核实 | ⏸=数据非 GEO（GitHub/CELLxGENE）| ❓=待查

## A. GEO 已确认（5 个，可直接入队）

| GSE | 组织 | 分组 | 备注 |
|-----|------|------|------|
| GSE173193 | 胎盘 | 2 PE + 2 对照 + 2 GDM + 2 高龄 | ✅ 已处理完毕（Phase 0 试金石） |
| GSE282038 | 胎盘滋养层 | EOPE+FGR | 亚型分析关键 |
| GSE298119 | 胎盘 | 足月 PE vs 血压正常 | LOPE |
| GSE240922 | 蜕膜 + 绒毛 | （PE 相关，需确认分组） | 检索确认 |
| GSE266604 | 胎盘 Visium FFPE 空间转录组 | PE | 空间验证用（非 scRNA） |

## B. 研究确认、GSE 号待核实（6 个）

| 研究 | 组织 | 分组 | 状态 |
|------|------|------|------|
| Zhang T. 2021, Mol Genet Genomic Med (PMID 34212522, PMC8404237) | 胎盘（剖宫产定点取样） | PE vs 健康 | ❓ GSE 未索引，查 PMC 数据声明；可能未公开 |
| Yang J. 2023, Front Immunol 1142273 | 胎盘 6 + 蜕膜 7 | LOPE 3+3 vs 正常 | ⏸ **数据在 GitHub**（JustMoveOnnn/preeclampsia，预处理矩阵，需确认是否 raw count） |
| Luo F. 2023, Front Immunol 1234577 (PMID 37854606) | 胎盘蜕膜 | PE，HLA-F 免疫耐受 | ❓ 待查 GSE |
| Wei 2025 | 外周血 + 胎盘 | EOPE/LOPE | ❓ 亚型分析关键，待查 |
| Zhou 2024 (Nature 子刊, PBMC 免疫变异, 36 引) | PBMC | PE | ❓ 待查 |
| Yu 2025 (PMC) | 胎盘免疫微环境 | PE | ❓ 待查 |

## C. 特殊用途（非 GEO 主线）

| 来源 | 用途 |
|------|------|
| Admati 2023 Med（两分子亚型，Zeisel） | 亚型对话；数据声明待核（Cell/ScienceDirect） |
| Science Advances severe PE 图谱（sc+spatial） | 空间验证 |
| Vento-Tormo 2018 早孕参考（CELLxGENE） | 参考注释/映射 |
| Gomez-Lopez 预印本（蜕膜+外周血） | 备选 |
| Campbell 2023 bulk+sc | bulk 反卷积验证 |

## 核实方法（实例开机后第一件事）

NCBI eutils（实例上 NCBI 快）：
```
esearch.fcgi?db=gds&term=preeclampsia+AND+single+cell+RNA-seq&retmax=200
esummary.fcgi?db=gds&id=...   → 标题/样本数/平台，系统性补全 B 类 GSE 号
```
入队标准（方案 §4.2）：≥2 PE + ≥2 对照；10x/Smart-seq2；原始 count 可下载；孕周/深度元数据；PE 分型可判定。
