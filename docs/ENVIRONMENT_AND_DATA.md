# 项目环境与数据说明（pe_signature）

更新：2026-09-27

## 目录约定

- 代码：`D:\zcode\z0927\pe_signature\`（scripts / results / notes）
- 数据：`F:\pe_signature\`（raw / processed / external），venv 也在 `F:\pe_signature\env`
- 本机：Ryzen 5 5600G（6C12T）+ 32GB 内存 + 无独立 N 卡

## Python 环境

- venv：`F:\pe_signature\env`（Python 3.12.14，uv 创建）
- 激活：`F:\pe_signature\env\Scripts\activate`
- pip 镜像：清华镜像对 uv 返回 403，**用阿里云**：`--index-url https://mirrors.aliyun.com/pypi/simple/`
- 已装：scanpy 1.12.4 / anndata 0.13.4 / numpy 2.5.3 / scipy / matplotlib / h5py / pyarrow / leidenalg / igraph
- 待装验证：`sc-signature`（SIGnature 官方包）；其 GPU 归因部分本机无 N 卡跑不了，仅验证安装与 CPU 数据接口

## GSE173193（Phase 0 试金石）

- 8 样本 × 10x MTX 三件套，RAW.tar = 503,316,480 字节
- 分组（GEO SOFT 元数据）：
  | 样本 | 分组 | 年龄 |
  |------|------|------|
  | C1/C2 | control | 29/28 |
  | G1/G2 | GDM | 28/33 |
  | **P1/P2** | **PE** | **30/34** |
  | GL1/GL2 | elderly（高龄） | 39 |
- 核心对比：P1+P2 vs C1+C2；GDM/elderly 组留作敏感性分析（PE 特异 vs 共享程序）
- 下载：NCBI 国际带宽慢，用 8 线程 Range 分段并行 + cat 合并（脚本见 notes 下方命令）
- 下载后解压到 `F:\pe_signature\raw\extracted\<sample>\`，预处理脚本 `scripts/01_preprocess_gse173193.py`
  - QC 阈值：genes 200–6000、pct_mt ≤ 20%、gene min_cells=3
  - 输出：`GSE173193_qc.h5ad`（counts，供归因）、`GSE173193_normalized.h5ad`（lognorm+HVG3k+PCA/UMAP/leiden）、QC 报告 + marker sanity CSV

## Zenodo 预计算覆盖度（重要阴性发现）

- SIGnature 预计算归因分按细胞类型分 6 包（T/NK、神经元胶质、免疫、Other A-E/F-H/I-Z）+ 模型辅助文件（record 17903196）+ 教程数据（17905668）
- **I-Z 包 22 个文件中无 placenta/trophoblast/decidua/pregnancy 任何条目**
- 结论：官方 2200 万细胞图谱未覆盖胎盘 → ① Phase 1 归因必须自算（租 GPU）；② 证实基础模型对妊娠组织欠表征的担忧（写进论文的 gap 论证）；③ 新颖性进一步成立

## 下载实战经验（2026-09-27 实测，GSE173193 480MB）

| 路线 | 实测速度 | 结论 |
|------|---------|------|
| 家用宽带单线程 → NCBI | ~5KB/s/连接 | 不可用（100+ 分钟） |
| 家用宽带 8 线程 Range 并行 | ~2.5-10MB/min 聚合 | 兜底可用，约 1-3 小时 |
| **AutoDL 无卡模式** → NCBI | ~200KB/s 单连接 | **限速，不比家用快，别用** |
| **AutoDL 带卡实例** 单连接 | ~200KB/s | 同样受 NCBI 国际路由限制 |
| **AutoDL 带卡实例 8 线程** | **~2.1MB/s 聚合** | **✅ 实测 480MB 约 7 分钟拉完** |
| AutoDL → 本机 SFTP 回传（6 线程） | ~1.4MB/s | 480MB 约 7.5 分钟 |

**结论**：NCBI→中国 everywhere 都慢，无卡模式限速雪上加霜；带卡实例 + 8 线程并行是实测最优解（总耗时 ~15 分钟，卡费 ~¥0.6）。**以后批量下载数据集时直接在带卡实例上开多线程，再把数据/结果传回。** SFTP 回传 6 线程各 ~230KB/s，中继不是瓶颈。

**踩坑记录**：
1. `pkill -9 -f 'curl.*ncbi'` 会误杀自己的 exec shell（bash -c 命令行含该模式）——用 `pkill -9 curl`（按进程名）；
2. paramiko `exec_command` 的 test 探测命令失败是正常分支，需单独 probe() 不查退出码；
3. 下载校验必须 gate 后续步骤（echo 大小不算校验）；
4. 分段下载的目标路径与 SFTP 回传路径曾撞车覆写——回传走独立子目录；
5. AutoDL 数据盘基线占用 ~5GB（镜像自带 Miniconda），清完自己的文件后 df 显示 10% 是正常的。

## AutoDL 实例历史

- 北京B区 689机（无卡）：已弃用（限速），实例上无遗留文件
- 重庆A区 689机（RTX 4090）：当前实例，`connect.cqa1.seetacloud.com:33510`，临时文件已清空，**Phase 0 前可关机省钱（数据盘保留）**
