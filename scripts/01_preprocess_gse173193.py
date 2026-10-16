# -*- coding: utf-8 -*-
"""
GSE173193 预处理（Phase 0 试金石数据集）
样本分组（来自 GEO SOFT 元数据）：
  C1/C2 = control | G1/G2 = GDM | P1/P2 = PE | GL1/GL2 = elderly(高龄)
核心对比：PE (P1+P2) vs control (C1+C2)
输入：F:/pe_signature/raw/extracted/<sample>/  (10x MTX 三件套)
输出：
  F:/pe_signature/processed/GSE173193_qc.h5ad          （过滤后原始 counts）
  F:/pe_signature/processed/GSE173193_normalized.h5ad  （lognorm + HVG + PCA/UMAP）
  results/GSE173193_qc_report.txt / _per_sample.csv / _marker_sanity.csv
"""
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

RAW = Path("F:/pe_signature/raw/extracted")
OUT = Path("F:/pe_signature/processed")
RES = Path("D:/zcode/z0927/pe_signature/results")
OUT.mkdir(parents=True, exist_ok=True)
RES.mkdir(parents=True, exist_ok=True)

GROUPS = {"C": "control", "G": "GDM", "P": "PE", "GL": "elderly"}
# QC 阈值（胎盘组织常用，Phase 1 做敏感性分析时再调）
MIN_GENES, MAX_GENES, MAX_PCT_MT = 200, 6000, 20.0
MIN_CELLS_PER_GENE = 3
N_HVG, N_PCS = 3000, 30

sc.settings.verbosity = 1
sc.settings.figdir = RES


def group_of(sample_dir: str) -> str:
    m = re.match(r"[A-Z]+", sample_dir.split("_")[-1])
    return GROUPS[m.group()]


def main() -> None:
    sample_dirs = sorted(d for d in RAW.iterdir() if d.is_dir())
    if not sample_dirs:
        raise SystemExit(f"未找到样本目录，请先解压 RAW.tar 到 {RAW}")
    print(f"发现 {len(sample_dirs)} 个样本: {[d.name for d in sample_dirs]}")

    parts = []
    for d in sample_dirs:
        a = sc.read_10x_mtx(d, var_names="gene_symbols", make_unique=True, cache=True)
        a.var_names_make_unique()
        a.obs["sample"] = d.name
        a.obs["group"] = group_of(d.name)
        parts.append(a)
        print(f"  {d.name}: {a.n_obs} 细胞")

    adata = ad.concat(parts, index_unique="-")
    del parts
    print(f"合并后: {adata.n_obs} 细胞 × {adata.n_vars} 基因")

    # ---- QC ----
    adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL", "rps", "rpl"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True)

    per_sample = adata.obs.groupby("sample", observed=True).agg(
        group=("group", "first"),
        n_cells=("n_genes_by_counts", "size"),
        median_genes=("n_genes_by_counts", "median"),
        median_counts=("total_counts", "median"),
        median_pct_mt=("pct_counts_mt", "median"),
    )
    per_sample.to_csv(RES / "GSE173193_per_sample.csv")

    n0 = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=MIN_GENES)
    sc.pp.filter_cells(adata, max_genes=MAX_GENES)
    adata = adata[adata.obs["pct_counts_mt"] <= MAX_PCT_MT, :].copy()
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS_PER_GENE)
    print(f"QC: {n0} -> {adata.n_obs} 细胞, {adata.n_vars} 基因保留")

    adata.write(OUT / "GSE173193_qc.h5ad")  # 原始 counts，供归因阶段使用

    # ---- 归一化 + 降维（Phase 0 注释/聚类底图） ----
    adata_n = adata.copy()
    sc.pp.normalize_total(adata_n, target_sum=1e4)
    sc.pp.log1p(adata_n)
    adata_n.raw = adata_n
    sc.pp.highly_variable_genes(adata_n, n_top_genes=N_HVG, flavor="seurat", batch_key="sample")
    sc.pp.pca(adata_n, n_comps=N_PCS, use_highly_variable=True)
    sc.pp.neighbors(adata_n, n_neighbors=15, n_pcs=N_PCS)
    sc.tl.umap(adata_n)
    sc.tl.leiden(adata_n, resolution=1.0, key_added="leiden")
    adata_n.write(OUT / "GSE173193_normalized.h5ad")

    # ---- marker sanity（未注释前的粗检验） ----
    markers = ["KRT7", "CGA", "CGB3", "HLA-G", "PEG10", "PECAM1", "CDH5", "ACTA2",
               "NCAM1", "PTPRC", "LYZ", "CD3D", "COL1A1", "IGF1"]
    found = [m for m in markers if m in adata_n.raw.var_names]
    X = adata_n.raw[:, found].X
    expr = pd.DataFrame(np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X),
                        columns=found)
    expr["group"] = adata_n.obs["group"].to_numpy()
    sanity = expr.groupby("group", observed=True).agg(["mean", lambda s: (s > 0).mean()])
    sanity.columns = [f"{a}|{b}" for a, b in sanity.columns]
    sanity.to_csv(RES / "GSE173193_marker_sanity.csv")

    # ---- 报告 ----
    report = (
        f"样本数: {len(sample_dirs)}  合并: {n0} 细胞\n"
        f"QC 阈值: genes {MIN_GENES}-{MAX_GENES}, pct_mt<={MAX_PCT_MT}, gene min_cells={MIN_CELLS_PER_GENE}\n"
        f"QC 后: {adata.n_obs} 细胞 × {adata.n_vars} 基因\n\n"
        + per_sample.to_string()
        + "\n\n各组细胞数:\n"
        + adata_n.obs["group"].value_counts().to_string()
        + "\n\nMarker sanity (各组平均表达|阳性率):\n"
        + sanity.round(3).to_string()
    )
    (RES / "GSE173193_qc_report.txt").write_text(report, encoding="utf-8")

    sc.pl.umap(adata_n, color=["group", "leiden"], ncols=1, save="_gse173193_overview.png")
    print("完成。输出见", OUT, "与", RES)


if __name__ == "__main__":
    main()
