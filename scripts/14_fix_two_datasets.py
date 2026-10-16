# -*- coding: utf-8 -*-
"""
修复两个特殊数据集并补齐 QC h5ad：
A. GSE240922：10x h5 格式（filtered_feature_bc_matrix.h5），30 样本 = 早孕参考（reference_FT）
B. GSE290578：mtx 目录，逐样本体检（找出 barcodes/矩阵不匹配的坏样本并跳过）
"""
import json
from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc

BASE = Path("/root/pe_signature")
PROC = BASE / "data" / "processed"
MIN_GENES, MAX_GENES, MAX_PCT_MT, MIN_CELLS = 200, 6000, 20.0, 3

H5_OUT = PROC / "GSE240922_qc.h5ad"
MTX_OUT = PROC / "GSE290578_qc.h5ad"


def qc(adata):
    adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    sc.pp.filter_cells(adata, min_genes=MIN_GENES)
    sc.pp.filter_cells(adata, max_genes=MAX_GENES)
    adata = adata[adata.obs["pct_counts_mt"] <= MAX_PCT_MT, :].copy()
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS)
    return adata


def part_a_240922():
    ext = BASE / "data/gse/GSE240922/extracted"
    h5s = sorted(ext.glob("*_filtered_feature_bc_matrix.h5"))
    if not h5s:
        raise RuntimeError("no h5 files")
    parts = []
    for h in h5s:
        name = h.name[: -len("_filtered_feature_bc_matrix.h5")]
        a = sc.read_10x_h5(h)
        a.var_names_make_unique()
        a.obs_names = [name + "-" + str(x) for x in range(a.n_obs)]
        a.obs["sample"] = name
        a.obs["group"] = "reference_FT"
        a.obs["region"] = "decidua" if "_Dec_" in name else ("villi" if "_Vil_" in name else "NA")
        parts.append(a)
        print("  h5:", name, a.n_obs, "细胞", flush=True)
    adx = ad.concat(parts, index_unique="-")
    adx = qc(adx)
    adx.write_h5ad(H5_OUT)
    print("  A 完成: %d 细胞 × %d 基因 -> %s" % (adx.n_obs, adx.n_vars, H5_OUT), flush=True)


def part_b_290578():
    ext = BASE / "data/gse/GSE290578/extracted"
    parts, bad = [], []
    for sdir in sorted(p for p in ext.iterdir() if p.is_dir()):
        try:
            a = sc.read_10x_mtx(sdir, var_names="gene_symbols", make_unique=True, cache=True)
            a.var_names_make_unique()
            a.obs_names_make_unique()
            if a.n_obs != len(a.obs_names):
                raise ValueError("shape mismatch %d vs %d" % (a.n_obs, len(a.obs_names)))
        except Exception as e:
            bad.append((sdir.name, "%s: %s" % (type(e).__name__, e)))
            print("  坏样本跳过:", sdir.name, e, flush=True)
            continue
        g = "control" if "_Norm_" in sdir.name else ("PE_preterm" if "_Pt_" in sdir.name else "NA")
        region = "MFI" if "_MF" in sdir.name else ("deep" if "_P" in sdir.name else "NA")
        a.obs["sample"] = sdir.name
        a.obs["group"] = g
        a.obs["region"] = region
        parts.append(a)
        print("  ok:", sdir.name, a.n_obs, "细胞 group=%s" % g, flush=True)
    if len(parts) < 10:
        raise RuntimeError("可用样本过少: %d" % len(parts))
    adx = ad.concat(parts, index_unique="-")
    adx = qc(adx)
    adx.write_h5ad(MTX_OUT)
    print("  B 完成: %d 细胞 × %d 基因 -> %s" % (adx.n_obs, adx.n_vars, MTX_OUT), flush=True)
    print("  坏样本清单:", json.dumps(bad), flush=True)


if __name__ == "__main__":
    summary = {}
    for name, fn in [("A_240922", part_a_240922), ("B_290578", part_b_290578)]:
        try:
            fn()
            summary[name] = "ok"
        except Exception as e:
            summary[name] = "FAIL: %s: %s" % (type(e).__name__, e)
            print(summary[name], flush=True)
    print("FIX14-DONE", json.dumps(summary), flush=True)
