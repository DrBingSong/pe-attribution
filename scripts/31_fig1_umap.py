# -*- coding: utf-8 -*-
"""
F1 注释 UMAP（实例）：7 数据集，每数据集抽样 ≤25k 细胞
normalize → HVG → PCA → neighbors → UMAP，按 coarse 着色 + group 分面展示
输出：/root/pe_signature/results/figures/fig1_annotation_umap.png + 每数据集 png
"""
import random
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

BASE = Path("/root/pe_signature")
PROC = BASE / "data" / "processed"
ANN = BASE / "results" / "annotation"
OUT = BASE / "results" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

GSES = ["GSE173193", "GSE282038", "GSE298119", "GSE298602", "GSE240922", "GSE290578", "GSE329173"]
MAX_CELLS = 25000
COLORS = {"Trophoblast": "#b2182b", "Immune": "#2166ac", "Endothelial": "#1b7837",
          "Stromal": "#762a83", "Erythroid": "#e08214"}


def umap_one(gse):
    adata = sc.read_h5ad(PROC / ("%s_qc.h5ad" % gse))
    ann = pd.read_csv(ANN / ("%s_obs.csv" % gse))
    amap = dict(zip(ann.cell.astype(str), ann.coarse))
    adata.obs["coarse"] = [amap.get(str(x), "NA") for x in adata.obs_names]
    if adata.n_obs > MAX_CELLS:
        rng = np.random.default_rng(0)
        keep = rng.choice(adata.n_obs, MAX_CELLS, replace=False)
        adata = adata[keep].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2500, flavor="seurat")
    sc.pp.pca(adata, n_comps=30, mask_var="highly_variable")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata)
    adata.write(PROC / ("%s_umap25k.h5ad" % gse))
    sc.pl.umap(adata, color="coarse", save="_%s.png" % gse, size=8, palette=COLORS, title=gse)
    print("[%s] UMAP done (%d cells)" % (gse, adata.n_obs), flush=True)
    return gse


def main():
    import seaborn
    sc.settings.figdir = OUT
    for gse in GSES:
        try:
            umap_one(gse)
        except Exception as e:
            print("[%s] FAIL %s: %s" % (gse, type(e).__name__, e), flush=True)
    # 拼图
    fig, axes = plt.subplots(2, 4, figsize=(18, 9))
    for ax, gse in zip(axes.flat, GSES + [None]):
        if gse is None:
            ax.axis("off")
            import matplotlib.patches as mpatches
            ax.legend(handles=[mpatches.Patch(color=c, label=k) for k, c in COLORS.items()],
                      loc="center", fontsize=13, frameon=False)
            continue
        f = OUT / ("umap_%s.png" % gse)
        if f.exists():
            ax.imshow(plt.imread(f))
            ax.set_title(gse, fontsize=11)
        ax.axis("off")
    fig.suptitle("F1  Cell-type annotation across 7 PE cohorts (≤25k cells/dataset)", fontsize=14)
    fig.tight_layout()
    fig.savefig(OUT / "fig1_annotation_umap.png", dpi=200, bbox_inches="tight")
    print("F1-DONE", flush=True)


if __name__ == "__main__":
    main()
