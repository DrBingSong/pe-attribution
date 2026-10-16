# -*- coding: utf-8 -*-
"""
Phase 3 v1（全细胞水平）：跨队列 PE 归因共识挖掘
- 差分队列（PE vs control 可配对）：GSE173193 / GSE298119 / GSE298602 / GSE290578
- 投影队列（无内部对照）：GSE282038(EOPE) / GSE329173(severe_PE) / GSE240922(早孕参考)
输出：
  results/phase3/consensus_pe_up.csv / consensus_pe_down.csv（≥3/4 队列一致的基因）
  results/phase3/nmf_programs.csv（每程序 top40 基因 × 4）
  results/phase3/dataset_correlation.csv + heatmap.png
  results/phase3/projection_profiles.csv
  stdout: SUMMARY-JSON
注意：v1 为全细胞水平（含细胞组成混淆），类型分层版待正式注释后（实例）进行。
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

BASE = Path("/f/pe_signature") if Path("/f/pe_signature").exists() else Path("F:/pe_signature")
P1 = BASE / "results" / "phase1"
P0 = BASE / "results" / "phase0"
P3 = Path("D:/zcode/z0927/pe_signature/results/phase3")
P3.mkdir(parents=True, exist_ok=True)
OUT_NPZ = {
    "GSE173193": P0 / "atts_scimilarity_full.npz",
    "GSE282038": P1 / "atts_GSE282038.npz",
    "GSE298119": P1 / "atts_GSE298119.npz",
    "GSE298602": P1 / "atts_GSE298602.npz",
    "GSE240922": P1 / "atts_GSE240922.npz",
    "GSE290578": P1 / "atts_GSE290578.npz",
    "GSE329173": P1 / "atts_GSE329173.npz",
}
META = {
    "GSE173193": P0 / "atts_scimilarity_full_meta.npz",
    "GSE282038": P1 / "meta_GSE282038.npz",
    "GSE298119": P1 / "meta_GSE298119.npz",
    "GSE298602": P1 / "meta_GSE298602.npz",
    "GSE240922": P1 / "meta_GSE240922.npz",
    "GSE290578": P1 / "meta_GSE290578.npz",
    "GSE329173": P1 / "meta_GSE329173.npz",
}
DIFF_DATASETS = ["GSE173193", "GSE298119", "GSE298602", "GSE290578"]
PROJ_DATASETS = ["GSE282038", "GSE329173", "GSE240922"]
N_PROGRAMS = 4
TOP_N = 40


def as2d(a):
    import scipy.sparse as sp
    if sp.issparse(a):
        return np.asarray(a.todense(), dtype=np.float32)
    return np.asarray(a, dtype=np.float32)


def mean_att(gse):
    atts = sp.load_npz(OUT_NPZ[gse])
    meta = np.load(META[gse], allow_pickle=True)
    genes = [str(g) for g in meta["genes"]]
    groups = np.array([str(g) for g in meta["groups"]])
    return atts, groups, genes


def zrow(x):
    s = x.std()
    return (x - x.mean()) / (s if s > 0 else 1.0)


def main():
    diffs, profiles, gene_ref = {}, {}, None

    # ---- 差分数据集：PE mean - control mean（全细胞） ----
    for gse in DIFF_DATASETS:
        try:
            atts, groups, genes = mean_att(gse)
        except Exception as e:
            print("[%s] 跳过（文件损坏/缺失）: %s" % (gse, e), flush=True)
            continue
        if gene_ref is None:
            gene_ref = np.array(genes)
        pe_m = np.asarray(atts[groups == "PE"].mean(axis=0)).ravel()
        ct_m = np.asarray(atts[groups == "control"].mean(axis=0)).ravel()
        if (groups == "PE").sum() == 0 or (groups == "control").sum() == 0:
            print("[%s] 跳过（无 PE/control 分组）" % gse, flush=True)
            continue
        diffs[gse] = zrow(pe_m - ct_m)
        print("[%s] diff 完成" % gse, flush=True)
    diff_used = list(diffs.keys())

    D = np.vstack([diffs[g] for g in diff_used])
    df_D = pd.DataFrame(D, index=diff_used, columns=gene_ref)

    # ---- 跨队列相关 ----
    corr = df_D.T.corr(method="spearman")
    corr.to_csv(P3 / "dataset_correlation.csv")
    print("跨队列 PE 差异签名 Spearman 相关：")
    print(corr.round(3), flush=True)

    # ---- 一致性 Up / Down 基因集 ----
    n_pos = (D > 0).sum(axis=0)
    n_neg = (D < 0).sum(axis=0)
    strength = np.abs(D).mean(axis=0)
    need = max(2, int(np.ceil(0.67 * D.shape[0])))
    up_genes = gene_ref[(n_pos >= need) & (strength > np.quantile(strength, 0.75))]
    down_genes = gene_ref[(n_neg >= need) & (strength > np.quantile(strength, 0.75))]
    pd.DataFrame({"gene": up_genes}).to_csv(P3 / "consensus_pe_up.csv", index=False)
    pd.DataFrame({"gene": down_genes}).to_csv(P3 / "consensus_pe_down.csv", index=False)
    print("consensus PE-up: %d 基因 | PE-down: %d 基因" % (len(up_genes), len(down_genes)), flush=True)
    print("PE-up top30:", list(up_genes[:30]), flush=True)
    print("PE-down top30:", list(down_genes[:30]), flush=True)

    # ---- NMF 程序分解（移到非负：X - min） ----
    from sklearn.decomposition import NMF
    X = D - D.min()
    n_comp = min(N_PROGRAMS, D.shape[0])
    nmf = NMF(n_components=n_comp, init="nndsvd", random_state=0, max_iter=500)
    W = nmf.fit_transform(X)   # datasets × programs
    H = nmf.components_        # programs × genes
    prog_rows = []
    for p in range(n_comp):
        order = np.argsort(-H[p])[:TOP_N]
        prog_rows.append({"program": "P%d" % (p + 1),
                          "top_genes": ", ".join(gene_ref[order][:TOP_N])})
        pd.DataFrame({"gene": gene_ref[order], "weight": H[p][order]}).to_csv(
            P3 / ("nmf_P%d_genes.csv" % (p + 1)), index=False)
    pd.DataFrame(prog_rows).to_csv(P3 / "nmf_programs.csv", index=False)
    pd.DataFrame(W, index=diff_used,
                 columns=["P%d" % (p + 1) for p in range(n_comp)]).to_csv(
        P3 / "nmf_dataset_loadings.csv")
    print("NMF 数据集载荷：")
    print(pd.DataFrame(W, index=diff_used,
                       columns=["P%d" % (p + 1) for p in range(n_comp)]).round(2), flush=True)

    # ---- 热图（差异签名 z-score，top 变异基因） ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    var_order = np.argsort(-df_D.var(axis=0).values)[:150]
    fig, ax = plt.subplots(figsize=(14, 4))
    im = ax.imshow(df_D.iloc[:, var_order].values, aspect="auto", cmap="RdBu_r",
                   vmin=-2, vmax=2)
    ax.set_yticks(range(len(DIFF_DATASETS)))
    ax.set_yticklabels(DIFF_DATASETS)
    ax.set_title("PE vs control attribution differential (top150 var genes)")
    fig.colorbar(im, ax=ax, shrink=0.8)
    fig.tight_layout()
    fig.savefig(P3 / "dataset_heatmap.png", dpi=150)

    # ---- 投影队列：全剖面与 consensus 程序/基因集的相关 ----
    up_set = set(up_genes.tolist())
    down_set = set(down_genes.tolist())
    proj_rows = []
    for gse in PROJ_DATASETS:
        try:
            atts, groups, genes = mean_att(gse)
        except Exception as e:
            print("[投影跳过]", gse, e, flush=True)
            continue
        ai = {g: i for i, g in enumerate(genes)}
        up_i = [ai[g] for g in up_genes if g in ai]
        dn_i = [ai[g] for g in down_genes if g in ai]
        prof = np.asarray(atts.mean(axis=0)).ravel()
        zs = zrow(prof)
        up_z = float(np.mean([zs[ai[g]] for g in up_genes if g in ai])) if up_i else np.nan
        dn_z = float(np.mean([zs[ai[g]] for g in down_genes if g in ai])) if dn_i else np.nan
        proj_rows.append({"dataset": gse, "n_cells": atts.shape[0],
                          "PEup_mean_z": round(up_z, 3), "PEdown_mean_z": round(dn_z, 3),
                          "groups": ",".join(sorted(set(groups)))})
        print("[投影] %s groups=%s PEup_z=%.3f PEdown_z=%.3f" %
              (gse, ",".join(sorted(set(groups))), up_z, dn_z), flush=True)
    pd.DataFrame(proj_rows).to_csv(P3 / "projection_profiles.csv", index=False)

    summary = {"n_consensus_up": len(up_genes), "n_consensus_down": len(down_genes),
               "diff_datasets": diff_used, "proj_datasets": PROJ_DATASETS}
    print("SUMMARY-JSON-BEGIN")
    print(json.dumps(summary, indent=1))
    print("SUMMARY-JSON-END")


if __name__ == "__main__":
    main()
