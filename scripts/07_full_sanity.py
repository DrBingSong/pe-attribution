# -*- coding: utf-8 -*-
"""
Phase 0 全量 sanity：稀疏全量归因 + 全部细胞的类型门控分析与分组预览。
关键升级：组间比较在"同一细胞类型内部"做，剥离细胞组成混淆。
输出到 stdout（SANITY-JSON）。
"""
import json
import numpy as np
import scanpy as sc
import scipy.sparse as sp

BASE = "/root/pe_signature"
QC = BASE + "/data/processed/GSE173193_qc.h5ad"
ATT = BASE + "/results/phase0/atts_scimilarity_full.npz"
META = BASE + "/results/phase0/atts_scimilarity_full_meta.npz"

TYPE_GATES = {
    "trophoblast": (["KRT7"], ["PTPRC", "PECAM1"]),
    "immune": (["PTPRC"], ["KRT7"]),
    "endothelial": (["PECAM1"], ["KRT7", "PTPRC"]),
    "stromal": (["COL1A1"], ["KRT7", "PTPRC", "PECAM1"]),
}
CANON = {
    "trophoblast": ["KRT7", "CGA", "CGB3", "PEG10", "TFAP2C", "GATA3", "KRT18", "HLA-G"],
    "immune": ["PTPRC", "LYZ", "CD3D", "NKG7", "CD68", "CST3"],
    "endothelial": ["PECAM1", "CDH5", "VWF", "CLDN5"],
    "stromal": ["COL1A1", "DCN", "LUM"],
}
TF_LIST = ["TFAP2C", "GATA2", "GATA3", "TEAD4", "NR2F2", "HAND1", "MSX2", "DLX3", "ELF5", "KLF4"]
PREVIEW_SETS = {
    "inflammatory": ["S100A8", "S100A9", "S100A12", "IL1B", "TNF"],
    "hypoxia": ["VEGFA", "HK2", "PDK1", "LDHA", "NDRG1"],
    "angiogenic_imbalance": ["FLT1", "PGF", "ENG", "ENDOG"],
    "complement": ["C1QA", "C1QB", "C1QC", "C3", "CFH"],
}


def rank0(x):
    order = np.argsort(-x, kind="stable")
    r = np.empty_like(order)
    r[order] = np.arange(len(x))
    return r


def main():
    adata = sc.read_h5ad(QC)
    Xd = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    obs_names = adata.obs_names.astype(str)
    agene = {g: i for i, g in enumerate(adata.var_names)}
    obs = adata.obs

    atts = sp.load_npz(ATT)
    meta = np.load(META, allow_pickle=True)
    atts_genes = [str(g) for g in meta["genes"]]
    cells = [str(c) for c in meta["cells"]]
    groups = np.array([str(g) for g in meta["groups"]])
    samples = np.array([str(s) for s in meta["samples"]])
    aa = {g: i for i, g in enumerate(atts_genes)}
    order = {c: i for i, c in enumerate(obs_names)}
    row = np.array([order[c] for c in cells])
    sub = obs.iloc[row]
    print("cells=%d att_genes=%d" % (atts.shape[0], atts.shape[1]), flush=True)

    res = {"n_cells": int(atts.shape[0]), "n_genes": int(atts.shape[1])}

    # A. 类型内 marker 归因排名（全细胞）
    type_marks = {}
    for tname, (pos_m, neg_m) in TYPE_GATES.items():
        pm = [g for g in pos_m if g in agene]
        nm = [g for g in neg_m if g in agene]
        mask = (Xd[np.ix_(row, [agene[g] for g in pm])] > 0).any(axis=1)
        if nm:
            mask &= ~(Xd[np.ix_(row, [agene[g] for g in nm])] > 0).any(axis=1)
        n_sel = int(mask.sum())
        rk = rank0(np.asarray(atts[mask].mean(axis=0)).ravel())
        canon = {g: (int(rk[aa[g]]) if g in aa else "not_in_atts") for g in CANON.get(tname, [])}
        type_marks[tname] = {"n_cells": n_sel, "canon_att_rank": canon}
    res["A_marker_attribution_rank"] = type_marks

    # B. 低表达 TF 捞回
    tf_res = {}
    rk_att = rank0(np.asarray(atts.mean(axis=0)).ravel())
    rk_expr = rank0(Xd[row].mean(axis=0))
    for tf in TF_LIST:
        if tf in aa and tf in agene:
            ra, re = int(rk_att[aa[tf]]), int(rk_expr[agene[tf]])
            tf_res[tf] = {"att_rank": ra, "expr_rank": re, "gain_x": round((re + 1) / (ra + 1), 1)}
    res["B_lowexpr_TF_rescue"] = tf_res

    # C. 同类型内 PE vs control 的基因集归因（关键：剥离组成混淆）
    type_group = {}
    for tname in TYPE_GATES:
        pos_m, neg_m = TYPE_GATES[tname]
        pm = [g for g in pos_m if g in agene]
        nm = [g for g in neg_m if g in agene]
        mask = (Xd[np.ix_(row, [agene[g] for g in pm])] > 0).any(axis=1)
        if nm:
            mask &= ~(Xd[np.ix_(row, [agene[g] for g in nm])] > 0).any(axis=1)
        entry = {}
        for sname, gs in PREVIEW_SETS.items():
            ok = [g for g in gs if g in aa]
            if len(ok) < 2:
                continue
            gi = np.array([aa[g] for g in ok])
            by_group = {}
            for grp in ["PE", "control", "GDM", "elderly"]:
                m = (mask & (groups == grp))
                if m.sum() > 20:
                    by_group[grp] = (round(float(atts[m][:, gi].mean()), 4), int(m.sum()))
            entry[sname] = by_group
        type_group[tname] = entry
    res["C_within_type_group_preview"] = type_group
    res["sample_counts"] = {s: int((samples == s).sum()) for s in sorted(set(samples))}

    print("SANITY-JSON-BEGIN")
    print(json.dumps(res, indent=1, ensure_ascii=False))
    print("SANITY-JSON-END")


if __name__ == "__main__":
    main()
