# -*- coding: utf-8 -*-
"""
Phase 1 全数据集 sanity：对每个 (atts, meta, h5ad) 输出
A. 类型门控内标志基因归因排名；B. 低表达 TF 捞回；C. 类型内分组基因集预览。
全部打印到 stdout。
"""
import json
from pathlib import Path

import numpy as np
import scanpy as sc
import scipy.sparse as sp

BASE = Path("/root/pe_signature")
P1 = BASE / "results/phase1"
PROC = BASE / "data/processed"

DATASETS = ["GSE282038", "GSE298119", "GSE298602", "GSE240922", "GSE290578", "GSE329173"]

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
TF_LIST = ["TFAP2C", "GATA2", "GATA3", "TEAD4", "NR2F2", "HAND1", "KLF4", "ELF5"]
PREVIEW_SETS = {
    "inflammatory": ["S100A8", "S100A9", "S100A12", "IL1B"],
    "angiogenic_imbalance": ["FLT1", "PGF", "ENG", "ENDOG"],
    "complement": ["C1QA", "C1QB", "C1QC", "C3", "CFH"],
}


def rank0(x):
    order = np.argsort(-x, kind="stable")
    r = np.empty_like(order)
    r[order] = np.arange(len(x))
    return r


def analyze(gse):
    atts = sp.load_npz(P1 / ("atts_%s.npz" % gse))
    meta = np.load(P1 / ("meta_%s.npz" % gse), allow_pickle=True)
    atts_genes = [str(g) for g in meta["genes"]]
    cells = [str(c) for c in meta["cells"]]
    groups = np.array([str(g) for g in meta["groups"]])
    aa = {g: i for i, g in enumerate(atts_genes)}
    adata = sc.read_h5ad(PROC / ("%s_qc.h5ad" % gse))
    Xd = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    obs_names = adata.obs_names.astype(str)
    oidx = {c: i for i, c in enumerate(obs_names)}
    row = np.array([oidx[c] for c in cells])
    agene = {g: i for i, g in enumerate(adata.var_names)}
    res = {"n_cells": int(atts.shape[0]), "n_genes": int(atts.shape[1])}

    type_marks = {}
    for tname, (pos_m, neg_m) in TYPE_GATES.items():
        pm = [g for g in pos_m if g in agene]
        nm = [g for g in neg_m if g in agene]
        if not pm:
            continue
        mask = (Xd[np.ix_(row, [agene[g] for g in pm])] > 0).any(axis=1)
        if nm:
            mask &= ~(Xd[np.ix_(row, [agene[g] for g in nm])] > 0).any(axis=1)
        n_sel = int(mask.sum())
        if n_sel < 30:
            type_marks[tname] = {"n_cells": n_sel, "note": "too few"}
            continue
        rk = rank0(np.asarray(atts[mask].mean(axis=0)).ravel())
        canon = {g: (int(rk[aa[g]]) if g in aa else "not_in_atts") for g in CANON.get(tname, [])}
        type_marks[tname] = {"n_cells": n_sel, "canon_att_rank": canon}
    res["A_marker_attribution_rank"] = type_marks

    rk_att = rank0(np.asarray(atts.mean(axis=0)).ravel())
    rk_expr = rank0(Xd[row].mean(axis=0))
    tf_res = {}
    for tf in TF_LIST:
        if tf in aa and tf in agene:
            ra, re = int(rk_att[aa[tf]]), int(rk_expr[agene[tf]])
            tf_res[tf] = {"att_rank": ra, "expr_rank": re, "gain_x": round((re + 1) / (ra + 1), 1)}
    res["B_lowexpr_TF_rescue"] = tf_res

    prev = {}
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
            for grp in sorted(set(groups)):
                m = mask & (groups == grp)
                if m.sum() > 20:
                    by_group[grp] = round(float(atts[m][:, gi].mean()), 4)
            if by_group:
                entry[sname] = by_group
        if entry:
            prev[tname] = entry
    res["C_within_type_group_preview"] = prev
    return res


def main():
    report = {}
    for gse in DATASETS:
        try:
            report[gse] = analyze(gse)
            print("[%s] ok" % gse, flush=True)
        except Exception as e:
            report[gse] = {"error": "%s: %s" % (type(e).__name__, e)}
            print("[%s] FAIL %s" % (gse, e), flush=True)
    print("SANITY-ALL-BEGIN")
    print(json.dumps(report, indent=1, ensure_ascii=False))
    print("SANITY-ALL-END")


if __name__ == "__main__":
    main()
