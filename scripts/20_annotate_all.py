# -*- coding: utf-8 -*-
"""
正式细胞注释（marker 评分法，粗+细两级）：
- 粗：Trophoblast / Immune / Endothelial / Stromal / Erythroid
- 细：Trophoblast→CTB/STB/EVT；Immune→NK/T/Myeloid
normalize+log 后 sc.tl.score_genes → argmax。
输出：/root/pe_signature/results/annotation/<GSE>_obs.csv + all_obs_annotation.csv
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

BASE = Path("/root/pe_signature")
PROC = BASE / "data" / "processed"
OUT = BASE / "results" / "annotation"
OUT.mkdir(parents=True, exist_ok=True)

GSES = ["GSE173193", "GSE282038", "GSE298119", "GSE298602", "GSE240922", "GSE290578", "GSE329173"]
H5 = {g: PROC / ("%s_qc.h5ad" % g) for g in GSES}
OBS_CSV = {g: OUT / ("%s_obs.csv" % g) for g in GSES}

COARSE = {
    "Trophoblast": ["KRT7", "KRT18", "PEG10", "TFAP2C", "GATA3", "HSD3B1"],
    "Immune": ["PTPRC", "LYZ", "CD74"],
    "Endothelial": ["PECAM1", "CDH5", "CLDN5", "VWF"],
    "Stromal": ["COL1A1", "DCN", "LUM", "ACTA2"],
    "Erythroid": ["HBB", "HBA1", "ALAS2"],
}
FINE = {
    "Trophoblast": {"CTB": ["TP63", "EGFR", "ITGA6", "KRT15"],
                    "STB": ["CGA", "CGB3", "CSH1", "GH2"],
                    "EVT": ["HLA-G", "ITGA5", "MMP2", "NOTCH1"]},
    "Immune": {"NK": ["NCAM1", "NKG7", "GNLY", "KLRD1"],
               "T": ["CD3D", "CD3E", "TRAC"],
               "Myeloid": ["CD68", "LST1", "CD14", "FCGR3A"]},
}


def np_const(n, v):
    return np.full(n, v, dtype=float)


def score_sets(tmp, sets, prefix):
    scores = {}
    for key, genes in sets.items():
        ok = [g for g in genes if g in tmp.var_names]
        if len(ok) < 2:
            scores[key] = np_const(tmp.n_obs, -10.0)
            continue
        name = prefix + key
        sc.tl.score_genes(tmp, ok, score_name=name, random_state=0)
        scores[key] = tmp.obs[name].to_numpy()
    return scores


def annotate(adata):
    tmp = adata.copy()
    sc.pp.normalize_total(tmp, target_sum=1e4)
    sc.pp.log1p(tmp)

    coarse_scores = score_sets(tmp, COARSE, "sC_")
    keys = list(coarse_scores)
    mat = np.vstack([coarse_scores[k] for k in keys]).T
    coarse = [keys[i] for i in mat.argmax(axis=1)]

    fine = ["NA"] * len(coarse)
    for lineage, sub in FINE.items():
        idx = [i for i, c in enumerate(coarse) if c == lineage]
        if not idx:
            continue
        fine_scores = score_sets(tmp, sub, "sF_" + lineage + "_")
        fkeys = list(fine_scores)
        fmat = np.vstack([fine_scores[k] for k in fkeys]).T
        for i in idx:
            fine[i] = fkeys[int(fmat[i].argmax())]
    return coarse, fine


def main():
    summary = {}
    all_parts = []
    for gse in GSES:
        h5 = H5[gse]
        if not h5.exists():
            summary[gse] = "no-h5ad"
            print(summary[gse], gse, flush=True)
            continue
        try:
            adata = sc.read_h5ad(h5)
            coarse, fine = annotate(adata)
            df = pd.DataFrame({"cell": [str(x) for x in adata.obs_names],
                               "coarse": coarse, "fine": fine})
            df.to_csv(OBS_CSV[gse], index=False)
            all_parts.append(df)
            summary[gse] = {"n": len(df), "coarse": df["coarse"].value_counts().to_dict()}
            print("[%s]" % gse, json.dumps(summary[gse]), flush=True)
        except Exception as e:
            summary[gse] = "FAIL: %s: %s" % (type(e).__name__, e)
            print(summary[gse], flush=True)
    if all_parts:
        pd.concat(all_parts).to_csv(OUT / "all_obs_annotation.csv", index=False)
    print("ANNOTATE-DONE", json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
