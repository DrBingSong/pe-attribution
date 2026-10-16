# -*- coding: utf-8 -*-
"""
Phase 3 v2：类型分层跨队列共识挖掘（在正式注释后运行）。
- 差分队列（PE vs control）：GSE173193 / GSE298119 / GSE298602 / GSE290578
- 每个粗类型：类型内 PE mean - control mean → 数据集 × 基因差分矩阵
- 输出：每类型跨队列 Spearman（对比 v1 全细胞 ≈0）、共识基因（≥3/4 队列）、NMF 程序
- 滋养层/免疫再做 fine 层（CTB/STB/EVT；NK/T/Myeloid）
输出：/root/pe_signature/results/phase3_v2/ 下 CSV/JSON + stdout V2-SUMMARY
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.decomposition import NMF

BASE = Path("/root/pe_signature")
P1 = BASE / "results/phase1"
ANN = BASE / "results/annotation/all_obs_annotation.csv"
OUT = BASE / "results/phase3_v2"
OUT.mkdir(parents=True, exist_ok=True)

DIFF = ["GSE173193", "GSE298119", "GSE298602", "GSE290578"]  # v2.1: 4 队列
P0 = BASE / "results/phase0"
ATT = {g: P1 / ("atts_%s.npz" % g) for g in DIFF}
ATT["GSE173193"] = P0 / "atts_scimilarity_full.npz"
META = {g: P1 / ("meta_%s.npz" % g) for g in DIFF}
META["GSE173193"] = P0 / "atts_scimilarity_full_meta.npz"
TOP_N = 40
COARSE_TYPES = ["Trophoblast", "Immune", "Endothelial", "Stromal"]
FINE_TYPES = {"Trophoblast": ["CTB", "STB", "EVT"], "Immune": ["NK", "T", "Myeloid"]}


def rank0(x):
    order = np.argsort(-x, kind="stable")
    r = np.empty_like(order)
    r[order] = np.arange(len(x))
    return r


def spear(a, b):
    ra, rb = rank0(a), rank0(b)
    ra, rb = ra - ra.mean(), rb - rb.mean()
    den = np.sqrt((ra**2).sum() * (rb**2).sum())
    return float((ra * rb).sum() / den) if den > 0 else 0.0


def mine_type(type_label, per_ds_diff, gene_ref):
    """per_ds_diff: {gse: diff vector}; 返回相关矩阵/共识基因/NMF。"""
    names = list(per_ds_diff)
    D = np.vstack([per_ds_diff[g] for g in names])
    corr = {}
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            corr["%s_vs_%s" % (names[i], names[j])] = round(spear(D[i], D[j]), 3)
    need = max(2, int(np.ceil(0.67 * len(names))))
    n_pos = (D > 0).sum(axis=0)
    n_neg = (D < 0).sum(axis=0)
    strength = np.abs(D).mean(axis=0)
    up = gene_ref[(n_pos >= need) & (strength > np.quantile(strength, 0.75))]
    down = gene_ref[(n_neg >= need) & (strength > np.quantile(strength, 0.75))]
    n_comp = min(3, len(names))
    X = D - D.min()
    nmf = NMF(n_components=n_comp, init="nndsvd", random_state=0, max_iter=500)
    W = nmf.fit_transform(X)
    H = nmf.components_
    progs = {}
    for p in range(n_comp):
        order = np.argsort(-H[p])[:TOP_N]
        progs["P%d" % (p + 1)] = [str(g) for g in gene_ref[order]]
        pd.DataFrame({"gene": gene_ref[order], "weight": H[p][order]}).to_csv(
            OUT / ("%s_P%d_genes.csv" % (type_label, p + 1)), index=False)
    pd.DataFrame(W, index=names, columns=list(progs)).to_csv(
        OUT / ("%s_dataset_loadings.csv" % type_label))
    pd.DataFrame({"gene": up}).to_csv(OUT / ("%s_consensus_up.csv" % type_label), index=False)
    pd.DataFrame({"gene": down}).to_csv(OUT / ("%s_consensus_down.csv" % type_label), index=False)
    return {"datasets": names, "pairwise_spearman": corr,
            "n_consensus_up": int(len(up)), "n_consensus_down": int(len(down)),
            "mean_abs_spearman": round(float(np.mean(np.abs(list(corr.values())))), 3),
            "nmf_loadings": {n: {k: round(float(v), 1) for k, v in zip(progs, W[i])}
                             for i, n in enumerate(names)}}


def main():
    ann = pd.read_csv(ANN)
    ann_map = {str(c): (str(a), str(f)) for c, a, f in
               zip(ann["cell"], ann["coarse"], ann["fine"])}
    out = {}
    for level, types in [("coarse", COARSE_TYPES)] + [("fine", FINE_TYPES[k]) for k in FINE_TYPES]:
        for t in types:
            per_ds = {}
            gene_ref = None
            for gse in DIFF:
                try:
                    atts = sp.load_npz(ATT[gse])
                    meta = np.load(META[gse], allow_pickle=True)
                except Exception as e:
                    print("[%s|%s] 跳过 %s: %s" % (level, t, gse, e), flush=True)
                    continue
                cells = [str(c) for c in meta["cells"]]
                groups = np.array([str(g) for g in meta["groups"]])
                groups = np.where(groups == "PE_preterm", "PE", groups)  # v2.1: 纳入 GSE290578 早发 PE
                types_col = [ann_map.get(c, ("NA", "NA"))[0 if level == "coarse" else 1] for c in cells]
                tarr = np.array(types_col)
                m_pe = np.where((tarr == t) & (groups == "PE"))[0]
                m_ct = np.where((tarr == t) & (groups == "control"))[0]
                if m_pe.size < 50 or m_ct.size < 50:
                    print("[%s|%s] 跳过 %s（PE %d / CT %d 细胞不足）" % (level, t, gse, m_pe.size, m_ct.size), flush=True)
                    continue
                try:
                    pe_m = np.asarray(atts[m_pe].mean(axis=0)).ravel()
                    ct_m = np.asarray(atts[m_ct].mean(axis=0)).ravel()
                except ZeroDivisionError:
                    print("[%s|%s] 跳过 %s（稀疏均值除零）" % (level, t, gse), flush=True)
                    continue
                sd = pe_m.std()
                per_ds[gse] = (pe_m - ct_m) / (sd if sd > 0 else 1.0)
                if gene_ref is None:
                    gene_ref = np.array([str(g) for g in meta["genes"]])
                print("[%s|%s] %s: PE %d / CT %d" % (level, t, gse, m_pe.size, m_ct.size), flush=True)
            if len(per_ds) >= 3 and gene_ref is not None:
                res = mine_type("%s_%s" % (level, t), per_ds, gene_ref)
                out["%s_%s" % (level, t)] = res
            else:
                print("[%s|%s] 数据集不足（%d），跳过挖掘" % (level, t, len(per_ds)), flush=True)

    print("V2-SUMMARY-BEGIN")
    print(json.dumps(out, indent=1, ensure_ascii=False))
    print("V2-SUMMARY-END")


if __name__ == "__main__":
    main()
