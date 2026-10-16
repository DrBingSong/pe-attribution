# -*- coding: utf-8 -*-
"""
25 表达量 vs 归因头对头（论文核心图数据）：
对每个差分队列 × 每个粗细胞类型：
  A. 归因法差异签名：类型内 PE mean(atts) - CT mean(atts)，取 |z| top200 基因
  B. 表达法差异签名：类型内 PE mean(expr) - CT mean(expr)（lognorm），t 统计量 top200
  C. 跨队列两两重叠率（Jaccard）：归因法 vs 表达法
输出：results/headtohead/overlap_summary.csv + per-type 交集基因 CSV；stdout HEAD-JSON
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

BASE = Path("/root/pe_signature")
P1 = BASE / "results/phase1"
P0 = BASE / "results/phase0"
ANN = BASE / "results/annotation/all_obs_annotation.csv"
OUT = BASE / "results" / "headtohead"
OUT.mkdir(parents=True, exist_ok=True)

DIFF = ["GSE173193", "GSE298119", "GSE298602", "GSE290578"]
ATT = {g: P1 / ("atts_%s.npz" % g) for g in DIFF}
ATT["GSE173193"] = P0 / "atts_scimilarity_full.npz"
META = {g: P1 / ("meta_%s.npz" % g) for g in DIFF}
META["GSE173193"] = P0 / "atts_scimilarity_full_meta.npz"
H5 = {g: BASE / "data/processed" / ("%s_qc.h5ad" % g) for g in DIFF}
TYPES = ["Trophoblast", "Immune", "Endothelial", "Stromal"]
TOP = 200


def welch_t(x1, x2):
    m1, m2 = x1.mean(axis=0), x2.mean(axis=0)
    v1, v2 = x1.var(axis=0, ddof=1), x2.var(axis=0, ddof=1)
    n1, n2 = x1.shape[0], x2.shape[0]
    se = np.sqrt(v1 / n1 + v2 / n2)
    se[se == 0] = 1e-9
    return (m1 - m2) / se


def main():
    ann = pd.read_csv(ANN)
    ann_map = {str(c): str(a) for c, a in zip(ann["cell"], ann["coarse"])}

    sigs = {}  # sigs[method][type][gse] = set(top genes)
    for gse in DIFF:
        atts = sp.load_npz(ATT[gse])
        meta = np.load(META[gse], allow_pickle=True)
        cells = [str(c) for c in meta["cells"]]
        genes = [str(g) for g in meta["genes"]]
        groups = np.array([str(g) for g in meta["groups"]])
        groups = np.where(groups == "PE_preterm", "PE", groups)
        tarr = np.array([ann_map.get(c, "NA") for c in cells])

        adata = sc.read_h5ad(H5[gse])
        adata.obs["coarse"] = [ann_map.get(str(x), "NA") for x in adata.obs_names]
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        edense = adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X)
        egene = {g: i for i, g in enumerate(adata.var_names)}

        for t in TYPES:
            m_pe = np.where((tarr == t) & (groups == "PE"))[0]
            m_ct = np.where((tarr == t) & (groups == "control"))[0]
            if m_pe.size < 50 or m_ct.size < 50:
                continue
            # 归因签名（atts 行序 = meta cells 序）
            pe_a = np.asarray(atts[m_pe].mean(axis=0)).ravel()
            ct_a = np.asarray(atts[m_ct].mean(axis=0)).ravel()
            d_att = pe_a - ct_a
            top_att = np.argsort(-np.abs(d_att))[:TOP]
            # 表达签名：行索引映射到 h5ad obs 行序（与 atts 行序可能不同）
            obs_pos = {c: i for i, c in enumerate(adata.obs_names.astype(str))}
            h_pe = np.array([obs_pos[cells[i]] for i in m_pe])
            h_ct = np.array([obs_pos[cells[i]] for i in m_ct])
            col_idx = np.array([egene.get(g, -1) for g in genes])
            ok_col = col_idx >= 0
            pe_e = edense[np.ix_(h_pe, col_idx[ok_col])]
            ct_e = edense[np.ix_(h_ct, col_idx[ok_col])]
            tv = welch_t(pe_e, ct_e)
            tvals = np.full(len(genes), 0.0)
            tvals[ok_col] = tv
            order = np.argsort(-np.abs(tvals))[:TOP]
            for meth, idxs in [("attribution", set(np.array(genes)[top_att])),
                               ("expression", set(np.array(genes)[order]))]:
                sigs.setdefault(meth, {}).setdefault(t, {})[gse] = idxs
        print("[%s] 签名完成" % gse, flush=True)

    rows = []
    for t in TYPES:
        ds_attr = sigs.get("attribution", {}).get(t, {})
        ds_expr = sigs.get("expression", {}).get(t, {})
        common = sorted(set(ds_attr) & set(ds_expr))
        for i in range(len(common)):
            for j in range(i + 1, len(common)):
                a, b = common[i], common[j]
                ja = len(ds_attr[a] & ds_attr[b]) / len(ds_attr[a] | ds_attr[b])
                je = len(ds_expr[a] & ds_expr[b]) / len(ds_expr[a] | ds_expr[b])
                rows.append({"type": t, "pair": "%s|%s" % (a, b),
                             "attr_jaccard": round(ja, 4), "expr_jaccard": round(je, 4),
                             "ratio": round(ja / je, 2) if je > 0 else None})
        if len(common) >= 2:
            inter = set.intersection(*[ds_attr[g] for g in common])
            pd.DataFrame({"gene": sorted(inter)}).to_csv(
                OUT / ("%s_attr_consensus_top200_intersect.csv" % t), index=False)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "overlap_summary.csv", index=False)
    print(df.to_string(index=False), flush=True)
    agg = df.groupby("type")[["attr_jaccard", "expr_jaccard"]].mean().round(4)
    print("TYPE-MEAN-BEGIN")
    print(agg.to_string())
    print("TYPE-MEAN-END")
    print("HEAD-JSON-BEGIN")
    print(json.dumps({t: {"attr": float(v["attr_jaccard"]), "expr": float(v["expr_jaccard"])}
                      for t, v in agg.iterrows()}, indent=1))
    print("HEAD-JSON-END")


if __name__ == "__main__":
    main()
