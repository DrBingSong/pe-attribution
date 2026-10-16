# -*- coding: utf-8 -*-
"""
Phase 0 sanity 分析 v2：归因排名用 atts 自带基因表；表达/门控从 h5ad 取。
输出全部到 stdout（SANITY-JSON），不写文件。
"""
import json
import numpy as np
import scanpy as sc

QC = "/root/pe_signature/data/processed/GSE173193_qc.h5ad"
PHASE0 = "/root/pe_signature/results/phase0"
MODELS = ["scimilarity", "scfoundation"]

TYPE_GATES = {
    "trophoblast": (["KRT7"], ["PTPRC", "PECAM1"]),
    "immune": (["PTPRC"], ["KRT7"]),
    "endothelial": (["PECAM1"], ["KRT7", "PTPRC"]),
    "stromal": (["COL1A1"], ["KRT7", "PTPRC", "PECAM1"]),
}
CANON = {
    "trophoblast": ["KRT7", "CGA", "CGB3", "PEG10", "TFAP2C", "GATA3", "KRT18"],
    "immune": ["PTPRC", "LYZ", "CD3D", "NKG7", "CD68"],
    "endothelial": ["PECAM1", "CDH5", "VWF", "CLDN5"],
    "stromal": ["COL1A1", "DCN", "LUM"],
}
TF_LIST = ["TFAP2C", "GATA2", "GATA3", "TEAD4", "NR2F2", "HAND1", "MSX2", "DLX3"]
PREVIEW_SETS = {
    "inflammatory": ["S100A8", "S100A9", "S100A12", "IL1B", "TNF"],
    "hypoxia": ["VEGFA", "HK2", "PDK1", "LDHA", "NDRG1"],
    "angiogenic_imbalance": ["FLT1", "PGF", "ENG", "ENDOG"],
}


def rank0(x):  # 返回每个元素的排名（0 = 最大）
    order = np.argsort(-x, kind="stable")
    r = np.empty_like(order)
    r[order] = np.arange(len(x))
    return r


def main():
    adata = sc.read_h5ad(QC)
    Xd = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    obs_names = adata.obs_names.astype(str)
    agene_idx = {g: i for i, g in enumerate(adata.var_names)}
    obs = adata.obs

    report = {}
    for m in MODELS:
        try:
            z = np.load("%s/atts_%s.npz" % (PHASE0, m), allow_pickle=True)
        except Exception as e:
            report[m] = {"error": str(e)}
            continue
        atts = np.asarray(z["atts"], dtype=np.float32)
        atts_genes = [str(g) for g in z["genes"]]
        cells = [str(c) for c in z["cells"]]
        aa = {g: i for i, g in enumerate(atts_genes)}
        row = np.array([dict((c, i) for i, c in enumerate(obs_names))[c] for c in cells])
        sub = adata.obs.iloc[row]
        res = {"n_cells": atts.shape[0], "n_att_genes": atts.shape[1]}

        # A. 类型内 marker 归因排名
        type_marks = {}
        for tname, (pos_m, neg_m) in TYPE_GATES.items():
            pm = [g for g in pos_m if g in agene_idx]
            nm = [g for g in neg_m if g in agene_idx]
            if not pm:
                continue
            mask = (Xd[np.ix_(row, [agene_idx[g] for g in pm])] > 0).any(axis=1)
            if nm:
                mask &= ~(Xd[np.ix_(row, [agene_idx[g] for g in nm])] > 0).any(axis=1)
            n_sel = int(mask.sum())
            if n_sel < 30:
                type_marks[tname] = {"n_cells": n_sel, "note": "too few"}
                continue
            rk = rank0(atts[mask].mean(axis=0))
            canon = {}
            for g in CANON.get(tname, []):
                canon[g] = (int(rk[aa[g]]) if g in aa else "not_in_atts")
            type_marks[tname] = {"n_cells": n_sel, "canon_att_rank": canon}
        res["A_marker_attribution_rank"] = type_marks

        # B. 低表达 TF 捞回（全体细胞：归因排名 vs 表达排名）
        tf_res = {}
        rk_att = rank0(atts.mean(axis=0))
        expr_mean = Xd[row].mean(axis=0)
        rk_expr_full = rank0(expr_mean)
        for tf in TF_LIST:
            in_a, in_e = tf in aa, tf in agene_idx
            if in_a and in_e:
                ra, re = int(rk_att[aa[tf]]), int(rk_expr_full[agene_idx[tf]])
                tf_res[tf] = {"att_rank": ra, "expr_rank": re,
                              "gain_x": round((re + 1) / (ra + 1), 1)}
            else:
                tf_res[tf] = "in_atts=%s in_h5ad=%s" % (in_a, in_e)
        res["B_lowexpr_TF_rescue"] = tf_res

        # C. 分组基因集归因预览（PE vs 对照）
        prev = {}
        for sname, gs in PREVIEW_SETS.items():
            ok = [g for g in gs if g in aa]
            if not ok:
                continue
            gi = np.array([aa[g] for g in ok])
            by_group = {}
            for grp in ["PE", "control", "GDM", "elderly"]:
                mrow = (sub["group"].astype(str) == grp).to_numpy()
                if mrow.sum() > 0:
                    by_group[grp] = round(float(atts[mrow][:, gi].mean()), 4)
            prev[sname] = by_group
        res["C_group_preview"] = prev
        report[m] = res

    print("SANITY-JSON-BEGIN")
    print(json.dumps(report, indent=1, ensure_ascii=False))
    print("SANITY-JSON-END")


if __name__ == "__main__":
    main()
