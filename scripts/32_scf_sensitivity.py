# -*- coding: utf-8 -*-
"""
scFoundation 敏感性分析（实例 GPU）：GSE298602 抽 2000 细胞
→ scFoundation IG 归因 → 免疫门控内 PE vs control 签名 → 与 SCimilarity 共识免疫签名的 Spearman
输出：results/sensitivity/scf_immune_sig.csv + stdout SENS-JSON
"""
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

BASE = Path("/root/pe_signature")
PROC = BASE / "data" / "processed"
P1 = BASE / "results" / "phase1"
ANN = BASE / "results" / "annotation"
OUT = BASE / "results" / "sensitivity"
OUT.mkdir(parents=True, exist_ok=True)
GSE = "GSE298602"


def as2d(a):
    if sp.issparse(a):
        return np.asarray(a.todense(), dtype=np.float32)
    return np.asarray(a, dtype=np.float32)


def rank0(x):
    order = np.argsort(-x, kind="stable")
    r = np.empty_like(order)
    r[order] = np.arange(len(x))
    return r


def spear(a, b):
    ra, rb = rank0(a) - rank0(a).mean(), rank0(b) - rank0(b).mean()
    den = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / den) if den else 0.0


def main():
    adata = sc.read_h5ad(PROC / ("%s_qc.h5ad" % GSE))
    ann = pd.read_csv(ANN / ("%s_obs.csv" % GSE))
    amap = dict(zip(ann.cell.astype(str), zip(ann.coarse.astype(str), ann.fine.astype(str))))
    # 分层抽样：免疫 PE/control 各 1000
    immune = [c for c in adata.obs_names.astype(str) if amap.get(c, ("NA",))[0] == "Immune"]
    groups = adata.obs["group"].astype(str).to_numpy()
    idx_all = {str(c): i for i, c in enumerate(adata.obs_names)}
    pe_i = [idx_all[c] for c in immune if groups[idx_all[c]] == "PE"]
    ct_i = [idx_all[c] for c in immune if groups[idx_all[c]] == "control"]
    rng = np.random.default_rng(0)
    sel = np.concatenate([rng.choice(pe_i, min(1000, len(pe_i)), replace=False),
                          rng.choice(ct_i, min(1000, len(ct_i)), replace=False)])
    sub = adata[sel].copy()
    sub.layers["counts"] = sub.X.copy()
    print("子样本: %d（免疫 PE %d / CT %d）" % (sub.n_obs, len(pe_i), len(ct_i)), flush=True)

    from SIGnature.models.scfoundation import SCFoundationWrapper
    scf = SCFoundationWrapper(model_path=str(BASE / "data/model_files/scfoundation"), use_gpu=True)
    tdata = scf.preprocess_adata(sub.copy())
    t0 = time.time()
    atts = scf.calculate_attributions(tdata.X, batch_size=1, method="ig")
    atts = sp.csr_matrix(as2d(atts))
    print("scF 归因完成 %.0fs, shape=%s" % (time.time() - t0, atts.shape), flush=True)

    g2 = sub.obs["group"].astype(str).to_numpy()
    pe_sig = np.asarray(atts[g2 == "PE"].mean(axis=0)).ravel()
    ct_sig = np.asarray(atts[g2 == "control"].mean(axis=0)).ravel()
    d = pe_sig - ct_sig
    scf_genes = [str(x) for x in tdata.var_names]

    # SCimilarity 同数据集免疫签名（全量）
    att = sp.load_npz(P1 / ("atts_%s.npz" % GSE))
    meta = np.load(P1 / ("meta_%s.npz" % GSE), allow_pickle=True)
    cells = [str(c) for c in meta["cells"]]
    mgroups = np.array([str(g) for g in meta["groups"]])
    tarr = np.array([amap.get(c, ("NA",))[0] for c in cells])
    pe_s = np.asarray(att[(tarr == "Immune") & (mgroups == "PE")].mean(axis=0)).ravel()
    ct_s = np.asarray(att[(tarr == "Immune") & (mgroups == "control")].mean(axis=0)).ravel()
    d_sci = pe_s - ct_s
    sci_genes = [str(g) for g in meta["genes"]]

    # 共同基因空间上比较方向一致性
    common = [g for g in scf_genes if g in set(sci_genes)]
    si = {g: i for i, g in enumerate(sci_genes)}
    fi = {g: i for i, g in enumerate(scf_genes)}
    v1 = np.array([d_sci[si[g]] for g in common])
    v2 = np.array([d[fi[g]] for g in common])
    rho = spear(v1, v2)
    top200_sci = set(np.array(common)[np.argsort(-np.abs(v1))[:200]])
    top200_scf = set(np.array(common)[np.argsort(-np.abs(v2))[:200]])
    jac = len(top200_sci & top200_scf) / len(top200_sci | top200_scf)

    pd.DataFrame({"gene": common, "scfoundation_diff": v2, "scimilarity_diff": v1}).to_csv(
        OUT / "scf_immune_sig.csv", index=False)
    print("SENS-JSON: {\"spearman_scf_vs_sci\": %.3f, \"top200_jaccard\": %.3f, \"n_common\": %d}" % (rho, jac, len(common)), flush=True)
    print("SCF-SENS-DONE", flush=True)


if __name__ == "__main__":
    main()
