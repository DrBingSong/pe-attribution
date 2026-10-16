# -*- coding: utf-8 -*-
"""仅补跑 SCimilarity 归因（带分批返回值 vstack 修复），输出单一 npz。"""
import time
from os.path import join
from pathlib import Path

import numpy as np
import scanpy as sc

BASE = "/root/pe_signature"
QC = join(BASE, "data", "processed", "GSE173193_qc.h5ad")
OUT_NPZ = Path("/root/pe_signature/results/phase0/atts_scimilarity.npz")


def as2d(a):
    if isinstance(a, (list, tuple)):
        return np.vstack([np.asarray(x, dtype=np.float32) for x in a])
    arr = np.asarray(a, dtype=np.float32)
    if arr.dtype == object:
        return np.vstack([np.asarray(x, dtype=np.float32) for x in arr])
    return arr


def main():
    adata = sc.read_h5ad(QC)
    idx = (adata.obs.groupby("sample", observed=True)
           .sample(n=250, random_state=0).index)
    sub = adata[idx].copy()
    sub.layers["counts"] = sub.X.copy()
    print("子样本: %d 细胞" % sub.n_obs, flush=True)

    from SIGnature.models.scimilarity import SCimilarityWrapper
    scim = SCimilarityWrapper(model_path=join(BASE, "data/model_files/scimilarity"), use_gpu=True)
    tadata = scim.preprocess_adata(sub.copy())
    t0 = time.time()
    atts = as2d(scim.calculate_attributions(tadata.X, batch_size=100, method="ig"))
    print("归因完成 %.0fs, shape=%s" % (time.time() - t0, atts.shape), flush=True)

    OUT_NPZ.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        "/root/pe_signature/results/phase0/atts_scimilarity.npz",
        atts=atts,
        genes=np.array([str(x) for x in tadata.var_names], dtype=object),
        cells=np.array([str(x) for x in tadata.obs_names], dtype=object),
        samples=np.array([str(x) for x in tadata.obs["sample"]], dtype=object),
        groups=np.array([str(x) for x in tadata.obs["group"]], dtype=object))
    print("SCIM-ONLY-DONE", flush=True)


if __name__ == "__main__":
    main()
