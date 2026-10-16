# -*- coding: utf-8 -*-
"""全量 SCimilarity 归因（52,358 细胞），稀疏保存，避免稠密爆内存。"""
import time
from os.path import join

import numpy as np
import scanpy as sc
import scipy.sparse as sp

BASE = "/root/pe_signature"
QC = join(BASE, "data", "processed", "GSE173193_qc.h5ad")
OUT = "/root/pe_signature/results/phase0/atts_scimilarity_full.npz"
OUT_META = "/root/pe_signature/results/phase0/atts_scimilarity_full_meta.npz"


def main():
    adata = sc.read_h5ad(QC)
    adata.layers["counts"] = adata.X.copy()
    print("全量: %d 细胞 × %d 基因" % (adata.n_obs, adata.n_vars), flush=True)

    from SIGnature.models.scimilarity import SCimilarityWrapper
    scim = SCimilarityWrapper(model_path=join(BASE, "data/model_files/scimilarity"), use_gpu=True)
    tadata = scim.preprocess_adata(adata.copy())
    t0 = time.time()
    atts = scim.calculate_attributions(tadata.X, batch_size=100, method="ig")
    print("归因完成 %.0fs, type=%s" % (time.time() - t0, type(atts).__name__), flush=True)
    if not sp.issparse(atts):
        atts = sp.csr_matrix(np.asarray(atts, dtype=np.float32))
    else:
        atts = atts.astype(np.float32)

    sp.save_npz(OUT, atts)
    np.savez(OUT_META,
             genes=np.array([str(x) for x in tadata.var_names], dtype=object),
             cells=np.array([str(x) for x in tadata.obs_names], dtype=object),
             samples=np.array([str(x) for x in tadata.obs["sample"]], dtype=object),
             groups=np.array([str(x) for x in tadata.obs["group"]], dtype=object))
    print("FULL-DONE shape=%s" % (atts.shape,), flush=True)


if __name__ == "__main__":
    main()
