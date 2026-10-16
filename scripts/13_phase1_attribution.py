# -*- coding: utf-8 -*-
"""
Phase 1 归因：对每个已预处理数据集跑 SCimilarity 全量归因（稀疏保存）。
输入：/root/pe_signature/data/processed/<GSE>_qc.h5ad
输出：/root/pe_signature/results/phase1/atts_<GSE>.npz + <GSE>_meta.npz
"""
import json
import time
from os.path import join
from pathlib import Path

import numpy as np
import scanpy as sc
import scipy.sparse as sp

BASE = "/root/pe_signature"
MODEL_PATH = join(BASE, "data/model_files/scimilarity")
GSES = ["GSE282038", "GSE298119", "GSE298602", "GSE240922", "GSE290578", "GSE329173"]
H5ADS = {
    "GSE282038": Path(BASE + "/data/processed/GSE282038_qc.h5ad"),
    "GSE298119": Path(BASE + "/data/processed/GSE298119_qc.h5ad"),
    "GSE298602": Path(BASE + "/data/processed/GSE298602_qc.h5ad"),
    "GSE240922": Path(BASE + "/data/processed/GSE240922_qc.h5ad"),
    "GSE290578": Path(BASE + "/data/processed/GSE290578_qc.h5ad"),
    "GSE329173": Path(BASE + "/data/processed/GSE329173_qc.h5ad"),
}
OUT_DIR = Path(BASE + "/results/phase1")
OUT_NPZ = {
    "GSE282038": OUT_DIR / "atts_GSE282038.npz",
    "GSE298119": OUT_DIR / "atts_GSE298119.npz",
    "GSE298602": OUT_DIR / "atts_GSE298602.npz",
    "GSE240922": OUT_DIR / "atts_GSE240922.npz",
    "GSE290578": OUT_DIR / "atts_GSE290578.npz",
    "GSE329173": OUT_DIR / "atts_GSE329173.npz",
}
OUT_META = {
    "GSE282038": OUT_DIR / "meta_GSE282038.npz",
    "GSE298119": OUT_DIR / "meta_GSE298119.npz",
    "GSE298602": OUT_DIR / "meta_GSE298602.npz",
    "GSE240922": OUT_DIR / "meta_GSE240922.npz",
    "GSE290578": OUT_DIR / "meta_GSE290578.npz",
    "GSE329173": OUT_DIR / "meta_GSE329173.npz",
}


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    from SIGnature.models.scimilarity import SCimilarityWrapper
    scim = SCimilarityWrapper(model_path=MODEL_PATH, use_gpu=True)
    status = {}
    for gse in GSES:
        h5 = H5ADS[gse]
        if not h5.exists():
            status[gse] = "no-h5ad"
            print(status[gse], gse, flush=True)
            continue
        if OUT_NPZ[gse].exists():
            status[gse] = "already-done"
            print(status[gse], gse, flush=True)
            continue
        try:
            adata = sc.read_h5ad(h5)
            adata.layers["counts"] = adata.X.copy()
            print("[%s] %d 细胞 归因中 ..." % (gse, adata.n_obs), flush=True)
            t0 = time.time()
            tadata = scim.preprocess_adata(adata.copy())
            atts = scim.calculate_attributions(tadata.X, batch_size=100, method="ig")
            if sp.issparse(atts):
                atts = atts.astype(np.float32)
            else:
                atts = sp.csr_matrix(np.asarray(atts, dtype=np.float32))
            sp.save_npz(OUT_NPZ[gse], atts)
            np.savez(OUT_META[gse],
                     genes=np.array([str(x) for x in tadata.var_names], dtype=object),
                     cells=np.array([str(x) for x in tadata.obs_names], dtype=object),
                     samples=np.array([str(x) for x in tadata.obs["sample"]], dtype=object),
                     groups=np.array([str(x) for x in tadata.obs["group"]], dtype=object))
            status[gse] = "ok %ds shape=%s" % (round(time.time() - t0), atts.shape)
            print(status[gse], gse, flush=True)
        except Exception as e:
            status[gse] = "FAIL: %s: %s" % (type(e).__name__, e)
            print(status[gse], flush=True)
    print("PHASE1-DONE", json.dumps(status), flush=True)


if __name__ == "__main__":
    main()
