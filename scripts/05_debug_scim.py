# -*- coding: utf-8 -*-
"""定位 SCimilarity 保存失败的精确来源：逐字段构建并打印完整 traceback。"""
import time
import traceback
from os.path import join

import numpy as np
import scanpy as sc

BASE = "/root/pe_signature"
QC = join(BASE, "data", "processed", "GSE173193_qc.h5ad")


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

    from SIGnature.models.scimilarity import SCimilarityWrapper
    scim = SCimilarityWrapper(model_path=join(BASE, "data/model_files/scimilarity"), use_gpu=True)
    tadata = scim.preprocess_adata(sub.copy())
    print("var_names type:", type(tadata.var_names), "| 前两个:", list(tadata.var_names[:2]), flush=True)
    print("obs_names type:", type(tadata.obs_names), "| 前两个:", list(tadata.obs_names[:2]), flush=True)
    atts = scim.calculate_attributions(tadata.X, batch_size=100, method="ig")
    print("atts python type:", type(atts), flush=True)
    try:
        print("atts len:", len(atts), "elem0 shape:", np.asarray(atts[0]).shape, flush=True)
    except Exception:
        print("atts 无 len 或 elem 访问失败", flush=True)

    arr = as2d(atts)
    print("as2d shape:", arr.shape, arr.dtype, flush=True)

    fields = {}
    for key, val in [
        ("atts", arr),
        ("genes", np.array([str(x) for x in tadata.var_names], dtype=object)),
        ("cells", np.array([str(x) for x in tadata.obs_names], dtype=object)),
        ("samples", np.array([str(x) for x in tadata.obs["sample"]], dtype=object)),
        ("groups", np.array([str(x) for x in tadata.obs["group"]], dtype=object)),
    ]:
        try:
            np.savez_compressed("/root/pe_signature/results/phase0/dbg_%s.npz" % key, x=val)
            fields[key] = "ok"
        except Exception:
            fields[key] = "FAIL:\n" + traceback.format_exc(limit=3)
    print(json.dumps(fields, indent=1) if (json := __import__("json")) else fields)
    print("DEBUG-DONE", flush=True)


if __name__ == "__main__":
    main()
