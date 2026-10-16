# -*- coding: utf-8 -*-
"""
Phase 0 小样本试跑：在 GSE173193 QC 数据的分层子样本上跑 SIGnature 归因
- 每样本抽 250 细胞（共 ~2000），先跑权重就绪的 SCimilarity / SSL(RandomMask) / scFoundation
- scVI 需要 CZI census 权重与 ensembl 映射，单独阶段再跑
输出：results/phase0/ 下各模型 atts npz；状态汇总见日志 PHASE0-DONE 行
"""
import json
import time
from os.path import join
from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc

BASE = "/root/pe_signature"
MF = join(BASE, "data", "model_files")
QC = join(BASE, "data", "processed", "GSE173193_qc.h5ad")
N_PER_SAMPLE = 250
OUT_DIR = Path("/root/pe_signature/results/phase0")


def load_sub():
    adata = sc.read_h5ad(QC)
    idx = (adata.obs.groupby("sample", observed=True)
           .sample(n=min(N_PER_SAMPLE, 250), random_state=0).index)
    sub = adata[idx].copy()
    sub.layers["counts"] = sub.X.copy()  # SCimilarity 要求原始计数放在 layers['counts']
    print("子样本: %d 细胞 × %d 基因" % (sub.n_obs, sub.n_vars), flush=True)
    return sub


def timed(name, fn):
    t0 = time.time()
    out = fn()
    print("[%s] 完成 %.0fs" % (name, time.time() - t0), flush=True)
    return out


def as2d(a):
    """兼容稀疏矩阵 / 分批列表 / 数组的返回值，统一转稠密 float32。"""
    import scipy.sparse as sp
    if sp.issparse(a):
        return np.asarray(a.todense(), dtype=np.float32)
    if isinstance(a, (list, tuple)):
        flat = [x.todense() if sp.issparse(x) else x for x in a]
        return np.vstack([np.asarray(x, dtype=np.float32) for x in flat])
    arr = np.asarray(a, dtype=np.float32)
    if arr.dtype == object:
        return np.vstack([np.asarray(x, dtype=np.float32) for x in arr])
    return arr


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    sub = load_sub()
    base = MF  # tar 已解压到 /root/pe_signature/data/model_files/，内含 scimilarity/ scfoundation/ ssl/
    status = {}

    # ---- SCimilarity（权重含在 Zenodo 包内） ----
    try:
        from SIGnature.models.scimilarity import SCimilarityWrapper
        scim = SCimilarityWrapper(model_path=join(base, "scimilarity"), use_gpu=True)
        tadata = timed("scimilarity-prep", lambda: scim.preprocess_adata(sub.copy()))
        atts = timed("scimilarity-attr", lambda: scim.calculate_attributions(
            tadata.X, batch_size=100, method="ig"))
        np.savez_compressed(
            "/root/pe_signature/results/phase0/atts_scimilarity.npz",
            atts=as2d(atts),
            genes=np.array([str(x) for x in tadata.var_names], dtype=object),
            cells=np.array([str(x) for x in tadata.obs_names], dtype=object),
            samples=np.array([str(x) for x in tadata.obs["sample"]], dtype=object),
            groups=np.array([str(x) for x in tadata.obs["group"]], dtype=object))
        status["scimilarity"] = "ok"
    except Exception as e:
        status["scimilarity"] = "FAIL: %s: %s" % (type(e).__name__, e)
        print(status["scimilarity"], flush=True)

    # ---- SSL（scTab RandomMask，权重含在包内则可跑） ----
    try:
        from SIGnature.models.ssl import SSLWrapper
        ssl = SSLWrapper(model_path=join(base, "ssl"),
                         model_filename="RandomMask.ckpt", use_gpu=True)
        tadata = timed("ssl-prep", lambda: ssl.preprocess_adata(sub.copy()))
        atts = timed("ssl-attr", lambda: ssl.calculate_attributions(
            tadata.X, batch_size=100, method="ig"))
        np.savez_compressed(
            "/root/pe_signature/results/phase0/atts_ssl_randommask.npz",
            atts=as2d(atts),
            genes=np.array(tadata.var_names, dtype=object),
            cells=np.array(tadata.obs_names, dtype=object),
            samples=np.array(tadata.obs["sample"], dtype=object),
            groups=np.array(tadata.obs["group"], dtype=object))
        status["ssl_randommask"] = "ok"
    except Exception as e:
        status["ssl_randommask"] = "FAIL: %s: %s" % (type(e).__name__, e)
        print(status["ssl_randommask"], flush=True)

    # ---- scFoundation（全转录组 batch_size=1，最慢，放最后） ----
    try:
        from SIGnature.models.scfoundation import SCFoundationWrapper
        scf = SCFoundationWrapper(model_path=join(base, "scfoundation"), use_gpu=True)
        tadata = timed("scfoundation-prep", lambda: scf.preprocess_adata(sub.copy()))
        atts = timed("scfoundation-attr", lambda: scf.calculate_attributions(
            tadata.X, batch_size=1, method="ig"))
        np.savez_compressed(
            "/root/pe_signature/results/phase0/atts_scfoundation.npz",
            atts=as2d(atts),
            genes=np.array(tadata.var_names, dtype=object),
            cells=np.array(tadata.obs_names, dtype=object),
            samples=np.array(tadata.obs["sample"], dtype=object),
            groups=np.array(tadata.obs["group"], dtype=object))
        status["scfoundation"] = "ok"
    except Exception as e:
        status["scfoundation"] = "FAIL: %s: %s" % (type(e).__name__, e)
        print(status["scfoundation"], flush=True)

    print("PHASE0-DONE", json.dumps(status), flush=True)


if __name__ == "__main__":
    main()
