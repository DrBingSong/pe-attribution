# -*- coding: utf-8 -*-
"""
Phase 1 批量预处理：多 GSE 数据集 → 解压/重组 → 分组映射 → QC → h5ad
分组规则按数据集内置于 RULES；单数据集失败不阻塞其他。
安全：SOFT 仅允许 https+白名单主机且解析 IP 非私网；所有写出路径强制 contained。
输出：/root/pe_signature/data/processed/<GSE>_qc.h5ad
"""
import ipaddress
import json
import re
import socket
import tarfile
import urllib.parse
import urllib.request
from pathlib import Path

import anndata as ad
import scanpy as sc

BASE = Path("/root/pe_signature")
GSE_DIR = BASE / "data" / "gse"
PROC = BASE / "data" / "processed"
PROC.mkdir(parents=True, exist_ok=True)
SOFT_HOST = "www.ncbi.nlm.nih.gov"

GSES = ["GSE282038", "GSE298119", "GSE298602", "GSE240922", "GSE290578", "GSE329173"]
OUTS = {
    "GSE282038": PROC / "GSE282038_qc.h5ad",
    "GSE298119": PROC / "GSE298119_qc.h5ad",
    "GSE298602": PROC / "GSE298602_qc.h5ad",
    "GSE240922": PROC / "GSE240922_qc.h5ad",
    "GSE290578": PROC / "GSE290578_qc.h5ad",
    "GSE329173": PROC / "GSE329173_qc.h5ad",
}
MIN_GENES, MAX_GENES, MAX_PCT_MT, MIN_CELLS = 200, 6000, 20.0, 3


def safe_child(parent, *parts):
    """规范化并强制子路径不逃逸出 parent。"""
    p = Path(parent).joinpath(*parts).resolve()
    if Path(parent).resolve() not in p.parents:
        raise ValueError("path escape: %s" % p)
    return p


def fetch_soft(gse):
    if gse not in GSES:
        raise ValueError("dataset not in allowlist")
    url = "https://%s/geo/query/acc.cgi?acc=%s&targ=gsm&form=text&view=brief" % (SOFT_HOST, gse)
    parsed = urllib.parse.urlparse(url)
    assert parsed.scheme == "https" and parsed.hostname == SOFT_HOST
    ip = socket.getaddrinfo(parsed.hostname, 443)[0][4][0]
    assert not ipaddress.ip_address(ip).is_private, "resolved to private address"
    assert not ipaddress.ip_address(ip).is_loopback
    with urllib.request.urlopen(url, timeout=90) as r:
        return r.read().decode(errors="replace")


def parse_soft_groups(gse, soft):
    blocks = re.split(r"\^SAMPLE = ", soft)[1:]
    out = {}
    for b in blocks:
        m = re.search(r"!Sample_title = (.+)", b)
        if not m:
            continue
        title = m.group(1).strip()
        text = title + "\n" + b
        g = "NA"
        if gse == "GSE298119":
            g = "PE" if "preeclamptic" in text else ("control" if "normotensive" in text else "NA")
        elif gse == "GSE298602":
            g = "PE" if "PreE_SF" in text else ("ghtn" if "gHTN" in text else ("control" if "Control" in text else "NA"))
        elif gse == "GSE290578":
            g = "PE_preterm" if "preeclampsia" in text else ("control" if "normal" in text else "NA")
        elif gse == "GSE240922":
            g = "reference_FT"
        elif gse == "GSE329173":
            g = ("severe_PE" if re.search(r"severe eclampsia|severe preeclampsia", text, re.I)
                 else ("control" if re.search(r"normal|control|healthy", text, re.I) else "NA"))
        elif gse == "GSE282038":
            g = "EOPE" if "early onset" in text else "NA"
        out[title] = g
    return out


def group_of(soft_map, fname):
    low = fname.lower()
    best, bestlen = "NA", 0
    for title, g in soft_map.items():
        t = title.lower()
        if len(t) > bestlen and t in low:
            best, bestlen = g, len(t)
    return best


def extract_and_organize(gse):
    d = GSE_DIR / gse
    ext = d / "extracted"
    marker = safe_child(ext, ".organized")
    if marker.exists():
        return ext
    tars = list(d.glob("*.tar"))
    if not tars:
        raise RuntimeError("no RAW.tar")
    ext.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tars[0]) as tf:
        tf.extractall(ext)
    for bc in list(ext.rglob("*_barcodes.tsv.gz")):
        if ".." in bc.name or "/" in bc.stem:
            continue
        base_dir = bc.parent
        name = bc.name[: -len("_barcodes.tsv.gz")]
        tgt = safe_child(ext, name)
        tgt.mkdir(exist_ok=True)
        for suffix in ("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"):
            src = safe_child(base_dir, name + "_" + suffix)
            if src.exists():
                src.rename(safe_child(tgt, suffix))
        if base_dir != ext:
            final = safe_child(ext, name)
            if tgt != final:
                tgt.rename(final)
    marker.write_text("1")
    return ext


def load_dataset(gse):
    ext = extract_and_organize(gse)
    soft = fetch_soft(gse)
    smap = parse_soft_groups(gse, soft)
    parts, meta = [], []
    for sdir in sorted(p for p in ext.iterdir() if p.is_dir()):
        tri = [sdir / x for x in ("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz")]
        if not all(t.exists() for t in tri):
            continue
        a = sc.read_10x_mtx(sdir, var_names="gene_symbols", make_unique=True, cache=True)
        a.var_names_make_unique()
        a.obs_names_make_unique()
        g = group_of(smap, sdir.name.lower())
        a.obs["sample"] = sdir.name
        a.obs["group"] = g
        parts.append(a)
        meta.append((sdir.name, g, int(a.n_obs)))
    if not parts:
        raise RuntimeError("no 10x dirs recognized")
    adx = ad.concat(parts, index_unique="-")
    print("  样本映射:", meta, flush=True)
    return adx, meta


def qc_save(gse, adata):
    adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    sc.pp.filter_cells(adata, min_genes=MIN_GENES)
    sc.pp.filter_cells(adata, max_genes=MAX_GENES)
    adata = adata[adata.obs["pct_counts_mt"] <= MAX_PCT_MT, :].copy()
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS)
    adata.write(OUTS[gse])
    print("  %s QC后: %d 细胞 × %d 基因" % (gse, adata.n_obs, adata.n_vars), flush=True)
    print("  分组计数:", adata.obs["group"].value_counts().to_dict(), flush=True)


def main():
    summary = {}
    for gse in GSES:
        marker = GSE_DIR / gse / ".done"
        if not marker.exists():
            print("[%s] 下载未完成，跳过" % gse, flush=True)
            summary[gse] = "not-downloaded"
            continue
        if OUTS[gse].exists():
            print("[%s] 已处理，跳过" % gse, flush=True)
            summary[gse] = "already-done"
            continue
        try:
            print("[%s] 处理中 ..." % gse, flush=True)
            adx, _ = load_dataset(gse)
            qc_save(gse, adx)
            summary[gse] = "ok"
        except Exception as e:
            summary[gse] = "FAIL: %s: %s" % (type(e).__name__, e)
            print(summary[gse], flush=True)
    print("BATCH12-DONE", json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
