# -*- coding: utf-8 -*-
"""
28 亚型对比（本地）：EOPE(GSE290578) vs LOPE(GSE298119) vs severe(GSE298602) vs mixed(GSE173193)
每队列免疫+滋养层的类型内 PE 签名（z）：
  A. 与共识补体-单核模块（免疫交集基因）的相关：各亚型炎症强度排序
  B. EOPE vs LOPE 签名直接对比：差异最大/一致基因
输出 stdout + results/subtype/subtype_report.txt
"""
import numpy as np
import pandas as pd
import scipy.sparse as sp
from pathlib import Path

P1 = Path("F:/pe_signature/results/phase1")
P0 = Path("F:/pe_signature/results/phase0")
HH = Path(r"D:/zcode/z0927/pe_signature/results/headtohead")
OUT = Path(r"D:/zcode/z0927/pe_signature/results/subtype")
OUT.mkdir(parents=True, exist_ok=True)

GSES = {
    "GSE173193": {"att": P0 / "atts_scimilarity_full.npz", "meta": P0 / "atts_scimilarity_full_meta.npz", "subtype": "mixed"},
    "GSE298119": {"att": P1 / "atts_GSE298119.npz", "meta": P1 / "meta_GSE298119.npz", "subtype": "LOPE"},
    "GSE298602": {"att": P1 / "atts_GSE298602.npz", "meta": P1 / "meta_GSE298602.npz", "subtype": "severe"},
    "GSE290578": {"att": P1 / "atts_GSE290578.npz", "meta": P1 / "meta_GSE290578.npz", "subtype": "EOPE"},
}
ANN = pd.read_csv(r"F:/pe_signature/results/phase1/../phase1/all_annotation.csv") if False else None
# 注释 CSV 在实例上生成后已回传 headtohead 所在项目目录之外——本地用 marker 门控重建（免疫/滋养层）
IMM_MARKERS = ["PTPRC", "LYZ", "CD74"]
TROPH_MARKERS = ["KRT7", "KRT18", "PEG10"]


def gate(Xd, genes_idx, markers):
    mi = [genes_idx[m] for m in markers if m in genes_idx]
    if not mi:
        return np.zeros(Xd.shape[0], dtype=bool)
    return (Xd[:, mi] > 0).any(axis=1)


def zscore(v):
    s = v.std()
    return (v - v.mean()) / (s if s > 0 else 1.0)


def main():
    lines = []
    sigs_imm, sigs_troph = {}, {}
    import gc
    for gse, cfg in GSES.items():
        atts = sp.load_npz(cfg["att"])
        meta = np.load(cfg["meta"], allow_pickle=True)
        cells = [str(c) for c in meta["cells"]]
        genes = [str(g) for g in meta["genes"]]
        groups = np.array([str(g) for g in meta["groups"]])
        groups = np.where(groups == "PE_preterm", "PE", groups)
        gi = {g: i for i, g in enumerate(genes)}
        # 归因矩阵列上的基因表达不可得（本地无 h5ad），改用归因门控近似：
        # 免疫门控 = PTPRC/LYZ/CD74 归因高；滋养层 = KRT7/KRT18/PEG10 归因高
        A = atts[:2000] if atts.shape[0] > 2000 else atts  # 抽样估计门控阈值
        thr = {}
        for m in IMM_MARKERS + TROPH_MARKERS:
            if m in gi:
                thr[m] = np.quantile(np.asarray(A[:, gi[m]].todense()).ravel(), 0.75)
        imm_mask = np.zeros(atts.shape[0], dtype=bool)
        troph_mask = np.zeros(atts.shape[0], dtype=bool)
        dens_cols_imm = {m: np.asarray(atts[:, gi[m]].todense()).ravel() for m in IMM_MARKERS if m in gi}
        dens_cols_tro = {m: np.asarray(atts[:, gi[m]].todense()).ravel() for m in TROPH_MARKERS if m in gi}
        for m, col in dens_cols_imm.items():
            imm_mask |= col > thr[m]
        for m, col in dens_cols_tro.items():
            troph_mask |= col > thr[m]
        for tname, mask, store in [("Immune", imm_mask, sigs_imm), ("Trophoblast", troph_mask, sigs_troph)]:
            pe = np.asarray(atts[mask & (groups == "PE")].mean(axis=0)).ravel()
            ct = np.asarray(atts[mask & (groups == "control")].mean(axis=0)).ravel()
            store[cfg["subtype"]] = zscore(pe - ct)
            lines.append("%s[%s] PE %d / CT %d 细胞（归因门控）" %
                         (gse, tname, (mask & (groups == "PE")).sum(), (mask & (groups == "control")).sum()))
        del atts, dens_cols_imm, dens_cols_tro, imm_mask, troph_mask
        gc.collect()

    # A. 与共识补体-单核模块的相关
    imm_consensus = set(pd.read_csv(HH / "Immune_attr_consensus_top200_intersect.csv").gene)
    troph_consensus = set(pd.read_csv(HH / "Trophoblast_attr_consensus_top200_intersect.csv").gene)

    def corr_with(sig, genes, cset):
        hit = np.array([g in cset for g in genes])
        if hit.sum() < 3:
            return np.nan
        return float(np.corrcoef(sig, hit.astype(float))[0, 1])

    lines.append("\n== A. 与共识模块的相关（点二列）==")
    meta0 = np.load(GSES["GSE298119"]["meta"], allow_pickle=True)
    genes0 = [str(g) for g in meta0["genes"]]
    # 各队列基因空间一致（SCimilarity 28231），直接用 genes0
    for tname, store, cset in [("免疫", sigs_imm, imm_consensus), ("滋养层", sigs_troph, troph_consensus)]:
        for sub, sig in sorted(store.items()):
            r = corr_with(sig, genes0, cset)
            lines.append("  %s | %-8s vs 共识: r=%.3f" % (tname, sub, r))

    # B. EOPE vs LOPE 免疫签名对比
    if "EOPE" in sigs_imm and "LOPE" in sigs_imm:
        d = sigs_imm["EOPE"] - sigs_imm["LOPE"]
        order = np.argsort(-np.abs(d))
        lines.append("\n== B. 免疫签名 EOPE−LOPE 差异 top15 ==")
        lines.append("  " + ", ".join("%s(%+.2f)" % (genes0[i], d[i]) for i in order[:15]))
        agree = (np.sign(sigs_imm["EOPE"]) == np.sign(sigs_imm["LOPE"])).mean()
        lines.append("  方向一致率: %.1f%%" % (agree * 100))

    report = "\n".join(str(x) for x in lines)
    print(report)
    (OUT / "subtype_report.txt").write_text(report, encoding="utf-8")
    print("SUBTYPE-DONE")


if __name__ == "__main__":
    main()
