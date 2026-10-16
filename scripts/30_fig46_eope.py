# -*- coding: utf-8 -*-
"""
30 论文图 F4+F6 + EOPE 投影：
F4: 免疫/内皮/滋养层 top 通路富集条形图（-log10 adj.p）
F6: 低表达 TF 捞回 gain 倍数 × 7 数据集 点图
EOPE 投影: GSE282038(EOPE) / GSE329173(severe) 在 v1 共识 up/down 签名上的 z 值
数据均本地：results/enrichr、results/phase1_sanity_all.txt、F:/pe_signature
"""
import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp

R = Path(r"D:/zcode/z0927/pe_signature/results")
FIG = R / "figures"
FB = Path("F:/pe_signature/results")

# ---------------- F4 ----------------
fig, axes = plt.subplots(1, 3, figsize=(13, 4.8))
for ax, t, lib in zip(axes, ["Immune", "Endothelial", "Trophoblast"],
                      ["Reactome_2022", "GO_Biological_Process_2023", "GO_Biological_Process_2023"]):
    f = R / "enrichr" / ("%s_%s.csv" % (t, lib))
    if not f.exists():
        ax.set_title("%s (no data)" % t)
        continue
    df = pd.read_csv(f).head(8).iloc[::-1]
    names = [x.split(" R-HSA")[0].title()[:44] for x in df.term]
    vals = -np.log10(df.adj_pval)
    ax.barh(range(len(df)), vals, color="#2166ac", height=0.62)
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(names, fontsize=7.5)
    ax.set_xlabel("-log10(adj.p)", fontsize=9)
    ax.set_title("%s consensus module" % t, fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
fig.suptitle("F4  Functional annotation of cross-cohort consensus programs", fontsize=12, x=0.02, ha="left")
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(FIG / "fig4_enrichment.png", dpi=300, bbox_inches="tight")
print("F4 saved")

# ---------------- F6 ----------------
raw = (R / "phase1_sanity_all.txt").read_text(encoding="utf-8", errors="replace")
m = re.search(r"SANITY-ALL-BEGIN\n(\{.*?\n\})\nSANITY-ALL-END", raw, re.S)
san = json.loads(m.group(1))
# GSE173193 从 phase0 的 full sanity（夜间报告数字）手动并入
P0_GAINS = {"GATA2": 4.7, "NR2F2": 2.7, "KLF4": 2.0, "TEAD4": 1.9, "TFAP2C": 1.8,
            "GATA3": 0.5, "MSX2": 0.4, "HAND1": 1.1, "ELF5": 1.3, "DLX3": 0.8}
TFS = ["GATA2", "NR2F2", "KLF4", "TEAD4", "TFAP2C", "ELF5", "HAND1", "GATA3"]
dsets = [k for k in san if "B_lowexpr_TF_rescue" in san.get(k, {})]
fig, ax = plt.subplots(figsize=(8.5, 4.6))
rng = np.random.default_rng(1)
for j, tf in enumerate(TFS):
    xs, ys = [], []
    for i, ds in enumerate(dsets):
        v = san[ds]["B_lowexpr_TF_rescue"].get(tf)
        if isinstance(v, dict):
            xs.append(i + rng.uniform(-0.18, 0.18))
            ys.append(v["gain_x"])
    xs.append(len(dsets) + rng.uniform(-0.18, 0.18))
    ys.append(P0_GAINS.get(tf, np.nan))
    color = "#b2182b" if np.nanmean(ys) >= 1.5 else "#7f7f7f"
    ax.scatter(xs, ys, s=34, color=color, alpha=0.8, zorder=3)
ax.axhline(1.0, color="k", lw=0.8, ls="--")
ax.axhline(1.5, color="#b2182b", lw=0.8, ls=":")
labels = [d.replace("GSE", "") for d in dsets] + ["173193"]
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels, fontsize=9)
ax.set_ylabel("Attribution rank gain vs expression (×)", fontsize=10)
ax.set_title("F6  Low-expression TF rescue across 7 cohorts (red = mean ≥1.5×)", fontsize=11, loc="left")
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(FIG / "fig6_tf_rescue.png", dpi=300, bbox_inches="tight")
print("F6 saved")

# ---------------- EOPE 投影 ----------------
up = set(pd.read_csv(R / "phase3" / "consensus_pe_up.csv").gene)
dn = set(pd.read_csv(R / "phase3" / "consensus_pe_down.csv").gene)
rows = []
for gse, sub in [("GSE282038", "EOPE"), ("GSE329173", "severe_PE"), ("GSE240922", "reference_FT")]:
    atts = sp.load_npz(FB / "phase1" / ("atts_%s.npz" % gse))
    meta = np.load(FB / "phase1" / ("meta_%s.npz" % gse), allow_pickle=True)
    genes = [str(g) for g in meta["genes"]]
    ai = {g: i for i, g in enumerate(genes)}
    prof = np.asarray(atts.mean(axis=0)).ravel()
    z = (prof - prof.mean()) / (prof.std() or 1.0)
    upz = float(np.mean([z[ai[g]] for g in up if g in ai]))
    dnz = float(np.mean([z[ai[g]] for g in dn if g in ai]))
    rows.append({"dataset": gse, "subtype": sub, "n_cells": atts.shape[0],
                 "PEup_z": round(upz, 3), "PEdown_z": round(dnz, 3)})
    del atts
df = pd.DataFrame(rows)
df.to_csv(R / "subtype" / "eope_projection.csv", index=False)
print(df.to_string(index=False))
print("ALL30-DONE")
