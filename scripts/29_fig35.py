# -*- coding: utf-8 -*-
"""
29 论文图 F3+F5：
F3 共识图谱：四类型共识 up/down 基因数 + NMF 队列载荷热图（v2.1 日志）
F5 异质性证据链：v1 全细胞 ρ ≈ 0 → v2 类型分层 → v2.1 加 EOPE 队列（每类型 |ρ| 分布）
数据：results/phase3（v1）、results/phase3_v2（v2 3队列）、phase1/logs_v21.log（v2.1 4队列）
"""
import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

R = Path(r"D:/zcode/z0927/pe_signature/results")
FIG = R / "figures"
FIG.mkdir(parents=True, exist_ok=True)
TYPES = ["Trophoblast", "Immune", "Endothelial", "Stromal"]

# ---------- 数据收集 ----------
# v2（3 队列，本地 CSV）
v2_corr, v2_counts = {}, {}
for t in TYPES:
    f = R / "phase3_v2" / ("%s_dataset_loadings.csv" % ("coarse_" + t))
    up = pd.read_csv(R / "phase3_v2" / ("coarse_%s_consensus_up.csv" % t))
    dn = pd.read_csv(R / "phase3_v2" / ("coarse_%s_consensus_down.csv" % t))
    v2_counts[t] = (len(up), len(dn))

# v2.1（4 队列，从日志 JSON 解析）
v21 = {}
txt = Path("F:/pe_signature/results/phase1/logs_v21.log").read_text(encoding="utf-8", errors="replace")
m = re.search(r"V2-SUMMARY-BEGIN\n(\{.*?\n\})\nV2-SUMMARY-END", txt, re.S)
if m:
    v21 = json.loads(m.group(1))

# v1（全细胞，本地 CSV）
v1_corr = pd.read_csv(R / "phase3" / "dataset_correlation.csv", index_col=0)

# ---------- F3 ----------
fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw={"width_ratios": [1, 1.3]})

ax = axes[0]
x = np.arange(len(TYPES))
ups = [v2_counts[t][0] for t in TYPES]
dns = [v2_counts[t][1] for t in TYPES]
ax.bar(x, ups, 0.55, label="PE-up consensus", color="#b2182b")
ax.bar(x, [-d for d in dns], 0.55, label="PE-down consensus", color="#2166ac")
ax.axhline(0, color="k", lw=0.8)
for i, (u, d) in enumerate(zip(ups, dns)):
    ax.text(i, u + 12, str(u), ha="center", fontsize=9)
    ax.text(i, -d - 30, str(d), ha="center", fontsize=9)
ax.set_xticks(x)
ax.set_xticklabels(TYPES, fontsize=9)
ax.set_ylabel("Consensus genes (≥3/4 cohorts)", fontsize=10)
ax.set_title("A  Cell-type consensus PE programs", fontsize=11, loc="left")
ax.legend(fontsize=8, frameon=False, loc="lower right")
ax.spines[["top", "right"]].set_visible(False)

ax = axes[1]
labels = []
plot_types = [("Immune", "Immune"), ("Trophoblast", "Trophoblast")]
data_rows, data_cols = [], []
for t in TYPES:
    key = "coarse_" + t
    if key in v21 and "nmf_loadings" in v21[key]:
        load = v21[key]["nmf_loadings"]
        for g in load:
            data_rows.append([t[:4]] )
            data_cols.append(list(load[g].values()))
        labels += list(load.keys())
M = np.array(data_cols)
im = ax.imshow(M, aspect="auto", cmap="YlGnBu")
ax.set_xticks(range(M.shape[1]))
ax.set_xticklabels(["P1", "P2", "P3"], fontsize=9)
ax.set_yticks(range(len(labels)))
ax.set_yticklabels(labels, fontsize=8)
# 类型分隔线
bounds, cur = [], TYPES[0]
for i, lb in enumerate(labels):
    if not lb.startswith("GSE"):
        pass
sep = 0
prev_type = None
for i, t in enumerate([r[0] for r in data_rows]):
    if prev_type is not None and t != prev_type:
        ax.axhline(i - 0.5, color="k", lw=1.2)
    prev_type = t
ax.set_title("B  NMF loadings by cohort (v2.1, 4 cohorts)", fontsize=11, loc="left")
fig.colorbar(im, ax=ax, shrink=0.8)
fig.tight_layout()
fig.savefig(FIG / "fig3_consensus_map.png", dpi=300, bbox_inches="tight")
print("F3 saved")

# ---------- F5 ----------
fig, ax = plt.subplots(figsize=(8.5, 4.6))
stages = []
stages.append(("All-cell\n(v1, 3 cohorts)", [abs(v1_corr.iloc[i, j]) for i in range(3) for j in range(i + 1, 3)]))
if m:
    for t in TYPES:
        key = "coarse_" + t
        if key in v21:
            vals = [abs(v) for v in v21[key]["pairwise_spearman"].values()]
            stages.append(("%s\n(v2.1, 4 cohorts)" % t, vals))
    # v2（3 队列）从 phase3_v2 无法直接取 pairwise（当时只在 stdout）——用日志 v2.log 若有
    p2 = Path("F:/pe_signature/results/phase1/logs_v2.log")
    if p2.exists():
        txt2 = p2.read_text(encoding="utf-8", errors="replace")
        m2 = re.search(r"V2-SUMMARY-BEGIN\n(\{.*?\n\})\nV2-SUMMARY-END", txt2, re.S)
        if m2:
            for t in ["Trophoblast", "Immune"]:
                key = "coarse_" + t
                if key in json.loads(m2.group(1)):
                    vals = [abs(v) for v in json.loads(m2.group(1))[key]["pairwise_spearman"].values()]
                    stages.insert(1, ("%s\n(v2, 3 cohorts)" % t, vals))

xs, positions = [], []
pos = 0
for i, (name, vals) in enumerate(stages):
    ax.scatter([pos] * len(vals), vals, s=42, color="#b2182b", alpha=0.75, zorder=3)
    ax.scatter([pos], [np.mean(vals)], s=260, color="#2166ac", zorder=4, marker="_", linewidths=3)
    positions.append(pos)
    pos += 1
ax.set_xticks(positions)
ax.set_xticklabels([s[0] for s in stages], fontsize=8)
ax.set_ylabel("|Spearman ρ| between cohort signatures", fontsize=10)
ax.set_title("Evidence chain: composition confounding → stratification gains → cohort heterogeneity",
             fontsize=11, loc="left")
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(FIG / "fig5_heterogeneity.png", dpi=300, bbox_inches="tight")
print("F5 saved")
