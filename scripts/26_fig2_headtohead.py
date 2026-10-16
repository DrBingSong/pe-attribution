# -*- coding: utf-8 -*-
"""头对头核心图（论文 Fig.2 素材）：
A: 四细胞类型 × 两方法的跨队列 Jaccard（均值柱 + 24 配对散点）
B: 24 个配对的 attr vs expr 对角线散点 + Wilcoxon 符号秩检验
输出 results/figures/fig2_headtohead.png + 统计打印
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

SRC = Path(r"D:/zcode/z0927/pe_signature/results/headtohead/overlap_summary.csv")
OUT = Path(r"D:/zcode/z0927/pe_signature/results/figures")
OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(SRC)
types = ["Trophoblast", "Immune", "Endothelial", "Stromal"]
labels = ["Trophoblast", "Immune", "Endothelial", "Stromal"]

stat, p = wilcoxon(df["attr_jaccard"], df["expr_jaccard"])
print("Wilcoxon signed-rank: W=%.1f, p=%.2e (n=%d pairs)" % (stat, p, len(df)))
print("全局均值: attr=%.4f expr=%.4f ratio=%.1fx" %
      (df.attr_jaccard.mean(), df.expr_jaccard.mean(),
       df.attr_jaccard.mean() / df.expr_jaccard.mean()))
print("胜出配对: %d/%d" % ((df.attr_jaccard > df.expr_jaccard).sum(), len(df)))

fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw={"width_ratios": [1.4, 1]})

# Panel A
ax = axes[0]
x = np.arange(len(types))
w = 0.36
am = [df[df.type == t].attr_jaccard.mean() for t in types]
em = [df[df.type == t].expr_jaccard.mean() for t in types]
ax.bar(x - w / 2, am, w, label="Attribution (SIGnature)", color="#2166ac")
ax.bar(x + w / 2, em, w, label="Expression (t-test top200)", color="#d6604d")
rng = np.random.default_rng(0)
for i, t in enumerate(types):
    sub = df[df.type == t]
    ax.scatter(np.full(len(sub), i - w / 2) + rng.uniform(-0.07, 0.07, len(sub)),
               sub.attr_jaccard, s=12, color="k", alpha=0.45, zorder=3)
    ax.scatter(np.full(len(sub), i + w / 2) + rng.uniform(-0.07, 0.07, len(sub)),
               sub.expr_jaccard, s=12, color="k", alpha=0.45, zorder=3)
    ax.text(i, max(am[i], em[i]) + 0.032, "%.1fx" % (am[i] / em[i]),
            ha="center", fontsize=10, fontweight="bold", color="#333")
ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=10)
ax.set_ylabel("Cross-cohort gene-set overlap (Jaccard)", fontsize=10)
ax.set_title("A  PE differential signatures: cross-cohort reproducibility", fontsize=11, loc="left")
ax.legend(fontsize=9, frameon=False, loc="upper right")
ax.spines[["top", "right"]].set_visible(False)
ax.set_ylim(0, 0.44)

# Panel B
ax = axes[1]
lim = (0, 0.42)
ax.plot(lim, lim, "-", color="#999", lw=1, zorder=1)
colors = {"Trophoblast": "#4393c3", "Immune": "#f4a582", "Endothelial": "#92c5de", "Stromal": "#ca0020"}
for t in types:
    sub = df[df.type == t]
    ax.scatter(sub.expr_jaccard, sub.attr_jaccard, s=34, alpha=0.85,
               color=colors[t], label=t, edgecolor="white", linewidth=0.5, zorder=3)
ax.set_xlim(lim)
ax.set_ylim(lim)
ax.set_xlabel("Expression-based overlap", fontsize=10)
ax.set_ylabel("Attribution-based overlap", fontsize=10)
ax.set_title("B  %d/%d pairs above diagonal  (Wilcoxon p=%.1e)"
             % ((df.attr_jaccard > df.expr_jaccard).sum(), len(df), p),
             fontsize=11, loc="left")
ax.legend(fontsize=8, frameon=False, loc="lower right")
ax.spines[["top", "right"]].set_visible(False)

fig.tight_layout()
fig.savefig(OUT / "fig2_headtohead.png", dpi=300, bbox_inches="tight")
print("saved:", OUT / "fig2_headtohead.png")
