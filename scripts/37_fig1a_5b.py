# -*- coding: utf-8 -*-
"""
37 补绘图：Fig.1a（七队列总览：细胞数 × 分组构成堆叠条形）+ Fig.5b（亚型 × 共识模块相关条形）。
数据全部实测（QC 后细胞数 / 分组构成 / subtype_report 相关值）。
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

FIG = r"D:/zcode/z0927/pe_signature/results/figures"

# ---- Fig.1a 数据（QC 后）----
cohorts = [
    # (GSE, 标签, {组: 细胞数}, 亚型)
    ("GSE173193", "mixed PE", {"PE": 14750, "control": 11354, "GDM": 13292, "elderly": 12962}, "contrast"),
    ("GSE298119", "term PE", {"PE": 24406, "control": 17842}, "contrast"),
    ("GSE298602", "severe PE", {"PE": 26486, "control": 24692, "gHTN": 7745}, "contrast"),
    ("GSE290578", "preterm PE", {"PE": 62818, "control": 72391}, "contrast"),
    ("GSE282038", "EOPE", {"EOPE": 9305}, "profile"),
    ("GSE329173", "severe eclampsia", {"severe_PE": 33861}, "profile"),
    ("GSE240922", "1st-trim reference", {"reference_FT": 132823}, "reference"),
]
GROUP_COLORS = {"PE": "#b2182b", "EOPE": "#d6604d", "severe_PE": "#f4a582", "control": "#2166ac",
                "gHTN": "#92c5de", "GDM": "#4393c3", "elderly": "#67a9cf", "reference_FT": "#8c510a"}

fig, ax = plt.subplots(figsize=(9, 4.6))
names = [c[0] for c in cohorts]
x = np.arange(len(cohorts))
bottom = np.zeros(len(cohorts))
order = ["PE", "EOPE", "severe_PE", "gHTN", "control", "GDM", "elderly", "reference_FT"]
for g in order:
    vals = np.array([c[2].get(g, 0) for c in cohorts])
    if vals.sum() == 0:
        continue
    ax.bar(x, vals, 0.62, bottom=bottom, label=g, color=GROUP_COLORS[g])
    bottom += vals
for i, c in enumerate(cohorts):
    total = sum(c[2].values())
    role = {"contrast": "▪", "profile": "○", "reference": "△"}[c[3]]
    ax.text(i, total + 3000, "%s\n%d" % (role, total), ha="center", fontsize=8.5)
ax.set_xticks(x)
ax.set_xticklabels(names, rotation=28, ha="right", fontsize=9)
ax.set_ylabel("QC-passed cells", fontsize=10)
ax.set_title("Fig.1a  Seven preeclampsia scRNA-seq cohorts (464,727 cells;  ▪ internal contrast  ○ disease profile  △ reference)",
             fontsize=10.5, loc="left")
ax.legend(fontsize=8, ncol=4, frameon=False, loc="upper left")
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(FIG + "/fig1a_cohort_overview.png", dpi=300, bbox_inches="tight")
print("fig1a saved")

# ---- Fig.5b 数据（subtype_report 实测）----
subs = ["EOPE\n(preterm)", "severe", "LOPE\n(term)", "mixed"]
immune = [0.249, 0.209, -0.306, -0.320]
troph = [0.096, -0.006, -0.159, 0.272]

fig, ax = plt.subplots(figsize=(8.2, 4.4))
x = np.arange(len(subs))
ax.bar(x - 0.19, immune, 0.36, label="Immune complement module", color="#b2182b")
ax.bar(x + 0.19, troph, 0.36, label="Trophoblast consensus module", color="#7f7f7f")
for xi, v in zip(x - 0.19, immune):
    ax.text(xi, v + (0.015 if v > 0 else -0.04), "%+.2f" % v, ha="center", fontsize=9, color="#b2182b")
for xi, v in zip(x + 0.19, troph):
    ax.text(xi, v + (0.015 if v > 0 else -0.04), "%+.2f" % v, ha="center", fontsize=9, color="#555")
ax.axhline(0, color="k", lw=0.8)
ax.set_xticks(x)
ax.set_xticklabels(subs, fontsize=9.5)
ax.set_ylabel("Correlation with consensus module\n(point-biserial)", fontsize=10)
ax.set_title("Fig.5b  Two-face subtype structure: immune complement module\n(filled red) dissociates from trophoblast module (grey)", fontsize=10.5, loc="left")
ax.legend(fontsize=8.5, frameon=False, loc="lower left")
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(FIG + "/fig5b_subtype_module_corr.png", dpi=300, bbox_inches="tight")
print("fig5b saved")
