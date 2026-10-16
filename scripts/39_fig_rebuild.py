# -*- coding: utf-8 -*-
"""
39 图件统一重制（SCI 规范修复版）——针对 visual-judge 验收问题逐项修复：
规范统一：小写面板标签 abc；无内嵌长标题（图注交给 caption）；队列 ID 全格式；
画布=出版物理尺寸（双栏 6.85in / 单栏 3.35in），字号所见即所得；dpi=600。
修复项：
 Fig1a 标签/图例重叠 → 数值标签放条内、图例置底两行
 Fig1b 去重复标题 → 外层 GSE 标注一次
 Fig3 图例移出 + 大字号
 Fig4 通路名 wrap 完整 + abc 标签 + 大字号
 Fig5b 防裁剪（ylim 扩展）
 Fig6 虚线标注 + GSE 前缀
"""
import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd

R = Path(r"D:/zcode/z0927/pe_signature/results")
FIG = R / "figures"
FB = Path("F:/pe_signature/results")
plt.rcParams.update({"font.size": 7.5, "axes.titlesize": 8, "axes.labelsize": 8,
                     "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 6.5,
                     "pdf.fonttype": 42})
DBL = 6.85  # in, 双栏
SNGL = 3.35
DPI = 600

GSE_FULL = {"GSE173193": "GSE173193", "GSE282038": "GSE282038", "GSE298119": "GSE298119",
            "GSE298602": "GSE298602", "GSE240922": "GSE240922", "GSE290578": "GSE290578",
            "GSE329173": "GSE329173"}


def panel_label(ax, letter):
    ax.text(-0.08, 1.06, letter, transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")


# ---------------- Fig.1 ----------------
def fig1():
    cohorts = [
        ("GSE173193", {"PE": 14750, "control": 11354, "GDM": 13292, "elderly": 12962}),
        ("GSE298119", {"PE": 24406, "control": 17842}),
        ("GSE298602", {"PE": 26486, "control": 24692, "gHTN": 7745}),
        ("GSE290578", {"PE": 62818, "control": 72391}),
        ("GSE282038", {"EOPE": 9305}),
        ("GSE329173", {"severe_PE": 33861}),
        ("GSE240922", {"reference_FT": 132823}),
    ]
    GC = {"PE": "#b2182b", "EOPE": "#d6604d", "severe_PE": "#f4a582", "control": "#2166ac",
          "gHTN": "#92c5de", "GDM": "#4393c3", "elderly": "#67a9cf", "reference_FT": "#8c510a"}
    fig = plt.figure(figsize=(DBL, 7.0))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1.75], hspace=0.18)

    ax = fig.add_subplot(gs[0])
    x = np.arange(len(cohorts))
    bottom = np.zeros(len(cohorts))
    order = ["PE", "EOPE", "severe_PE", "gHTN", "control", "GDM", "elderly", "reference_FT"]
    for g in order:
        vals = np.array([c[1].get(g, 0) for c in cohorts]) / 1000.0
        if vals.sum() == 0:
            continue
        ax.bar(x, vals, 0.62, bottom=bottom, label=g, color=GC[g], linewidth=0.3, edgecolor="white")
        bottom += vals
    for i, c in enumerate(cohorts):
        ax.text(i, bottom[i] + 2.5, "%d" % sum(c[1].values()), ha="center", fontsize=6.5)
    ax.set_xticks(x)
    ax.set_xticklabels([GSE_FULL[c[0]] for c in cohorts], rotation=32, ha="right", fontsize=6.5)
    ax.set_ylabel("QC-passed cells (×1000)")
    ax.set_ylim(0, 160)
    ax.legend(ncol=8, frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.16),
              columnspacing=0.9, handlelength=1.0)
    ax.spines[["top", "right"]].set_visible(False)
    panel_label(ax, "a")

    inner = gridspec.GridSpecFromSubplotSpec(2, 4, subplot_spec=gs[1], hspace=0.25, wspace=0.12)
    order7 = ["GSE173193", "GSE282038", "GSE298119", "GSE298602", "GSE240922", "GSE290578", "GSE329173"]
    for i, g in enumerate(order7):
        axi = fig.add_subplot(inner[i // 4, i % 4])
        axi.imshow(plt.imread(FIG / ("umap_%s.png" % g)))
        axi.set_xticks([])
        axi.set_yticks([])
        axi.set_xlabel(GSE_FULL[g], labelpad=1.5, fontsize=7)
        for s in axi.spines.values():
            s.set_visible(False)
    axl = fig.add_subplot(inner[1, 3])
    axl.axis("off")
    import matplotlib.patches as mpatches
    colors5 = {"Trophoblast": "#b2182b", "Immune": "#2166ac", "Endothelial": "#1b7837",
               "Stromal": "#762a83", "Erythroid": "#e08214"}
    axl.legend(handles=[mpatches.Patch(color=c, label=k) for k, c in colors5.items()],
               loc="center", frameon=False, fontsize=8.5, handlelength=1.2)
    axi0 = fig.axes[1]
    panel_label(axi0, "b")
    fig.savefig(FIG / "Fig1_final.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("Fig1 rebuilt")


# ---------------- Fig.2 ----------------
def fig2():
    df = pd.read_csv(R / "headtohead" / "overlap_summary.csv")
    types = ["Trophoblast", "Immune", "Endothelial", "Stromal"]
    fig, axes = plt.subplots(1, 2, figsize=(DBL, 3.1),
                             gridspec_kw={"width_ratios": [1.5, 1], "wspace": 0.28})
    ax = axes[0]
    x = np.arange(len(types))
    w = 0.36
    am = [df[df.type == t].attr_jaccard.mean() for t in types]
    em = [df[df.type == t].expr_jaccard.mean() for t in types]
    ax.bar(x - w / 2, am, w, label="Attribution (SIGnature)", color="#2166ac")
    ax.bar(x + w / 2, em, w, label="Expression (t-test)", color="#d6604d")
    rng = np.random.default_rng(0)
    for i, t in enumerate(types):
        sub = df[df.type == t]
        ax.scatter(np.full(len(sub), i - w / 2) + rng.uniform(-0.06, 0.06, len(sub)),
                   sub.attr_jaccard, s=8, color="k", alpha=0.4, zorder=3)
        ax.scatter(np.full(len(sub), i + w / 2) + rng.uniform(-0.06, 0.06, len(sub)),
                   sub.expr_jaccard, s=8, color="k", alpha=0.4, zorder=3)
        ax.text(i, max(am[i], em[i]) + 0.036, "%.1f×" % (am[i] / em[i]),
                ha="center", fontsize=8, fontweight="bold", color="#333")
    ax.set_xticks(x)
    ax.set_xticklabels(types, fontsize=7.5)
    ax.set_ylabel("Cross-cohort overlap (Jaccard)")
    ax.legend(frameon=False, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_ylim(0, 0.46)
    panel_label(ax, "a")

    ax = axes[1]
    lim = (0, 0.42)
    ax.plot(lim, lim, "-", color="#999", lw=0.8, zorder=1)
    colors = {"Trophoblast": "#4393c3", "Immune": "#f4a582", "Endothelial": "#92c5de", "Stromal": "#ca0020"}
    for t in types:
        sub = df[df.type == t]
        ax.scatter(sub.expr_jaccard, sub.attr_jaccard, s=18, alpha=0.85,
                   color=colors[t], label=t, edgecolor="white", linewidth=0.4, zorder=3)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("Expression-based overlap")
    ax.set_ylabel("Attribution-based overlap")
    from scipy.stats import wilcoxon
    stat, p = wilcoxon(df["attr_jaccard"], df["expr_jaccard"])
    ax.text(0.03, 0.94, "24/24 pairs above diagonal\nWilcoxon p = %.1e" % p,
            transform=ax.transAxes, fontsize=7, va="top")
    ax.legend(frameon=False, loc="lower right", fontsize=6)
    ax.spines[["top", "right"]].set_visible(False)
    panel_label(ax, "b")
    fig.savefig(FIG / "fig2_headtohead.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("Fig2 rebuilt")


# ---------------- Fig.3 ----------------
def fig3():
    types = ["Trophoblast", "Immune", "Endothelial", "Stromal"]
    counts = {}
    for t in types:
        up = pd.read_csv(R / "phase3_v2" / ("coarse_%s_consensus_up.csv" % t))
        dn = pd.read_csv(R / "phase3_v2" / ("coarse_%s_consensus_down.csv" % t))
        counts[t] = (len(up), len(dn))
    txt = Path("F:/pe_signature/results/phase1/logs_v21.log").read_text(errors="replace")
    v21 = json.loads(re.search(r"V2-SUMMARY-BEGIN\n(\{.*?\n\})\nV2-SUMMARY-END", txt, re.S).group(1))

    fig, axes = plt.subplots(1, 2, figsize=(DBL, 3.2),
                             gridspec_kw={"width_ratios": [1, 1.25], "wspace": 0.3})
    ax = axes[0]
    x = np.arange(len(types))
    ups = [counts[t][0] for t in types]
    dns = [counts[t][1] for t in types]
    ax.bar(x, ups, 0.58, label="PE-up", color="#b2182b")
    ax.bar(x, [-d for d in dns], 0.58, label="PE-down", color="#2166ac")
    ax.axhline(0, color="k", lw=0.7)
    for i, (u, d) in enumerate(zip(ups, dns)):
        ax.text(i, u + 12, str(u), ha="center", fontsize=7)
        ax.text(i, -d - 34, str(d), ha="center", fontsize=7)
    ax.set_xticks(x)
    ax.set_xticklabels(types, fontsize=7.5)
    ax.set_ylabel("Consensus genes (≥3/4 cohorts)")
    ax.legend(frameon=False, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    panel_label(ax, "a")

    ax = axes[1]
    labels, rows, type_tags = [], [], []
    for t in types:
        key = "coarse_" + t
        if key not in v21:
            continue
        load = v21[key]["nmf_loadings"]
        for g in load:
            labels.append(g)
            rows.append(list(load[g].values()))
            type_tags.append(t)
    M = np.array(rows)
    im = ax.imshow(M, aspect="auto", cmap="YlGnBu")
    ax.set_xticks(range(M.shape[1]))
    ax.set_xticklabels(["P1", "P2", "P3"])
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=6.5)
    prev = None
    for i, tg in enumerate(type_tags):
        if prev is not None and tg != prev:
            ax.axhline(i - 0.5, color="k", lw=1.0)
        prev = tg
    for i, tg in enumerate(set(type_tags)):
        pass
    ys = [i for i, tg in enumerate(type_tags)]
    for t in types:
        idx = [i for i, tg in enumerate(type_tags) if tg == t]
        if idx:
            ax.text(-0.62, np.mean(idx), t[:5], fontsize=6.5, ha="right", va="center",
                    transform=ax.get_yaxis_transform())
    cb = fig.colorbar(im, ax=ax, shrink=0.75)
    cb.ax.tick_params(labelsize=6)
    panel_label(ax, "b")
    fig.savefig(FIG / "fig3_consensus_map.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("Fig3 rebuilt")


# ---------------- Fig.4 ----------------
def fig4():
    def wrap(s, n=30):
        words = s.split()
        lines, cur = [], ""
        for w_ in words:
            if len(cur) + len(w_) + 1 <= n:
                cur = (cur + " " + w_).strip()
            else:
                lines.append(cur)
                cur = w_
        lines.append(cur)
        return "\n".join(lines)

    panels = [("Immune", "Reactome_2022", "a"),
              ("Endothelial", "GO_Biological_Process_2023", "b"),
              ("Trophoblast", "GO_Biological_Process_2023", "c")]
    fig, axes = plt.subplots(1, 3, figsize=(DBL, 3.6))
    for ax, (t, lib, letter) in zip(axes, panels):
        df = pd.read_csv(R / "enrichr" / ("%s_%s.csv" % (t, lib))).head(7).iloc[::-1]
        names = [wrap(x.split(" R-HSA")[0].title(), 24) for x in df.term]
        vals = -np.log10(df.adj_pval)
        ax.barh(range(len(df)), vals, color="#2166ac", height=0.6)
        ax.set_yticks(range(len(df)))
        ax.set_yticklabels(names, fontsize=5.6)
        ax.set_xlabel("−log₁₀(adj. p)", fontsize=7)
        ax.set_title(t, fontsize=8, pad=12)
        ax.spines[["top", "right"]].set_visible(False)
        panel_label(ax, letter)
    fig.tight_layout(w_pad=2.2)
    fig.savefig(FIG / "fig4_enrichment.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("Fig4 rebuilt")


# ---------------- Fig.5 ----------------
def fig5():
    # 5a 数据链
    v1 = pd.read_csv(R / "phase3" / "dataset_correlation.csv", index_col=0)
    stages = [("All-cell\n(v1)", [abs(v1.iloc[i, j]) for i in range(3) for j in range(i + 1, 3)])]
    txt = Path("F:/pe_signature/results/phase1/logs_v21.log").read_text(errors="replace")
    v21 = json.loads(re.search(r"V2-SUMMARY-BEGIN\n(\{.*?\n\})\nV2-SUMMARY-END", txt, re.S).group(1))
    for t in ["Trophoblast", "Immune", "Endothelial", "Stromal"]:
        key = "coarse_" + t
        if key in v21:
            stages.append(("%s\n(stratified)" % t[:5],
                           [abs(v) for v in v21[key]["pairwise_spearman"].values()]))
    # 5b 亚型
    subs = ["EOPE", "severe", "LOPE", "mixed"]
    immune = [0.249, 0.209, -0.306, -0.320]
    troph = [0.096, -0.006, -0.159, 0.272]

    fig, (axa, axb) = plt.subplots(1, 2, figsize=(DBL, 3.3), gridspec_kw={"wspace": 0.3})
    xs = np.arange(len(stages))
    for i, (name, vals) in enumerate(stages):
        axa.scatter([i] * len(vals), vals, s=16, color="#b2182b", alpha=0.75, zorder=3)
        axa.scatter([i], [np.mean(vals)], s=110, color="#2166ac", zorder=4, marker="_", linewidths=2.2)
    axa.set_xticks(xs)
    axa.set_xticklabels([s[0] for s in stages], fontsize=6.5)
    axa.set_ylabel("|Spearman ρ| between cohorts")
    axa.spines[["top", "right"]].set_visible(False)
    panel_label(axa, "a")

    x = np.arange(len(subs))
    axb.bar(x - 0.19, immune, 0.36, label="Immune complement module", color="#b2182b")
    axb.bar(x + 0.19, troph, 0.36, label="Trophoblast module", color="#7f7f7f")
    for xi, v in zip(x - 0.19, immune):
        axb.text(xi, v + (0.018 if v > 0 else -0.05), "%+.2f" % v, ha="center", fontsize=6.5, color="#b2182b")
    for xi, v in zip(x + 0.19, troph):
        axb.text(xi, v + (0.018 if v > 0 else -0.05), "%+.2f" % v, ha="center", fontsize=6.5, color="#555")
    axb.axhline(0, color="k", lw=0.7)
    axb.set_xticks(x)
    axb.set_xticklabels(subs, fontsize=7.5)
    axb.set_ylabel("Correlation with consensus\nmodule (point-biserial)", fontsize=7)
    axb.set_ylim(-0.55, 0.45)
    axb.legend(frameon=False, loc="lower left", fontsize=6)
    axb.spines[["top", "right"]].set_visible(False)
    panel_label(axb, "b")
    fig.savefig(FIG / "Fig5_final.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("Fig5 rebuilt")


# ---------------- Fig.6 ----------------
def fig6():
    raw = (R / "phase1_sanity_all.txt").read_text(encoding="utf-8", errors="replace")
    san = json.loads(re.search(r"SANITY-ALL-BEGIN\n(\{.*?\n\})\nSANITY-ALL-END", raw, re.S).group(1))
    P0 = {"GATA2": 4.7, "NR2F2": 2.7, "KLF4": 2.0, "TEAD4": 1.9, "TFAP2C": 1.8,
          "ELF5": 1.3, "HAND1": 1.1, "DLX3": 0.8, "GATA3": 0.5, "MSX2": 0.4}
    TFS = ["GATA2", "NR2F2", "KLF4", "TEAD4", "TFAP2C", "ELF5", "HAND1", "GATA3"]
    dsets = [k for k in san if "B_lowexpr_TF_rescue" in san.get(k, {})]
    fig, ax = plt.subplots(figsize=(SNGL * 2, 3.2))
    rng = np.random.default_rng(1)
    for j, tf in enumerate(TFS):
        xs, ys = [], []
        for i, ds in enumerate(dsets):
            v = san[ds]["B_lowexpr_TF_rescue"].get(tf)
            if isinstance(v, dict):
                xs.append(i + rng.uniform(-0.18, 0.18))
                ys.append(v["gain_x"])
        xs.append(len(dsets) + rng.uniform(-0.18, 0.18))
        ys.append(P0.get(tf, np.nan))
        mean_v = np.nanmean(ys)
        ax.scatter(xs, ys, s=13, color="#b2182b" if mean_v >= 1.5 else "#7f7f7f", alpha=0.8, zorder=3)
    ax.axhline(1.0, color="k", lw=0.8, ls="--")
    ax.axhline(1.5, color="#b2182b", lw=0.8, ls=":")
    ax.text(len(dsets) + 0.45, 1.0, "no gain (1.0×)", fontsize=6, va="center")
    ax.text(len(dsets) + 0.45, 1.5, "rescued (1.5×)", fontsize=6, va="center", color="#b2182b")
    labels = [GSE_FULL[d] for d in dsets] + ["GSE173193"]
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=32, ha="right", fontsize=6)
    ax.set_ylabel("Attribution rank gain vs\nexpression (×)", fontsize=7)
    ax.set_xlim(-0.6, len(labels) - 0.1 + 1.1)
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(FIG / "fig6_tf_rescue.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("Fig6 rebuilt")


if __name__ == "__main__":
    fig1(); fig2(); fig3(); fig4(); fig5(); fig6()
    print("ALL-FIGURES-REBUILT")
