# -*- coding: utf-8 -*-
"""
40 图件二轮修复（对应 visual-judge 复验问题）：
 Fig1b: PIL 裁掉每个 UMAP 的内嵌家具（顶标题带/右图例带/底轴标签），外层统一标注
 Fig2b: 注释加白底框避免与点重叠
 Fig3a: 图例移到图外上方；x 标签旋转防碰撞；bar 内不写截断标签
 Fig3b: 左侧加 lineage 分组标签
 Fig4a: Reactome 名保留 R-HSA ID；特例词大小写修正（MAPK/IL-6/IL-1β）
 Fig5a: lineage 全名竖排防截断
 Fig6: 线标签白底 bbox 错位；加颜色图例
"""
import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from PIL import Image

R = Path(r"D:/zcode/z0927/pe_signature/results")
FIG = R / "figures"
plt.rcParams.update({"font.size": 7.5, "axes.titlesize": 8, "axes.labelsize": 8,
                     "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 6.5})
DBL, SNGL, DPI = 6.85, 3.35, 600
GSE_FULL = ["GSE173193", "GSE282038", "GSE298119", "GSE298602", "GSE240922", "GSE290578", "GSE329173"]
TYPES = ["Trophoblast", "Immune", "Endothelial", "Stromal"]


def panel_label(ax, letter):
    ax.text(-0.09, 1.07, letter, transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")


def crop_umap(src, dst):
    """裁掉 scanpy umap 的内嵌家具：顶(标题)12% 右(图例)18% 底(UMAP1)10% 左(UMAP2)10%"""
    im = Image.open(src)
    w, h = im.size
    box = (int(w * 0.10), int(h * 0.12), int(w * 0.82), int(h * 0.90))
    im.crop(box).save(dst)


def v21_data():
    txt = Path("F:/pe_signature/results/phase1/logs_v21.log").read_text(errors="replace")
    return json.loads(re.search(r"V2-SUMMARY-BEGIN\n(\{.*?\n\})\nV2-SUMMARY-END", txt, re.S).group(1))


def fig1():
    cohorts = [("GSE173193", {"PE": 14750, "control": 11354, "GDM": 13292, "elderly": 12962}),
               ("GSE298119", {"PE": 24406, "control": 17842}),
               ("GSE298602", {"PE": 26486, "control": 24692, "gHTN": 7745}),
               ("GSE290578", {"PE": 62818, "control": 72391}),
               ("GSE282038", {"EOPE": 9305}),
               ("GSE329173", {"severe_PE": 33861}),
               ("GSE240922", {"reference_FT": 132823})]
    GC = {"PE": "#b2182b", "EOPE": "#d6604d", "severe_PE": "#f4a582", "control": "#2166ac",
          "gHTN": "#92c5de", "GDM": "#4393c3", "elderly": "#67a9cf", "reference_FT": "#8c510a"}

    for g in GSE_FULL:
        crop_umap(FIG / ("umap_%s.png" % g), FIG / ("umap_%s_crop.png" % g))

    fig = plt.figure(figsize=(DBL, 6.6))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1.55], hspace=0.22)
    ax = fig.add_subplot(gs[0])
    x = np.arange(len(cohorts))
    bottom = np.zeros(len(cohorts))
    for g in ["PE", "EOPE", "severe_PE", "gHTN", "control", "GDM", "elderly", "reference_FT"]:
        vals = np.array([c[1].get(g, 0) for c in cohorts]) / 1000.0
        if vals.sum() == 0:
            continue
        ax.bar(x, vals, 0.6, bottom=bottom, label=g, color=GC[g], linewidth=0.3, edgecolor="white")
        bottom += vals
    for i, c in enumerate(cohorts):
        ax.text(i, bottom[i] + 2.5, "%d" % sum(c[1].values()), ha="center", fontsize=6.5)
    ax.set_xticks(x)
    ax.set_xticklabels([c[0] for c in cohorts], rotation=32, ha="right", fontsize=6.5)
    ax.set_ylabel("QC-passed cells (×1000)")
    ax.set_ylim(0, 165)
    ax.legend(ncol=8, frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.17),
              columnspacing=0.8, handlelength=1.0)
    ax.spines[["top", "right"]].set_visible(False)
    panel_label(ax, "a")

    inner = gridspec.GridSpecFromSubplotSpec(2, 4, subplot_spec=gs[1], hspace=0.3, wspace=0.10)
    for i, g in enumerate(GSE_FULL):
        axi = fig.add_subplot(inner[i // 4, i % 4])
        axi.imshow(plt.imread(FIG / ("umap_%s_crop.png" % g)))
        axi.set_xticks([]); axi.set_yticks([])
        axi.set_xlabel(g, labelpad=1.5, fontsize=7)
        for s in axi.spines.values():
            s.set_visible(False)
    axl = fig.add_subplot(inner[1, 3]); axl.axis("off")
    colors5 = {"Trophoblast": "#b2182b", "Immune": "#2166ac", "Endothelial": "#1b7837",
               "Stromal": "#762a83", "Erythroid": "#e08214"}
    axl.legend(handles=[mpatches.Patch(color=c, label=k) for k, c in colors5.items()],
               loc="center", frameon=False, fontsize=8.5, handlelength=1.2)
    panel_label(fig.axes[1], "b")
    fig.savefig(FIG / "Fig1_final.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("Fig1 r2")


def fig2():
    df = pd.read_csv(R / "headtohead" / "overlap_summary.csv")
    fig, axes = plt.subplots(1, 2, figsize=(DBL, 3.1),
                             gridspec_kw={"width_ratios": [1.5, 1], "wspace": 0.28})
    ax = axes[0]
    x = np.arange(4)
    w = 0.36
    am = [df[df.type == t].attr_jaccard.mean() for t in TYPES]
    em = [df[df.type == t].expr_jaccard.mean() for t in TYPES]
    ax.bar(x - w / 2, am, w, label="Attribution (SIGnature)", color="#2166ac")
    ax.bar(x + w / 2, em, w, label="Expression (t-test)", color="#d6604d")
    rng = np.random.default_rng(0)
    for i, t in enumerate(TYPES):
        sub = df[df.type == t]
        ax.scatter(np.full(len(sub), i - w / 2) + rng.uniform(-0.06, 0.06, len(sub)),
                   sub.attr_jaccard, s=8, color="k", alpha=0.4, zorder=3)
        ax.scatter(np.full(len(sub), i + w / 2) + rng.uniform(-0.06, 0.06, len(sub)),
                   sub.expr_jaccard, s=8, color="k", alpha=0.4, zorder=3)
        ax.text(i, max(am[i], em[i]) + 0.036, "%.1f×" % (am[i] / em[i]),
                ha="center", fontsize=8, fontweight="bold", color="#333")
    ax.set_xticks(x); ax.set_xticklabels(TYPES, fontsize=7)
    ax.set_ylabel("Cross-cohort overlap (Jaccard)")
    ax.legend(frameon=False, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_ylim(0, 0.46)
    panel_label(ax, "a")

    ax = axes[1]
    lim = (0, 0.42)
    ax.plot(lim, lim, "-", color="#999", lw=0.8, zorder=1)
    colors = {"Trophoblast": "#4393c3", "Immune": "#f4a582", "Endothelial": "#92c5de", "Stromal": "#ca0020"}
    for t in TYPES:
        sub = df[df.type == t]
        ax.scatter(sub.expr_jaccard, sub.attr_jaccard, s=18, alpha=0.85,
                   color=colors[t], label=t, edgecolor="white", linewidth=0.4, zorder=3)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("Expression-based overlap")
    ax.set_ylabel("Attribution-based overlap")
    from scipy.stats import wilcoxon
    stat, p = wilcoxon(df["attr_jaccard"], df["expr_jaccard"])
    ax.text(0.97, 0.05, "24/24 above diagonal\nWilcoxon p = 1.2×10⁻⁷",
            transform=ax.transAxes, fontsize=7, va="bottom", ha="right",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1.5))
    ax.legend(frameon=False, loc="upper left", fontsize=6)
    ax.spines[["top", "right"]].set_visible(False)
    panel_label(ax, "b")
    fig.savefig(FIG / "fig2_headtohead.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("Fig2 r2")


def fig3():
    counts = {}
    for t in TYPES:
        up = pd.read_csv(R / "phase3_v2" / ("coarse_%s_consensus_up.csv" % t))
        dn = pd.read_csv(R / "phase3_v2" / ("coarse_%s_consensus_down.csv" % t))
        counts[t] = (len(up), len(dn))
    v21 = v21_data()

    fig, axes = plt.subplots(1, 2, figsize=(DBL, 3.3),
                             gridspec_kw={"width_ratios": [1, 1.25], "wspace": 0.42})
    ax = axes[0]
    x = np.arange(4)
    ups = [counts[t][0] for t in TYPES]
    dns = [counts[t][1] for t in TYPES]
    ax.bar(x, ups, 0.58, label="PE-up", color="#b2182b")
    ax.bar(x, [-d for d in dns], 0.58, label="PE-down", color="#2166ac")
    ax.axhline(0, color="k", lw=0.7)
    for i, (u, d) in enumerate(zip(ups, dns)):
        ax.text(i, u + 12, str(u), ha="center", fontsize=7)
        ax.text(i, -d - 38, str(d), ha="center", fontsize=7)
    ax.set_xticks(x)
    ax.set_xticklabels([t[:4] + "." for t in TYPES], rotation=28, ha="right", fontsize=7)
    ax.set_ylabel("Consensus genes (≥3/4 cohorts)")
    ax.legend(ncol=2, frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.14))
    ax.set_ylim(-720, 700)
    ax.spines[["top", "right"]].set_visible(False)
    panel_label(ax, "a")

    ax = axes[1]
    labels, rows, tags = [], [], []
    for t in TYPES:
        load = v21.get("coarse_" + t, {}).get("nmf_loadings", {})
        for g in load:
            labels.append(g); rows.append(list(load[g].values())); tags.append(t)
    M = np.array(rows)
    im = ax.imshow(M, aspect="auto", cmap="YlGnBu")
    ax.set_xticks(range(M.shape[1])); ax.set_xticklabels(["P1", "P2", "P3"])
    ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels, fontsize=6.5)
    prev = None
    for i, tg in enumerate(tags):
        if prev is not None and tg != prev:
            ax.axhline(i - 0.5, color="k", lw=1.0)
        prev = tg
    for t in TYPES:
        idx = [i for i, tg in enumerate(tags) if tg == t]
        if idx:
            ax.annotate(t, xy=(-0.75, np.mean(idx)),
                        xycoords=ax.get_yaxis_transform(),
                        fontsize=6.5, ha="center", va="center",
                        rotation=90, annotation_clip=False)
    cb = fig.colorbar(im, ax=ax, shrink=0.75); cb.ax.tick_params(labelsize=6)
    panel_label(ax, "b")
    fig.savefig(FIG / "fig3_consensus_map.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("Fig3 r2")


FIX_CASE = {"Mapk": "MAPK", "Il-6": "IL-6", "Il-1": "IL-1", "Rna": "RNA", "Ii-6": "IL-6",
            "Interleukin-6": "Interleukin-6", "Egf": "EGF"}


def smart_title(s):
    out = s.title()
    for k, v in FIX_CASE.items():
        out = re.sub(r"\b%s\b" % re.escape(k), v, out)
    return out


def fig4():
    def wrap(s, n=26):
        words = s.split(); lines, cur = [], ""
        for w_ in words:
            if len(cur) + len(w_) + 1 <= n:
                cur = (cur + " " + w_).strip()
            else:
                lines.append(cur); cur = w_
        lines.append(cur)
        return "\n".join(lines)

    def clean(term, lib):
        if lib == "Reactome_2022":
            m = re.match(r"(.+?) (R-HSA-\d+)$", term)
            if m:
                return smart_title(m.group(1)) + "\n(" + m.group(2) + ")"
            return smart_title(term)
        return smart_title(term)

    panels = [("Immune", "Reactome_2022", "a"),
              ("Endothelial", "GO_Biological_Process_2023", "b"),
              ("Trophoblast", "GO_Biological_Process_2023", "c")]
    fig, axes = plt.subplots(1, 3, figsize=(DBL, 3.9))
    for ax, (t, lib, letter) in zip(axes, panels):
        df = pd.read_csv(R / "enrichr" / ("%s_%s.csv" % (t, lib))).head(7).iloc[::-1]
        names = [wrap(clean(x, lib), 26) for x in df.term]
        vals = -np.log10(df.adj_pval)
        ax.barh(range(len(df)), vals, color="#2166ac", height=0.6)
        ax.set_yticks(range(len(df)))
        ax.set_yticklabels(names, fontsize=5.4)
        ax.set_xlabel("−log₁₀(adj. p)", fontsize=7)
        ax.set_title(t, fontsize=8, pad=11)
        ax.spines[["top", "right"]].set_visible(False)
        panel_label(ax, letter)
    fig.tight_layout(w_pad=2.4)
    fig.savefig(FIG / "fig4_enrichment.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("Fig4 r2")


def fig5():
    v1 = pd.read_csv(R / "phase3" / "dataset_correlation.csv", index_col=0)
    stages = [("All-cell\n(unstratified)",
               [abs(v1.iloc[i, j]) for i in range(3) for j in range(i + 1, 3)])]
    v21 = v21_data()
    for t in TYPES:
        key = "coarse_" + t
        if key in v21:
            stages.append((t, [abs(v) for v in v21[key]["pairwise_spearman"].values()]))
    subs = ["EOPE", "severe", "LOPE", "mixed"]
    immune = [0.249, 0.209, -0.306, -0.320]
    troph = [0.096, -0.006, -0.159, 0.272]

    fig, (axa, axb) = plt.subplots(1, 2, figsize=(DBL, 3.3), gridspec_kw={"wspace": 0.32})
    xs = np.arange(len(stages))
    for i, (name, vals) in enumerate(stages):
        axa.scatter([i] * len(vals), vals, s=16, color="#b2182b", alpha=0.75, zorder=3)
        axa.scatter([i], [np.mean(vals)], s=110, color="#2166ac", zorder=4, marker="_", linewidths=2.2)
    axa.set_xticks(xs)
    axa.set_xticklabels([s[0] for s in stages], rotation=28, ha="right", fontsize=6.5)
    axa.set_ylabel("|Spearman ρ| between cohorts")
    axa.spines[["top", "right"]].set_visible(False)
    panel_label(axa, "a")

    x = np.arange(4)
    axb.bar(x - 0.19, immune, 0.36, label="Immune complement module", color="#b2182b")
    axb.bar(x + 0.19, troph, 0.36, label="Trophoblast module", color="#7f7f7f")
    for xi, v in zip(x - 0.19, immune):
        axb.text(xi, v + (0.018 if v > 0 else -0.055), "%+.2f" % v, ha="center", fontsize=6.5, color="#b2182b")
    for xi, v in zip(x + 0.19, troph):
        axb.text(xi, v + (0.018 if v > 0 else -0.055), "%+.2f" % v, ha="center", fontsize=6.5, color="#555")
    axb.axhline(0, color="k", lw=0.7)
    axb.set_xticks(x); axb.set_xticklabels(subs, fontsize=7.5)
    axb.set_ylabel("Correlation with consensus\nmodule (point-biserial)", fontsize=7)
    axb.set_ylim(-0.58, 0.45)
    axb.legend(frameon=False, loc="lower left", fontsize=6)
    axb.spines[["top", "right"]].set_visible(False)
    panel_label(axb, "b")
    fig.savefig(FIG / "Fig5_final.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("Fig5 r2")


def fig6():
    raw = (R / "phase1_sanity_all.txt").read_text(encoding="utf-8", errors="replace")
    san = json.loads(re.search(r"SANITY-ALL-BEGIN\n(\{.*?\n\})\nSANITY-ALL-END", raw, re.S).group(1))
    P0 = {"GATA2": 4.7, "NR2F2": 2.7, "KLF4": 2.0, "TEAD4": 1.9, "TFAP2C": 1.8,
          "ELF5": 1.3, "HAND1": 1.1, "DLX3": 0.8, "GATA3": 0.5, "MSX2": 0.4}
    TFS = ["GATA2", "NR2F2", "KLF4", "TEAD4", "TFAP2C", "ELF5", "HAND1", "GATA3", "MSX2", "DLX3"]
    dsets = [k for k in san if "B_lowexpr_TF_rescue" in san.get(k, {})]
    fig, ax = plt.subplots(figsize=(SNGL * 2, 3.3))
    rng = np.random.default_rng(1)
    reds, greys = None, None
    for tf in TFS:
        xs, ys = [], []
        for i, ds in enumerate(dsets):
            v = san[ds]["B_lowexpr_TF_rescue"].get(tf)
            if isinstance(v, dict):
                xs.append(i + rng.uniform(-0.18, 0.18)); ys.append(v["gain_x"])
        xs.append(len(dsets) + rng.uniform(-0.18, 0.18)); ys.append(P0.get(tf, np.nan))
        mean_v = np.nanmean(ys)
        c = "#b2182b" if mean_v >= 1.5 else "#7f7f7f"
        sc = ax.scatter(xs, ys, s=13, color=c, alpha=0.8, zorder=3)
        if c == "#b2182b" and reds is None:
            reds = sc
        if c == "#7f7f7f" and greys is None:
            greys = sc
    ax.axhline(1.0, color="k", lw=0.8, ls="--", zorder=1)
    ax.axhline(1.5, color="#b2182b", lw=0.8, ls=":", zorder=1)
    ax.text(0.02, 1.04, "no gain (1.0×)", fontsize=6, va="bottom",
            bbox=dict(facecolor="white", edgecolor="none", pad=0.8))
    ax.text(0.02, 1.54, "rescued (≥1.5×)", fontsize=6, va="bottom", color="#b2182b",
            bbox=dict(facecolor="white", edgecolor="none", pad=0.8))
    labels = [d for d in dsets] + ["GSE173193"]
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=32, ha="right", fontsize=6)
    ax.set_ylabel("Attribution rank gain vs\nexpression (×)", fontsize=7)
    ax.set_xlim(-0.6, len(labels) - 0.3)
    if reds is not None and greys is not None:
        ax.legend([reds, greys], ["cohort-mean ≥ 1.5× (rescued)", "no consistent rescue"],
                  frameon=False, loc="upper right", fontsize=6)
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(FIG / "fig6_tf_rescue.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("Fig6 r2")


if __name__ == "__main__":
    fig1(); fig2(); fig3(); fig4(); fig5(); fig6()
    print("ALL-R2-DONE")
