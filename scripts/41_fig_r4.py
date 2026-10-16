# -*- coding: utf-8 -*-
"""41 图件 R4 精准修复（仅 4 图）：
Fig1: UMAP 右侧裁 18%→30%（图例点列裁净）
Fig2b: 注释放左下（对角线下方无数据），图例右下
Fig3a: 去 ylim 截条；删竖标签（x 轴已表意）；Fig3b: lineage 竖标签贴近热图（x=-0.32）
Fig6: 参考线标签移到数据区右侧空白
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
from PIL import Image

R = Path(r"D:/zcode/z0927/pe_signature/results")
FIG = R / "figures"
plt.rcParams.update({"font.size": 7.5, "axes.labelsize": 8,
                     "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 6.5})
DBL, SNGL, DPI = 6.85, 3.35, 600
GSE_FULL = ["GSE173193", "GSE282038", "GSE298119", "GSE298602", "GSE240922", "GSE290578", "GSE329173"]
TYPES = ["Trophoblast", "Immune", "Endothelial", "Stromal"]


def panel_label(ax, letter):
    ax.text(-0.09, 1.07, letter, transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")


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
        im = Image.open(FIG / ("umap_%s.png" % g))
        w, h = im.size
        im.crop((int(w * 0.10), int(h * 0.12), int(w * 0.70), int(h * 0.90))).save(
            FIG / ("umap_%s_crop.png" % g))

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
    import matplotlib.patches as mpatches
    axl = fig.add_subplot(inner[1, 3]); axl.axis("off")
    colors5 = {"Trophoblast": "#b2182b", "Immune": "#2166ac", "Endothelial": "#1b7837",
               "Stromal": "#762a83", "Erythroid": "#e08214"}
    axl.legend(handles=[mpatches.Patch(color=c, label=k) for k, c in colors5.items()],
               loc="center", frameon=False, fontsize=8.5, handlelength=1.2)
    panel_label(fig.axes[1], "b")
    fig.savefig(FIG / "Fig1_final.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("Fig1 r4")


def fig2():
    df = pd.read_csv(R / "headtohead" / "overlap_summary.csv")
    fig, axes = plt.subplots(1, 2, figsize=(DBL, 3.1),
                             gridspec_kw={"width_ratios": [1.5, 1], "wspace": 0.30})
    ax = axes[0]
    x = np.arange(4); w = 0.36
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
    ax.text(0.04, 0.06, "24/24 above diagonal\nWilcoxon p = 1.2×10⁻⁷",
            transform=ax.transAxes, fontsize=7, va="bottom", ha="left",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.9, pad=1.5))
    ax.legend(frameon=False, loc="upper right", fontsize=6,
              bbox_to_anchor=(1.02, 1.02))
    ax.spines[["top", "right"]].set_visible(False)
    panel_label(ax, "b")
    fig.savefig(FIG / "fig2_headtohead.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("Fig2 r4")


def fig3():
    counts = {}
    for t in TYPES:
        up = pd.read_csv(R / "phase3_v2" / ("coarse_%s_consensus_up.csv" % t))
        dn = pd.read_csv(R / "phase3_v2" / ("coarse_%s_consensus_down.csv" % t))
        counts[t] = (len(up), len(dn))
    v21 = v21_data()

    fig, axes = plt.subplots(1, 2, figsize=(DBL, 3.3),
                             gridspec_kw={"width_ratios": [1, 1.25], "wspace": 0.50})
    ax = axes[0]
    x = np.arange(4)
    ups = [counts[t][0] for t in TYPES]
    dns = [counts[t][1] for t in TYPES]
    ax.bar(x, ups, 0.58, label="PE-up", color="#b2182b")
    ax.bar(x, [-d for d in dns], 0.58, label="PE-down", color="#2166ac")
    ax.axhline(0, color="k", lw=0.7)
    for i, (u, d) in enumerate(zip(ups, dns)):
        ax.text(i, u + 25, str(u), ha="center", fontsize=7)
        ax.text(i, -d - 45, str(d), ha="center", fontsize=7)
    ax.set_xticks(x)
    ax.set_xticklabels([t[:4] + "." for t in TYPES], rotation=28, ha="right", fontsize=7)
    ax.set_ylabel("Consensus genes (≥3/4 cohorts)")
    ax.legend(ncol=2, frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.14))
    ax.set_ylim(-820, 1250)
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
            ax.annotate(t[:4] + ".", xy=(-0.32, np.mean(idx)),
                        xycoords=ax.get_yaxis_transform(),
                        fontsize=6.5, ha="center", va="center", rotation=90,
                        annotation_clip=False)
    cb = fig.colorbar(im, ax=ax, shrink=0.75); cb.ax.tick_params(labelsize=6)
    panel_label(ax, "b")
    fig.savefig(FIG / "fig3_consensus_map.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("Fig3 r4")


def fig6():
    raw = (R / "phase1_sanity_all.txt").read_text(encoding="utf-8", errors="replace")
    san = json.loads(re.search(r"SANITY-ALL-BEGIN\n(\{.*?\n\})\nSANITY-ALL-END", raw, re.S).group(1))
    P0 = {"GATA2": 4.7, "NR2F2": 2.7, "KLF4": 2.0, "TEAD4": 1.9, "TFAP2C": 1.8,
          "ELF5": 1.3, "HAND1": 1.1, "DLX3": 0.8, "GATA3": 0.5, "MSX2": 0.4}
    TFS = ["GATA2", "NR2F2", "KLF4", "TEAD4", "TFAP2C", "ELF5", "HAND1", "GATA3"]
    dsets = [k for k in san if "B_lowexpr_TF_rescue" in san.get(k, {})]
    n_total = len(dsets) + 1
    fig, ax = plt.subplots(figsize=(SNGL * 2 + 0.9, 3.3))
    rng = np.random.default_rng(1)
    reds = greys = None
    for tf in TFS:
        xs, ys = [], []
        for i, ds in enumerate(dsets):
            v = san[ds]["B_lowexpr_TF_rescue"].get(tf)
            if isinstance(v, dict):
                xs.append(i + rng.uniform(-0.18, 0.18)); ys.append(v["gain_x"])
        xs.append(len(dsets) + rng.uniform(-0.18, 0.18)); ys.append(P0.get(tf, np.nan))
        c = "#b2182b" if np.nanmean(ys) >= 1.5 else "#7f7f7f"
        sc = ax.scatter(xs, ys, s=13, color=c, alpha=0.8, zorder=3)
        if c == "#b2182b" and reds is None:
            reds = sc
        if c == "#7f7f7f" and greys is None:
            greys = sc
    ax.axhline(1.0, color="k", lw=0.8, ls="--", zorder=1)
    ax.axhline(1.5, color="#b2182b", lw=0.8, ls=":", zorder=1)
    ax.text(n_total - 0.35, 1.0, "no gain\n(1.0×)", fontsize=6, va="center", ha="left")
    ax.text(n_total - 0.35, 1.5, "rescued\n(≥1.5×)", fontsize=6, va="center", ha="left", color="#b2182b")
    labels = [d for d in dsets] + ["GSE173193"]
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=32, ha="right", fontsize=6)
    ax.set_ylabel("Attribution rank gain vs\nexpression (×)", fontsize=7)
    ax.set_xlim(-0.6, n_total + 0.75)
    if reds is not None and greys is not None:
        ax.legend([reds, greys], ["cohort-mean ≥ 1.5× (rescued)", "no consistent rescue"],
                  frameon=False, loc="upper left", fontsize=6)
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(FIG / "fig6_tf_rescue.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("Fig6 r4")


if __name__ == "__main__":
    fig1(); fig2(); fig3(); fig6()
    print("R4-DONE")
