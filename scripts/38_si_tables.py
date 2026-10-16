# -*- coding: utf-8 -*-
"""38 SI Tables 导出：results/ 各 CSV → SI_Tables.xlsx（多 sheet）。"""
from pathlib import Path

import pandas as pd

R = Path(r"D:/zcode/z0927/pe_signature/results")
OUT = R / "SI_Tables.xlsx"

# S1 队列元数据（实测）
S1 = pd.DataFrame([
    ["GSE173193", "placenta", "PE 2 / control 2 / GDM 2 / elderly 2", 52358, "contrast (mixed)"],
    ["GSE298119", "placenta", "term PE 2 / control 1", 42248, "contrast (LOPE)"],
    ["GSE298602", "placenta+membranes", "severe-PE 3+ / control 3+ / gHTN 3+", 58923, "contrast (severe)"],
    ["GSE290578", "placenta MFI+deep", "preterm PE 8 / term control 6*", 135209, "contrast (EOPE proxy)"],
    ["GSE282038", "placenta", "EOPE 1 (pooled)", 9305, "profile"],
    ["GSE329173", "placenta", "severe eclampsia 3", 33861, "profile"],
    ["GSE240922", "decidua+villi 7-10w", "reference 30", 132823, "reference"],
], columns=["GSE", "tissue", "groups", "QC_cells", "role"])
# *two author-corrupt samples excluded

S8 = pd.DataFrame([
    ["GSE25906", "Illumina GPL6102 mixed", 23, 37, -1.13, "2.64e-01", "PE slightly low"],
    ["GSE306864", "RNA-seq LOPE-dominant", 13, 45, -0.22, "8.26e-01", "PE slightly low"],
], columns=["cohort", "assay/subtype", "n_PE", "n_CT", "module_t", "module_p", "direction"])

S6 = pd.DataFrame([
    ["Immune", "EOPE", 0.249], ["Immune", "severe", 0.209], ["Immune", "LOPE", -0.306], ["Immune", "mixed", -0.320],
    ["Trophoblast", "EOPE", 0.096], ["Trophoblast", "severe", -0.006], ["Trophoblast", "LOPE", -0.159], ["Trophoblast", "mixed", 0.272],
], columns=["module", "subtype", "point_biserial_r"])

S7 = pd.DataFrame([  # 7 队列均值（GSE173193 来自 phase0，一位小数）；详见 Note
    ["GATA2", 5.3, 7, 7], ["NR2F2", 4.6, 7, 7], ["KLF4", 1.6, 7, 5], ["TEAD4", 1.6, 7, 5],
    ["TFAP2C", 1.6, 7, 6], ["ELF5", 1.4, 7, 3], ["HAND1", 1.3, 7, 2], ["GATA3", 0.8, 7, 1],
    ["DLX3", 0.8, 1, 0], ["MSX2", 0.4, 1, 0],
], columns=["TF", "mean_gain_x", "n_cohorts_measurable", "cohorts_with_gain_ge_1.5x"])
S7_note = ("mean_gain_x = mean rescue factor across measurable cohorts (GSE173193 gains from the phase-0 "
           "full-cohort run, one-decimal reporting). Rescue factor = (expression rank + 1) / (attribution rank + 1). "
           "DLX3 and MSX2 were evaluable in GSE173193 only.")
S7.loc[len(S7)] = [None, None, None, None]
S7.loc[len(S7)] = ["Note", S7_note, None, None]

sheets = {
    "S1_cohorts": S1,
    "S3_headtohead": pd.read_csv(R / "headtohead" / "overlap_summary.csv"),
    "S5_immune_enrichr": pd.read_csv(R / "enrichr" / "Immune_KEGG_2021_Human.csv"),
    "S5b_immune_reactome": pd.read_csv(R / "enrichr" / "Immune_Reactome_2022.csv"),
    "S5c_troph_GO": pd.read_csv(R / "enrichr" / "Trophoblast_GO_Biological_Process_2023.csv"),
    "S6_subtype": S6,
    "S7_TF_rescue": S7,
    "S8_bulk": S8,
    "S9_scf_sens": pd.read_csv(R / "sensitivity" / "scf_immune_sig.csv"),
}
with pd.ExcelWriter(OUT, engine="openpyxl") as w:
    for name, df in sheets.items():
        df.to_excel(w, sheet_name=name[:31], index=False)
print("SI_Tables.xlsx saved:", OUT)
print("sheets:", list(sheets))
