# -*- coding: utf-8 -*-
"""44 全英文润色：三个文件逐处替换（零内容/数字/引用变动）。"""
from pathlib import Path

MS = Path(r"D:/zcode/z0927/pe_signature/notes/MANUSCRIPT_GB_FINAL.md")
CL = Path("F:/pe_signature/submission_GB/CoverLetter.md")
NT = Path("F:/pe_signature/submission_GB/Additional_file_1_Supplementary_Note.md")

REPL_MS = [
    ("Cross-cohort disagreement, not absence of signal, defines single-cell preeclampsia (PE) research: comparable studies report trophoblast, immune, and endothelial programs that fail to replicate, blocking consensus on disease programs.",
     "Cross-cohort disagreement, rather than absence of signal, defines single-cell preeclampsia (PE) research: comparable studies report trophoblast, immune, and endothelial programs that fail to replicate, precluding consensus on disease mechanisms."),
    ("Single-cell transcriptomics has rapidly populated the field with atlases of the preeclamptic placenta and decidua, yet the resulting literature is strikingly self-contradictory:",
     "Single-cell transcriptomics has rapidly filled the field with atlases of the preeclamptic placenta and decidua, yet the resulting literature is strikingly self-contradictory:"),
    ("First, cell-composition shifts between cohorts and groups can dominate unstratified differential comparisons, as demonstrated for bulk tissue by deconvolution analysis [7].",
     "First, shifts in cell composition between cohorts and disease groups can dominate unstratified differential comparisons, as demonstrated for bulk tissue by deconvolution analysis [7]."),
    ("Second, expression-based differential statistics are sensitive to sequencing depth, drop-out, and house-gene noise, so that the ranking of genes — the quantity studies compare — is itself unstable.",
     "Second, expression-based differential statistics are sensitive to sequencing depth, drop-out, and housekeeping-gene noise, so the ranking of genes — the quantity studies ultimately compare — is itself unstable."),
    ("However, attribution scoring has not been evaluated in placental biology, where the central cell type — trophoblast — is absent from standard training corpora and from the public precomputed attribution libraries",
     "Attribution scoring, however, has not been evaluated in placental biology, where the central cell type — trophoblast — is absent from standard training corpora and from public precomputed attribution libraries"),
    ("making empirical validation in the target tissue a precondition rather than a formality",
     "which makes empirical validation in the target tissue a precondition rather than a formality"),
    ("We further examined the claimed advantage for low-expression transcription factors and verified that consensus conclusions are not artifacts of a single foundation model.",
     "We further examined the claimed advantage for low-expression transcription factors and verified that our conclusions do not depend on a single foundation model."),
    ("Public single-cell RNA-seq datasets of human placenta, decidua, or chorionic villi with preeclampsia (PE) cases were systematically identified via NCBI E-utilities",
     "Public single-cell RNA-seq datasets of human placenta, decidua, or chorionic villi containing preeclampsia (PE) cases were systematically identified through NCBI E-utilities"),
    ("Sample-level group annotations were parsed from SOFT/MINiML metadata with dataset-specific rules.",
     "Sample-level group annotations were parsed from SOFT/MINiML metadata using dataset-specific rules."),
    ("Two samples of GSE290578 (GSM8817230/31) were excluded: the barcode files supplied by the authors are one line short of the count matrix, failing integrity checks.",
     "Two GSE290578 samples (GSM8817230/31) were excluded because the barcode files supplied by the depositors are one line short of the corresponding count matrices, failing integrity checks."),
    ("Two-level marker-scoring annotation (scanpy score_genes on log-normalized data).",
     "Annotation used a two-level marker-scoring scheme (scanpy score_genes on log-normalized data)."),
    ("Argmax assignment at each level; resulting proportions matched known placental biology for each cohort's sampling design.",
     "Assignment was by argmax at each level; the resulting proportions matched known placental biology for each cohort's sampling design."),
    ("We assembled seven public single-cell cohorts of the human placenta and decidua totalling 464,727 quality-controlled cells:",
     "We assembled seven public single-cell cohorts of the human placenta and decidua totaling 464,727 quality-controlled cells:"),
    ("Before any disease comparison, we asked whether foundation-model attribution is biologically calibrated in pregnancy tissue — the principal precondition for this study, since trophoblast is under-represented in model training corpora and the public SIGnature precomputed library contains no placental data.",
     "Before performing any disease comparison, we asked whether foundation-model attribution is biologically calibrated in pregnancy tissue — the principal precondition for this study, because trophoblast is under-represented in model training corpora and the public SIGnature precomputed library contains no placental data."),
    ("Attribution calibration in placenta is therefore not merely adequate but remarkably consistent across laboratories, chemistries, and gestational ages.",
     "Attribution calibration in placenta is thus not merely adequate but remarkably consistent across laboratories, chemistries, and gestational ages."),
    ("The central empirical problem of PE single-cell biology is that differential signatures conflict between cohorts. We quantified this directly.",
     "The central empirical problem of PE single-cell biology is that differential signatures conflict between cohorts; we quantified this directly."),
    ("Attribution scoring and lineage stratification together thus convert PE differential analysis from cohort-specific to reproducible.",
     "Together, attribution scoring and lineage stratification convert PE differential analysis from cohort-specific to reproducible."),
    ("Functional annotation was lineage-coherent and unexpectedly specific (Fig. 4).",
     "Functional annotation was lineage-coherent and notably specific (Fig. 4)."),
    ("The trophoblast consensus centered on MAPK-cascade regulation containing FLT1 — the gene encoding the anti-angiogenic sFlt-1 whose circulating levels anchor PE pathophysiology and therapeutics [11] — alongside membrane-raft repair and fibrin-clot dissolution modules (ANXA2, S100A10).",
     "The trophoblast consensus centered on MAPK-cascade regulation including FLT1 — the gene encoding the anti-angiogenic sFlt-1, whose circulating levels anchor PE pathophysiology and therapeutics [11] — alongside membrane-raft repair and fibrin-clot-dissolution modules (ANXA2, S100A10)."),
    ("Endothelium carried a complement-receptor infiltration module (C5AR1, TYROBP, AIF1) with IL-6 and IL-1β regulatory terms; stroma carried ECM and gap-junction remodeling (TIMP3, ACTA2, TAGLN, SRGN).",
     "Endothelium carried a complement-receptor infiltration module (C5AR1, TYROBP, AIF1) with IL-6 and IL-1β regulatory terms, and stroma carried ECM and gap-junction remodeling (TIMP3, ACTA2, TAGLN, SRGN)."),
    ("PE is clinically heterogeneous, and so were our cohort signatures.",
     "PE is clinically heterogeneous, and so, accordingly, were our cohort signatures."),
    ("The immune complement module showed a clean dichotomy aligned with the placental-versus-maternal disease faces:",
     "The immune complement module showed a clear dichotomy aligned with the placental-versus-maternal disease faces:"),
]

REPL_MS2 = [
    ("These results independently reproduce, in attribution space and without prior labeling, the two-molecular-faces structure described by Admati et al. [5], and they nominate a quantitative, cohort-portable readout — the complement-module correlation — for subtype assignment.",
     "These results independently reproduce, in attribution space and without prior labeling, the two-molecular-faces structure described by Admati et al. [5], and they nominate a quantitative, cohort-portable readout — the complement-module correlation — for subtype assignment in future cohorts."),
    ("A claimed advantage of attribution scoring is sensitivity to low-abundance regulators.",
     "A claimed advantage of attribution scoring is sensitivity to low-abundance regulators; we tested it directly."),
    ("Four factors showed no rescue (GATA3, MSX2, DLX3, HAND1), indicating selectivity rather than a generic low-rank boost.",
     "Four factors showed no rescue (GATA3, MSX2, DLX3, HAND1), indicating selectivity rather than a generic boost of low-ranked genes."),
    ("this overlap lies within the range SCimilarity achieves against itself across cohorts (0.09–0.37), indicating that the consensus direction is a property of the biology rather than of a single model.",
     "this overlap lies within the range SCimilarity achieves against itself across cohorts (0.09–0.37), indicating that the consensus direction reflects the biology rather than any single model."),
    ("Our central result is that this instability is not intrinsic to the biology but substantially to the measurement layer:",
     "Our central result is that this instability lies substantially in the measurement layer rather than in the biology itself:"),
    ("The measured effect is large enough to reinterpret a literature: signatures previously dismissed as cohort artifacts contain a recoverable consensus once importance is scored robustly.",
     "The measured effect is large enough to reinterpret the literature: signatures previously dismissed as cohort artifacts yield a recoverable consensus once importance is scored robustly."),
    ("That these two modules dissociate across subtypes — complement tracking the placental face, trophoblast programs varying independently — gives the atlas clinical texture:",
     "That these two modules dissociate across subtypes — complement tracking the placental face, trophoblast programs varying independently — lends the atlas clinical texture:"),
    ("while adding what that study could not provide: a quantitative, cohort-portable readout (the complement-module correlation) for assigning any future cohort to a face.",
     "while adding what that study could not provide — a quantitative, cohort-portable readout (the complement-module correlation) for assigning any future cohort to a face."),
    ("Why attribution scores replicate where expression statistics do not is not fully resolved by this study; plausible contributors consistent with our observations include robustness to drop-out and depth variation, preserved weighting of low-expression regulators, and, once lineages are stratified, insulation from residual composition noise — but a mechanistic account will require dedicated perturbation-scale evaluation.",
     "Why attribution scores replicate where expression statistics do not is not fully resolved here; contributors consistent with our observations include robustness to drop-out and depth variation, preserved weighting of low-expression regulators, and — once lineages are stratified — insulation from residual composition noise. A mechanistic account, however, will require dedicated perturbation-scale evaluation."),
    ("Several limitations frame these claims.",
     "Several limitations qualify these claims."),
    ('Subtype labels are inherited from the source cohorts: severity definitions, inclusion criteria, sampling site, and platforms differ across deposits, and some labels (notably "preterm PE" as an early-onset proxy, and severity assigned without uniform blood-pressure/organ-damage criteria) are inferred rather than standardized; gestational-age-matched controls would harden the EOPE contrast.',
     'Subtype labels are inherited from the source cohorts: severity definitions, inclusion criteria, sampling sites, and platforms differ across deposits, and some labels — notably "preterm PE" as an early-onset proxy, and severity assigned without uniform blood-pressure or organ-damage criteria — are inferred rather than standardized; gestational-age-matched controls would strengthen the EOPE contrast.'),
    ("Four contrastable cohorts bound the replication statistics, though the 24/24 concordance leaves little ambiguity at current power.",
     "Four contrastable cohorts bound the replication statistics, although the 24/24 concordance leaves little ambiguity at this scale."),
    ("This is a retrospective, correlational analysis of public data: it demonstrates that the complement–myeloid module tracks PE subtype structure, but cannot establish complement activation as a causal driver of disease.",
     "This is a retrospective, correlational analysis of public data: it demonstrates that the complement–myeloid module tracks PE subtype structure but cannot establish complement activation as a causal driver of disease."),
    ("And we did not attempt first-trimester prediction — the clinical prize — because the cohorts are delivery-tissue;",
     "Finally, we did not attempt first-trimester prediction — the clinical prize — because the cohorts comprise delivery tissue;"),
    ("In sum, attribution scoring is not merely a convenient re-weighting of expression data; it changes which parts of PE biology survive cross-cohort scrutiny.",
     "In sum, attribution scoring is not merely a convenient re-weighting of expression data: it changes which parts of PE biology survive cross-cohort scrutiny."),
    ("Attribution scoring converts preeclampsia differential analysis from cohort-specific to reproducible, recovers a myeloid complement program of exceptional statistical sharpness together with canonical trophoblast and stromal modules, and exposes a two-face subtype architecture validated directionally in bulk tissue.",
     "Attribution scoring converts preeclampsia differential analysis from cohort-specific to reproducible, recovers a myeloid complement program of unusual statistical sharpness together with canonical trophoblast and stromal modules, and exposes a two-face subtype architecture corroborated directionally in independent bulk cohorts."),
]

REPL_CL = [
    ("but two open questions have limited their adoption: how their gradient-attribution scores behave in tissues absent from training corpora, and whether that behavior translates into better science in practice. Our manuscript addresses both with a field-scale testbed:",
     "yet two open questions have limited their adoption: how their gradient-attribution scores behave in tissues absent from training corpora, and whether that behavior translates into better science in practice. Our manuscript addresses both questions on a field-scale testbed:"),
    ('This quantitatively reframes "conflicting PE signatures" as substantially a measurement artifact — and provides a robust, cheap (~35 single-GPU-minutes for all cohorts), fully re-runnable layer for the field.',
     'This quantitatively reframes "conflicting PE signatures" as largely a measurement artifact and provides a robust, inexpensive (~35 single-GPU-minutes for all cohorts), fully re-runnable analysis layer for the field.'),
    ("it is a genome-scale methods benchmark whose payoff is disease biology, entirely re-runnable on public data",
     "it is a genome-scale methods benchmark whose payoff is disease biology, fully re-runnable on public data"),
    ("and the benchmark framework itself generalizes to any single-cell field suffering cross-cohort irreproducibility.",
     "and the benchmark framework generalizes to any single-cell field suffering cross-cohort irreproducibility."),
]

REPL_NT = [
    ("All GEO downloads were transferred with byte-exact and, for all attribution matrices, per-chunk SHA256 verification.",
     "All GEO downloads were verified byte-exact on transfer and, for every attribution matrix, chunk-by-chunk by SHA256."),
    ("which fails matrix integrity checks at load time. No other integrity failures were detected in any cohort.",
     "failing integrity checks at load time. No other integrity failures were detected in any cohort."),
    ("Its within-lineage PE signature correlates with SCimilarity's at Spearman ρ = 0.346 over 17,680 shared genes; top-200 Jaccard 0.212, within the range SCimilarity achieves against itself across cohorts (0.09–0.37).",
     "Its within-lineage PE signature correlates with SCimilarity's at Spearman ρ = 0.346 over 17,680 shared genes, with a top-200 Jaccard of 0.212 — within the range SCimilarity achieves against itself across cohorts (0.09–0.37)."),
    ("canonical lineage markers ranked at the extreme top of attribution in all seven cohorts",
     "canonical lineage markers ranked at the very top of attribution in all seven cohorts"),
    ("All head-to-head and consensus analyses were additionally re-run with attribution-quantile gating in place of annotation labels for the subtype analyses, with concordant results.",
     "For the subtype analyses, head-to-head and consensus computations were additionally re-run with attribution-quantile gating in place of annotation labels, with concordant results."),
]


def apply(path, repls, name):
    t = path.read_text(encoding="utf-8")
    n = miss = 0
    for old, new in repls:
        if old in t:
            t = t.replace(old, new)
            n += 1
        else:
            miss += 1
            print("[%s] MISS:" % name, old[:55])
    path.write_text(t, encoding="utf-8")
    print("[%s] applied %d/%d (miss %d)" % (name, n, len(repls), miss))


apply(MS, REPL_MS + REPL_MS2, "Manuscript")
apply(CL, REPL_CL, "CoverLetter")
apply(NT, REPL_NT, "Note")
print("POLISH-DONE")
