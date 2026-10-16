# pe-attribution

Reproducibility package for:

> **Attribution scores from single-cell foundation models reveal reproducible, subtype-specific preeclampsia programs across cohorts**

Cross-cohort benchmark of foundation-model attribution signatures vs expression-based differential signatures in preeclampsia (PE) single-cell transcriptomics: 7 public cohorts, 464,727 placental/decidual cells, 24/24 cohort-pair replication advantage (mean top-200 Jaccard 0.205 vs 0.056, Wilcoxon p = 1.2e-07), consensus program mining, two-face subtype architecture, and independent bulk validation.

## Repository layout

```
pe-attribution/
├── scripts/          # All analysis code (49 scripts, numbered by pipeline stage)
├── results/          # Figure source data, SI tables, per-analysis outputs (read-only record)
│   ├── figures/      # Final publication figures (Fig1-Fig6)
│   ├── SI_Tables.xlsx
│   ├── headtohead/   # Attribution vs expression head-to-head benchmark (Fig 2)
│   ├── phase3_v2/    # Consensus program mining (Fig 3-4)
│   ├── subtype/      # Subtype analyses (Fig 5)
│   ├── enrichr/      # Functional annotation
│   ├── bulk_validation/  # GSE25906 + GSE306864 directional validation
│   └── sensitivity/  # scFoundation sensitivity analysis
└── docs/
    ├── INTEGRITY_AUDIT_LOG.md    # Complete integrity audit: Phase 0 debug + fix chain
    ├── ENVIRONMENT_AND_DATA.md   # Runtime environment and data acquisition log
    └── DATASET_INVENTORY.md      # Cohort screening inventory
```

## Data availability

All input datasets are public (NCBI GEO):
GSE173193, GSE290578, GSE298602, GSE298119, GSE282038, GSE329173, GSE240922 (single-cell);
GSE25906, GSE306864 (bulk validation).
Attribution model helper files: Zenodo record 17903196; scFoundation weights: HuggingFace `genbio-ai/scFoundation`.

## Environment

- Python 3.12, scanpy 1.12.4, scikit-learn 1.9, SIGnature 1.0.0 (`pip install sc-signature`)
- Attribution compute: single NVIDIA RTX 4090, ~35 GPU-minutes total for all seven cohorts
- Downstream analyses run on CPU
- Reference setup: `scripts/recreate_env.sh`

## Reproduction

1. Download raw GEO archives (see `scripts/10_harvest_gds.sh`, `11_batch_download.sh`).
2. Preprocess per cohort (`01`, `12`, `14`, `15`).
3. Attribution scoring (`02`, `13`; sensitivity `32`).
4. Benchmark & mining (`25` head-to-head, `18` consensus, `27` enrichr, `28` subtype, `35`-`36` bulk).
5. Figures (`26`, `29`, `30`, `31`, `37`, `39`-`41`); SI tables (`38`).

Every run records git commit / dependency versions / seed / GPU model in the results JSON `meta` field. Results files are append-only.

## License

MIT
