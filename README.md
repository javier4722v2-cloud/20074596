# TLR7→BAFF→IFN Loop in Lupus — Pan-Variant Transcriptomic Validation

Analysis code and paper draft for:

> **The TLR7→BAFF Self-Sustaining Loop in Lupus: A Pan-Variant Transcriptomic Validation in Blood Across SLE, CLE, LN, and Primary APS**

## Overview

Transcriptomic meta-analysis of 10 publicly available GEO datasets (>1,000 patients, 5 lupus variants) testing whether the TLR7→MYD88→IRF7→BAFF→IFN amplification loop is the common blood-circuit across all forms of lupus.

## Repository Structure

```
├── paper_draft.md          # Paper manuscript
├── paper_draft.docx        # Paper manuscript (DOCX)
├── scripts/
│   ├── data/               # Dataset-specific processing
│   │   ├── process_gse65391.py
│   │   ├── comprehensive_reanalysis.py
│   │   ├── gse81622/
│   │   ├── gse11909/
│   │   ├── gse121239/
│   │   └── belimumab/
│   └── analysis/           # Core analyses
│       ├── correlation_analysis.py
│       ├── healthy_ref.py
│       ├── analyze_by_organ.py
│       ├── analyze_by_sex.py
│       ├── cle_data/       # CLE, LN subtype analyses
│       └── loop_analysis/  # SLE loop modelling
└── figures/
    └── generate_figures.py
```

## Datasets Used

| Dataset | Tissue | Condition |
|---------|--------|-----------|
| GSE65391 | Whole blood | SLE paediatric, all SLEDAI (n=924) |
| GSE81622 | PBMC | SLE adult (n=28) |
| GSE11909 | Whole blood | SLE paediatric (n=22) |
| GSE50772 | PBMC | SLE adult, active (n=34) |
| GSE121239 | PBMC | SLE adult (n=36) |
| GSE167923 | Whole blood | CLE (n=62) |
| GSE72747 | Whole blood | LN active renal (n=30) |
| GSE109248 | Skin biopsy | CLE (n=25) |
| GSE112943 | Skin biopsy | CLE (n=16) |
| GSE205465 | Whole blood | Primary APS (n=62) |

## Dependencies

- Python 3.8+
- numpy, scipy, pandas, matplotlib, seaborn
- GEOparse (optional, for raw data processing)

## Reproducibility

All data is publicly available from NCBI GEO. Scripts process GEO2R-preprocessed tables where available. Contact the author for pre-processed expression matrices.

## License

MIT
