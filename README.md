# TLR7→BAFF→IFN Loop in SLE — Transcriptomic Validation Across 979 Patients

Analysis code and paper for:

> **The TLR7→BAFF Self-Sustaining Loop in Systemic Lupus Erythematosus: Transcriptomic Validation Across 979 Patients and All Organ Domains**

## Overview

Transcriptomic meta-analysis of 11 publicly available GEO datasets (979 patients, 501 healthy controls) testing whether the TLR7→MYD88→IRF7→BAFF→IFN amplification loop is the common blood-circuit across SLE (all organ domains), LN (SLE renal), and CLE.

**v2 correction:** APS removed — primary antiphospholipid syndrome is a distinct thrombotic autoimmune disorder, not a lupus variant.

## Repository Structure

```
├── paper_lupus.md          # Paper manuscript (Markdown)
├── paper_lupus.docx        # Paper manuscript (DOCX)
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
├── figures/
│   └── generate_figures.py
└── README.md
```

## Datasets Analysed

| Dataset | Type | n |
|---------|------|:--:|
| GSE65391 | SLE paediatric whole blood | 611 |
| GSE11909 | SLE paediatric PBMC | 22 |
| GSE81622 | SLE adult PBMC | 28 |
| GSE50772 | SLE adult PBMC | 34 |
| GSE121239 | SLE adult PBMC | 65 |
| GSE72326 | LN whole blood | 157 |
| GSE167923 | CLE whole blood | 62 |
| GSE72747 | LN blood (correlations only) | 30 |
| GSE109248 | CLE skin biopsies | 20 |
| GSE112943 | CLE skin biopsies | 21 |
| **Total** | | **979 patients + 501 controls** |

## Dependencies

- Python 3.8+
- numpy, pandas, scipy, scikit-learn
- matplotlib, seaborn
- docx (python-docx)
- GEOparse (for GEO data access)

## Usage

```bash
# Download and process individual datasets
python scripts/data/process_gse65391.py
python scripts/data/gse81622/process_gse81622.py

# Run core analyses
python scripts/analysis/correlation_analysis.py
python scripts/analysis/analyze_by_organ.py

# Generate figures
python figures/generate_figures.py
```

## Contact

javier4722v2@gmail.com
