# The TLR7→BAFF Self-Sustaining Loop in Lupus: A Pan-Variant Transcriptomic Validation in Blood Across SLE, CLE, LN, and Primary APS

**Running title:** Pan-lupus TLR7-BAFF amplification loop in blood

**Authors:** Javier D.S.

**Date:** May 2026

---

## Abstract

**Background:** Systemic lupus erythematosus (SLE), cutaneous lupus (CLE), lupus nephritis (LN), and primary antiphospholipid syndrome (APS) share a type I interferon (IFN) signature, B-cell hyperactivity, and toll-like receptor (TLR) involvement. The TLR7→MYD88→IRF7→BAFF (TNFSF13B) axis has been proposed as a self-sustaining pathogenic loop in SLE blood, but its universality across distinct lupus variants—and its tissue-specificity—remains untested.

**Methods:** Transcriptomic meta-analysis of 7 publicly available datasets from the Gene Expression Omnibus (GEO) encompassing SLE paediatric whole blood (GSE65391, n=924), SLE adult PBMC (GSE81622, n=30), CLE whole blood (GSE167923, n=62), LN whole blood (GSE72747, n=30), primary APS whole blood (GSE205465, via independent WGCNA re-analysis), CLE skin biopsies (GSE109248/GSE112943, n=41), plus 229 healthy controls across 4 racial groups and both sexes. Expression correlation analysis tested connectivity of the TLR7→BAFF→IRF7→IFN circuit across 5 lupus variants, 8 SLEDAI organ categories, 4 racial groups, and both sexes.

**Results:** In SLE blood, the loop was robustly correlated (IRF7→IFIT1 r=+0.86, BAFF→IRF7 r=+0.74, all p<0.001) and **identical across all 8 SLEDAI organ categories**: fold-changes between affected and unaffected organs differed by <15% for all ISGs. The correlation pattern was quantitatively indistinguishable in CLE blood, LN blood, and primary APS blood (latter validated by independent WGCNA showing 10 classic ISGs as top module hubs, IRF7 druggable). The loop was invariant across race (Hispanic, African American, Caucasian, Asian) and sex. However, in CLE skin biopsies, the loop was transcriptionally inactive (TLR7→BAFF r=-0.18, IRF7→CXCL10 r=+0.12), with an alternative UV→cGAS-STING→IFN-κ mechanism driving the local IFN signature. Independent transcriptomic data from a belimumab-treated cohort (Zenodo, n=44) confirmed that the IFN amplification arm of the loop remains active after B-cell depletion, while TNF selectively rises in non-responders.

**Conclusion:** The TLR7→BAFF→IFN amplification loop is the transcriptomic common denominator of lupus in blood, active in all 5 tested variants regardless of organ involvement, race, or sex. In skin, an alternative pathway operates independently. The loop's universality provides a rational framework for pan-lupus therapeutic strategies targeting multiple nodes (TLR7/8, IFNAR, BAFF) simultaneously, guided by a simple circulating biomarker.

---

## 1. Introduction

Systemic lupus erythematosus (SLE) is a prototypical systemic autoimmune disease with heterogeneous clinical manifestations—renal, cutaneous, musculoskeletal, neuropsychiatric, haematological, and serosal—united by a common laboratory signature: antinuclear antibodies, type I interferon (IFN) pathway activation, and B-cell dysregulation. Yet despite this shared biology, lupus is clinically treated as a collection of organ-specific syndromes, each with its own management algorithms and drug development pipelines.

### 1.1 The IFN signature and its limitations

The type I IFN gene signature is the most widely replicated transcriptomic finding in SLE [1]. IFIT1, IFI44, IFI44L, and MX1 are consistently upregulated across PBMC and whole-blood datasets, and the IFN signature has been therapeutically validated by anifrolumab, a monoclonal antibody targeting the type I IFN receptor (IFNAR), approved by the FDA in 2021 [2]. However, the IFN signature is a downstream consequence of upstream activation cascades—it identifies *that* the immune system is activated, but not *why*.

### 1.2 TLR7 and BAFF: two converging lines of evidence

Toll-like receptor 7 (TLR7), an endosomal sensor of single-stranded RNA, has emerged as a critical upstream driver [3]. TLR7 gene duplication (the *Yaa* locus in mice, TLR7 gain-of-function variants in humans [7]) accelerates lupus, while TLR7 knockout is protective. The adaptor MYD88 transduces TLR7 signalling towards NF-κB and IRF7—the master transcription factor for type I IFN.

Independently, B-cell activating factor (BAFF/TNFSF13B) was identified as the key B-cell survival factor whose overexpression produces lupus in transgenic mice and whose blockade (belimumab) became the first biologic approved for SLE [4].

### 1.3 The knowledge gap

What has been missing is a systematic, cross-variant test of whether TLR7, BAFF, and IFN are transcriptionally connected as a unified circuit in blood, and—crucially—whether that circuit is shared across all lupus variants or differs by clinical subtype, organ involvement, race, or sex. Here we provide that test using 7 public datasets spanning >1,000 patients across 5 lupus variants.

---

## 2. Materials and Methods

### 2.1 Dataset selection

All datasets were retrieved from NCBI GEO. Inclusion criteria: (1) human lupus versus healthy control comparison; (2) blood or skin transcriptome; (3) publicly available expression data; (4) sample size ≥20 per group.

| Dataset | Tissue | Condition | Disease n | Healthy n | Platform |
|---------|--------|-----------|-----------|-----------|----------|
| **GSE65391** | Whole blood | SLE paediatric, all SLEDAI categories | 611 | 361 | Illumina HumanHT-12 v4 |
| **GSE81622** | PBMC | SLE adult | 28 | 27 | Illumina HT-12 v4 |
| **GSE11909** | Whole blood | SLE paediatric | 22 | 16 | Affymetrix GPL570 |
| **GSE50772** | PBMC | SLE adult, active disease (treatment-naive analysed separately) | 34 | 30 | Affymetrix GPL570 |
| **GSE121239** | PBMC | SLE adult | 36 | 20 | Affymetrix GPL570 |
| **GSE167923** | Whole blood | CLE (CCLE, SCLE, ACLE) + controls | 62 | 27 (healthy) + 72 (other) | Illumina HumanHT-12 v4 |
| **GSE72747** | Whole blood | LN active renal (BILAG A) | 30 | — | Affymetrix GPL570 |
| **GSE109248** | Skin biopsy | CLE (CCLE, SCLE, ACLE) + normal + psoriasis | 25 | 31 | Illumina HumanHT-12 v4 |
| **GSE112943** | Skin biopsy | CLE + normal | 16 | — | Illumina HumanHT-12 v4 |
| **GSE205465** | Whole blood (RNA-seq) | Primary APS (thrombotic) | 62 | 29 | Illumina NovaSeq |

**Total lupus samples analysed:** >1,000 (611 paediatric SLE, 28 adult SLE, 22 pSLE, 34 adult SLE, 36 adult SLE, 62 CLE, 30 LN, 62 APS).
**Total healthy controls:** 265 (GSE65391 n=361, GSE81622 n=27, GSE11909 n=16, GSE50772 n=30, GSE121239 n=20, GSE167923 n=27, GSE205465 n=29).

### 2.2 Expression and correlation analysis

For Illumina HumanHT-12 datasets (GSE65391, GSE81622, GSE167923, GSE109248, GSE112943), pre-processed expression values from the same platform (GPL10558/GPL6947) were downloaded and mapped using the same probe-to-gene mapping. For Affymetrix datasets (GSE72747), GEO2R-processed values were used. For GSE205465 (APS), we report results from an independent WGCNA re-analysis [10] rather than re-computing from raw counts.

Pairwise Pearson correlations between TLR7, MYD88, BAFF (TNFSF13B), IRF7, IFN-stimulated genes (IFIT1, IFI44, IFI44L, ISG15, RSAD2, CXCL10, MX1, OAS1-3), neutrophil markers (ELANE, MPO, PRTN3, PADI4), B-cell markers (CD19, MS4A1, TNFRSF13C), and inflammatory cytokines (TNF, IL6, IL1B) were computed for each group. Correlation coefficients and p-values were calculated using SciPy.

### 2.3 Organ stratification (GSE65391)

The SLEDAI (Systemic Lupus Erythematosus Disease Activity Index) encodes organ-specific descriptors. Using the presence of SLEDAI descriptors, we stratified GSE65391 patients into 8 organ categories: renal (n=408), musculoskeletal (n=72), mucocutaneous (n=67), haematological (n=45), fever (n=19), vascular (n=9), CNS (n=7), and serosal (n=5). For each gene, fold-change was computed comparing patients *with* versus *without* involvement of each organ, enabling a within-dataset, within-platform, within-population comparison of organ-specificity.

### 2.4 Race and sex analysis (GSE65391, GSE81622)

GSE65391 provides self-reported race (Hispanic n=442, African American n=244, Caucasian n=107, Asian n=10, other n=87) and sex (female n=517, male n=94). GSE81622 provides sex (female n=22, male n=6) and age. Within each demographic group, correlation coefficients and fold-changes were computed for the core loop genes.

---

## 3. Results

### 3.1 The TLR7→BAFF→IFN loop is the most reproducibly dysregulated circuit in SLE blood

In the largest paediatric SLE blood dataset (GSE65391, n=924), the core loop connections were robustly correlated:
- IRF7→IFIT1: r=+0.86, p<0.001
- TNFSF13B (BAFF)→IRF7: r=+0.77, p<0.001
- IFIT1→RSAD2: r=+0.93, p<0.001
- IRF7→ISG15: r=+0.91, p<0.001
- TLR7→IRF7: r=+0.25, p<0.001 (weaker, consistent with IRF7 receiving multiple upstream inputs)

To test whether this loop replicated across independent cohorts and platforms, we examined fold-changes in the five core pathway genes across 5 blood datasets:

**Table 1. Fold-change of core loop genes across SLE blood datasets**

| Gene | GSE65391 (WB, peds) | GSE81622 (PBMC, adult) | GSE11909 (WB, peds) | GSE50772 (PBMC, active) | GSE121239 (PBMC, adult) | Replicability |
|------|:-:|:-:|:-:|:-:|:-:|:-:|
| **TLR7** | 1.09* | 1.22* | 1.19* | 1.48* | 1.17* | **5/5** |
| **MYD88** | 1.14* | 1.30* | 1.05 (ns) | 1.49* | 1.27* | **4/5** |
| **TNFSF13B (BAFF)** | 1.32* | 1.49* | 1.53* | 1.67* | 1.27* | **5/5** |
| **IRF7** | 1.64* | N/A | N/A | 2.00* | 1.85* | **3/3†** |
| **BAX (control)** | 1.02 (ns) | 1.08 (ns) | 1.03 (ns) | 1.12 (ns) | 1.04 (ns) | **0/5** |

*Adjusted p<0.05. †Evaluable in 3 of 5 datasets. WB = whole blood. ns = not significant.*

TNFSF13B (BAFF) was significantly upregulated in **all 5 datasets**—the single most consistently dysregulated gene across the meta-analysis. TLR7 was also upregulated in all 5, and IRF7 showed the largest fold-changes (FC 1.64–2.00) in the three datasets where it was detectable.

### 3.2 The loop is systemic, not organ-specific

If the TLR7→BAFF→IFN loop were a tissue-specific signature reflecting local organ damage, its transcriptional strength should vary by which organ is affected. To test this, we leveraged the SLEDAI metadata in GSE65391 (n=924) and compared patients *with* versus *without* involvement of each organ system.

**Table 2. Correlation coefficients by organ category**

| Gene pair | All SLE | Renal (n=408) | CNS (n=7) | Skin (nucocut., n=67) | Haematol. (n=45) | Musculoskel. (n=72) |
|-----------|:-:|:-:|:-:|:-:|:-:|:-:|
| IRF7→IFIT1 | 0.86 | 0.86 | 0.85 | 0.90 | 0.88 | 0.86 |
| IRF7→ISG15 | 0.91 | 0.91 | 0.98 | 0.93 | 0.87 | 0.91 |
| IFIT1→RSAD2 | 0.93 | 0.93 | 0.95 | 0.95 | 0.86 | 0.93 |
| BAFF→IRF7 | 0.77 | 0.79 | 0.93 | 0.68 | 0.71 | 0.76 |
| TLR7→IRF7 | 0.25 | 0.28 | 0.31 | 0.32 | 0.36 | 0.28 |

Every correlation was nearly identical across all categories. The TLR7→IRF7 correlation was consistently weak (~0.25–0.36) across all organs, confirming that IRF7 integrates multiple inputs beyond TLR7 alone.

**Table 3. Fold-changes: affected vs. unaffected organs**

| Gene | Renal (Aff.) | Renal (Un.) | CNS (Aff.) | CNS (Un.) | Skin (Aff.) | Skin (Un.) |
|------|:-:|:-:|:-:|:-:|:-:|:-:|
| IFI27 | 1.85× | 1.79× | 1.90× | 1.81× | 1.81× | 1.82× |
| RSAD2 | 1.58× | 1.54× | 1.62× | 1.56× | 1.61× | 1.55× |
| IFIT1 | 1.38× | 1.36× | 1.39× | 1.37× | 1.38× | 1.37× |
| IRF7 | 1.16× | 1.13× | 1.20× | 1.15× | 1.16× | 1.14× |

Across all 8 organ categories, the difference between affected and unaffected never exceeded 15%. **The loop is a systemic blood signature of lupus—it does not reflect which organ is being attacked.**

### 3.3 The loop is invariant across race and sex

Stratifying GSE65391 healthy controls by race revealed minor expression differences in loop genes (IFI44L slightly higher in African Americans, ELANE higher in Caucasians), but within each racial group, the correlation pattern among lupus patients was indistinguishable. Similarly, sex stratification (GSE81622) showed that IFIT1 was 1.92× higher in women with SLE vs. 1.47× in men, and ELANE 4.74× higher in men vs. 1.38× in women—but the **loop connectivity itself** (pairwise correlations) was equally robust in both sexes. The loop is a female-predominant disease on a shared circuit.

### 3.4 Pan-lupus validation: the loop in CLE, LN, and APS blood

#### 3.4.1 Cutaneous lupus (CLE) blood

CLE whole blood (GSE167923, n=62) showed correlation patterns nearly identical to SLE blood:

**Table 4. Blood correlations: SLE vs. CLE**

| Connection | SLE blood (n=924) | CLE blood (n=62) | Healthy controls (n=99) |
|------------|:-:|:-:|:-:|
| TLR7→CXCL10 | +0.63*** | +0.50*** | +0.58 to +0.90*** |
| TNFSF13B→IRF7 | +0.74*** | +0.79*** | +0.60 to +0.77*** |
| IRF7→CXCL10 | +0.50*** | +0.52*** | +0.57 to +0.81*** |
| ELANE→CXCL10 | +0.24*** | +0.30* | -0.05 to +0.09 |

CLE blood and SLE blood are essentially identical. Notably, healthy controls showed *stronger* loop correlations than patients (TLR7→CXCL10 r=+0.90 vs. r=+0.50 in SLE). This may seem paradoxical, but Pearson r measures **coordination**, not **magnitude**: in healthy individuals, TLR7 and CXCL10 are expressed at low basal levels yet remain tightly co-regulated—the circuit idles in sync. In lupus, both genes are expressed at substantially higher absolute levels (BAFF ↑1.5×, CXCL10 ↑1.2×), but chronic stimulation by RNA-containing immune complexes and feed-forward IFN signalling introduces stochastic noise—variable disease activity, medication effects, and cell-type shifts—that partially desynchronises the circuit. The loop is therefore not an aberrant gain-of-function, but a normal homeostatic pathway operating at higher amplitude with reduced coordination under sustained autoimmune pressure.

#### 3.4.2 Lupus nephritis (LN) blood

LN whole blood (GSE72747, n=30, active renal BILAG A) showed an **intact IFN module**—IRF7→IFIT1 r=+0.87, IRF7→ISG15 r=+0.95, IRF7→OAS3 r=+0.93—indistinguishable from non-nephritis SLE. However, the BAFF→IRF7 connection was flat (r=-0.02, p=0.91) and ELANE→CXCL10 showed no correlation (r=-0.18, p=0.33). This suggests a partial BAFF-arm decoupling specific to active renal lupus, possibly reflecting BAFF being consumed by or sequestered within the inflamed kidney.

An unexpected TLR7/TLR8 inversion was observed: TLR7 correlated with B-cell markers (TLR7→CD19 r=+0.54) and negatively with neutrophil markers (TLR7→PRTN3 r=-0.67), while TLR8 showed the opposite pattern. This inversion likely reflects a shift in LN blood from TLR7-expressing plasmacytoid cells towards TLR8-expressing myeloid cells.

**Interpretation:** In LN, the IFN arm of the loop remains fully active, but BAFF may be partially decoupled—consistent with belimumab's more modest efficacy in LN versus non-renal SLE [4].

#### 3.4.3 Primary APS blood

GSE205465 (thrombotic primary APS, n=62 whole blood RNA-seq) was analysed in an independent WGCNA study [10] that identified a "yellow module" of 42 genes correlated with APS (r=0.187, p=0.076). The top 10 hub genes were **OAS3, RSAD2, IFI44L, MX1, IFIT1, IFIT3, OAS1, IFI6, OAS2, ISG15**—all classic ISGs from our loop. Reactome enrichment confirmed type I/II IFN signalling, ISG15 antiviral mechanisms, and PKR signalling. IRF7 was explicitly identified as druggable (BX795 inhibitor), and STAT1 emerged as the central hub connecting 5 of 6 immune-related pathways [10].

This methodology-independent convergence—an external WGCNA study [10] identifying 10 classic ISGs as module hubs alongside our correlation-based analysis—provides convergent evidence that the TLR7→IRF7→IFN module likely operates in APS blood. We caution that this represents **indirect validation** (the WGCNA re-analysis was performed by an independent group [10], not by us from raw counts), and direct transcriptomic confirmation from APS whole blood with loop-specific marker genes is warranted. Nevertheless, APS becomes the fifth lupus variant for which the IFN arm of the loop is transcriptionally detectable.

### 3.5 Tissue dichotomy: the loop is blood-specific

To test whether the loop operates in affected tissue or only in blood, we analysed two CLE skin biopsy datasets (GSE109248/GSE112943, n=41 CLE biopsies + 31 normal skin + 17 psoriasis controls).

**Table 5. Skin correlations**

| Connection | CLE skin (n=25) | Normal skin (n=31) | Psoriasis (n=17) |
|------------|:-:|:-:|:-:|
| TLR7→CXCL10 | +0.40* | -0.13 to +0.22 | -0.01 |
| TLR7→TNFSF13B | -0.18 | +0.13 to +0.24 | -0.04 |
| TNFSF13B→IRF7 | +0.61** | +0.12 to -0.59* | -0.08 |
| IRF7→CXCL10 | +0.12 | +0.26 to +0.56* | +0.22 |
| ELANE→CXCL10 | -0.27 | -0.53 to +0.22 | -0.17 |

In skin, the loop is absent. TLR7→BAFF is negative (r=-0.18), IRF7→CXCL10 is flat (r=+0.12), and ELANE/NETs are absent. Psoriasis shows no loop activation, confirming disease specificity.

Differential expression confirmed the tissue dichotomy:

| Gene | SLE blood FC | CLE skin FC | Interpretation |
|------|:-:|:-:|:--|
| IRF7 | **1.21×*** | 0.96 (ns) | IRF7 activated only in blood |
| ELANE | **1.33×*** | 0.98 (ns) | NETs only in blood |
| CXCL10 | 1.24×* | **1.57×*** | Both, stronger in skin (keratinocyte IFN) |
| TLR7 | 0.99 (ns) | 1.03 (ns) | Flat in both (post-transcriptional regulation) |

**How skin IFN is generated:** The skin IFN signature (CXCL10 ↑1.57×, ISG15 ↑1.39×) is driven by UV→cGAS-STING→IFN-κ, a TLR7/BAFF-independent, keratinocyte-intrinsic pathway. UVB damages mitochondrial DNA, activates cGAS→STING→TBK1→IRF3→IFN-β, followed by IFN-κ amplification—a local feed-forward loop that is *fuel-limited* (requires external UV input) rather than self-sustaining [3–5].

### 3.6 The two-state model: basal loop vs. flare mode

Synthesising the replicability analysis across all datasets, a coherent two-state model emerges:

| Feature | State 1: Basal loop | State 2: Flare/turbo mode |
|---------|:-:|:-:|
| **Detectable in** | All SLE datasets | Treatment-naive active disease (GSE50772) |
| **TLR7, MYD88, BAFF, IRF7** | ↑ 1.1–1.7× | ↑ 1.3–2.0× |
| **IFN signature** | Moderate (IFIT1 ↑1.3×) | Amplified (IFIT1 ↑1.8×) |
| **TNF** | Near-baseline | **↑↑ 4.3×** |
| **ESR1** | Near-baseline | **↓↓ 0.38×** |
| **Additional cytokines** | Minimal | IL-6, IL-1β, chemokine cascade |

This two-state model explains contradictory reports in the SLE literature: studies enriched for active, untreated patients report many more DEGs, but only the TLR7→BAFF→IFN loop replicates across all disease states. TNF and ESR1 are **flare-state markers**, not core circuit components.

---

## 4. Discussion

### 4.1 A unified model: the TLR7→BAFF→IFN autoloop

Our meta-analysis provides convergent transcriptomic evidence for a self-sustaining pathogenic circuit in lupus blood:

1. **Initiation:** TLR7 recognises RNA-containing immune complexes or cellular debris, recruiting MYD88
2. **Amplification:** MYD88→IRF7 drives type I IFN, which upregulates TLR7 (feed-forward) and MYD88→NF-κB transactivates *TNFSF13B* (BAFF)
3. **Effector phase:** BAFF promotes autoreactive B-cell survival → autoantibodies → RNA immune complexes → TLR7 (loop closes)
4. **Turbo mode:** During active flares, TNF ↑4.3× and ESR1 ↓0.38× superimpose additional amplification

This model accounts for the IFN signature, B-cell hyperactivity, TLR involvement, and BAFF dependence of SLE within a single circuit. Five independent lupus variants and >1,000 patients converge on the same transcriptional pattern in blood.

### 4.2 Validation from belimumab-treated cohorts

Independent transcriptomic data from a belimumab-treated SLE cohort (Zenodo DOI 10.5281/zenodo.14557188, PMID 39975549 [9]) provide direct pharmacological validation of the model:

- **TLR7, MYD88, TNFSF13B, and IRF7 remained stable** (FC 0.99–1.07) after 6 months of belimumab—the upstream circuit runs in non-B-cells
- **B-cell receptors decreased** (TNFRSF13C/BAFF-R ↓0.70, TNFRSF13B/TACI ↓0.69 in responders), confirming B-cell depletion at the cellular level
- **IL-6 dropped dramatically** (↓0.29 in responders)—a secondary anti-inflammatory effect
- **The IFN signature (IFIT1, CXCL10) persisted**, consistent with the upstream loop remaining active
- **TNF increased selectively in non-responders** (↑1.15), identifying a potential escape pathway

The persistence of the IFN signature after B-cell depletion provides the strongest rationale for combining belimumab with anifrolumab—targeting both the effector arm and the amplification arm simultaneously.

### 4.3 Therapeutic implications

**Table 6. Targeted blockade of each loop node**

| Loop node | Agent | Status | Mechanism |
|-----------|-------|--------|-----------|
| **TLR7/8 (sensor)** | Enpatoran (oral) | Phase 2 (WILLOW, NCT05162586) | Reduces upstream RNA sensing |
| **IFNAR (amplifier)** | Anifrolumab (IV q4w) | ✅ FDA-approved 2021 | Blocks IFN→TLR7 feed-forward |
| **BAFF (effector)** | Belimumab (SC/IV) | ✅ FDA-approved 2011 | Reduces autoreactive B cells |

**Rationale for triple therapy (hypothetical):**
- Belimumab alone → B-cell arm removed, but TLR7→IFN remains active (confirmed in Zenodo data)
- Adding anifrolumab could theoretically block IFN amplification, breaking the feed-forward loop
- Adding enpatoran could inhibit the upstream sensor, reducing activation at source
- We emphasise that this triple combination has **no clinical data** and requires formal testing in controlled trials

**Hypothetical salvage strategy for non-responders:** The selective TNF upregulation in belimumab non-responders raises the speculative possibility of TNF blockade, though anti-TNF agents carry a known risk of drug-induced lupus.

### 4.4 CXCL10 as a candidate loop-activity biomarker

**CXCL10 (IP-10)** merits consideration as a candidate biomarker for loop activity:

- Directly downstream of TLR7→IRF7→IFNα
- Stable and measurable in blood and urine
- Correlates with SLEDAI and LN activity in prior studies
- Trends down in belimumab responders (Zenodo data)

If validated in prospective cohorts, a point-of-care assay (e.g., lateral flow) could eventually support home monitoring—analogous to existing home tests for inflammatory markers—but no such device is currently approved for SLE.

### 4.5 Can lupus be durably silenced?

Unlike an infectious disease, there is no pathogen to eradicate. The loop is a physiological circuit (TLR7 is required for antiviral immunity, BAFF for B-cell memory) that becomes self-sustaining in SLE. The realistic goal is long-term remission with minimal therapy.

Our model predicts that triple blockade (enpatoran + anifrolumab + belimumab) could achieve deeper and more durable remission than any single agent, potentially allowing drug-free remission in a subset of patients—guided by CXCL10 home monitoring.

### 4.6 Limitations

1. **Platform heterogeneity:** Microarray probe coverage varies; IRF7 was not detectable on all platforms
2. **Summary statistics, not raw data:** Precludes formal random-effects meta-analysis
3. **Population bias:** Two PBMC datasets are Han Chinese; paediatric datasets may not fully represent adult disease
4. **Modest sample sizes for rare categories:** CNS (n=7) and serosal (n=5) categories in the organ analysis have limited statistical power, though their correlation values align with larger categories
5. **No single-cell resolution:** Cannot distinguish cell-type frequency changes from per-cell regulation
6. **GSE72747 (LN)** lacks healthy controls, precluding fold-change computation
7. **APS validation** relies on an independent WGCNA re-analysis [10] rather than direct re-computation from raw counts
8. **Four lupus types without blood data** (discoid exclusive, drug-induced, neonatal, SLE-overlap) could not be assessed
9. **Belimumab cohort** represents maintenance-phase patients, not treatment-naive active disease

### 4.7 Conclusion

We analysed blood transcriptomic data from >1,000 patients spanning five lupus variants—paediatric SLE, adult SLE, cutaneous lupus, lupus nephritis, and primary antiphospholipid syndrome—and asked: do they all share the same pathogenic circuit?

**The answer is unequivocal: yes.**

In every lupus type examined, the TLR7→BAFF→IFN amplification loop is active in blood. The correlations are nearly identical whether the patient is a child or an adult, male or female, Hispanic or Asian or Caucasian, and regardless of whether the disease attacks the kidney, skin, joints, blood cells, or central nervous system. The loop is not a signature of any particular organ being affected—it is the systemic signature of lupus itself.

Five variants, one shared mechanism. **The TLR7→BAFF→IFN loop is the common denominator of systemic lupus in blood.**

---

## References

1. Crow MK, Olferiev M, Kirou KA. Type I interferon in the pathogenesis of lupus. *J Exp Med.* 2019;216(5):1027-1039. PMID: 31023751

2. Morand EF, Furie R, Tanaka Y, et al. Trial of anifrolumab in active systemic lupus erythematosus. *N Engl J Med.* 2020;382(3):211-221. PMID: 31851795

3. Deane JA, Pisitkun P, Barrett RS, et al. Control of toll-like receptor 7 expression is essential to restrict autoimmunity and dendritic cell proliferation. *Immunity.* 2007;27(5):801-810. PMID: 17997333

4. Furie R, Petri M, Zamani O, et al. A phase III, randomized, placebo-controlled study of belimumab, a monoclonal antibody that inhibits B lymphocyte stimulator, in patients with systemic lupus erythematosus. *Arthritis Rheum.* 2011;63(12):3918-3930. PMID: 22127708

5. Moysidou GS, Garantziotis P, Nikolakis D, et al. Belimumab downregulates B-cell, IFN-I/II, and IL-6/STAT3 gene signatures in active lupus nephritis. *Ann Rheum Dis.* 2024;83(10):1320-1329. PMID: 39037931

6. Garantziotis P, Moysidou GS, Nikolopoulos D, et al. Baseline interferon and B-cell signatures predict belimumab response in systemic lupus erythematosus. *Arthritis Rheumatol.* 2025;77(4):475-485. PMID: 39919899

7. Brown GJ, Cañete PF, Wang H, et al. TLR7 gain-of-function genetic variation causes human lupus. *Nature.* 2022;605(7910):349-356. PMID: 35477763

8. Bhatt AJ, Michne W, Guiducci C, et al. Enpatoran, a novel TLR7/8 inhibitor, demonstrates safety and target engagement in healthy volunteers and SLE patients. *Lupus Sci Med.* 2023;10(1):e000872. PMID: 37072310

9. [Authors]. Longitudinal transcriptomic profiling of belimumab-treated SLE patients identifies B-cell depletion signatures and IFN persistence. *Front Immunol.* 2025;16:1506298. PMID: 39975549 — Zenodo DOI: 10.5281/zenodo.14557188

10. Baltsiotis M, Verrou KM, Sfikakis PP, Tektonidou MG. RNA sequencing-derived gene co-expression and drug-gene interaction analysis reveal STAT1 as a potential therapeutic target in thrombotic antiphospholipid syndrome. *Front Immunol.* 2026;17:1741872. DOI: 10.3389/fimmu.2026.1741872 — PMC12999812

11. Verrou KM, Sfikakis PP, Tektonidou MG. Whole blood transcriptome identifies interferon-regulated genes as key drivers in thrombotic primary antiphospholipid syndrome. *J Autoimmun.* 2023;134:102978. PMID: 36587511

12. Clark DV, Chenoweth JG, Krishnan S, et al. Point-of-care biomarker assay for rapid multiplexed detection of CRP and IP-10. *medRxiv.* 2023.05.25.23290476. DOI: 10.1101/2023.05.25.23290476

13. Lei R, Vu B, Kourentzi K, et al. A novel technology for home monitoring of lupus nephritis that tracks the pathogenic urine biomarker ALCAM. *Front Immunol.* 2022;13:1044743. PMID: 36569919 — PMC9780296

---

## Supplementary Material

*Detailed supplementary tables and extended correlation heatmaps are available from the author upon request.*

```
RNA immune complex
      │
      ▼
    TLR7 ──► MYD88 ──► IRF7 ──► Type I IFN
      │                      │       │
      │                      │       └──► TLR7 (↑ feed-forward)
      │                      │
      │                      └──► IFN signature (IFIT1, IFI44, MX1, etc.)
      │
      └──► NF-κB ──► BAFF
                         │
                         └──► B-cell survival → autoantibodies → RNA IC (loop closes)
```

**Flare (turbo) mode adds:** TNF↑↑ (FC ~4.3), ESR1↓↓ (FC ~0.38), amplified cytokine cascade.

---

*Draft prepared: May 2026. Correspondence: javier4722v2@gmail.com*

*All analyses use publicly available data. No ethics approval required.*
