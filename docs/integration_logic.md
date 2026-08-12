# Integration Logic: FANTOM5 + FANTOM6 Mast Cell Data

## Overview

This document details the complete data integration logic and statistical modeling strategy for combining CAGE-seq mast cell gene expression datasets from two independent consortia projects (**FANTOM5** and **FANTOM6**) into a unified, batch-corrected analysis pipeline.

The primary objective is to evaluate tissue-specific transcriptomic differences between **Breast-derived Mast Cells (BsMC)** and **Foreskin-derived Mast Cells (FsMC)** while accounting for activation state (Native vs. IgE-stimulated) and platform batch effects.

---

## 1. Input Data & Sample Cohort (N=21)

### 1.1 FANTOM6 Dataset (Lenti-CAGE, 2022)
- **Source Files:**
  - `data/20220504/BR1256789_Native_Stimulated_gene_counts.txt` (Raw CAGE counts)
  - `data/20220504/BR1256789_Native_Stimulated_sampleinfo.txt` (Sample metadata)
  - `data/20220504/F6_CAT.gene.info_shorted.tsv` (CAT gene ID to HGNC symbol mapping)
- **Selection:** Restricted strictly to **Negative Control ASO (`aso_NC`)** samples (excluding targeted ASO knockdowns `aso_M1` and `aso_M2`).
- **Technical Replicate Handling:** FANTOM6 CAGE libraries contain 2–3 technical replicates per biological sample (e.g. `BR1_aso_NC`, `BR1_aso_NC.1`, `BR1_aso_NC.2`). Raw counts were **summed across technical replicates** for each biological sample prior to low-count filtering and TMM normalization.
- **Samples (14 total):** 7 donors (BR1, BR2, BR5, BR6, BR7, BR8, BR9).
  - **Breast (female):** BR5, BR8 (both Native & Stimulated = 4 samples)
  - **Foreskin (male):** BR1, BR2 (Native-only = 2 samples); BR6, BR7, BR9 (both Native & Stimulated = 6 samples)

### 1.2 FANTOM5 Dataset (HeliScope CAGE, hg38-remapped)
- **Source File:** `data/fantom5_hg38/fantom5_hg38_mast_cell_gene_counts.tsv`
- **Samples (7 total):** All derived from **Breast tissue** across 6 donors (F5_D1, F5_D2, F5_D3, F5_D4, F5_D5, F5_D8).
  - **Native state:** F5_D1, F5_D2, F5_D3, F5_D4, F5_D5 (expanded), F5_D8 (expanded)
  - **Stimulated state:** F5_D1 (stimulated), F5_D5 (expanded and stimulated), F5_D8 (expanded and stimulated)

### 1.3 Unified Cohort Breakdown

| Dataset | Donors | Tissue | Native Samples | Stimulated Samples | Total Samples |
|---------|--------|--------|----------------|--------------------|---------------|
| **FANTOM6** | BR5, BR8 | Breast | 2 | 2 | 4 |
| **FANTOM6** | BR1, BR2, BR6, BR7, BR9 | Foreskin | 5 | 3 | 8 |
| **FANTOM5** | F5_D1, D2, D3, D4, D5, D8 | Breast | 6 | 3 | 9 |
| **Total** | **13 unique donors** | **Breast (13), Foreskin (8)** | **13** | **8** | **21 biological samples** |

---

## 2. Integration Pipeline Workflow

```
FANTOM6 Raw Counts (CAT IDs)          FANTOM5 Raw Counts (HGNC)
       │                                     │
 [Filter: aso_NC only]                [All 9 MC samples]
       │                                     │
 [Map CAT → HGNC symbols]             [Parse Donor & State]
       │                                     │
 [Sum Technical Replicates]           [Map expanded → Native]
       │                                     │
 14 F6 biological samples              9 F5 biological samples
       │                                     │
       └────────────── INNER JOIN ───────────┘
                         │
                  10,047 shared genes
                  21 biological samples
                         │
             [edgeR::filterByExpr()]
                         │
             [edgeR::calcNormFactors(TMM)]
                         │
             [limma::voom()]
                         │
     [duplicateCorrelation(block = Donor)]
                         │
    [lmFit: ~ Tissue + State + Dataset]
```

---

## 3. Statistical Modeling Strategy

### 3.1 Model Formula
To estimate the tissue-specific expression differences while controlling for activation state and platform batch effects, we fit the following linear model:

$$\text{Expression} \sim \text{Tissue} + \text{State} + \text{Dataset}$$

- **`Tissue`**: Factor with levels `c("Foreskin", "Breast")`. The estimated coefficient `TissueBreast` represents the contrast **Breast-derived Mast Cells vs. Foreskin-derived Mast Cells (BsMC vs. FsMC)**.
  - Positive $\text{log}_2\text{FC} = \text{higher expression in Breast (BsMC)}$
  - Negative $\text{log}_2\text{FC} = \text{higher expression in Foreskin (FsMC)}$
- **`State`**: Factor with levels `c("Native", "Stimulated")`. Adjusts for IgE-mediated activation response.
- **`Dataset`**: Factor with levels `c("FANTOM6", "FANTOM5")`. Adjusts for platform batch effects between HeliScope CAGE and Lenti-CAGE.

### 3.2 Intra-Donor Correlation Handling (`duplicateCorrelation`)
The dataset exhibits a **mixed design**: 8 donors are paired (providing both Native and Stimulated samples), while 5 donors are unpaired (Native only). 

To maximize statistical power without violating independence assumptions:
1. `limma::duplicateCorrelation(v, design, block = meta$Donor)` computes the consensus intra-donor correlation ($\rho \approx 0.433$).
2. `limma::lmFit()` incorporates this correlation structure into generalized least squares (GLS) estimation.

This approach recovers the full power of a paired analysis for the activation effect (3,456 DE genes) while utilizing all N=21 samples to estimate the baseline tissue effect.

---

## 4. Key Results Summary

| Contrast | Method | Sig. Genes ($p_{\text{adj}} < 0.05$) | Key Markers |
|----------|--------|--------------------------------------|-------------|
| **Tissue (BsMC vs FsMC)** | `limma-voom` + `dupCor` | **11** | `TPSD1` ($\delta$-tryptase), `RPS4Y1`, `EIF1AY`, `DDX3Y`, `KDM5D`, `UTY`, `XIST` |
| **Activation (Stimulated vs Native)** | `limma-voom` + `dupCor` | **3,456** | `CCL1`, `CXCL8`, `TNF`, `FOS`, `JUN`, `IL6`, `IL3` |
| **Interaction (Tissue × State)** | `limma-voom` + `dupCor` | **0** | Activation response is perfectly conserved across tissue origins |

---

## 5. Publication Concordance Plots

To validate the integration results against the independent single-dataset analysis from the manuscript, we performed cross-pipeline concordance analyses against the exported Excel tables (`Supplement Tables.xlsx`).

### 5.1 Tissue Effect Concordance (Breast vs. Foreskin)
![Tissue Concordance](../results/Concordance_Tissue_Excel_vs_Limma.png)

- **Comparison:** Manuscript DREAM (`Supplement Tables.xlsx`, Sheet `S3a DGE_TissueBreast_gene_Dream`) vs. Integration `limma-voom` (`DE_BsMC_vs_FsMC.tsv`).
- **Concordance Metric:** Pearson correlation **$r = 0.953$** across 9,869 matched genes.
- **Median Absolute Difference:** $\text{Median } |\Delta \text{log}_2\text{FC}| = 0.23$ across the 23 key genes.

### 5.2 State Effect Concordance (Stimulated vs. Native)
![State Concordance](../results/Concordance_State_Excel_vs_Limma.png)

- **Comparison:** Manuscript DREAM (`Supplement Tables.xlsx`, Sheet `S2a_Dream_DEG_Sti.vs.base.`) vs. Integration `limma-voom` (`DE_Stimulated_vs_Native.tsv`).
- **Concordance Metric:** Pearson correlation **$r = 0.892$** across 9,131 matched genes.

---

## 6. Batch Correction for Visualization

For PCA and heatmap visualizations, the `Dataset` platform effect is removed using `limma::removeBatchEffect()`:

```r
expr_corrected <- removeBatchEffect(v$E, batch=meta$Dataset,
                                     design=model.matrix(~ Tissue + State, data=meta))
```

*Note: Batch correction via `removeBatchEffect` is used exclusively for plotting. The linear model operates directly on unadjusted logCPM values while treating `Dataset` as a fixed covariate.*
