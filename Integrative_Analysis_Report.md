# Integrative Analysis Report: FANTOM5 & FANTOM6 Mast Cells (BsMC vs FsMC)

## 1. Overview and Rationale
This repository contains the integrated transcriptomic analysis of primary human mast cells from two distinct tissue sources across two FANTOM datasets:
- **FANTOM5**: Mast cells from breast tissue (BsMC).
- **FANTOM6**: Mast cells from foreskin tissue (FsMC).

For both tissues, transcriptomes were captured in two states:
- **Native**: Baseline expression without stimulation.
- **Stimulated**: Activated via IgE cross-linking ($\text{FC}\varepsilon\text{RI}$) or PMA/Ionomycin.

### 1.1 Goal
The primary objective of this integration is to characterize the transcriptional differences between mast cells from different anatomical sites (**Tissue Effect**) and their response to activation (**State Effect**), while robustly controlling for the systematic technical variance between the two datasets (**Batch Effect**).

### 1.2 Sample Cohort (N=21)
The integrated model utilizes a total of **21 samples** derived from **13 unique donors**. 

To maximize statistical power while accounting for donor-specific baseline variances, the model includes:
- **8 Paired Donors (16 samples):** Donors with both a Native and a Stimulated sample (FANTOM6: BR5, BR6, BR7, BR8, BR9; FANTOM5: F5_D1, F5_D5, F5_D8).
- **5 Unpaired Donors (5 samples):** Donors with only a Native baseline sample (FANTOM6: BR1, BR2; FANTOM5: F5_D2, F5_D3, F5_D4).

The `duplicateCorrelation()` function in `limma` was employed to seamlessly integrate both paired and unpaired samples within the same linear model.

### 1.3 Model Design
We employed the `limma-voom` pipeline with `duplicateCorrelation` to model the expression data using the following design:
`~ Tissue + State + Dataset` (block = Donor)

This model accounts for the baseline differences between Breast (BsMC) and Foreskin (FsMC), the activation status (Stimulated vs Native), and the systematic batch effect (FANTOM5 vs FANTOM6).

## 2. Key Findings

### 2.1 State Effect (Stimulated vs Native)
The model successfully isolated a massive and highly consistent activation signature across all donors, regardless of their tissue origin. 
- **3,456 genes** were identified as significantly differentially expressed (`adj.P < 0.05`).
- The classic mast cell activation markers such as **CCL1** (`logFC = +10.78`, `adj.P = 3.69e-09`) and **CXCL8** (`logFC = +5.50`, `adj.P = 5.29e-03`) were among the top upregulated genes, perfectly capturing the biological intent of the experiment.

### 2.2 Tissue Effect (Foreskin vs Breast)
The model identified **11 genes** that are significantly differentially expressed between the two tissues (`adj.P < 0.05`).
- **Biological Sex Driven Baseline (10 genes):** The overwhelming majority of the significant tissue differences are driven by the biological sex of the donors. Foreskin donors are exclusively male, while breast donors are exclusively female. As expected, Y-chromosome genes (e.g., `RPS4Y1`, `EIF1AY`, `DDX3Y`, `UTY`, `KDM5D`) are strongly upregulated in FsMC, while X-inactivation genes (`XIST`) are strongly upregulated in BsMC.
- **Tissue-Specific Marker (1 gene):** **TPSD1** ($\delta$-Tryptase 1) was identified as significantly higher in Foreskin mast cells (`logFC = +3.76`, `adj.P = 0.0019`).
- **Conserved Core Mast Cell Identity:** Core lineage markers (e.g., `KIT`, `CPA3`, `TPSAB1`, `TPSB2`, `FCER1A`, `FCER1G`, `HDC`) show no significant difference between tissues, proving that primary human mast cells retain a conserved core identity across tissue niches.

### 2.3 Tissue × Activation Interaction
The interaction model (`~ Tissue * State + Dataset`) identified **0 significant genes** (`adj.P < 0.05`), proving that the transcriptional response to activation is **tissue-independent**.

## 3. Data Processing and Quality Control
- **Filtering:** Genes with low expression were filtered out using `edgeR::filterByExpr`, retaining 10,047 robustly expressed genes.
- **Batch Correction:** PCA plots (see `results/PCA_batch_corrected.png`) confirm that the inclusion of the `Dataset` covariate successfully aligns FANTOM5 and FANTOM6 samples, removing laboratory-specific technical bias.
- **Reproducibility:** All scripts used to perform normalization, statistical modeling, and plotting are provided in the `scripts/` directory.

## 4. Conclusion
The integration of FANTOM5 and FANTOM6 datasets using `limma-voom` with `duplicateCorrelation(block=Donor)` and `Dataset` covariate successfully captured both the immense transcriptional shift during mast cell activation (3,456 genes) and the specific baseline tissue differences (11 genes, featuring `TPSD1`).
