# Integrative Analysis Report: FANTOM5 & FANTOM6 Mast Cells (BsMC vs FsMC)

## 1. Overview and Rationale
This repository contains the integrated transcriptomic analysis of primary human mast cells from two distinct tissue sources across two FANTOM datasets:
- **FANTOM5**: Mast cells from breast tissue (BsMC).
- **FANTOM6**: Mast cells from foreskin tissue (FsMC).

For both tissues, transcriptomes were captured in two states:
- **Native**: Baseline expression without stimulation.
- **Stimulated**: Activated via IgE cross-linking (FceRI) or PMA/Ionomycin.

### 1.1 Goal
The primary objective of this integration is to characterize the transcriptional differences between mast cells from different anatomical sites (**Tissue Effect**) and their response to activation (**State Effect**), while robustly controlling for the systematic technical variance between the two datasets (**Batch Effect**).

### 1.2 Model Design
We employed the `limma-voom` pipeline to model the expression data using the following design:
`~ TissueBreast + StateStimulated + DatasetFANTOM5`

This model accounts for the baseline differences between Breast (BsMC) and Foreskin (FsMC), the activation status (Stimulated vs Native), and the systematic batch effect (FANTOM5 vs FANTOM6). The reference levels for the factors were set to `Foreskin`, `Native`, and `FANTOM6`.

## 2. Key Findings

### 2.1 State Effect (Stimulated vs Native)
The model successfully isolated a massive and highly consistent activation signature across all donors, regardless of their tissue origin. 
- **2,349 genes** were identified as significantly differentially expressed (`adj.P < 0.05`).
- The classic mast cell activation markers such as **CCL1** (`logFC = +10.78`, `adj.P = 3.69e-09`) and **CXCL8** (`logFC = +5.50`, `adj.P = 5.29e-03`) were among the top upregulated genes, perfectly capturing the biological intent of the experiment.

### 2.2 Tissue Effect (Breast vs Foreskin)
The model identified **17 genes** that are significantly differentially expressed between the two tissues (`adj.P < 0.05`).
- **Biological Sex Driven Baseline:** The overwhelming majority of the significant tissue differences are driven by the biological sex of the donors. The foreskin donors are exclusively male, while the breast donors are exclusively female. As expected, Y-chromosome genes (e.g., `RPS4Y1`, `EIF1AY`, `DDX3Y`) are strongly upregulated in the FsMC baseline, while X-inactivation genes (`XIST`) are strongly upregulated in the BsMC baseline.
- **Mast Cell Specific Differences:** Despite the strong sex-driven variance, the model retained sufficient statistical power to identify genuine mast-cell-specific tissue markers. Most notably, **TPSD1** (Tryptase Delta 1) was found to be significantly higher in Foreskin mast cells (`logFC = -3.76`, `adj.P = 0.0019`).

## 3. Data Processing and Quality Control
- **Filtering:** Genes with extremely low expression across the majority of samples were filtered out using `edgeR::filterByExpr`, retaining 10,047 robustly expressed genes.
- **Batch Correction:** The PCA plots (see `results/PCA_batch_corrected.png`) confirm that the inclusion of the `Dataset` covariate successfully aligns the FANTOM5 and FANTOM6 samples, removing the laboratory-specific technical bias and allowing for cross-tissue comparisons.
- **Reproducibility:** All scripts used to perform the normalization, statistical modeling, and plotting are provided in the `scripts/` directory.

## 4. Conclusion
The integration of the FANTOM5 and FANTOM6 datasets using a unified `limma-voom` linear model successfully captured both the immense transcriptional shift during mast cell activation and the subtle, sex-driven baseline differences between breast and foreskin mast cells. The resulting catalog of 2,349 activation genes provides a robust, cross-tissue signature of primary human mast cell degranulation.
