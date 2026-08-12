
# Supplementary Materials

## Integrative Transcriptomic Analysis of Breast-derived (BsMC) vs. Foreskin-derived (FsMC) Mast Cells Across FANTOM5 and FANTOM6 CAGE-seq Platforms

---

## Supplementary Tables

### Table S1. Complete Sample Metadata (N=21 Biological Samples)

Complete metadata for all 21 biological samples included in the integrative analysis across FANTOM5 and FANTOM6 datasets. Technical replicates from FANTOM6 (2–3 per biological sample) were summed at the raw count level prior to normalization.

| Sample ID | Donor | Tissue | State | Dataset | Original Tech Reps |
|:---|:---|:---|:---|:---|:---|
| BR1_aso_NC | BR1 | Foreskin | Native | FANTOM6 | 3 |
| BR2_aso_NC | BR2 | Foreskin | Native | FANTOM6 | 3 |
| BR5_aso_NC | BR5 | Breast | Native | FANTOM6 | 3 |
| BR5_AER_aso_NC | BR5 | Breast | Stimulated | FANTOM6 | 3 |
| BR6_aso_NC | BR6 | Foreskin | Native | FANTOM6 | 2 |
| BR6_AER_aso_NC | BR6 | Foreskin | Stimulated | FANTOM6 | 2 |
| BR7_aso_NC | BR7 | Foreskin | Native | FANTOM6 | 2 |
| BR7_AER_aso_NC | BR7 | Foreskin | Stimulated | FANTOM6 | 2 |
| BR8_aso_NC | BR8 | Breast | Native | FANTOM6 | 3 |
| BR8_AER_aso_NC | BR8 | Breast | Stimulated | FANTOM6 | 3 |
| BR9_aso_NC | BR9 | Foreskin | Native | FANTOM6 | 2 |
| BR9_AER_aso_NC | BR9 | Foreskin | Stimulated | FANTOM6 | 2 |
| Mast cell, donor 1 | F5_D1 | Breast | Native | FANTOM5 | 1 |
| Mast cell – stimulated, donor 1 | F5_D1 | Breast | Stimulated | FANTOM5 | 1 |
| Mast cell, donor 2 | F5_D2 | Breast | Native | FANTOM5 | 1 |
| Mast cell, donor 3 | F5_D3 | Breast | Native | FANTOM5 | 1 |
| Mast cell, donor 4 | F5_D4 | Breast | Native | FANTOM5 | 1 |
| Mast cell, expanded, donor 5 | F5_D5 | Breast | Native | FANTOM5 | 1 |
| Mast cell, expanded and stimulated, donor 5 | F5_D5 | Breast | Stimulated | FANTOM5 | 1 |
| Mast cell, expanded, donor 8 | F5_D8 | Breast | Native | FANTOM5 | 1 |
| Mast cell, expanded and stimulated, donor 8 | F5_D8 | Breast | Stimulated | FANTOM5 | 1 |

---

### Table S2. Differentially Expressed Genes — Tissue Contrast (BsMC vs FsMC)

All 11 significant genes (Benjamini–Hochberg adjusted P < 0.05) from the tissue contrast (Breast vs Foreskin) using limma-voom with model `~ Tissue + State + Dataset` and `duplicateCorrelation(block=Donor)`. Positive log2FC indicates higher expression in Breast MCs; negative log2FC indicates higher expression in Foreskin MCs.

| Rank | Gene Symbol | log2FC (BsMC / FsMC) | AveExpr | t-statistic | P-value | adj. P-value | Genomic Origin / Category |
|:---|:---|:---|:---|:---|:---|:---|:---|
| 1 | RPS4Y1 | -9.61 | 1.61 | -16.03 | 1.07e-12 | 1.07e-08 | Y-chromosome |
| 2 | EIF1AY | -8.22 | 0.85 | -13.53 | 2.65e-11 | 1.33e-07 | Y-chromosome |
| 3 | TXLNGY | -7.42 | 0.57 | -13.06 | 4.88e-11 | 1.63e-07 | Y-chromosome |
| 4 | DDX3Y | -7.72 | 0.91 | -11.66 | 3.32e-10 | 8.33e-07 | Y-chromosome |
| 5 | KDM5D | -7.11 | 0.46 | -10.66 | 1.48e-09 | 2.97e-06 | Y-chromosome |
| 6 | UTY | -5.94 | -0.03 | -10.32 | 2.50e-09 | 4.19e-06 | Y-chromosome |
| 7 | TTTY14 | -5.28 | -0.17 | -10.16 | 3.24e-09 | 4.65e-06 | Y-chromosome |
| 8 | ZFY | -5.24 | -0.05 | -8.76 | 3.88e-08 | 4.87e-05 | Y-chromosome |
| 9 | PPIAP29 | -5.02 | -0.13 | -8.66 | 4.70e-08 | 5.25e-05 | Y-chromosome |
| 10 | TPSD1 | -3.86 | 11.57 | -6.21 | 3.93e-06 | 3.95e-03 | Tissue-specific (Delta-tryptase) |
| 11 | XIST | +4.41 | 3.30 | 5.56 | 1.83e-05 | 1.67e-02 | X-inactivation |

---

### Table S3. Canonical Mast Cell Markers — Conserved Across Tissues

Expression of canonical mast cell markers showing no significant tissue-dependent differences (adj.P > 0.05), confirming a conserved core mast cell identity.

| Gene | Name / Product | Function | log2FC (BsMC / FsMC) | adj. P-value | Status |
|:---|:---|:---|:---|:---|:---|
| KIT | KIT Proto-Oncogene | SCF receptor | +0.04 | 1.00 | Conserved |
| CPA3 | Carboxypeptidase A3 | Secretory granule protease | +0.13 | 1.00 | Conserved |
| TPSAB1 | Tryptase Alpha/Beta 1 | Major MC tryptase | +0.61 | 1.00 | Conserved |
| TPSB2 | Tryptase Beta 2 | Major MC tryptase | -0.25 | 1.00 | Conserved |
| TPSG1 | Tryptase Gamma 1 | Membrane-bound tryptase | +0.63 | 1.00 | Conserved |
| FCER1A | Fc Epsilon Receptor Ia | High-affinity IgE receptor alpha | +0.48 | 1.00 | Conserved |
| FCER1G | Fc Epsilon Receptor Ig | High-affinity IgE receptor gamma | +0.51 | 1.00 | Conserved |
| HDC | Histidine Decarboxylase | Histamine synthesis | -0.24 | 1.00 | Conserved |
| CMA1 | Chymase 1 | Serine protease | -1.39 | 0.93 | Conserved |

---

### Table S4. Top Activation DE Genes (Stimulated vs Native)

Top 10 upregulated and top 10 downregulated genes upon IgE-mediated stimulation across N=21 samples (Total: 3,456 DE genes, adj.P < 0.05).

#### Top Upregulated Genes
| Rank | Gene | log2FC | P-value | adj. P-value | Function |
|:---|:---|:---|:---|:---|:---|
| 1 | CCL1 | +10.78 | 3.67e-13 | 3.69e-09 | T-cell chemoattractant |
| 2 | AQP2 | +8.26 | 7.34e-13 | 3.69e-09 | Water channel |
| 3 | XIRP1 | +10.06 | 3.64e-11 | 1.22e-07 | Actin-binding |
| 4 | GBP2 | +4.61 | 9.99e-11 | 2.51e-07 | Guanylate-binding protein |
| 5 | GZMB | +8.27 | 2.06e-10 | 4.15e-07 | Cytotoxic serine protease |
| 6 | NRP2 | +2.31 | 2.07e-10 | 4.15e-07 | VEGF co-receptor |
| 7 | IL3 | +8.59 | 4.12e-09 | 5.92e-06 | Interleukin-3 |
| 8 | CCL2 | +6.22 | 4.20e-09 | 6.03e-06 | Monocyte chemoattractant |
| 9 | FASLG | +5.80 | 4.20e-09 | 6.03e-06 | Fas ligand |
| 10 | CXCL8 | +5.50 | 5.26e-06 | 5.29e-03 | Neutrophil chemoattractant (IL-8) |

#### Top Downregulated Genes
| Rank | Gene | log2FC | P-value | adj. P-value | Function |
|:---|:---|:---|:---|:---|:---|
| 1 | FLI1 | -1.59 | 2.23e-09 | 3.74e-06 | Transcription factor |
| 2 | STX3 | -1.51 | 2.61e-09 | 3.74e-06 | Vesicle fusion syntaxin |
| 3 | ARRB1 | -1.34 | 1.89e-08 | 2.11e-05 | Beta-arrestin 1 |
| 4 | RASSF1 | -1.81 | 2.50e-08 | 2.51e-05 | Tumor suppressor |
| 5 | LAIR1 | -1.60 | 7.62e-08 | 6.38e-05 | Inhibitory immune receptor |

---

### Table S5. Tissue x State Interaction Analysis

No significant interaction genes were detected (0 genes with adj.P < 0.05), confirming that the IgE-mediated activation response is tissue-independent.

| Gene | Interaction log2FC | adj. P-value | Interpretation |
|:---|:---|:---|:---|
| GZMB | -4.94 | 0.99 | Non-significant |
| TNFRSF4 | +3.53 | 0.99 | Non-significant |
| CCL1 | -3.58 | 0.99 | Non-significant |
| FASLG | -3.46 | 0.99 | Non-significant |
| TPSD1 | +0.32 | 0.99 | Non-significant |
| CXCL8 | -1.27 | 0.99 | Non-significant |

---

## Supplementary Figures

### Figure S1. Batch-Corrected PCA Plot
![Batch-Corrected PCA](c:/Users/guerkan.bal/OneDrive - Charité - Universitätsmedizin Berlin/DEG_Fantom6/data/fantom5_hg38/limma_results/PCA_batch_corrected.png)

### Figure S2. Volcano Plot — Tissue Contrast (BsMC vs FsMC)
![Volcano Plot](c:/Users/guerkan.bal/OneDrive - Charité - Universitätsmedizin Berlin/DEG_Fantom6/data/fantom5_hg38/limma_results/Volcano_BsMC_vs_FsMC.png)

### Figure S3. Expression of Key Mast Cell & Tissue Genes
![Key Genes Panel](c:/Users/guerkan.bal/OneDrive - Charité - Universitätsmedizin Berlin/DEG_Fantom6/data/fantom5_hg38/limma_results/Boxplots_KeyGenes_Panel.png)

### Figure S4. Heatmap — Top 30 Tissue DE Genes
![Heatmap Tissue](c:/Users/guerkan.bal/OneDrive - Charité - Universitätsmedizin Berlin/DEG_Fantom6/data/fantom5_hg38/limma_results/Heatmap_Top30_Tissue_DE.png)

### Figure S5. Heatmap — Top 40 Activation DE Genes
![Heatmap Activation](c:/Users/guerkan.bal/OneDrive - Charité - Universitätsmedizin Berlin/DEG_Fantom6/data/fantom5_hg38/limma_results/Heatmap_Top40_Activation_DE.png)

### Figure S6. Publication Concordance: Tissue Effect (Breast vs Foreskin)
![Tissue Concordance](c:/Users/guerkan.bal/OneDrive - Charité - Universitätsmedizin Berlin/DEG_Fantom6/data/fantom5_hg38/limma_results/Concordance_Tissue_Excel_vs_Limma.png)

### Figure S7. Publication Concordance: State Effect (Stimulated vs Native)
![State Concordance](c:/Users/guerkan.bal/OneDrive - Charité - Universitätsmedizin Berlin/DEG_Fantom6/data/fantom5_hg38/limma_results/Concordance_State_Excel_vs_Limma.png)

---

## Supplementary Methods

### Data Processing Pipeline & Statistical Rationale
1. **Count Matrix Construction**: FANTOM6 CAT gene IDs were mapped to HGNC symbols. FANTOM6 technical replicates (2–3 per sample) were summed at raw count level. FANTOM5 hg38-remapped counts were merged by exact HGNC symbol (10,047 shared genes across N=21 samples).
2. **Filtering & Normalization**: Lowly expressed genes were removed using `edgeR::filterByExpr()`. Count matrices were normalized via TMM (`edgeR::calcNormFactors()`) and transformed to logCPM using `limma::voom()`.
3. **Linear Mixed Modeling**: Model formula `~ Tissue + State + Dataset` with `limma::duplicateCorrelation(block=meta$Donor)` to estimate intra-donor consensus correlation (r = 0.433).

