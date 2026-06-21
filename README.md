# Fantom Integration F5 & F6

## Integrative Transcriptomic Analysis of Foreskin-derived vs Breast-derived Mast Cells

This repository contains the analysis code, results, and supplementary materials for the integrative differential gene expression analysis comparing **Foreskin-derived Mast Cells (FsMC)** and **Breast-derived Mast Cells (BsMC)** across FANTOM5 and FANTOM6 CAGE-seq platforms.

## Key Findings

- **TPSD1** (δ-Tryptase) is significantly upregulated in foreskin mast cells (log₂FC = +3.76, adj.P = 0.002)
- **17 tissue DE genes** (adj.P < 0.05), 12 of which are sex-linked
- **2,349 activation DE genes** upon IgE-mediated stimulation
- **0 interaction genes** — the activation response is tissue-independent
- Core mast cell markers (KIT, CPA3, TPSAB1, FCER1G) show no tissue differences

## Statistical Method

**Unpaired, fixed-effects limma-voom** with model:

```
~ Tissue + State + Dataset
```

- **Tissue**: Breast vs Foreskin (primary contrast)
- **State**: Native vs Stimulated (covariate)
- **Dataset**: FANTOM5 vs FANTOM6 (batch covariate)

No paired/random effects design was used because donors are completely nested within tissue (no donor provides both foreskin and breast tissue).

## Repository Structure

```
├── scripts/
│   ├── run_limma_integration.R      # Main DE analysis
│   ├── run_extended_analysis.R      # Boxplots, heatmaps, interaction model
│   └── run_dream_analysis.R         # Exploratory mixed-effects (not final)
├── results/
│   ├── DE_FsMC_vs_BsMC.tsv          # Full tissue DE table (10,047 genes)
│   ├── DE_Stimulated_vs_Native.tsv  # Full activation DE table
│   ├── DE_Interaction_Tissue_x_State.tsv  # Interaction results
│   ├── limma_sampleinfo.tsv         # Sample metadata
│   ├── PCA_batch_corrected.png      # PCA plot
│   ├── Volcano_FsMC_vs_BsMC.png     # Volcano plot
│   ├── Heatmap_Top40_Activation_DE.png  # Activation heatmap
│   └── Boxplots_KeyGenes_Panel.png  # Key gene boxplots
├── Integrative_Analysis_Report_FsMC_vs_BsMC.docx  # Full report
├── Supplementary_Materials.docx     # Supplementary tables and figures
└── README.md
```

## Data Sources

| Dataset | Platform | Samples |
|---------|----------|---------|
| FANTOM6 | CAGE-seq (lentiviral) | 12 biological (NC-only, tech reps summed) |
| FANTOM5 | CAGE-seq (HeliScope) | 9 (remapped to hg38) |

## Requirements

- R ≥ 4.4
- limma, edgeR, ggplot2, data.table, pheatmap, gridExtra, variancePartition

## Contact

Gürkan BAL — guerkan.bal@charite.de  
Charité – Universitätsmedizin Berlin
