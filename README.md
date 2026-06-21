# Fantom Integration F5 & F6

## Integrative Transcriptomic Analysis of Foreskin-derived vs Breast-derived Mast Cells

This repository contains the complete analysis pipeline — input data, scripts, results, and documentation — for the integrative differential gene expression analysis comparing **Foreskin-derived Mast Cells (FsMC)** and **Breast-derived Mast Cells (BsMC)** across FANTOM5 and FANTOM6 CAGE-seq platforms.

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

No paired/random effects design was used because donors are completely nested within tissue (no donor provides both foreskin and breast tissue). See [`docs/integration_logic.md`](docs/integration_logic.md) for a detailed description of the integration pipeline, data processing steps, and design rationale.

## Repository Structure

```
├── data/
│   ├── fantom6/
│   │   ├── BR1256789_Native_Stimulated_gene_counts.txt   # F6 raw CAGE counts (all ASO conditions)
│   │   ├── BR1256789_Native_Stimulated_sampleinfo.txt    # F6 sample metadata (tissue, state, batch)
│   │   └── F6_CAT.gene.info_shorted.tsv                  # CAT gene ID → HGNC symbol mapping
│   └── fantom5_hg38/
│       └── fantom5_hg38_mast_cell_gene_counts.tsv        # F5 raw counts (remapped to hg38)
├── scripts/
│   ├── run_limma_integration.R      # Main DE analysis (limma-voom, unpaired)
│   ├── run_extended_analysis.R      # Boxplots, heatmaps, interaction model
│   └── run_dream_analysis.R         # Exploratory mixed-effects (not used for final results)
├── results/
│   ├── DE_FsMC_vs_BsMC.tsv          # Full tissue DE table (10,047 genes)
│   ├── DE_Stimulated_vs_Native.tsv  # Full activation DE table
│   ├── DE_Interaction_Tissue_x_State.tsv  # Interaction results
│   ├── limma_sampleinfo.tsv         # Final sample metadata (21 samples)
│   ├── PCA_batch_corrected.png      # PCA plot (batch-corrected for visualization)
│   ├── Volcano_FsMC_vs_BsMC.png     # Volcano plot (tissue contrast)
│   ├── Heatmap_Top40_Activation_DE.png  # Activation heatmap
│   └── Boxplots_KeyGenes_Panel.png  # Key gene boxplots (12 genes)
├── docs/
│   └── integration_logic.md         # Detailed integration pipeline documentation
├── Integrative_Analysis_Report_FsMC_vs_BsMC.docx  # Full analysis report
├── Supplementary_Materials.docx     # Supplementary tables S1–S5 and figures S1–S5
└── README.md
```

## Input Data

| File | Dataset | Description | Size |
|------|---------|-------------|------|
| `BR1256789_Native_Stimulated_gene_counts.txt` | FANTOM6 | Raw CAGE counts; rows = CAT gene IDs, columns = CAGE libraries (NC, M1, M2 × Native/Stimulated × 2–3 tech reps) | ~24 MB |
| `BR1256789_Native_Stimulated_sampleinfo.txt` | FANTOM6 | Library-level metadata: donor, tissue, state, batch | 6 KB |
| `F6_CAT.gene.info_shorted.tsv` | FANTOM6 | Gene annotation: CAT gene ID → HGNC symbol mapping | ~7 MB |
| `fantom5_hg38_mast_cell_gene_counts.tsv` | FANTOM5 | Raw CAGE counts remapped to hg38; 9 mast cell samples, all breast tissue | ~780 KB |

## Integration Pipeline

```
FANTOM6 Raw Counts                    FANTOM5 Raw Counts (hg38)
       |                                       |
  [Select NC-only]                     [All 9 MC samples]
       |                                       |
  [CAT → HGNC mapping]                [Already HGNC symbols]
       |                                       |
  [Sum technical replicates]           [Map expanded → Native]
       |                                       |
  12 biological samples                  9 samples
       |                                       |
       +───────────── INNER JOIN ──────────────+
                         |
                  10,047 shared genes
                  21 biological samples
                         |
                  [filterByExpr → TMM → voom]
                         |
              [limma: ~ Tissue + State + Dataset]
                         |
              ┌──────────┼──────────┐
              ↓          ↓          ↓
         Tissue DE   Activation  Interaction
         (17 genes)  DE (2,349)  (0 genes)
```

For full details on each processing step, see [`docs/integration_logic.md`](docs/integration_logic.md).

## Sample Composition

| Dataset | Tissue | Native | Stimulated | Total |
|---------|--------|--------|------------|-------|
| FANTOM6 | Foreskin | 5 | 3 | 8 |
| FANTOM6 | Breast | 2 | 2 | 4 |
| FANTOM5 | Breast | 6 | 3 | 9 |
| **Total** | | **13** | **8** | **21** |

## Requirements

- R ≥ 4.4
- limma, edgeR, ggplot2, data.table, pheatmap, gridExtra, variancePartition

## How to Run

```r
# 1. Main DE analysis
source("scripts/run_limma_integration.R")

# 2. Extended analyses (boxplots, heatmaps, interaction)
source("scripts/run_extended_analysis.R")
```

> **Note:** The scripts contain hard-coded paths. Update `setwd()` and file paths to match your local directory structure before running.

## Contact

Gürkan BAL — guerkan.bal@charite.de  
Charité – Universitätsmedizin Berlin
