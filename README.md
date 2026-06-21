# Fantom Integration F5 & F6

## Integrative Transcriptomic Analysis of Foreskin-derived vs Breast-derived Mast Cells

[![DOI](https://img.shields.io/badge/Platform-CAGE--seq-blue)]()
[![License](https://img.shields.io/badge/License-MIT-green)]()
[![R](https://img.shields.io/badge/R-≥4.4-blue)]()

This repository contains the **complete, reproducible analysis pipeline** — input data, scripts, results, and documentation — for the integrative differential gene expression analysis comparing **Foreskin-derived Mast Cells (FsMC)** and **Breast-derived Mast Cells (BsMC)** across FANTOM5 and FANTOM6 CAGE-seq platforms.

---

## Key Findings

| Finding | Detail |
|---------|--------|
| 🧬 **TPSD1 validated** | δ-Tryptase significantly upregulated in FsMC (log₂FC = +3.76, adj.P = 0.002) |
| 📊 **17 tissue DE genes** | adj.P < 0.05; 12 sex-linked, 5 tissue-specific (incl. TPSD1, GPX1, CLDN23, MPDZ, C2) |
| 🔥 **2,349 activation DE genes** | IgE-mediated stimulation: 1,224 up / 1,125 down |
| 🔄 **0 interaction genes** | Activation response is **tissue-independent** |
| ✅ **Conserved MC identity** | KIT, CPA3, TPSAB1, TPSB2, FCER1A, FCER1G, HDC — no tissue difference |

---

## Statistical Method

**Unpaired, fixed-effects limma-voom** with model:

```r
design <- model.matrix(~ Tissue + State + Dataset, data = meta)
```

| Term | Levels | Role |
|------|--------|------|
| **Tissue** | Breast / Foreskin | Primary contrast |
| **State** | Native / Stimulated | Covariate (adjusts for activation) |
| **Dataset** | FANTOM5 / FANTOM6 | Covariate (adjusts for platform batch) |

### Why unpaired?

Donors are **completely nested within tissue** — no single donor contributes both foreskin and breast tissue. A mixed-effects model with `(1|Donor)` was evaluated but rejected: it cannot reduce variance for the cross-tissue contrast and artificially constrained degrees of freedom (yielding only 14 vs 17 DE genes). See [`docs/integration_logic.md`](docs/integration_logic.md) for the full rationale.

---

## Repository Structure

```
Fantom-Integration-F5-F6/
│
├── data/                              # ── Input Data ──
│   ├── fantom6/
│   │   ├── BR1256789_Native_Stimulated_gene_counts.txt   # F6 raw CAGE counts
│   │   ├── BR1256789_Native_Stimulated_sampleinfo.txt    # F6 sample metadata
│   │   └── F6_CAT.gene.info_shorted.tsv                  # CAT gene ID → HGNC mapping
│   └── fantom5_hg38/
│       └── fantom5_hg38_mast_cell_gene_counts.tsv        # F5 raw counts (hg38)
│
├── scripts/                           # ── Analysis Scripts ──
│   ├── run_limma_integration.R        # Main DE analysis (limma-voom, unpaired)
│   ├── run_extended_analysis.R        # Boxplots, heatmaps, interaction model
│   └── run_dream_analysis.R           # Exploratory mixed-effects (rejected)
│
├── results/                           # ── Output Files ──
│   ├── DE_FsMC_vs_BsMC.tsv            # Tissue DE table (10,047 genes)
│   ├── DE_Stimulated_vs_Native.tsv    # Activation DE table (10,047 genes)
│   ├── DE_Interaction_Tissue_x_State.tsv  # Interaction term results
│   ├── limma_sampleinfo.tsv           # Final metadata (21 samples)
│   ├── PCA_batch_corrected.png        # PCA (batch-corrected for visualization)
│   ├── Volcano_FsMC_vs_BsMC.png       # Volcano plot (tissue contrast)
│   ├── Heatmap_Top40_Activation_DE.png    # Activation heatmap (top 40 genes)
│   └── Boxplots_KeyGenes_Panel.png    # Key gene boxplots (12 genes)
│
├── docs/                              # ── Documentation ──
│   └── integration_logic.md           # Detailed pipeline & design rationale
│
├── Integrative_Analysis_Report_FsMC_vs_BsMC.docx   # Full report
├── Supplementary_Materials.docx                     # Tables S1–S5, Figures S1–S5
└── README.md
```

---

## Input Data

| File | Dataset | Format | Description |
|------|---------|--------|-------------|
| `BR1256789_Native_Stimulated_gene_counts.txt` | FANTOM6 | TSV, ~24 MB | Raw CAGE counts. Rows = 63,276 CAT gene IDs. Columns = 90 CAGE libraries (7 donors × 3 ASO conditions × Native/Stimulated × 2–3 tech reps). |
| `BR1256789_Native_Stimulated_sampleinfo.txt` | FANTOM6 | TSV, 6 KB | Library-level metadata: `library_name`, `SampleName`, `library.size`, `Rep` (donor), `State`, `Tissue`, `Batch`. |
| `F6_CAT.gene.info_shorted.tsv` | FANTOM6 | TSV, ~7 MB | Gene annotation file mapping FANTOM6 CAT gene identifiers to HGNC symbols. Entries with `__na` have no HGNC mapping. |
| `fantom5_hg38_mast_cell_gene_counts.tsv` | FANTOM5 | TSV, ~780 KB | Raw CAGE counts remapped to hg38 (STAR aligner). Rows = HGNC gene symbols. 9 breast-derived mast cell samples from 5 donors. |

---

## Integration Pipeline

```
FANTOM6 Raw Counts                      FANTOM5 Raw Counts (hg38)
  (63,276 CAT IDs × 90 libs)             (HGNC symbols × 9 samples)
         │                                         │
    ① Select NC-only columns               ⑤ All 9 MC samples
         │                                         │
    ② CAT → HGNC symbol mapping            ⑥ Map "expanded" → Native
         │                                    "expanded+stimulated" → Stimulated
    ③ Aggregate multi-mapping genes                │
         │                                         │
    ④ Sum technical replicates                     │
         │                                         │
    12 biological samples                    9 biological samples
         │                                         │
         └───────────── ⑦ INNER JOIN ──────────────┘
                              │
                    10,047 shared genes
                    21 biological samples
                              │
                    ⑧ filterByExpr() — remove low-count genes
                              │
                    ⑨ TMM normalization (edgeR)
                              │
                    ⑩ voom transformation (limma)
                              │
                    ⑪ lmFit + eBayes: ~ Tissue + State + Dataset
                              │
                 ┌────────────┼────────────┐
                 ↓            ↓            ↓
           Tissue DE     Activation    Interaction
           17 genes      2,349 genes   0 genes
```

For step-by-step details with code snippets, see [`docs/integration_logic.md`](docs/integration_logic.md).

---

## Sample Composition

### FANTOM6 — Tissue Assignment from Sampleinfo

| Donor | Tissue | Sex | Native | Stimulated |
|-------|--------|-----|--------|------------|
| BR1 | Foreskin | ♂ | ✅ | — |
| BR2 | Foreskin | ♂ | ✅ | — |
| BR6 | Foreskin | ♂ | ✅ | ✅ |
| BR7 | Foreskin | ♂ | ✅ | ✅ |
| BR9 | Foreskin | ♂ | ✅ | ✅ |
| BR5 | Breast | ♀ | ✅ | ✅ |
| BR8 | Breast | ♀ | ✅ | ✅ |

### FANTOM5 — All Breast-derived

| Donor | State Mapping | Original Label |
|-------|--------------|----------------|
| F5_D1 | Native | Mast cell, donor 1 |
| F5_D1 | Stimulated | Mast cell – stimulated, donor 1 |
| F5_D2 | Native | Mast cell, donor 2 |
| F5_D3 | Native | Mast cell, donor 3 |
| F5_D4 | Native | Mast cell, donor 4 |
| F5_D5 | Native | Mast cell, **expanded**, donor 5 |
| F5_D5 | Stimulated | Mast cell, **expanded and stimulated**, donor 5 |
| F5_D8 | Native | Mast cell, **expanded**, donor 8 |
| F5_D8 | Stimulated | Mast cell, **expanded and stimulated**, donor 8 |

### Summary

| Dataset | Tissue | Native | Stimulated | Total |
|---------|--------|--------|------------|-------|
| FANTOM6 | Foreskin | 5 | 3 | 8 |
| FANTOM6 | Breast | 2 | 2 | 4 |
| FANTOM5 | Breast | 6 | 3 | 9 |
| **Total** | | **13** | **8** | **21** |

---

## Scripts

### `scripts/run_limma_integration.R` — Main Analysis

The primary analysis script performing the complete pipeline:
1. Loads F6 raw counts and selects NC-only columns
2. Maps CAT gene IDs → HGNC symbols
3. Sums technical replicates per biological sample
4. Loads F5 hg38-remapped counts
5. Merges both datasets by HGNC symbol (inner join)
6. Builds sample metadata (Tissue, State, Dataset, Donor)
7. Filters low-count genes (`filterByExpr`), applies TMM normalization
8. Runs limma-voom with `~ Tissue + State + Dataset`
9. Extracts tissue and activation contrasts
10. Generates PCA and volcano plots

### `scripts/run_extended_analysis.R` — Extended Visualisation & Interaction

- Boxplots of 12 key genes (TPSD1, CXCL8, CCL1, RPS4Y1, XIST, GPX1, tryptases, KIT, FCER1G)
- Heatmaps of top tissue and activation DE genes (pheatmap, Ward's D2 clustering)
- Tissue × State interaction model (`~ Tissue * State + Dataset`)

### `scripts/run_dream_analysis.R` — Exploratory (Not Used)

Mixed-effects model with `(1|Donor)` random effect using `variancePartition::dream`. Rejected for final analysis due to loss of power from nested donor structure.

---

## Requirements

```r
# Core
install.packages(c("data.table", "ggplot2", "pheatmap", "gridExtra"))

# Bioconductor
BiocManager::install(c("limma", "edgeR", "variancePartition"))
```

- **R** ≥ 4.4
- **limma** ≥ 3.60
- **edgeR** ≥ 4.2
- **variancePartition** ≥ 1.34
- **ggplot2** ≥ 3.5
- **data.table** ≥ 1.15
- **pheatmap** ≥ 1.0
- **gridExtra** ≥ 2.3

---

## How to Run

```r
# Set your working directory to the cloned repo
setwd("/path/to/Fantom-Integration-F5-F6")

# 1. Main DE analysis — generates DE tables, PCA, volcano
source("scripts/run_limma_integration.R")

# 2. Extended analyses — boxplots, heatmaps, interaction model
source("scripts/run_extended_analysis.R")
```

> ⚠️ **Note:** The scripts contain hard-coded paths (`setwd()`). Update these to match your local directory structure before running.

---

## Contact

**Gürkan BAL**  
Charité – Universitätsmedizin Berlin  
📧 guerkan.bal@charite.de
