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
| 📊 **11 tissue DE genes** | adj.P < 0.05; 10 sex-linked, 1 tissue-specific (TPSD1) |
| 🔥 **3,456 activation DE genes** | IgE-mediated stimulation perfectly captured across both tissues |
| 🔄 **0 interaction genes** | Activation response is **tissue-independent** |
| ✅ **Conserved MC identity** | KIT, CPA3, TPSAB1, TPSB2, FCER1A, FCER1G, HDC — no tissue difference |

---

## Statistical Method

**Mixed-design limma-voom** with `duplicateCorrelation` to handle paired donors:

```r
design <- model.matrix(~ Tissue + State + Dataset, data = meta)
corfit <- duplicateCorrelation(v, design, block = meta$Donor)
```

| Term | Levels | Role |
|------|--------|------|
| **Tissue** | Breast / Foreskin | Primary contrast |
| **State** | Native / Stimulated | Covariate (adjusts for activation) |
| **Dataset** | FANTOM5 / FANTOM6 | Covariate (adjusts for platform batch) |

### Why duplicateCorrelation?

The dataset is a **mixed design**: 8 donors are paired (Native + Stimulated), while 5 donors are unpaired (Native only). By using `duplicateCorrelation(block=Donor)`, we calculate the intra-donor consensus correlation (r = 0.433) and penalize between-donor variance. This elegantly recovers the high statistical power of a paired design for the State effect (finding >1,000 additional activation genes) while utilizing all 21 samples to robustly estimate the Tissue effect.

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
│   ├── run_limma_integration.R        # FsMC vs BsMC (Foreskin / Breast)
│   ├── run_extended_analysis.R        # Boxplots, heatmaps, interaction
│   ├── run_limma_integration_BvF.R    # BsMC vs FsMC (Breast / Foreskin) ← REVERSED
│   ├── run_extended_analysis_BvF.R    # Extended analysis for BvF contrast
│   ├── run_concordance_plot.R         # Concordance plots (FvB↔BvF & DREAM↔limma)
│   └── run_dream_analysis.R           # Exploratory mixed-effects (rejected)
│
├── results/                           # ── FsMC vs BsMC Results ──
│   ├── DE_FsMC_vs_BsMC.tsv            # Tissue DE table (10,047 genes)
│   ├── DE_Stimulated_vs_Native.tsv    # Activation DE table
│   ├── DE_Interaction_Tissue_x_State.tsv
│   ├── limma_sampleinfo.tsv           # Final metadata (21 samples)
│   ├── PCA_batch_corrected.png
│   ├── Volcano_FsMC_vs_BsMC.png
│   ├── Heatmap_Top40_Activation_DE.png
│   ├── Boxplots_KeyGenes_Panel.png
│   ├── Concordance_FvB_vs_BvF.png     # Self-concordance (r = −1.000)
│   └── Concordance_DREAM_vs_Limma_BvF.png  # DREAM vs limma (r = −0.988)
│
├── results/limma_results_BvF/         # ── BsMC vs FsMC Results (REVERSED) ──
│   ├── DE_BsMC_vs_FsMC.tsv            # Tissue DE table (reversed logFC)
│   ├── DE_Stimulated_vs_Native.tsv
│   ├── DE_Interaction_Tissue_x_State.tsv
│   ├── Volcano_BsMC_vs_FsMC.png
│   ├── PCA_batch_corrected.png
│   ├── Heatmap_*.png
│   ├── Boxplots_KeyGenes_Panel.png
│   └── Concordance_*.png
│
├── docs/
│   └── integration_logic.md           # Pipeline & design rationale
│
├── Integrative_Analysis_Report_FsMC_vs_BsMC.docx
├── Supplementary_Materials.docx
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
                    ⑪ limma::duplicateCorrelation(block=Donor)
                              │
                    ⑫ lmFit + eBayes: ~ Tissue + State + Dataset
                              │
                 ┌────────────┼────────────┐
                 ↓            ↓            ↓
           Tissue DE     Activation    Interaction
           11 genes      3,456 genes   0 genes
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

### `scripts/run_limma_integration.R` — FsMC vs BsMC (primary)

The primary analysis script: positive logFC = higher in Foreskin.
1. Loads F6 raw counts, selects NC-only columns, maps CAT → HGNC, sums tech reps
2. Merges with F5 hg38 counts (inner join → 10,047 genes)
3. Runs limma-voom `~ Tissue + State + Dataset` (ref: Breast)
4. Extracts tissue, activation contrasts; generates PCA + volcano

### `scripts/run_limma_integration_BvF.R` — BsMC vs FsMC (reversed)

Identical pipeline with **reversed reference level**: positive logFC = higher in Breast.
Only the factor levels change: `Tissue = factor(..., levels=c("Foreskin", "Breast"))`. P-values are identical; logFC signs are mirrored.

### `scripts/run_extended_analysis.R` / `_BvF.R` — Extended Analyses

- Boxplots of 12 key genes (TPSD1, CXCL8, CCL1, RPS4Y1, XIST, GPX1, tryptases, KIT, FCER1G)
- Heatmaps of top tissue and activation DE genes (pheatmap, Ward's D2 clustering)
- Tissue × State interaction model (`~ Tissue * State + Dataset`)

### `scripts/run_concordance_plot.R` — Concordance Analysis

Generates two concordance scatter plots:
1. **Self-concordance** (FsMC-vs-BsMC ↔ BsMC-vs-FsMC): validates that reversing the contrast produces perfectly mirrored logFCs (Pearson r = −1.000)
2. **Method concordance** (Manuscript DREAM ↔ Integration limma-voom): compares the mixed-effects model from the manuscript to the fixed-effects integration model (Pearson r = −0.988)

Both include labeled key genes, significance categories, and zoomed insets for core MC genes.

### `scripts/run_dream_analysis.R` — Exploratory (Not Used)

Mixed-effects model with `(1|Donor)` using `variancePartition::dream`. Rejected due to nested donor structure.

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

# 1. Main DE analysis (FsMC vs BsMC) — DE tables, PCA, volcano
source("scripts/run_limma_integration.R")

# 2. Extended analyses — boxplots, heatmaps, interaction model
source("scripts/run_extended_analysis.R")

# 3. Reversed contrast (BsMC vs FsMC) — same analysis, mirrored logFC
source("scripts/run_limma_integration_BvF.R")
source("scripts/run_extended_analysis_BvF.R")

# 4. Concordance plots — compare FvB↔BvF and DREAM↔limma
source("scripts/run_concordance_plot.R")
```

> ⚠️ **Note:** The scripts contain hard-coded paths (`setwd()`). Update these to match your local directory structure before running.

---

## Contact

**Gürkan BAL**  
Charité – Universitätsmedizin Berlin  
📧 guerkan.bal@charite.de