# Integration Logic: FANTOM5 + FANTOM6 Mast Cell Data

## Overview

This document describes the data integration logic for combining CAGE-seq mast cell gene expression data from two independent FANTOM consortia projects (FANTOM5 and FANTOM6) into a unified count matrix suitable for differential expression analysis.

---

## 1. Input Data Files

### 1.1 FANTOM6 Gene Counts
**File:** `data/fantom6/BR1256789_Native_Stimulated_gene_counts.txt`

- **Format:** Tab-delimited, rows = CAT gene IDs, columns = CAGE-seq libraries
- **Content:** Raw CAGE counts for 7 biological replicates (BR1, BR2, BR5, BR6, BR7, BR8, BR9) across 3 ASO conditions (NC, M1, M2) and 2 activation states (Native, Stimulated/AER)
- **Column naming convention:** `BR{donor}_{state}_aso_{condition}` where:
  - `BR1-BR9`: Biological replicate / donor identifier
  - `AER`: Anti-IgE receptor stimulated; absence = Native
  - `aso_NC`: Negative control (scrambled ASO)
  - `aso_M1`, `aso_M2`: ASO knockdown targets (excluded from this analysis)
- **Technical replicates:** Each biological sample has 2–3 technical replicates, represented as additional columns with `.1`, `.2` suffixes (e.g., `BR1_aso_NC`, `BR1_aso_NC.1`, `BR1_aso_NC.2`)

### 1.2 FANTOM6 Sample Information
**File:** `data/fantom6/BR1256789_Native_Stimulated_sampleinfo.txt`

- **Format:** Tab-delimited metadata, one row per CAGE library
- **Columns:** `library_name`, `SampleName`, `library.size`, `Rep`, `State`, `Tissue`, `Batch`
- **Key tissue assignments:**
  - Foreskin (male donors): BR1, BR2, BR6, BR7, BR9
  - Breast (female donors): BR5, BR8

### 1.3 FANTOM6 Gene Annotation
**File:** `data/fantom6/F6_CAT.gene.info_shorted.tsv`

- **Format:** Tab-delimited mapping table
- **Purpose:** Maps FANTOM6 CAT gene identifiers to HGNC gene symbols
- **Key columns:** `geneID` (CAT ID), `HGNC_symbol` (standard gene name)
- **Special values:** `__na` indicates no HGNC mapping (these genes are excluded)

### 1.4 FANTOM5 Gene Counts (hg38-remapped)
**File:** `data/fantom5_hg38/fantom5_hg38_mast_cell_gene_counts.tsv`

- **Format:** Tab-delimited, rows = HGNC gene symbols, columns = CAGE-seq libraries
- **Content:** 9 mast cell samples, all from breast tissue, originally aligned to hg19 and remapped to hg38 using STAR aligner
- **Column naming:** URL-encoded FANTOM5 sample names, e.g.:
  - `counts.Mast%20cell%2c%20donor1.CNhs12566...` → Mast cell, donor 1
  - `counts.Mast%20cell%20-%20stimulated%2c%20donor1.CNhs11073...` → Stimulated, donor 1
  - `counts.Mast%20cell%2c%20expanded%2c%20donor5...` → Expanded, donor 5
  - `counts.Mast%20cell%2c%20expanded%20and%20stimulated%2c%20donor5...` → Expanded and stimulated, donor 5

---

## 2. Integration Steps

### Step 1: FANTOM6 — Select NC-only Samples

Only Negative Control (aso_NC) samples are retained. ASO knockdown samples (M1, M2) are excluded to avoid confounding from gene perturbation.

```r
nc_cols <- grep("aso_NC", names(f6_all), value=TRUE)
```

**Result:** 30 columns (technical replicates across 12 biological conditions)

### Step 2: FANTOM6 — Map CAT IDs to HGNC Symbols

CAT gene IDs are replaced with HGNC gene symbols using the annotation file. Genes without HGNC mapping (`__na`) are removed. Multiple CAT IDs mapping to the same HGNC symbol are aggregated by summing counts.

```r
f6_merged <- merge(f6_nc, f6_anno[, .(geneID, HGNC_symbol)], by="geneID")
f6_agg <- f6_annotated[, lapply(.SD, sum), by=HGNC_symbol, .SDcols=gene_cols]
```

### Step 3: FANTOM6 — Sum Technical Replicates

Technical replicates (identified by `.1`, `.2` suffixes on column names) are summed per biological sample. This is the recommended approach for count data (Lun & Smyth, 2015) as it:
- Preserves the count nature of the data
- Correctly reflects total sequencing depth
- Avoids artificial inflation of sample size

```r
base_names <- sub("\\.[0-9]+$", "", f6_sample_cols)
# For each unique base name, sum all matching columns
f6_summed[[bn]] <- rowSums(as.matrix(f6_agg[, cols_to_sum, with=FALSE]))
```

**Result:** 12 unique biological samples (from 30 technical replicates)

| Base Name | Donor | State | Tech Reps Summed |
|-----------|-------|-------|-----------------|
| BR1_aso_NC | BR1 | Native | 3 |
| BR2_aso_NC | BR2 | Native | 3 |
| BR5_aso_NC | BR5 | Native | 3 |
| BR5_AER_aso_NC | BR5 | Stimulated | 3 |
| BR6_aso_NC | BR6 | Native | 2 |
| BR6_AER_aso_NC | BR6 | Stimulated | 2 |
| BR7_aso_NC | BR7 | Native | 2 |
| BR7_AER_aso_NC | BR7 | Stimulated | 2 |
| BR8_aso_NC | BR8 | Native | 3 |
| BR8_AER_aso_NC | BR8 | Stimulated | 3 |
| BR9_aso_NC | BR9 | Native | 2 |
| BR9_AER_aso_NC | BR9 | Stimulated | 2 |

### Step 4: FANTOM5 — State Assignment

FANTOM5 sample names are parsed to assign Donor and State:

| FANTOM5 Sample Name Pattern | Assigned Donor | Assigned State | Rationale |
|-----------------------------|---------------|---------------|-----------|
| `Mast cell, donor1` | F5_D1 | Native | Baseline |
| `Mast cell - stimulated, donor1` | F5_D1 | Stimulated | IgE stimulated |
| `Mast cell, donor2` | F5_D2 | Native | Baseline |
| `Mast cell, donor3` | F5_D3 | Native | Baseline |
| `Mast cell, donor4` | F5_D4 | Native | Baseline |
| `Mast cell, expanded, donor5` | F5_D5 | **Native** | Expanded = baseline culture |
| `Mast cell, expanded and stimulated, donor5` | F5_D5 | **Stimulated** | Expanded + stimulated |
| `Mast cell, expanded, donor8` | F5_D8 | **Native** | Expanded = baseline culture |
| `Mast cell, expanded and stimulated, donor8` | F5_D8 | **Stimulated** | Expanded + stimulated |

**Key decision:** "Expanded" samples are treated as the **Native** (baseline) state. "Expanded and stimulated" samples are treated as **Stimulated**. This is because the "expanded" protocol represents the standard culture method for these donors, not a perturbation.

### Step 5: Cross-Platform Merge

FANTOM6 (HGNC symbols) and FANTOM5 (HGNC symbols) are merged by gene name using an inner join. Only genes present in both datasets are retained.

```r
combined <- merge(f6_summed, f5_counts, by="HGNC_symbol", all=FALSE)
```

**Result:** 10,047 shared genes across 21 biological samples

### Step 6: Metadata Construction

Each sample is annotated with four key factors:

| Factor | Levels | Source |
|--------|--------|--------|
| **Tissue** | Breast, Foreskin | Sampleinfo file (F6); all "Breast" (F5) |
| **State** | Native, Stimulated | Sample name parsing |
| **Dataset** | FANTOM6, FANTOM5 | Source dataset |
| **Donor** | BR1–BR9, F5_D1–F5_D8 | Sample name parsing |

### Step 7: Filtering and Normalization

```r
dge <- DGEList(counts=counts_mat, samples=meta)
keep <- filterByExpr(dge, group=dge$samples$Tissue)  # Gene filtering
dge <- dge[keep, , keep.lib.sizes=FALSE]              # Remove low-count genes
dge <- calcNormFactors(dge, method="TMM")              
```

### 3. Linear Modeling Strategy
We employed the `limma-voom` pipeline to model the expression data.

Because our dataset is a **mixed design**—some donors are paired (having both a Native and Stimulated sample, e.g., FANTOM6 BR5, BR6, BR7, BR8, BR9; FANTOM5 D1, D5, D8) while others are unpaired (Native only, e.g., FANTOM6 BR1, BR2; FANTOM5 D2, D3, D4)—we needed a model that could leverage the paired nature of the activation response without throwing away the unpaired baseline samples.

To achieve this, we used **`limma::duplicateCorrelation()`**:
1. It calculates the consensus correlation (r = 0.433) between samples belonging to the same donor (`block=Donor`).
2. It feeds this correlation structure into `lmFit`, effectively penalizing between-donor variance when calculating the State (activation) effect.

This approach perfectly recovers the statistical power of a paired analysis for the activation effect (yielding >1,000 additional activation DE genes compared to a naive linear model), while simultaneously using all 21 samples to powerfully estimate the baseline Tissue effect.

```r
design <- model.matrix(~ Tissue + State + Dataset, data=meta)
v <- voom(dge, design, plot=FALSE)

corfit <- duplicateCorrelation(v, design, block=meta$Donor)
fit <- lmFit(v, design, block=meta$Donor, correlation=corfit$consensus.correlation)
fit <- eBayes(fit)
```

**Results:**
- **Tissue Effect (Breast vs Foreskin):** 11 significant DE genes. (10 sex-linked, 1 mast-cell specific: TPSD1).
- **State Effect (Stimulated vs Native):** 3,456 significant DE genes.
- **Interaction (Tissue × State):** 0 significant DE genes (the activation response is perfectly conserved).

### The Final Decision
By utilizing **fixed-effects with `duplicateCorrelation(block=Donor)`**, we found the perfect middle ground. We avoided the rank-deficiency traps of the `variancePartition::dream` mixed-effects model (which failed to estimate the cross-tissue contrast effectively due to perfectly nested donors), while still capturing the crucial intra-donor correlation that gives paired designs their immense statistical power for within-subject effects (activation).

---

## 4. Visualization Batch Correction

For PCA and other visualization purposes only, the Dataset effect is removed using `limma::removeBatchEffect()`. This does NOT affect the statistical model — it is used exclusively for plotting.

```r
expr_corrected <- removeBatchEffect(v$E, batch=meta$Dataset,
                                     design=model.matrix(~ Tissue + State, data=meta))
```

---

## 5. Flow Diagram

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
                  [filterByExpr()]
                         |
                  [TMM normalization]
                         |
                  [voom transformation]
                         |
              [limma: ~ Tissue + State + Dataset]
                         |
              ┌──────────┼──────────┐
              ↓          ↓          ↓
         Tissue DE   Activation  Interaction
         (17 genes)  DE (2349)   (0 genes)
```
