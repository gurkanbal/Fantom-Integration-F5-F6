##############################################################################
# run_limma_integration.R
# Integrative BsMC vs FsMC analysis using limma-voom
# Model: ~ Tissue + State + Dataset (with duplicateCorrelation)
# Tissue Contrast: BsMC vs FsMC (Breast / Foreskin)
# Cohort: N=21 biological samples (13 unique donors)
# F6: NC-only samples, tech reps SUMMED at raw count level
##############################################################################

library(data.table)
library(edgeR)
library(limma)
library(ggplot2)

setwd("c:/Users/guerkan.bal/OneDrive - Charité - Universitätsmedizin Berlin/DEG_Fantom6")
out_dir <- "data/fantom5_hg38/limma_results"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# Helper function to find input data files (supporting gz and uncompressed)
find_data_file <- function(paths) {
  for (p in paths) {
    if (file.exists(p)) return(p)
  }
  stop("Input file not found in paths: ", paste(paths, collapse=", "))
}

f6_counts_path <- find_data_file(c(
  "data/fantom6/BR1256789_Native_Stimulated_gene_counts.txt.gz",
  "data/fantom6/BR1256789_Native_Stimulated_gene_counts.txt",
  "20220504/BR1256789_Native_Stimulated_gene_counts.txt.gz",
  "20220504/BR1256789_Native_Stimulated_gene_counts.txt"
))

f6_info_path   <- find_data_file(c(
  "data/fantom6/BR1256789_Native_Stimulated_sampleinfo.txt",
  "20220504/BR1256789_Native_Stimulated_sampleinfo.txt"
))

f6_anno_path   <- find_data_file(c(
  "data/fantom6/F6_CAT.gene.info_shorted.tsv.gz",
  "data/fantom6/F6_CAT.gene.info_shorted.tsv",
  "20220504/F6_CAT.gene.info_shorted.tsv.gz",
  "20220504/F6_CAT.gene.info_shorted.tsv"
))

f5_counts_path <- find_data_file(c(
  "data/fantom5_hg38/fantom5_hg38_mast_cell_gene_counts.tsv.gz",
  "data/fantom5_hg38/fantom5_hg38_mast_cell_gene_counts.tsv"
))

###########################################################################
# 1. LOAD F6 DATA (NC-only)
###########################################################################
message("Loading F6 raw counts (NC-only)...")
f6_all <- as.data.table(read.delim(f6_counts_path, check.names=TRUE, stringsAsFactors=FALSE))
setnames(f6_all, names(f6_all)[1], "geneID")

# Keep only NC (Negative Control ASO) columns
nc_cols <- grep("aso_NC", names(f6_all), value=TRUE)
f6_nc <- f6_all[, c("geneID", nc_cols), with=FALSE]

# Map CAT IDs to HGNC symbols
f6_anno <- as.data.table(read.delim(f6_anno_path, stringsAsFactors=FALSE))
f6_merged <- merge(f6_nc, f6_anno[, .(geneID, HGNC_symbol)], by="geneID", all.x=TRUE)
f6_annotated <- f6_merged[HGNC_symbol != "__na" & !is.na(HGNC_symbol) & HGNC_symbol != ""]
f6_gene_cols <- setdiff(names(f6_annotated), c("geneID", "HGNC_symbol"))
f6_agg <- f6_annotated[, lapply(.SD, sum, na.rm=TRUE), by=HGNC_symbol, .SDcols=f6_gene_cols]

###########################################################################
# 2. SUM TECHNICAL REPLICATES (F6)
###########################################################################
message("Summing technical replicates per biological sample...")
f6_sample_cols <- setdiff(names(f6_agg), "HGNC_symbol")
base_names <- sub("\\.[0-9]+$", "", f6_sample_cols)
unique_bases <- unique(base_names)

f6_summed <- data.table(HGNC_symbol = f6_agg$HGNC_symbol)
for (bn in unique_bases) {
  cols_to_sum <- f6_sample_cols[base_names == bn]
  if (length(cols_to_sum) == 1) {
    f6_summed[[bn]] <- f6_agg[[cols_to_sum]]
  } else {
    f6_summed[[bn]] <- rowSums(as.matrix(f6_agg[, cols_to_sum, with=FALSE]))
  }
}

###########################################################################
# 3. LOAD F5 DATA (all mast cell samples)
###########################################################################
message("Loading F5 mast cell raw counts...")
f5_counts <- fread(f5_counts_path)
setnames(f5_counts, names(f5_counts)[1], "HGNC_symbol")

###########################################################################
# 4. MERGE DATASETS BY EXACT HGNC SYMBOL
###########################################################################
message("Merging F5 and F6 count matrices...")
combined <- merge(f6_summed, f5_counts, by="HGNC_symbol", all=FALSE)

###########################################################################
# 5. BUILD METADATA COHORT (N=21)
###########################################################################
samples <- setdiff(names(combined), "HGNC_symbol")
meta <- data.frame(SampleID = samples, stringsAsFactors = FALSE)

f6_sampleinfo <- fread(f6_info_path)
tissue_map <- unique(f6_sampleinfo[, .(Rep, Tissue)])

for (i in seq_len(nrow(meta))) {
  sid <- meta$SampleID[i]
  
  if (grepl("^BR", sid)) {
    meta$Dataset[i] <- "FANTOM6"
    meta$Donor[i] <- sub("^(BR[0-9]+).*", "\\1", sid)
    meta$State[i] <- ifelse(grepl("AER", sid), "Stimulated", "Native")
    donor_tissue <- tissue_map$Tissue[tissue_map$Rep == meta$Donor[i]]
    meta$Tissue[i] <- donor_tissue[1]
  } else {
    meta$Dataset[i] <- "FANTOM5"
    meta$Tissue[i] <- "Breast"
    if (grepl("stimulated.*donor1", sid, ignore.case=TRUE)) {
      meta$Donor[i] <- "F5_D1"; meta$State[i] <- "Stimulated"
    } else if (grepl("donor1", sid, ignore.case=TRUE)) {
      meta$Donor[i] <- "F5_D1"; meta$State[i] <- "Native"
    } else if (grepl("donor2", sid, ignore.case=TRUE)) {
      meta$Donor[i] <- "F5_D2"; meta$State[i] <- "Native"
    } else if (grepl("donor3", sid, ignore.case=TRUE)) {
      meta$Donor[i] <- "F5_D3"; meta$State[i] <- "Native"
    } else if (grepl("donor4", sid, ignore.case=TRUE)) {
      meta$Donor[i] <- "F5_D4"; meta$State[i] <- "Native"
    } else if (grepl("expanded.*stimulated.*donor5", sid, ignore.case=TRUE)) {
      meta$Donor[i] <- "F5_D5"; meta$State[i] <- "Stimulated"
    } else if (grepl("expanded.*stimulated.*donor8", sid, ignore.case=TRUE)) {
      meta$Donor[i] <- "F5_D8"; meta$State[i] <- "Stimulated"
    } else if (grepl("expanded.*donor5", sid, ignore.case=TRUE)) {
      meta$Donor[i] <- "F5_D5"; meta$State[i] <- "Native"
    } else if (grepl("expanded.*donor8", sid, ignore.case=TRUE)) {
      meta$Donor[i] <- "F5_D8"; meta$State[i] <- "Native"
    }
  }
}

# Factor levels: Foreskin is reference for Tissue so TissueBreast coefficient = BsMC vs FsMC
meta$Tissue <- factor(meta$Tissue, levels=c("Foreskin", "Breast"))
meta$State <- factor(meta$State, levels=c("Native", "Stimulated"))
meta$Dataset <- factor(meta$Dataset, levels=c("FANTOM6", "FANTOM5"))
rownames(meta) <- meta$SampleID

write.table(meta, file=file.path(out_dir, "limma_sampleinfo.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

###########################################################################
# 6. DGEList, FILTERING & TMM NORMALIZATION
###########################################################################
message("Filtering expressed genes & TMM normalization...")
counts_mat <- as.matrix(combined[, -1, with=FALSE])
rownames(counts_mat) <- combined$HGNC_symbol
counts_mat <- counts_mat[, meta$SampleID]

dge <- DGEList(counts=counts_mat, samples=meta)
keep <- filterByExpr(dge, group=dge$samples$Tissue)
dge <- dge[keep, , keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge, method="TMM")

###########################################################################
# 7. VOOM + LIMMA (~ Tissue + State + Dataset, duplicateCorrelation)
###########################################################################
message("Running limma-voom with duplicateCorrelation(block=Donor)...")
design <- model.matrix(~ Tissue + State + Dataset, data=meta)

v <- voom(dge, design, plot=FALSE)
corfit <- duplicateCorrelation(v, design, block=meta$Donor)
message(paste("Consensus correlation:", round(corfit$consensus.correlation, 3)))

fit <- lmFit(v, design, block=meta$Donor, correlation=corfit$consensus.correlation)
fit <- eBayes(fit)

# Tissue contrast (Breast vs Foreskin -> BsMC vs FsMC)
tt_tissue <- topTable(fit, coef="TissueBreast", number=Inf)
write.table(tt_tissue, file=file.path(out_dir, "DE_BsMC_vs_FsMC.tsv"),
            sep="\t", quote=FALSE)

n_sig_tissue <- sum(tt_tissue$adj.P.Val < 0.05)
message(paste("Tissue DE genes (BsMC vs FsMC, adj.P < 0.05):", n_sig_tissue))

# State contrast (Stimulated vs Native)
tt_state <- topTable(fit, coef="StateStimulated", number=Inf)
write.table(tt_state, file=file.path(out_dir, "DE_Stimulated_vs_Native.tsv"),
            sep="\t", quote=FALSE)
message(paste("Activation DE genes (Stimulated vs Native, adj.P < 0.05):", sum(tt_state$adj.P.Val < 0.05)))

###########################################################################
# 8. PCA PLOT (Batch-corrected)
###########################################################################
message("Generating Batch-Corrected PCA plot...")
expr_corrected <- removeBatchEffect(v$E, batch=meta$Dataset,
                                     design=model.matrix(~ Tissue + State, data=meta))

pca <- prcomp(t(expr_corrected), scale.=TRUE)
pca_df <- as.data.frame(pca$x)
pca_df$SampleID <- rownames(pca_df)
pca_df <- merge(pca_df, meta, by="SampleID")
var_exp <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

p_pca <- ggplot(pca_df, aes(x=PC1, y=PC2, color=Tissue, shape=State)) +
  geom_point(size=4, alpha=0.8) +
  theme_bw(base_size=14) +
  labs(title="Batch-Corrected PCA (NC-only, N=21)",
       subtitle="~ Tissue + State + Dataset | limma-voom",
       x=paste0("PC1 (", var_exp[1], "%)"),
       y=paste0("PC2 (", var_exp[2], "%)")) +
  scale_color_manual(values=c("Breast"="#E74C3C", "Foreskin"="#3498DB"))

ggsave(file.path(out_dir, "PCA_batch_corrected.png"), plot=p_pca, width=9, height=6, dpi=300)

###########################################################################
# 9. VOLCANO PLOT - Tissue (BsMC vs FsMC)
###########################################################################
message("Generating Volcano plot...")
tt_tissue$gene <- rownames(tt_tissue)
tt_tissue$sig <- ifelse(tt_tissue$adj.P.Val < 0.05 & abs(tt_tissue$logFC) > 1, "Significant", "NS")
top_genes <- head(tt_tissue[order(tt_tissue$adj.P.Val), ], 20)

p_volcano <- ggplot(tt_tissue, aes(x=logFC, y=-log10(adj.P.Val), color=sig)) +
  geom_point(alpha=0.5, size=1.5) +
  scale_color_manual(values=c("Significant"="#E74C3C", "NS"="grey60")) +
  geom_text(data=top_genes, aes(label=gene), size=3.5, vjust=-0.8, color="black") +
  theme_bw(base_size=14) +
  labs(title="Volcano Plot: BsMC vs FsMC (limma-voom)",
       subtitle="~ Tissue + State + Dataset",
       x="log2 Fold Change (Breast / Foreskin)",
       y="-log10(adj. P-Value)") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey40") +
  geom_vline(xintercept=c(-1, 1), linetype="dashed", color="grey40") +
  theme(legend.position="none")

ggsave(file.path(out_dir, "Volcano_BsMC_vs_FsMC.png"), plot=p_volcano, width=9, height=7, dpi=300)

message("Integration analysis completed successfully!")
