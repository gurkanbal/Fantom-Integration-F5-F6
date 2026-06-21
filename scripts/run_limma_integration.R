##############################################################################
# run_limma_integration.R
# Integrative BsMC vs FsMC analysis using limma-voom
# Model: ~ Tissue + State + Dataset
# F6: NC-only samples, tech reps SUMMED | F5: all samples
##############################################################################

library(data.table)
library(edgeR)
library(limma)
library(ggplot2)

setwd("c:/Users/guerkan.bal/OneDrive - Charité - Universitätsmedizin Berlin/DEG_Fantom6")
out_dir <- "data/fantom5_hg38/limma_results"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

###########################################################################
# 1. LOAD F6 DATA (NC-only)
###########################################################################
message("Loading F6 counts (NC-only)...")
f6_all <- as.data.table(read.delim("20220504/BR1256789_Native_Stimulated_gene_counts.txt",
                                    check.names=TRUE, stringsAsFactors=FALSE))
setnames(f6_all, names(f6_all)[1], "geneID")

# Keep only NC columns
nc_cols <- grep("aso_NC", names(f6_all), value=TRUE)
f6_nc <- f6_all[, c("geneID", nc_cols), with=FALSE]
message(paste("F6 NC columns selected:", length(nc_cols)))

# Map CAT IDs to HGNC
f6_anno <- as.data.table(read.delim("20220504/F6_CAT.gene.info_shorted.tsv", stringsAsFactors=FALSE))
f6_merged <- merge(f6_nc, f6_anno[, .(geneID, HGNC_symbol)], by="geneID", all.x=TRUE)
f6_annotated <- f6_merged[HGNC_symbol != "__na" & !is.na(HGNC_symbol) & HGNC_symbol != ""]
f6_gene_cols <- setdiff(names(f6_annotated), c("geneID", "HGNC_symbol"))
f6_agg <- f6_annotated[, lapply(.SD, sum, na.rm=TRUE), by=HGNC_symbol, .SDcols=f6_gene_cols]

###########################################################################
# 2. SUM TECHNICAL REPLICATES (F6)
###########################################################################
message("Summing technical replicates...")
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
message(paste("F6 biological samples after summing:", length(unique_bases),
              "(from", length(f6_sample_cols), "tech reps)"))

###########################################################################
# 3. LOAD F5 DATA (all samples)
###########################################################################
message("Loading F5 counts...")
f5_counts <- fread("data/fantom5_hg38/fantom5_hg38_mast_cell_gene_counts.tsv")
setnames(f5_counts, names(f5_counts)[1], "HGNC_symbol")

###########################################################################
# 4. MERGE
###########################################################################
message("Merging F5 and F6...")
combined <- merge(f6_summed, f5_counts, by="HGNC_symbol", all=FALSE)
message(paste("Genes after merge:", nrow(combined)))

###########################################################################
# 5. BUILD METADATA
###########################################################################
message("Building metadata...")
samples <- setdiff(names(combined), "HGNC_symbol")
meta <- data.frame(SampleID = samples, stringsAsFactors = FALSE)

# Tissue lookup from official sampleinfo
f6_sampleinfo <- fread("20220504/BR1256789_Native_Stimulated_sampleinfo.txt")
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
    if (grepl("stimulated%2c%20donor1", sid)) {
      meta$Donor[i] <- "F5_D1"; meta$State[i] <- "Stimulated"
    } else if (grepl("donor1", sid)) {
      meta$Donor[i] <- "F5_D1"; meta$State[i] <- "Native"
    } else if (grepl("donor2", sid)) {
      meta$Donor[i] <- "F5_D2"; meta$State[i] <- "Native"
    } else if (grepl("donor3", sid)) {
      meta$Donor[i] <- "F5_D3"; meta$State[i] <- "Native"
    } else if (grepl("donor4", sid)) {
      meta$Donor[i] <- "F5_D4"; meta$State[i] <- "Native"
    } else if (grepl("expanded%20and%20stimulated.*donor5", sid)) {
      meta$Donor[i] <- "F5_D5"; meta$State[i] <- "Stimulated"
    } else if (grepl("expanded%20and%20stimulated.*donor8", sid)) {
      meta$Donor[i] <- "F5_D8"; meta$State[i] <- "Stimulated"
    } else if (grepl("expanded.*donor5", sid)) {
      meta$Donor[i] <- "F5_D5"; meta$State[i] <- "Native"
    } else if (grepl("expanded.*donor8", sid)) {
      meta$Donor[i] <- "F5_D8"; meta$State[i] <- "Native"
    }
  }
}

# Factorize
meta$Tissue <- factor(meta$Tissue, levels=c("Breast", "Foreskin"))
meta$State <- factor(meta$State, levels=c("Native", "Stimulated"))
meta$Dataset <- factor(meta$Dataset, levels=c("FANTOM6", "FANTOM5"))
rownames(meta) <- meta$SampleID

# Print summary
message("\n=== SAMPLE SUMMARY ===")
message(paste("Total biological samples:", nrow(meta)))
cat("\nTissue x State:\n")
print(table(meta$Tissue, meta$State))
cat("\nTissue x Dataset:\n")
print(table(meta$Tissue, meta$Dataset))
cat("\nAll samples:\n")
print(meta[, c("SampleID", "Donor", "Tissue", "State", "Dataset")])

write.table(meta, file=file.path(out_dir, "limma_sampleinfo.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

###########################################################################
# 6. DGEList, FILTER, NORMALIZE
###########################################################################
message("\nFiltering & Normalization...")
counts_mat <- as.matrix(combined[, -1, with=FALSE])
rownames(counts_mat) <- combined$HGNC_symbol
counts_mat <- counts_mat[, meta$SampleID]

dge <- DGEList(counts=counts_mat, samples=meta)
keep <- filterByExpr(dge, group=dge$samples$Tissue)
dge <- dge[keep, , keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge, method="TMM")
message(paste("Genes retained:", sum(keep)))

###########################################################################
# 7. VOOM + LIMMA (~ Tissue + State + Dataset)
###########################################################################
message("Running limma-voom...")
design <- model.matrix(~ Tissue + State + Dataset, data=meta)
message("Design matrix columns:")
print(colnames(design))

v <- voom(dge, design, plot=FALSE)
fit <- lmFit(v, design)
fit <- eBayes(fit)

# Tissue contrast (Foreskin vs Breast)
tt_tissue <- topTable(fit, coef="TissueForeskin", number=Inf)
write.table(tt_tissue, file=file.path(out_dir, "DE_FsMC_vs_BsMC.tsv"),
            sep="\t", quote=FALSE)

n_sig <- sum(tt_tissue$adj.P.Val < 0.05)
n_sig_fc <- sum(tt_tissue$adj.P.Val < 0.05 & abs(tt_tissue$logFC) > 1)
message(paste("DE genes (adj.P < 0.05):", n_sig))
message(paste("DE genes (adj.P < 0.05 & |logFC| > 1):", n_sig_fc))

# State contrast (Stimulated vs Native)
tt_state <- topTable(fit, coef="StateStimulated", number=Inf)
write.table(tt_state, file=file.path(out_dir, "DE_Stimulated_vs_Native.tsv"),
            sep="\t", quote=FALSE)
message(paste("Activation DE genes (adj.P < 0.05):", sum(tt_state$adj.P.Val < 0.05)))

###########################################################################
# 8. PCA PLOT (Batch-corrected)
###########################################################################
message("Generating PCA plot...")
expr_corrected <- removeBatchEffect(v$E, batch=meta$Dataset,
                                     design=model.matrix(~ Tissue + State, data=meta))

pca <- prcomp(t(expr_corrected), scale.=TRUE)
pca_df <- as.data.frame(pca$x)
pca_df$SampleID <- rownames(pca_df)
pca_df <- merge(pca_df, meta, by="SampleID")

var_exp <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

p_pca <- ggplot(pca_df, aes(x=PC1, y=PC2, color=Tissue, shape=State)) +
  geom_point(size=4, alpha=0.8) +
  theme_bw() +
  labs(title="Batch-Corrected PCA (NC-only, tech reps summed)",
       subtitle="~ Tissue + State + Dataset | limma-voom",
       x=paste0("PC1 (", var_exp[1], "%)"),
       y=paste0("PC2 (", var_exp[2], "%)")) +
  scale_color_manual(values=c("Breast"="#E74C3C", "Foreskin"="#3498DB")) +
  theme(text=element_text(size=14))

ggsave(file.path(out_dir, "PCA_batch_corrected.png"), plot=p_pca, width=9, height=6)

###########################################################################
# 9. VOLCANO PLOT - Tissue
###########################################################################
message("Generating Volcano plots...")
tt_tissue$gene <- rownames(tt_tissue)
tt_tissue$sig <- ifelse(tt_tissue$adj.P.Val < 0.05 & abs(tt_tissue$logFC) > 1, "Significant", "NS")
top_genes <- head(tt_tissue[order(tt_tissue$adj.P.Val), ], 20)

p_volcano <- ggplot(tt_tissue, aes(x=logFC, y=-log10(adj.P.Val), color=sig)) +
  geom_point(alpha=0.5, size=1.5) +
  scale_color_manual(values=c("Significant"="#E74C3C", "NS"="grey60")) +
  geom_text(data=top_genes, aes(label=gene), size=3, vjust=-0.8, color="black") +
  theme_bw() +
  labs(title="Volcano Plot: FsMC vs BsMC (limma-voom)",
       subtitle="~ Tissue + State + Dataset",
       x="log2 Fold Change (Foreskin / Breast)",
       y="-log10(adj. P-Value)") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey40") +
  geom_vline(xintercept=c(-1, 1), linetype="dashed", color="grey40") +
  theme(text=element_text(size=13), legend.position="none")

ggsave(file.path(out_dir, "Volcano_FsMC_vs_BsMC.png"), plot=p_volcano, width=9, height=7)

###########################################################################
# 10. KEY GENE CHECK
###########################################################################
message("\n=== KEY GENE RESULTS: Tissue (FsMC vs BsMC) ===")
key_genes <- c("TPSD1", "CXCL8", "RPS4Y1", "XIST", "CCL1", "KIT",
               "TPSAB1", "TPSB2", "TPSG1", "HDC", "FCER1A", "FCER1G",
               "CMA1", "CPA3", "GPX1")
for (g in key_genes) {
  if (g %in% rownames(tt_tissue)) {
    row <- tt_tissue[g, ]
    sig <- ifelse(row$adj.P.Val < 0.05, " ***", "")
    message(sprintf("  %s: logFC=%.2f, adj.P=%.2e%s", g, row$logFC, row$adj.P.Val, sig))
  }
}

message("\n=== KEY GENE RESULTS: Activation (Stimulated vs Native) ===")
act_genes <- c("CCL1", "CXCL8", "RILPL2", "ACSL1", "ATP13A3", "CLIC4")
for (g in act_genes) {
  if (g %in% rownames(tt_state)) {
    row <- tt_state[g, ]
    sig <- ifelse(row$adj.P.Val < 0.05, " ***", "")
    message(sprintf("  %s: logFC=%.2f, adj.P=%.2e%s", g, row$logFC, row$adj.P.Val, sig))
  }
}

message("\nDone! All results saved to: ", out_dir)
