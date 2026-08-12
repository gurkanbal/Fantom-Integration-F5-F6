##############################################################################
# run_extended_analysis.R
# Heatmaps, Boxplots, and Interaction Analysis for BsMC vs FsMC
##############################################################################

library(data.table)
library(edgeR)
library(limma)
library(ggplot2)
library(pheatmap)

setwd("c:/Users/guerkan.bal/OneDrive - Charité - Universitätsmedizin Berlin/DEG_Fantom6")
out_dir <- "data/fantom5_hg38/limma_results"

# 1. Load DE tables and metadata
meta <- read.delim(file.path(out_dir, "limma_sampleinfo.tsv"))
rownames(meta) <- meta$SampleID

tt_tissue <- read.delim(file.path(out_dir, "DE_BsMC_vs_FsMC.tsv"))
tt_state <- read.delim(file.path(out_dir, "DE_Stimulated_vs_Native.tsv"))

# Load normalized expression matrix
# Run voom again to retrieve matrix E
f6_all <- as.data.table(read.delim("20220504/BR1256789_Native_Stimulated_gene_counts.txt", check.names=TRUE, stringsAsFactors=FALSE))
setnames(f6_all, names(f6_all)[1], "geneID")
nc_cols <- grep("aso_NC", names(f6_all), value=TRUE)
f6_nc <- f6_all[, c("geneID", nc_cols), with=FALSE]

f6_anno <- as.data.table(read.delim("20220504/F6_CAT.gene.info_shorted.tsv", stringsAsFactors=FALSE))
f6_merged <- merge(f6_nc, f6_anno[, .(geneID, HGNC_symbol)], by="geneID", all.x=TRUE)
f6_annotated <- f6_merged[HGNC_symbol != "__na" & !is.na(HGNC_symbol) & HGNC_symbol != ""]
f6_gene_cols <- setdiff(names(f6_annotated), c("geneID", "HGNC_symbol"))
f6_agg <- f6_annotated[, lapply(.SD, sum, na.rm=TRUE), by=HGNC_symbol, .SDcols=f6_gene_cols]

f6_sample_cols <- setdiff(names(f6_agg), "HGNC_symbol")
base_names <- sub("\\.[0-9]+$", "", f6_sample_cols)
unique_bases <- unique(base_names)

f6_summed <- data.table(HGNC_symbol = f6_agg$HGNC_symbol)
for (bn in unique_bases) {
  cols_to_sum <- f6_sample_cols[base_names == bn]
  if (length(cols_to_sum) == 1) f6_summed[[bn]] <- f6_agg[[cols_to_sum]]
  else f6_summed[[bn]] <- rowSums(as.matrix(f6_agg[, cols_to_sum, with=FALSE]))
}

f5_counts <- fread("data/fantom5_hg38/fantom5_hg38_mast_cell_gene_counts.tsv")
setnames(f5_counts, names(f5_counts)[1], "HGNC_symbol")
combined <- merge(f6_summed, f5_counts, by="HGNC_symbol", all=FALSE)

counts_mat <- as.matrix(combined[, -1, with=FALSE])
rownames(counts_mat) <- combined$HGNC_symbol
counts_mat <- counts_mat[, meta$SampleID]

meta$Tissue <- factor(meta$Tissue, levels=c("Foreskin", "Breast"))
meta$State <- factor(meta$State, levels=c("Native", "Stimulated"))
meta$Dataset <- factor(meta$Dataset, levels=c("FANTOM6", "FANTOM5"))

dge <- DGEList(counts=counts_mat, samples=meta)
keep <- filterByExpr(dge, group=dge$samples$Tissue)
dge <- dge[keep, , keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge, method="TMM")

design <- model.matrix(~ Tissue + State + Dataset, data=meta)
v <- voom(dge, design, plot=FALSE)

expr_corrected <- removeBatchEffect(v$E, batch=meta$Dataset, design=model.matrix(~ Tissue + State, data=meta))

# -----------------------------------------------------------------------------
# 2. HEATMAP: Top 30 Tissue DE Genes
# -----------------------------------------------------------------------------
top_tissue_genes <- head(rownames(tt_tissue[order(tt_tissue$adj.P.Val), ]), 30)
mat_tissue <- expr_corrected[top_tissue_genes, ]
mat_tissue_z <- t(scale(t(mat_tissue)))

annotation_col <- meta[, c("Tissue", "State")]
annotation_colors <- list(
  Tissue = c("Breast"="#E74C3C", "Foreskin"="#3498DB"),
  State = c("Native"="#2ECC71", "Stimulated"="#9B59B6")
)

png(file.path(out_dir, "Heatmap_Top30_Tissue_DE.png"), width=900, height=800, res=100)
pheatmap(mat_tissue_z,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         main = "Top 30 Tissue DE Genes (z-score logCPM)",
         show_colnames = TRUE,
         fontsize_row = 10,
         clustering_distance_cols = "euclidean")
dev.off()

# -----------------------------------------------------------------------------
# 3. HEATMAP: Top 40 Activation DE Genes
# -----------------------------------------------------------------------------
top_state_genes <- head(rownames(tt_state[order(tt_state$adj.P.Val), ]), 40)
mat_state <- expr_corrected[top_state_genes, ]
mat_state_z <- t(scale(t(mat_state)))

png(file.path(out_dir, "Heatmap_Top40_Activation_DE.png"), width=900, height=900, res=100)
pheatmap(mat_state_z,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         main = "Top 40 Activation DE Genes (z-score logCPM)",
         show_colnames = TRUE,
         fontsize_row = 9,
         clustering_distance_cols = "euclidean")
dev.off()

# -----------------------------------------------------------------------------
# 4. KEY GENE BOXPLOTS PANEL
# -----------------------------------------------------------------------------
key_panel_genes <- c("TPSD1", "CXCL8", "CCL1", "RPS4Y1", "EIF1AY", "KIT", "CPA3", "TPSAB1", "TPSB2")
present_genes <- intersect(key_panel_genes, rownames(expr_corrected))

df_box_list <- list()
for (g in present_genes) {
  tmp <- data.frame(
    SampleID = colnames(expr_corrected),
    Expression = expr_corrected[g, ],
    Gene = g,
    Tissue = meta$Tissue,
    State = meta$State,
    stringsAsFactors = FALSE
  )
  df_box_list[[g]] <- tmp
}
df_box <- rbindlist(df_box_list)

p_box <- ggplot(df_box, aes(x=Tissue, y=Expression, fill=State)) +
  geom_boxplot(outlier.shape=NA, alpha=0.7) +
  geom_jitter(position=position_jitterdodge(jitter.width=0.2), size=1.5, alpha=0.8) +
  facet_wrap(~Gene, scales="free_y") +
  theme_bw(base_size=12) +
  scale_fill_manual(values=c("Native"="#2ECC71", "Stimulated"="#9B59B6")) +
  labs(title="Expression of Key Mast Cell & Tissue Genes",
       subtitle="Batch-corrected log2 CPM across N=21 samples",
       y="log2 CPM") +
  theme(legend.position="bottom")

ggsave(file.path(out_dir, "Boxplots_KeyGenes_Panel.png"), plot=p_box, width=10, height=8, dpi=300)

# -----------------------------------------------------------------------------
# 5. INTERACTION ANALYSIS (~ Tissue * State + Dataset)
# -----------------------------------------------------------------------------
design_int <- model.matrix(~ Tissue * State + Dataset, data=meta)
v_int <- voom(dge, design_int, plot=FALSE)
corfit_int <- duplicateCorrelation(v_int, design_int, block=meta$Donor)
fit_int <- lmFit(v_int, design_int, block=meta$Donor, correlation=corfit_int$consensus.correlation)
fit_int <- eBayes(fit_int)

tt_int <- topTable(fit_int, coef="TissueBreast:StateStimulated", number=Inf)
write.table(tt_int, file=file.path(out_dir, "DE_Interaction_Tissue_x_State.tsv"),
            sep="\t", quote=FALSE)

message("Extended analysis completed successfully!")
