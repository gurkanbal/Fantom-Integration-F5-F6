##############################################################################
# run_extended_analysis.R
# Extended analyses: Boxplots, Heatmap, Interaction model
##############################################################################

library(data.table)
library(edgeR)
library(limma)
library(ggplot2)
library(pheatmap)

if (!require("pheatmap", quietly=TRUE)) install.packages("pheatmap", repos="https://cran.rstudio.com/")
library(pheatmap)

setwd("c:/Users/guerkan.bal/OneDrive - Charité - Universitätsmedizin Berlin/DEG_Fantom6")
out_dir <- "data/fantom5_hg38/limma_results"

###########################################################################
# 1. RELOAD DATA (same pipeline as run_limma_integration.R)
###########################################################################
message("Loading and processing data...")
f6_all <- as.data.table(read.delim("20220504/BR1256789_Native_Stimulated_gene_counts.txt",
                                    check.names=TRUE, stringsAsFactors=FALSE))
setnames(f6_all, names(f6_all)[1], "geneID")
nc_cols <- grep("aso_NC", names(f6_all), value=TRUE)
f6_nc <- f6_all[, c("geneID", nc_cols), with=FALSE]

f6_anno <- as.data.table(read.delim("20220504/F6_CAT.gene.info_shorted.tsv", stringsAsFactors=FALSE))
f6_merged <- merge(f6_nc, f6_anno[, .(geneID, HGNC_symbol)], by="geneID", all.x=TRUE)
f6_annotated <- f6_merged[HGNC_symbol != "__na" & !is.na(HGNC_symbol) & HGNC_symbol != ""]
f6_gene_cols <- setdiff(names(f6_annotated), c("geneID", "HGNC_symbol"))
f6_agg <- f6_annotated[, lapply(.SD, sum, na.rm=TRUE), by=HGNC_symbol, .SDcols=f6_gene_cols]

# Sum tech reps
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

f5_counts <- fread("data/fantom5_hg38/fantom5_hg38_mast_cell_gene_counts.tsv")
setnames(f5_counts, names(f5_counts)[1], "HGNC_symbol")
combined <- merge(f6_summed, f5_counts, by="HGNC_symbol", all=FALSE)

# Metadata
samples <- setdiff(names(combined), "HGNC_symbol")
meta <- data.frame(SampleID = samples, stringsAsFactors = FALSE)
f6_sampleinfo <- fread("20220504/BR1256789_Native_Stimulated_sampleinfo.txt")
tissue_map <- unique(f6_sampleinfo[, .(Rep, Tissue)])

for (i in seq_len(nrow(meta))) {
  sid <- meta$SampleID[i]
  if (grepl("^BR", sid)) {
    meta$Dataset[i] <- "FANTOM6"
    meta$Donor[i] <- sub("^(BR[0-9]+).*", "\\1", sid)
    meta$State[i] <- ifelse(grepl("AER", sid), "Stimulated", "Native")
    meta$Tissue[i] <- tissue_map$Tissue[tissue_map$Rep == meta$Donor[i]][1]
  } else {
    meta$Dataset[i] <- "FANTOM5"
    meta$Tissue[i] <- "Breast"
    if (grepl("stimulated%2c%20donor1", sid)) { meta$Donor[i] <- "F5_D1"; meta$State[i] <- "Stimulated"
    } else if (grepl("donor1", sid)) { meta$Donor[i] <- "F5_D1"; meta$State[i] <- "Native"
    } else if (grepl("donor2", sid)) { meta$Donor[i] <- "F5_D2"; meta$State[i] <- "Native"
    } else if (grepl("donor3", sid)) { meta$Donor[i] <- "F5_D3"; meta$State[i] <- "Native"
    } else if (grepl("donor4", sid)) { meta$Donor[i] <- "F5_D4"; meta$State[i] <- "Native"
    } else if (grepl("expanded%20and%20stimulated.*donor5", sid)) { meta$Donor[i] <- "F5_D5"; meta$State[i] <- "Stimulated"
    } else if (grepl("expanded%20and%20stimulated.*donor8", sid)) { meta$Donor[i] <- "F5_D8"; meta$State[i] <- "Stimulated"
    } else if (grepl("expanded.*donor5", sid)) { meta$Donor[i] <- "F5_D5"; meta$State[i] <- "Native"
    } else if (grepl("expanded.*donor8", sid)) { meta$Donor[i] <- "F5_D8"; meta$State[i] <- "Native"
    }
  }
}

meta$Tissue <- factor(meta$Tissue, levels=c("Breast", "Foreskin"))
meta$State <- factor(meta$State, levels=c("Native", "Stimulated"))
meta$Dataset <- factor(meta$Dataset, levels=c("FANTOM6", "FANTOM5"))
meta$Group <- paste0(meta$Tissue, "_", meta$State)
rownames(meta) <- meta$SampleID

# DGEList
counts_mat <- as.matrix(combined[, -1, with=FALSE])
rownames(counts_mat) <- combined$HGNC_symbol
counts_mat <- counts_mat[, meta$SampleID]
dge <- DGEList(counts=counts_mat, samples=meta)
keep <- filterByExpr(dge, group=dge$samples$Tissue)
dge <- dge[keep, , keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge, method="TMM")

# Voom
design <- model.matrix(~ Tissue + State + Dataset, data=meta)
v <- voom(dge, design, plot=FALSE)

# Batch-corrected logCPM for plotting
expr_corrected <- removeBatchEffect(v$E, batch=meta$Dataset,
                                     design=model.matrix(~ Tissue + State, data=meta))

###########################################################################
# 2. BOXPLOTS OF KEY GENES
###########################################################################
message("Generating boxplots of key genes...")

key_genes <- c("TPSD1", "CXCL8", "CCL1", "RPS4Y1", "XIST", "GPX1",
               "TPSB2", "TPSG1", "TPSAB1", "CPA3", "KIT", "FCER1G")

plot_list <- list()
for (g in key_genes) {
  if (g %in% rownames(expr_corrected)) {
    df <- data.frame(
      Expression = expr_corrected[g, ],
      Tissue = meta$Tissue,
      State = meta$State,
      Dataset = meta$Dataset,
      Group = meta$Group
    )
    
    p <- ggplot(df, aes(x=Tissue, y=Expression, fill=State)) +
      geom_boxplot(alpha=0.7, outlier.shape=NA, width=0.6) +
      geom_jitter(aes(shape=Dataset), width=0.15, size=2.5, alpha=0.8) +
      theme_bw() +
      scale_fill_manual(values=c("Native"="#3498DB", "Stimulated"="#E74C3C")) +
      scale_shape_manual(values=c("FANTOM6"=16, "FANTOM5"=17)) +
      labs(title=g, y="Batch-corrected logCPM") +
      theme(text=element_text(size=13), plot.title=element_text(face="bold", size=16))
    
    plot_list[[g]] <- p
  }
}

# Save individual boxplots
for (g in names(plot_list)) {
  ggsave(file.path(out_dir, paste0("Boxplot_", g, ".png")),
         plot=plot_list[[g]], width=6, height=5)
}

# Combined panel of 12 genes
library(gridExtra)
if (length(plot_list) >= 12) {
  p_combined <- arrangeGrob(grobs=plot_list[1:12], ncol=4, nrow=3)
  ggsave(file.path(out_dir, "Boxplots_KeyGenes_Panel.png"),
         plot=p_combined, width=22, height=14)
}
message("Saved boxplots.")

###########################################################################
# 3. HEATMAP OF TOP DE GENES (TISSUE)
###########################################################################
message("Generating heatmap of top tissue DE genes...")
tt_tissue <- read.delim(file.path(out_dir, "DE_FsMC_vs_BsMC.tsv"), row.names=1)
# Top 30 tissue DE genes by p-value (exclude Y-chrom dominated list for diversity)
top_tissue <- head(tt_tissue[order(tt_tissue$P.Value), ], 30)
top_genes_tissue <- rownames(top_tissue)

# Heatmap data (batch-corrected, z-scored)
hm_data <- expr_corrected[top_genes_tissue, ]
hm_scaled <- t(scale(t(hm_data)))

# Annotation
ann_col <- data.frame(
  Tissue = meta$Tissue,
  State = meta$State,
  Dataset = meta$Dataset,
  row.names = meta$SampleID
)

ann_colors <- list(
  Tissue = c("Breast"="#E74C3C", "Foreskin"="#3498DB"),
  State = c("Native"="#2ECC71", "Stimulated"="#E67E22"),
  Dataset = c("FANTOM6"="#9B59B6", "FANTOM5"="#1ABC9C")
)

png(file.path(out_dir, "Heatmap_Top30_Tissue_DE.png"), width=1000, height=800)
pheatmap(hm_scaled,
         annotation_col = ann_col,
         annotation_colors = ann_colors,
         clustering_method = "ward.D2",
         show_colnames = FALSE,
         fontsize_row = 9,
         main = "Top 30 Tissue DE Genes (FsMC vs BsMC)\nZ-scored, batch-corrected logCPM")
dev.off()
message("Saved tissue heatmap.")

###########################################################################
# 4. HEATMAP OF TOP ACTIVATION GENES
###########################################################################
message("Generating heatmap of top activation DE genes...")
tt_state <- read.delim(file.path(out_dir, "DE_Stimulated_vs_Native.tsv"), row.names=1)
top_state <- head(tt_state[order(tt_state$P.Value), ], 40)
top_genes_state <- rownames(top_state)

hm_data_act <- expr_corrected[top_genes_state, ]
hm_scaled_act <- t(scale(t(hm_data_act)))

png(file.path(out_dir, "Heatmap_Top40_Activation_DE.png"), width=1000, height=900)
pheatmap(hm_scaled_act,
         annotation_col = ann_col,
         annotation_colors = ann_colors,
         clustering_method = "ward.D2",
         show_colnames = FALSE,
         fontsize_row = 9,
         main = "Top 40 Activation DE Genes (Stimulated vs Native)\nZ-scored, batch-corrected logCPM")
dev.off()
message("Saved activation heatmap.")

###########################################################################
# 5. INTERACTION MODEL: Tissue x State
###########################################################################
message("\nRunning interaction model: ~ Tissue * State + Dataset ...")
design_int <- model.matrix(~ Tissue * State + Dataset, data=meta)
colnames(design_int)
v_int <- voom(dge, design_int, plot=FALSE)
fit_int <- lmFit(v_int, design_int)
fit_int <- eBayes(fit_int)

# Interaction term: genes where activation DIFFERS between tissues
tt_interaction <- topTable(fit_int, coef="TissueForeskin:StateStimulated", number=Inf)
write.table(tt_interaction, file=file.path(out_dir, "DE_Interaction_Tissue_x_State.tsv"),
            sep="\t", quote=FALSE)

n_int <- sum(tt_interaction$adj.P.Val < 0.05)
message(paste("Interaction DE genes (adj.P < 0.05):", n_int))

# Top interaction genes
message("\n=== TOP INTERACTION GENES (different activation in FsMC vs BsMC) ===")
top_int <- head(tt_interaction[order(tt_interaction$P.Value), ], 15)
for (i in seq_len(nrow(top_int))) {
  g <- rownames(top_int)[i]
  message(sprintf("  %s: logFC=%.2f, adj.P=%.2e",
                  g, top_int$logFC[i], top_int$adj.P.Val[i]))
}

# Key genes in interaction
message("\n=== KEY GENES IN INTERACTION ===")
key_int <- c("TPSD1", "CXCL8", "CCL1", "RILPL2", "ACSL1", "GPX1")
for (g in key_int) {
  if (g %in% rownames(tt_interaction)) {
    row <- tt_interaction[g, ]
    sig <- ifelse(row$adj.P.Val < 0.05, " ***", "")
    message(sprintf("  %s: interaction logFC=%.2f, adj.P=%.2e%s",
                    g, row$logFC, row$adj.P.Val, sig))
  }
}

message("\nAll extended analyses done! Results in: ", out_dir)
