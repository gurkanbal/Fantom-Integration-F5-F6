##############################################################################
# run_concordance_plot.R
# Concordance scatter plot: Exported Excel DREAM Data vs limma-voom Integration
##############################################################################

library(readxl)
library(ggplot2)
library(ggrepel)
library(data.table)

setwd("c:/Users/guerkan.bal/OneDrive - Charité - Universitätsmedizin Berlin/DEG_Fantom6")

out_dir_bvf <- "data/fantom5_hg38/limma_results_BvF"
out_dir_fvb <- "data/fantom5_hg38/limma_results"
dir.create(out_dir_bvf, showWarnings=FALSE, recursive=TRUE)

# 1. Load Excel exported DREAM data
excel_file <- "C:/Users/guerkan.bal/OneDrive - Charité - Universitätsmedizin Berlin/Allergy_revision_antigravitiy/Supplement Tables.xlsx"

if (!file.exists(excel_file)) {
  stop("Excel file not found at: ", excel_file)
}

dream_df <- as.data.frame(read_excel(excel_file, sheet="S3a DGE_TissueBreast_gene_Dream"))
dream_df$Gene <- dream_df$HGNC_symbol

# 2. Load limma-voom integration results
file_bvf <- file.path(out_dir_bvf, "DE_BsMC_vs_FsMC.tsv")
file_fvb <- file.path(out_dir_fvb, "DE_FsMC_vs_BsMC.tsv")

key_genes <- c("TPSD1", "RPS4Y1", "EIF1AY", "XIST", "GPX1", "CLDN23", "C2", "MPDZ", 
               "KIT", "CPA3", "TPSAB1", "FCER1G", "HDC", "DDX3Y", "KDM5D", "UTY", 
               "TPSB2", "CCL1", "CXCL8")

# -----------------------------------------------------------------------------
# PLOT 1: BsMC vs FsMC (Both Breast vs Foreskin -> positive slope y = x)
# -----------------------------------------------------------------------------
if (file.exists(file_bvf)) {
  limma_bvf <- read.delim(file_bvf)
  limma_bvf$Gene <- rownames(limma_bvf)
  
  df_bvf <- merge(dream_df[, c("Gene", "logFC", "adj.P.Val")],
                  limma_bvf[, c("Gene", "logFC", "adj.P.Val")],
                  by="Gene", suffixes=c(".excel", ".limma"))
                  
  df_bvf$sig_excel <- !is.na(df_bvf$adj.P.Val.excel) & df_bvf$adj.P.Val.excel < 0.05
  df_bvf$sig_limma <- !is.na(df_bvf$adj.P.Val.limma) & df_bvf$adj.P.Val.limma < 0.05
  
  df_bvf$category <- "NS in both"
  df_bvf$category[df_bvf$sig_excel | df_bvf$sig_limma] <- "Sig in ONE"
  df_bvf$category[df_bvf$sig_excel & df_bvf$sig_limma] <- "Sig in BOTH"
  df_bvf$category <- factor(df_bvf$category, levels=c("Sig in BOTH", "Sig in ONE", "NS in both"))
  
  r_val <- cor(df_bvf$logFC.excel, df_bvf$logFC.limma, method="pearson", use="complete.obs")
  median_delta <- median(abs(df_bvf$logFC.excel - df_bvf$logFC.limma), na.rm=TRUE)
  
  p_main_bvf <- ggplot(df_bvf, aes(x=logFC.excel, y=logFC.limma, color=category)) +
    geom_point(alpha=0.5, size=1.5) +
    scale_color_manual(values=c("Sig in BOTH"="#2A76D2", "Sig in ONE"="#ED713A", "NS in both"="#A0A0A0")) +
    geom_abline(intercept=0, slope=1, linetype="dashed", color="black", linewidth=0.8) +
    geom_text_repel(data=df_bvf[df_bvf$Gene %in% key_genes, ], aes(label=Gene),
                    color="black", size=3.8, max.overlaps=Inf, seed=42) +
    theme_bw(base_size=14) +
    labs(title="Concordance: Excel Exported DREAM vs Integration (limma-voom)",
         subtitle=sprintf("Contrast: Breast vs Foreskin (BsMC vs FsMC) | Pearson r = %.3f | Median |Δ log2FC| = %.3f", r_val, median_delta),
         x="log2FC (DREAM - Excel S3a: Breast / Foreskin)",
         y="log2FC (limma-voom: Breast / Foreskin)",
         color="Significance") +
    theme(legend.position="bottom")
    
  p_inset_bvf <- ggplot(df_bvf, aes(x=logFC.excel, y=logFC.limma, color=category)) +
    geom_point(alpha=0.6, size=1.2) +
    scale_color_manual(values=c("Sig in BOTH"="#2A76D2", "Sig in ONE"="#ED713A", "NS in both"="#A0A0A0")) +
    geom_abline(intercept=0, slope=1, linetype="dashed", color="black") +
    coord_cartesian(xlim=c(-2, 2), ylim=c(-2, 2)) +
    theme_bw(base_size=10) +
    theme(legend.position="none", axis.title=element_blank(), panel.background=element_rect(fill="white"))
    
  p_final_bvf <- p_main_bvf + annotation_custom(ggplotGrob(p_inset_bvf), xmin=-11, xmax=-2, ymin=2, ymax=11)
  
  ggsave(file.path(out_dir_bvf, "Concordance_Excel_vs_Limma_BvF.png"), plot=p_final_bvf, width=9, height=9, dpi=300)
  ggsave(file.path(out_dir_fvb, "Concordance_Excel_vs_Limma_BvF.png"), plot=p_final_bvf, width=9, height=9, dpi=300)
  message("Saved Concordance_Excel_vs_Limma_BvF.png")
}

# -----------------------------------------------------------------------------
# PLOT 2: FsMC vs BsMC (DREAM Breast/Foreskin vs Limma Foreskin/Breast -> slope y = -x)
# -----------------------------------------------------------------------------
if (file.exists(file_fvb)) {
  limma_fvb <- read.delim(file_fvb)
  limma_fvb$Gene <- rownames(limma_fvb)
  
  df_fvb <- merge(dream_df[, c("Gene", "logFC", "adj.P.Val")],
                  limma_fvb[, c("Gene", "logFC", "adj.P.Val")],
                  by="Gene", suffixes=c(".excel", ".limma"))
                  
  df_fvb$sig_excel <- !is.na(df_fvb$adj.P.Val.excel) & df_fvb$adj.P.Val.excel < 0.05
  df_fvb$sig_limma <- !is.na(df_fvb$adj.P.Val.limma) & df_fvb$adj.P.Val.limma < 0.05
  
  df_fvb$category <- "NS in both"
  df_fvb$category[df_fvb$sig_excel | df_fvb$sig_limma] <- "Sig in ONE"
  df_fvb$category[df_fvb$sig_excel & df_fvb$sig_limma] <- "Sig in BOTH"
  df_fvb$category <- factor(df_fvb$category, levels=c("Sig in BOTH", "Sig in ONE", "NS in both"))
  
  r_val <- cor(df_fvb$logFC.excel, df_fvb$logFC.limma, method="pearson", use="complete.obs")
  median_delta <- median(abs(df_fvb$logFC.excel + df_fvb$logFC.limma), na.rm=TRUE)
  
  p_main_fvb <- ggplot(df_fvb, aes(x=logFC.excel, y=logFC.limma, color=category)) +
    geom_point(alpha=0.5, size=1.5) +
    scale_color_manual(values=c("Sig in BOTH"="#2A76D2", "Sig in ONE"="#ED713A", "NS in both"="#A0A0A0")) +
    geom_abline(intercept=0, slope=-1, linetype="dashed", color="black", linewidth=0.8) +
    geom_text_repel(data=df_fvb[df_fvb$Gene %in% key_genes, ], aes(label=Gene),
                    color="black", size=3.8, max.overlaps=Inf, seed=42) +
    theme_bw(base_size=14) +
    labs(title="Concordance: Excel Exported DREAM vs Integration (limma-voom)",
         subtitle=sprintf("DREAM (Breast/Foreskin) vs limma-voom (Foreskin/Breast) | Pearson r = %.3f | Median |Δ log2FC| = %.3f", r_val, median_delta),
         x="log2FC (DREAM - Excel S3a: Breast / Foreskin)",
         y="log2FC (limma-voom: Foreskin / Breast)",
         color="Significance") +
    theme(legend.position="bottom")
    
  p_inset_fvb <- ggplot(df_fvb, aes(x=logFC.excel, y=logFC.limma, color=category)) +
    geom_point(alpha=0.6, size=1.2) +
    scale_color_manual(values=c("Sig in BOTH"="#2A76D2", "Sig in ONE"="#ED713A", "NS in both"="#A0A0A0")) +
    geom_abline(intercept=0, slope=-1, linetype="dashed", color="black") +
    coord_cartesian(xlim=c(-2, 2), ylim=c(-2, 2)) +
    theme_bw(base_size=10) +
    theme(legend.position="none", axis.title=element_blank(), panel.background=element_rect(fill="white"))
    
  p_final_fvb <- p_main_fvb + annotation_custom(ggplotGrob(p_inset_fvb), xmin=-11, xmax=-2, ymin=-11, ymax=-2)
  
  ggsave(file.path(out_dir_fvb, "Concordance_Excel_vs_Limma_FvB.png"), plot=p_final_fvb, width=9, height=9, dpi=300)
  ggsave(file.path(out_dir_bvf, "Concordance_Excel_vs_Limma_FvB.png"), plot=p_final_fvb, width=9, height=9, dpi=300)
  message("Saved Concordance_Excel_vs_Limma_FvB.png")
}

message("Excel DREAM vs Limma concordance plots generated successfully!")
