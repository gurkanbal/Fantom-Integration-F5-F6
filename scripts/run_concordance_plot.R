##############################################################################
# run_concordance_plot.R
# Concordance scatter plot: FsMC vs BsMC vs BsMC vs FsMC
# Concordance scatter plot: DREAM vs LIMMA
##############################################################################

library(data.table)
library(ggplot2)

setwd("c:/Users/guerkan.bal/OneDrive - Charité - Universitätsmedizin Berlin/DEG_Fantom6")

out_dir_bvf <- "data/fantom5_hg38/limma_results_BvF"
out_dir_fvb <- "data/fantom5_hg38/limma_results"
dir.create(out_dir_bvf, showWarnings=FALSE, recursive=TRUE)

# 1. Concordance FvB vs BvF
old_file <- file.path(out_dir_fvb, "DE_FsMC_vs_BsMC.tsv")
new_file <- file.path(out_dir_bvf, "DE_BsMC_vs_FsMC.tsv")

if (file.exists(old_file) && file.exists(new_file)) {
    old_df <- read.delim(old_file)
    new_df <- read.delim(new_file)
    
    old_df$Gene <- rownames(old_df)
    new_df$Gene <- rownames(new_df)
    
    df <- merge(old_df[, c("Gene", "logFC", "adj.P.Val")], 
                new_df[, c("Gene", "logFC", "adj.P.Val")], 
                by="Gene", suffixes=c(".old", ".new"))
                
    # Classify significance
    df$sig_old <- df$adj.P.Val.old < 0.05
    df$sig_new <- df$adj.P.Val.new < 0.05
    
    df$category <- "NS"
    df$category[df$sig_old | df$sig_new] <- "Sig in ONE"
    df$category[df$sig_old & df$sig_new] <- "Sig in BOTH"
    
    # Calculate stats
    r_val <- cor(df$logFC.old, df$logFC.new, method="pearson")
    median_delta <- median(abs(df$logFC.old + df$logFC.new)) # Should be near 0
    
    key_genes <- c("TPSD1", "RPS4Y1", "EIF1AY", "XIST", "GPX1", "CLDN23", "C2", "MPDZ", 
                   "KIT", "CPA3", "TPSAB1", "FCER1G", "HDC", "DDX3Y", "KDM5D", "UTY", 
                   "TPSB2", "CCL1", "CXCL8")
    
    p_main <- ggplot(df, aes(x=logFC.old, y=logFC.new, color=category)) +
        geom_point(alpha=0.6, size=1.5) +
        scale_color_manual(values=c("NS"="grey80", "Sig in ONE"="orange", "Sig in BOTH"="blue")) +
        geom_abline(intercept=0, slope=-1, linetype="dashed", color="black", size=1) +
        geom_text(data=df[df$Gene %in% key_genes, ], aes(label=Gene), color="black", size=3.5, vjust=-0.8) +
        theme_bw() +
        labs(title="Concordance: FsMC-vs-BsMC vs BsMC-vs-FsMC",
             subtitle=sprintf("Pearson r = %.3f | Median |Δ log2FC| = %.3f", r_val, median_delta),
             x="log2FC (Foreskin / Breast) [Old]",
             y="log2FC (Breast / Foreskin) [New]") +
        theme(text=element_text(size=14))
        
    p_inset <- ggplot(df, aes(x=logFC.old, y=logFC.new, color=category)) +
        geom_point(alpha=0.6, size=1) +
        scale_color_manual(values=c("NS"="grey80", "Sig in ONE"="orange", "Sig in BOTH"="blue")) +
        geom_abline(intercept=0, slope=-1, linetype="dashed", color="black") +
        coord_cartesian(xlim=c(-2, 2), ylim=c(-2, 2)) +
        theme_bw() +
        theme(legend.position="none", 
              axis.title=element_blank(),
              panel.background = element_rect(fill = "white", colour = "black"))

    p_final <- p_main + annotation_custom(ggplotGrob(p_inset), xmin=min(df$logFC.old), xmax=0, ymin=0, ymax=max(df$logFC.new))
    
    ggsave(file.path(out_dir_bvf, "Concordance_FvB_vs_BvF.png"), plot=p_final, width=10, height=10)
    ggsave(file.path(out_dir_fvb, "Concordance_FvB_vs_BvF.png"), plot=p_final, width=10, height=10)
}

# 2. Concordance DREAM vs LIMMA BvF
dream_file <- "data/fantom5_hg38/dream_results/dream_DE_FsMC_vs_BsMC.tsv"

if (file.exists(dream_file) && file.exists(new_file)) {
    dream_df <- read.delim(dream_file)
    new_df <- read.delim(new_file)
    
    dream_df$Gene <- rownames(dream_df)
    new_df$Gene <- rownames(new_df)
    
    df <- merge(dream_df[, c("Gene", "logFC", "adj.P.Val")], 
                new_df[, c("Gene", "logFC", "adj.P.Val")], 
                by="Gene", suffixes=c(".dream", ".limma"))
                
    # Classify significance
    df$sig_dream <- df$adj.P.Val.dream < 0.05
    df$sig_limma <- df$adj.P.Val.limma < 0.05
    
    df$category <- "NS"
    df$category[df$sig_dream | df$sig_limma] <- "Sig in ONE"
    df$category[df$sig_dream & df$sig_limma] <- "Sig in BOTH"
    
    # Both are now Breast vs Foreskin, so they should be correlated positively
    r_val <- cor(df$logFC.dream, df$logFC.limma, method="pearson", use="complete.obs")
    median_delta <- median(abs(df$logFC.dream - df$logFC.limma), na.rm=TRUE)
    
    p_main <- ggplot(df, aes(x=logFC.dream, y=logFC.limma, color=category)) +
        geom_point(alpha=0.6, size=1.5) +
        scale_color_manual(values=c("NS"="grey80", "Sig in ONE"="orange", "Sig in BOTH"="blue")) +
        geom_abline(intercept=0, slope=1, linetype="dashed", color="black", size=1) +
        geom_text(data=df[df$Gene %in% key_genes, ], aes(label=Gene), color="black", size=3.5, vjust=-0.8) +
        theme_bw() +
        labs(title="Concordance: Manuscript (DREAM) vs Integration (limma-voom)",
             subtitle=sprintf("Both contrasts: Breast vs Foreskin | Pearson r = %.3f | Median |Δ log2FC| = %.3f", r_val, median_delta),
             x="log2FC (DREAM Mixed Model)",
             y="log2FC (limma-voom)") +
        theme(text=element_text(size=14))
        
    p_inset <- ggplot(df, aes(x=logFC.dream, y=logFC.limma, color=category)) +
        geom_point(alpha=0.6, size=1) +
        scale_color_manual(values=c("NS"="grey80", "Sig in ONE"="orange", "Sig in BOTH"="blue")) +
        geom_abline(intercept=0, slope=1, linetype="dashed", color="black") +
        coord_cartesian(xlim=c(-2, 2), ylim=c(-2, 2)) +
        theme_bw() +
        theme(legend.position="none", 
              axis.title=element_blank(),
              panel.background = element_rect(fill = "white", colour = "black"))

    p_final <- p_main + annotation_custom(ggplotGrob(p_inset), xmin=min(df$logFC.dream, na.rm=T), xmax=0, ymin=0, ymax=max(df$logFC.limma, na.rm=T))
    
    ggsave(file.path(out_dir_bvf, "Concordance_DREAM_vs_Limma_BvF.png"), plot=p_final, width=10, height=10)
    ggsave(file.path(out_dir_fvb, "Concordance_DREAM_vs_Limma_BvF.png"), plot=p_final, width=10, height=10)
}

message("Concordance plots generated successfully!")