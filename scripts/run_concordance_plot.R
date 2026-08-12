##############################################################################
# run_concordance_plot.R
# Publication Concordance Plots: Tissue Effect & State Effect
# F5-F6 Integration (limma-voom N=21) vs Manuscript DREAM (Excel Supplement)
##############################################################################

library(readxl)
library(ggplot2)
library(ggrepel)
library(data.table)

setwd("c:/Users/guerkan.bal/OneDrive - Charité - Universitätsmedizin Berlin/DEG_Fantom6")
out_dir <- "data/fantom5_hg38/limma_results"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

excel_file <- "C:/Users/guerkan.bal/OneDrive - Charité - Universitätsmedizin Berlin/Allergy_revision_antigravitiy/Supplement Tables.xlsx"

if (!file.exists(excel_file)) {
  stop("Excel file not found at: ", excel_file)
}

# =============================================================================
# 1. TISSUE EFFECT CONCORDANCE (Breast vs Foreskin / BsMC vs FsMC)
# =============================================================================
message("Generating Tissue Effect Concordance plot...")
dream_tissue <- as.data.frame(read_excel(excel_file, sheet="S3a DGE_TissueBreast_gene_Dream"))
dream_tissue$Gene <- dream_tissue$HGNC_symbol

limma_tissue <- read.delim(file.path(out_dir, "DE_BsMC_vs_FsMC.tsv"))
limma_tissue$Gene <- rownames(limma_tissue)

df_tissue_all <- merge(dream_tissue[, c("Gene", "logFC", "adj.P.Val")], 
                       limma_tissue[, c("Gene", "logFC", "adj.P.Val")], 
                       by="Gene", suffixes=c(".dream", ".limma"))

genes_23 <- c("RPS4Y1", "EIF1AY", "DDX3Y", "KDM5D", "MPDZ", "UTY", "TTTY14", "ZFY", "PRKY", "TPSD1", "C2", 
              "GPX1", "CLDN23", "XIST", "HDC", "TPSB2", "KIT", "CPA3", "FCER1A", "TPSAB1", "EIF1AX", "FCER1G", "TPSG1")

df_t <- df_tissue_all[df_tissue_all$Gene %in% genes_23, ]

df_t$sig_dream <- df_t$adj.P.Val.dream < 0.05
df_t$sig_limma <- df_t$adj.P.Val.limma < 0.05

df_t$category <- "not significant in either"
df_t$category[df_t$sig_dream | df_t$sig_limma] <- "significant in one only"
df_t$category[df_t$sig_dream & df_t$sig_limma] <- "significant in both"
df_t$category <- factor(df_t$category, levels=c("significant in both", "significant in one only", "not significant in either"))

median_delta_t <- median(abs(df_t$logFC.limma - df_t$logFC.dream), na.rm=TRUE)

title_t <- "Effect sizes agree closely; significance calls do not always"
sub_t <- sprintf("median |\u0394 log\u2082FC| = %.2f over 23 genes \u00b7 dashed line = exact agreement", median_delta_t)

col_sig_both <- "#2A76D2" # blue
col_sig_one <- "#ED713A"  # orange
col_sig_none <- "#8D8D8D" # grey

box_xmin <- -1.7; box_xmax <- 1.7
box_ymin <- -1.7; box_ymax <- 1.7

df_t$main_label <- df_t$Gene
in_box_t <- df_t$logFC.dream >= box_xmin & df_t$logFC.dream <= box_xmax & 
            df_t$logFC.limma >= box_ymin & df_t$logFC.limma <= box_ymax
df_t$main_label[in_box_t] <- ""

p_main_t <- ggplot(df_t, aes(x=logFC.dream, y=logFC.limma, color=category)) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="#8c8c8c", linewidth=0.8) +
  annotate("rect", xmin=box_xmin, xmax=box_xmax, ymin=box_ymin, ymax=box_ymax, 
           fill=NA, color="#8c8c8c", linetype="dashed", linewidth=0.6) +
  annotate("text", x=-7.5, y=-2.5, label="core MC genes + tryptases\n\u2192 see inset", 
           color="#595959", size=4, hjust=0) +
  geom_point(size=4) +
  scale_color_manual(values=c("significant in both"=col_sig_both, 
                              "significant in one only"=col_sig_one, 
                              "not significant in either"=col_sig_none)) +
  geom_text_repel(aes(label=main_label), color="#595959", size=4, 
                  point.padding=0.3, box.padding=0.4, max.overlaps=Inf,
                  segment.color="transparent", seed=42, min.segment.length=0) +
  theme_minimal(base_size=14) +
  labs(title=title_t,
       subtitle=sub_t,
       x="Manuscript (DREAM)  log\u2082FC  breast vs foreskin",
       y="F5-F6 integration (limma-voom)  log\u2082FC  breast vs foreskin") +
  theme(
    plot.title = element_text(face="plain", size=16, margin=margin(b=5)),
    plot.subtitle = element_text(size=14, margin=margin(b=20)),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size=12, color="#595959"),
    axis.title = element_text(color="#595959"),
    axis.text = element_text(color="#8c8c8c"),
    panel.grid.major = element_line(color="#f0f0f0"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill="white", color=NA),
    plot.background = element_rect(fill="white", color=NA)
  ) +
  scale_x_continuous(limits=c(-11, 11), breaks=seq(-10, 10, by=2.5)) +
  scale_y_continuous(limits=c(-11, 11), breaks=seq(-10, 10, by=2.5))

df_inset_t <- df_t[in_box_t, ]
p_inset_t <- ggplot(df_inset_t, aes(x=logFC.dream, y=logFC.limma, color=category)) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="#8c8c8c", linewidth=0.6) +
  geom_point(size=3) +
  scale_color_manual(values=c("significant in both"=col_sig_both, 
                              "significant in one only"=col_sig_one, 
                              "not significant in either"=col_sig_none)) +
  geom_text_repel(aes(label=Gene), color="#595959", size=3, max.overlaps=Inf,
                  point.padding=0.2, box.padding=0.3, segment.color="transparent", seed=42) +
  theme_minimal(base_size=10) +
  coord_cartesian(xlim=c(box_xmin, box_xmax), ylim=c(box_ymin, box_ymax)) +
  theme(
    legend.position = "none", axis.title = element_blank(),
    axis.text = element_text(color="#8c8c8c"),
    panel.grid.major = element_line(color="#f0f0f0"), panel.grid.minor = element_blank(),
    panel.border = element_rect(color="#f0f0f0", fill=NA, linewidth=1),
    panel.background = element_rect(fill="#fafafa", color=NA)
  )

p_final_t <- p_main_t + annotation_custom(ggplotGrob(p_inset_t), xmin=3, xmax=11, ymin=2, ymax=10)
ggsave(file.path(out_dir, "Concordance_Tissue_Excel_vs_Limma.png"), plot=p_final_t, width=9, height=9, dpi=300)

# =============================================================================
# 2. STATE EFFECT CONCORDANCE (Stimulated vs Native)
# =============================================================================
message("Generating State Effect Concordance plot...")
dream_state <- as.data.frame(read_excel(excel_file, sheet="S2a_Dream_DEG_Sti.vs.base."))
dream_state$Gene <- dream_state$HGNC_symbol

limma_state <- read.delim(file.path(out_dir, "DE_Stimulated_vs_Native.tsv"))
limma_state$Gene <- rownames(limma_state)

df_s <- merge(dream_state[, c("Gene", "logFC", "adj.P.Val")], 
              limma_state[, c("Gene", "logFC", "adj.P.Val")], 
              by="Gene", suffixes=c(".dream", ".limma"))

df_s$sig_dream <- !is.na(df_s$adj.P.Val.dream) & df_s$adj.P.Val.dream < 0.05
df_s$sig_limma <- !is.na(df_s$adj.P.Val.limma) & df_s$adj.P.Val.limma < 0.05

df_s$category <- "not significant in either"
df_s$category[df_s$sig_dream | df_s$sig_limma] <- "significant in one only"
df_s$category[df_s$sig_dream & df_s$sig_limma] <- "significant in both"
df_s$category <- factor(df_s$category, levels=c("significant in both", "significant in one only", "not significant in either"))

r_val_s <- cor(df_s$logFC.dream, df_s$logFC.limma, method="pearson", use="complete.obs")
median_delta_s <- median(abs(df_s$logFC.limma - df_s$logFC.dream), na.rm=TRUE)

top_act_genes <- c("CCL1", "CXCL8", "RILPL2", "ACSL1", "ATP13A3", "CLIC4", "EGR1", "FOS", "JUN", "TNF", "IL6")

p_state <- ggplot(df_s, aes(x=logFC.dream, y=logFC.limma, color=category)) +
  geom_point(alpha=0.4, size=1.5) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="#8c8c8c", linewidth=0.8) +
  scale_color_manual(values=c("significant in both"=col_sig_both, 
                              "significant in one only"=col_sig_one, 
                              "not significant in either"=col_sig_none)) +
  geom_text_repel(data=df_s[df_s$Gene %in% top_act_genes, ], aes(label=Gene),
                  color="black", size=4, max.overlaps=Inf, seed=42) +
  theme_minimal(base_size=14) +
  labs(title="Activation Concordance: Manuscript (DREAM) vs Integration (limma-voom)",
       subtitle=sprintf("Pearson r = %.3f \u00b7 median |\u0394 log\u2082FC| = %.2f across %d genes", r_val_s, median_delta_s, nrow(df_s)),
       x="Manuscript (DREAM)  log\u2082FC  stimulated vs native",
       y="F5-F6 integration (limma-voom)  log\u2082FC  stimulated vs native") +
  theme(
    plot.title = element_text(face="plain", size=16, margin=margin(b=5)),
    plot.subtitle = element_text(size=14, margin=margin(b=20)),
    legend.position = "bottom", legend.title = element_blank(),
    legend.text = element_text(size=12, color="#595959"),
    axis.title = element_text(color="#595959"), axis.text = element_text(color="#8c8c8c"),
    panel.grid.major = element_line(color="#f0f0f0"), panel.grid.minor = element_blank(),
    panel.background = element_rect(fill="white", color=NA),
    plot.background = element_rect(fill="white", color=NA)
  )

ggsave(file.path(out_dir, "Concordance_State_Excel_vs_Limma.png"), plot=p_state, width=9, height=9, dpi=300)

message("Concordance plots generated successfully!")
