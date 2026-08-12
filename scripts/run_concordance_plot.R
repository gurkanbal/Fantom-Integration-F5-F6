##############################################################################
# run_concordance_plot.R
# Concordance scatter plot: Manuscript (DREAM) vs Integration (limma-voom)
# Exactly matching the publication figure design
##############################################################################

library(readxl)
library(ggplot2)
library(ggrepel)
library(data.table)

setwd("c:/Users/guerkan.bal/OneDrive - Charité - Universitätsmedizin Berlin/DEG_Fantom6")

out_dir_fvb <- "data/fantom5_hg38/limma_results"
out_dir_bvf <- "data/fantom5_hg38/limma_results_BvF"
dir.create(out_dir_fvb, showWarnings=FALSE, recursive=TRUE)
dir.create(out_dir_bvf, showWarnings=FALSE, recursive=TRUE)

# 1. Read Manuscript (DREAM) data from Excel
excel_file <- "C:/Users/guerkan.bal/OneDrive - Charité - Universitätsmedizin Berlin/Allergy_revision_antigravitiy/Supplement Tables.xlsx"

if (!file.exists(excel_file)) {
  stop("Excel file not found at: ", excel_file)
}

dream_df <- as.data.frame(read_excel(excel_file, sheet="S3a DGE_TissueBreast_gene_Dream"))
dream_df$Gene <- dream_df$HGNC_symbol

# 2. Read Integration (limma-voom) FsMC vs BsMC data
limma_file <- file.path(out_dir_fvb, "DE_FsMC_vs_BsMC.tsv")
limma_df <- read.delim(limma_file)
limma_df$Gene <- rownames(limma_df)

# Merge
df_all <- merge(dream_df[, c("Gene", "logFC", "adj.P.Val")], 
                limma_df[, c("Gene", "logFC", "adj.P.Val")], 
                by="Gene", suffixes=c(".dream", ".limma"))

# The 23 key genes
genes_23 <- c("RPS4Y1", "EIF1AY", "DDX3Y", "KDM5D", "MPDZ", "UTY", "TTTY14", "ZFY", "PRKY", "TPSD1", "C2", 
              "GPX1", "CLDN23", "XIST", "HDC", "TPSB2", "KIT", "CPA3", "FCER1A", "TPSAB1", "EIF1AX", "FCER1G", "TPSG1")

df <- df_all[df_all$Gene %in% genes_23, ]

# Classify significance
df$sig_dream <- df$adj.P.Val.dream < 0.05
df$sig_limma <- df$adj.P.Val.limma < 0.05

df$category <- "not significant in either"
df$category[df$sig_dream | df$sig_limma] <- "significant in one only"
df$category[df$sig_dream & df$sig_limma] <- "significant in both"

df$category <- factor(df$category, levels=c("significant in both", "significant in one only", "not significant in either"))

# Calculate median delta
median_delta <- median(abs(df$logFC.limma + df$logFC.dream), na.rm=TRUE)

title_str <- "Effect sizes agree closely; significance calls do not always"
sub_str <- sprintf("median |\u0394 log\u2082FC| = %.2f over 23 genes \u00b7 dashed line = exact mirror agreement", median_delta)

# Colors matching publication figure
col_sig_both <- "#2A76D2" # blue
col_sig_one <- "#ED713A"  # orange
col_sig_none <- "#8D8D8D" # grey

# Prepare inset box coordinates
box_xmin <- -1.7
box_xmax <- 1.7
box_ymin <- -1.7
box_ymax <- 1.7

# Hide labels inside inset box for main plot
df$main_label <- df$Gene
in_box <- df$logFC.dream >= box_xmin & df$logFC.dream <= box_xmax & 
          df$logFC.limma >= box_ymin & df$logFC.limma <= box_ymax
df$main_label[in_box] <- ""

# Main plot
p_main <- ggplot(df, aes(x=logFC.dream, y=logFC.limma, color=category)) +
  geom_abline(intercept=0, slope=-1, linetype="dashed", color="#8c8c8c", linewidth=0.8) +
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
                  segment.color="transparent", seed=42, 
                  min.segment.length=0) +
  theme_minimal(base_size=14) +
  labs(title=title_str,
       subtitle=sub_str,
       x="Manuscript (DREAM)  log\u2082FC  breast vs foreskin",
       y="F5-F6 integration (limma-voom)  log\u2082FC  foreskin vs breast") +
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

# Inset plot
df_inset <- df[in_box, ]

p_inset <- ggplot(df_inset, aes(x=logFC.dream, y=logFC.limma, color=category)) +
  geom_abline(intercept=0, slope=-1, linetype="dashed", color="#8c8c8c", linewidth=0.6) +
  geom_point(size=3) +
  scale_color_manual(values=c("significant in both"=col_sig_both, 
                              "significant in one only"=col_sig_one, 
                              "not significant in either"=col_sig_none)) +
  geom_text_repel(aes(label=Gene), color="#595959", size=3, max.overlaps=Inf,
                  point.padding=0.2, box.padding=0.3,
                  segment.color="transparent", seed=42) +
  theme_minimal(base_size=10) +
  coord_cartesian(xlim=c(box_xmin, box_xmax), ylim=c(box_ymin, box_ymax)) +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    axis.text = element_text(color="#8c8c8c"),
    panel.grid.major = element_line(color="#f0f0f0"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color="#f0f0f0", fill=NA, linewidth=1),
    panel.background = element_rect(fill="#fafafa", color=NA)
  )

p_final <- p_main + annotation_custom(ggplotGrob(p_inset), 
                                      xmin=3, xmax=11, 
                                      ymin=2, ymax=10)

# Save exact plot to results
out_file1 <- file.path(out_dir_fvb, "Concordance_Excel_vs_Limma.png")
out_file2 <- file.path(out_dir_bvf, "Concordance_Excel_vs_Limma.png")

ggsave(out_file1, plot=p_final, width=9, height=9, dpi=300)
ggsave(out_file2, plot=p_final, width=9, height=9, dpi=300)

message("Concordance plot generated successfully and saved as Concordance_Excel_vs_Limma.png")
