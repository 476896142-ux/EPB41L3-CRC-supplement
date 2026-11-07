# 全基因 GSEA：独立重绘，不重新计算 GSEA ----
suppressPackageStartupMessages({library(enrichplot); library(ggplot2); library(stringr)})
output_dir <- "E:/codex/2026-10-31/dui/outputs/GSEA_all_genes"
gsea_go <- readRDS(file.path(output_dir, "GO_GSEA.rds"))
gsea_kegg <- readRDS(file.path(output_dir, "KEGG_GSEA.rds"))
save_figure <- function(p, name, width, height) {
  ggsave(file.path(output_dir, paste0(name, ".png")), p, width = width, height = height,
    dpi = 300, bg = "white", device = ragg::agg_png)
  ggsave(file.path(output_dir, paste0(name, ".pdf")), p, width = width, height = height,
    device = grDevices::cairo_pdf, bg = "white")
}
keep_significant <- function(object) {
  object@result <- object@result[is.finite(object@result$p.adjust) & object@result$p.adjust < 0.05, ]
  object@result <- object@result[order(object@result$p.adjust, -abs(object@result$NES), object@result$ID), ]
  object
}
go_sig <- keep_significant(gsea_go)
kegg_sig <- keep_significant(gsea_kegg)
if (nrow(kegg_sig@result)) save_figure(dotplot(kegg_sig, showCategory = 15) + theme_bw(base_size = 11), "KEGG_dotplot", 9, 10)
if (nrow(go_sig@result)) {
  p <- dotplot(go_sig, showCategory = 8, split = "ONTOLOGY", label_format = 42) +
    facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") + theme_bw(base_size = 11)
  save_figure(p, "GO_dotplot", 10, 12)
}
plot_curves <- function(object, label) {
  selected <- head(object@result, 5)
  if (!nrow(selected)) return(invisible(NULL))
  object@result$Description <- str_wrap(object@result$Description, 48)
  p <- gseaplot2(object, geneSetID = selected$ID, title = paste(label, "GSEA: all-gene ranking"),
    color = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F")[seq_len(nrow(selected))],
    base_size = 11, subplots = 1:3, rel_heights = c(2, 0.6, 0.8), pvalue_table = FALSE)
  save_figure(p, paste0(label, "_top5_curves"), 11, 9)
  write.csv(selected, file.path(output_dir, paste0(label, "_top5_curve_terms.csv")), row.names = FALSE)
}
plot_curves(kegg_sig, "KEGG")
plot_curves(go_sig, "GO")

# 联合棒棒糖图：GO BP/CC/MF + KEGG，每组 FDR 前五 ----
columns <- c("ID", "Description", "setSize", "Count", "GeneRatio", "NES", "p.adjust")
g <- go_sig@result[, columns, drop = FALSE]
g$Group <- paste("GO", go_sig@result$ONTOLOGY)
k <- kegg_sig@result[, columns, drop = FALSE]
k$Group <- "KEGG"
all_results <- rbind(g, k)
group_order <- c("GO BP", "GO CC", "GO MF", "KEGG")
palette <- c("GO BP" = "#56A45A", "GO CC" = "#369BC1", "GO MF" = "#A66BB0", "KEGG" = "#53689A")
plot_data <- do.call(rbind, lapply(group_order, function(group) {
  rows <- all_results[all_results$Group == group, ]
  rows <- rows[order(rows$p.adjust, -abs(rows$NES), rows$ID), ]
  head(rows, 5)
}))
if (nrow(plot_data)) {
  plot_data$Group <- factor(plot_data$Group, levels = group_order)
  plot_data$y <- rev(seq_len(nrow(plot_data))) + (4 - as.integer(plot_data$Group)) * 0.9
  plot_data$Label <- str_wrap(plot_data$Description, 43)
  plot_data$neg_log10_FDR <- -log10(pmax(plot_data$p.adjust, .Machine$double.xmin))
  ranges <- do.call(rbind, lapply(group_order, function(group) {
    y <- plot_data$y[plot_data$Group == group]
    if (!length(y)) return(NULL)
    data.frame(Group = group, ymin = min(y) - 0.42, ymax = max(y) + 0.42)
  }))
  p <- ggplot(plot_data, aes(y = y)) +
    geom_segment(aes(x = 0, xend = GeneRatio, yend = y, color = Group), linewidth = 0.65, show.legend = FALSE) +
    geom_segment(data = ranges, aes(x = 0, xend = 0, y = ymin, yend = ymax, color = Group),
      inherit.aes = FALSE, linewidth = 2.3) +
    geom_point(aes(x = GeneRatio, size = Count, fill = neg_log10_FDR), shape = 21, color = "#777777", stroke = 0.4) +
    geom_text(aes(x = -0.025, label = Label, color = Group), hjust = 1, size = 3.5,
      fontface = "bold", lineheight = 0.95, show.legend = FALSE) +
    scale_color_manual(values = palette, breaks = group_order, name = "Category") +
    scale_fill_gradientn(colors = c("#2596C1", "#A8B3AF", "#F3A000"), name = expression(-log[10](FDR))) +
    scale_size_continuous(name = "Core genes", range = c(3, 7.5)) +
    scale_x_continuous(breaks = seq(0, 1, 0.2), expand = expansion(mult = c(0, 0.02))) +
    scale_y_continuous(expand = expansion(add = 0.8)) + coord_cartesian(xlim = c(0, 1.02), clip = "off") +
    labs(x = "Core enrichment gene ratio", y = NULL, title = "GO and KEGG GSEA enrichment",
      subtitle = "All-gene logFC ranking | Top 5 terms per category | FDR < 0.05",
      caption = "Gene ratio = leading-edge genes / gene set size within the ranked list.") +
    guides(fill = guide_colorbar(order = 1), color = guide_legend(order = 2), size = guide_legend(order = 3)) +
    theme_classic(base_size = 11) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      axis.line.y = element_blank(), axis.title.x = element_text(margin = margin(t = 10)),
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(color = "#555555", margin = margin(b = 16)),
      plot.caption = element_text(size = 8, color = "#666666", margin = margin(t = 14)),
      legend.position = "right", legend.title = element_text(size = 10), legend.text = element_text(size = 9),
      plot.margin = margin(t = 16, r = 16, b = 12, l = 270))
  save_figure(p, "GO_KEGG_lollipop_all_genes", 12.5, max(8, nrow(plot_data) * 0.42 + 2))
  write.csv(plot_data, file.path(output_dir, "lollipop_selected_terms.csv"), row.names = FALSE)
}
message("全基因图表已输出。")
