# 统一导出与三个图件；只使用已经计算的 MR 结果，不在绘图时重复运行 MR。
save_eqtl_plot <- function(plot, name, output_dir, width, height) {
  ggplot2::ggsave(file.path(output_dir, paste0(name, ".png")), plot,
                 device = ragg::agg_png, width = width, height = height, dpi = 300, bg = "white")
  ggplot2::ggsave(file.path(output_dir, paste0(name, ".pdf")), plot,
                 device = grDevices::cairo_pdf, width = width, height = height, bg = "white")
}
make_eqtl_plots <- function(summary_table, mr, loo, harmonised, output_dir) {
  colors <- ggsci::pal_jama()(3)
  plot_theme <- theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          strip.background = element_rect(fill = "grey95"),
          strip.text = element_text(face = "bold"),
          plot.margin = margin(10, 20, 10, 10),
          plot.title = element_text(face = "bold", size = 13))
  forest_data <- summary_table |>
    arrange(OR) |>
    mutate(GeneSymbol = factor(GeneSymbol, levels = rev(GeneSymbol)),
           label = sprintf("%.3f (%.3f-%.3f)   P=%.3g", OR, OR_low, OR_high, p_value))
  forest <- ggplot(forest_data, aes(OR, GeneSymbol)) +
    geom_vline(xintercept = 1, linetype = 2, color = "grey50") +
    geom_errorbar(aes(xmin = OR_low, xmax = OR_high), orientation = "y", width = .18) +
    geom_point(size = 2.5, color = colors[1]) + scale_x_log10() +
    labs(x = "Odds ratio (95% CI; logarithmic scale)", y = NULL) + plot_theme
  labels <- ggplot(forest_data, aes(0, GeneSymbol, label = label)) +
    geom_text(hjust = 0, size = 3.2) +
    scale_x_continuous(limits = c(0, 1)) + theme_void() +
    labs(title = "OR (95% CI) and P value") +
    theme(plot.title = element_text(size = 10, face = "bold"))
  combined <- patchwork::wrap_plots(forest, labels, widths = c(1.1, 1.15)) +
    patchwork::plot_annotation(title = "eQTL Mendelian randomization: colorectal cancer")
  save_eqtl_plot(combined, "MR_forest", output_dir, 10, max(4.5, nrow(forest_data) * .34 + 1.2))
  if (nrow(loo)) {
    loo_data <- loo |>
      mutate(label = if_else(SNP == "All", "All SNPs", paste("Drop", SNP)),
             label = factor(label, levels = rev(c(sort(unique(label[label != "All SNPs"])), "All SNPs"))),
             full = SNP == "All")
    loo_plot <- ggplot(loo_data, aes(b, label, color = full)) +
      geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
      geom_errorbar(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se),
                    orientation = "y", width = .18) +
      geom_point(size = 2) +
      scale_color_manual(values = c(colors[2], colors[1]), guide = "none") +
      facet_wrap(~ exposure, scales = "free_y", ncol = 2) + plot_theme +
      labs(x = "IVW beta (95% CI)", y = NULL, title = "Leave-one-out sensitivity analysis",
           subtitle = "All SNPs: full estimate; Drop: one instrument omitted")
    save_eqtl_plot(loo_plot, "MR_leave_one_out", output_dir, 10, ceiling(length(unique(loo$exposure)) / 2) * 2.8)
  }
  funnel_genes <- summary_table$GeneSymbol[summary_table$n_snp >= 3]
  if (length(funnel_genes)) {
    funnel_data <- harmonised |>
      filter(GeneSymbol %in% funnel_genes) |>
      mutate(wald_beta = beta.outcome / beta.exposure,
             wald_se = se.outcome / abs(beta.exposure), precision = 1 / wald_se)
    reference <- mr |>
      filter(exposure %in% funnel_genes, method %in% c("Inverse variance weighted", "MR Egger")) |>
      mutate(GeneSymbol = exposure)
    funnel <- ggplot(funnel_data, aes(wald_beta, precision)) +
      geom_vline(data = reference, aes(xintercept = b, color = method, linetype = method)) +
      geom_point(size = 2.2) +
      ggrepel::geom_text_repel(aes(label = SNP), size = 2.7, seed = 123,
                               max.overlaps = Inf, min.segment.length = 0) +
      scale_color_manual(values = setNames(colors[1:2], c("Inverse variance weighted", "MR Egger"))) +
      scale_x_continuous(expand = expansion(mult = .25)) +
      scale_y_continuous(expand = expansion(mult = .20)) +
      facet_wrap(~ GeneSymbol, scales = "free", ncol = 2) + plot_theme +
      theme(legend.position = "bottom") +
      labs(x = "Single-SNP Wald ratio", y = "Precision (1 / SE of Wald ratio)",
           title = "Instrument precision plots", color = NULL, linetype = NULL,
           subtitle = "Approximate SE; few instruments per gene: descriptive inspection only")
    save_eqtl_plot(funnel, "MR_funnel", output_dir, 10, ceiling(length(funnel_genes) / 2) * 3)
  }
}
