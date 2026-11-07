rm(list=ls())
library(tidyverse)
# Set IEUGWASR_TOKEN (or OPENGWAS_JWT) in your local environment before running.
# Never place a token in this script or commit it to version control.
tok <- Sys.getenv("IEUGWASR_TOKEN", unset = Sys.getenv("OPENGWAS_JWT", unset = ""))
if (!nzchar(tok)) stop("Set IEUGWASR_TOKEN locally before running this script.")
Sys.setenv(IEUGWASR_TOKEN = tok, OPENGWAS_JWT = tok)

library(ieugwasr)
setwd("D:/Bioinfor/colorectal/MR/")
# ========== 参数配置 ==========
ld_ref <- "D:/Bioinfor/colorectal/MR/g1000_eur/g1000_eur"  # 面板前缀
plink_path <- "D:/Bioinfor/colorectal/MR/LD/plink_win64/plink.exe"

p_thresh <- 5e-8      # 顺式pQTL显著性阈值
r2_thresh <- 0.001    # LD r²阈值
kb_window <- 10000    # 10Mb窗口
f_thresh <- 10        # 弱工具F统计量阈值

out_rds  <- "pQTL_exp_dat_clump.rds"   # RDS 保存路径(保留所有列)
out_csv  <- "pQTL_exp_dat_clump.csv"   # CSV 保存路径(便于查看/共享)
result_dir <- file.path("D:/Bioinfor/colorectal/MR/pQTL")
# ========== 1. 读取并格式化数据 ==========
exp_dat <- read_csv("pQTL/cis.exposure_data.csv", show_col_types = FALSE) %>%
  mutate(
    pval = as.numeric(pval.exposure),
    beta = as.numeric(beta.exposure),
    se = as.numeric(se.exposure),
    F_stat = (beta / se)^2
  )

# ========== 2. 按蛋白分组做独立LD clump ==========
clump_one_protein <- function(df) {
  tryCatch({
    # 关键：重命名为ieugwasr要求的标准列名 rsid + pval
    snp_input <- df %>% select(rsid = SNP, pval)
    
    clump_res <- ieugwasr::ld_clump(
      dat = snp_input,
      clump_kb = kb_window,
      clump_r2 = r2_thresh,
      clump_p = p_thresh,
      bfile = ld_ref,
      plink_bin = plink_path
    )
    
    # 返回保留的独立SNP对应行
    df %>% filter(SNP %in% clump_res$rsid)
  }, error = function(e) {
    message("蛋白 ", unique(df$exposure), " 处理失败: ", e$message)
    return(tibble())
  })
}

# 先按P值初筛，再按蛋白分组clump
exp_dat_clump <- exp_dat %>%
  filter(pval < p_thresh) %>%       # 先初筛显著SNP，减少计算量
  group_by(id.exposure, exposure) %>%
  group_split() %>%
  lapply(clump_one_protein) %>%
  bind_rows()
# ---------- 保存 ----------
saveRDS(exp_dat_clump, file.path(result_dir,"exp_dat_clump.rds"))
write.csv(exp_dat_clump, file.path(result_dir,"exp_dat_clump.csv"))

# ========== 3. 弱工具筛选（F > 10） ==========
exp_dat_final <- exp_dat_clump %>%
  filter(F_stat > f_thresh)

# ========== 4. 结果汇总与导出 ==========
iv_summary <- exp_dat_final %>%
  group_by(id.exposure, exposure) %>%
  summarise(
    n_iv = n(),
    min_pval = min(pval),
    max_F = max(F_stat),
    chr = unique(chr.exposure),
    .groups = "drop"
  )

write_csv(exp_dat_final,
          file.path(result_dir, "IVs_selected_cis_pQTL.csv"))

write_csv(iv_summary,file.path(result_dir, "IV_summary_per_protein.csv"))

cat("筛选完成：共保留", nrow(iv_summary), "个蛋白，", nrow(exp_dat_final), "个独立工具变量\n")
##############----------------Mr分析------------------------###################
# ---- 读取已完成的工具变量筛选结果 ----
iv_raw <- readr::read_csv(file.path(result_dir,'IVs_selected_cis_pQTL.csv'), show_col_types = FALSE) |>
  transmute(
    id.exposure, exposure, SNP,
    beta.exposure = as.numeric(beta.exposure),
    se.exposure = as.numeric(se.exposure),
    effect_allele.exposure = toupper(effect_allele.exposure),
    other_allele.exposure = toupper(other_allele.exposure),
    eaf.exposure = as.numeric(eaf.exposure),
    pval.exposure = as.numeric(pval.exposure),
    samplesize.exposure = as.numeric(samplesize.exposure),
    chr = str_remove(as.character(chr.exposure), "^chr"),
    pos = as.integer(pos.exposure)
  ) |>
  filter(!is.na(SNP), !is.na(chr), !is.na(pos),
         !is.na(beta.exposure), !is.na(se.exposure)) |>
  distinct(SNP, .keep_all = TRUE)
message("IV rows: ", nrow(iv_raw), "; proteins: ", n_distinct(iv_raw$exposure))

# ---- 使用 VCF tabix 索引提取 IV 的 outcome 关联 ----
library(ieugwasr)

Sys.getenv("IEUGWASR_TOKEN")
gwasinfo("ebi-a-GCST90013862",
         opengwas_jwt = Sys.getenv("IEUGWASR_TOKEN"))
# 拿 ncase / ncontrol / samplesize / build
extract_outcome <- function(iv, vcf_file) {
  query <- GRanges(seqnames = iv$chr,
                   ranges = IRanges(start = iv$pos, end = iv$pos))
  param <- ScanVcfParam(which = query,
                        geno = c("ES", "SE", "LP", "AF"))
  vcf <- readVcf(vcf_file, genome = "hg19", param = param)
  if (!length(vcf)) return(tibble())
  tibble(
    SNP = as.character(rowRanges(vcf)$ID),
    effect_allele.outcome = toupper(as.character(unlist(alt(vcf)))),
    other_allele.outcome = toupper(as.character(ref(vcf))),
    beta.outcome = as.numeric(geno(vcf)$ES[, 1]),
    se.outcome = as.numeric(geno(vcf)$SE[, 1]),
    pval.outcome = 10^(-as.numeric(geno(vcf)$LP[, 1])),
    eaf.outcome = as.numeric(geno(vcf)$AF[, 1])
  ) |>
    filter(!is.na(SNP), !str_detect(effect_allele.outcome, ","),
           !str_detect(other_allele.outcome, ",")) |>
    distinct(SNP, .keep_all = TRUE)
}

library(dplyr)
library(TwoSampleMR)


outcome_raw <- extract_outcome_data(
  snps     = iv_raw$SNP,        # 字符向量，API 路径按 rsID 查，不受 build 影响
  outcomes = "ebi-a-GCST90013862",
  proxies  = TRUE
) |> 
  mutate(id.outcome = "ebi-a-GCST90013862",
         outcome = "Colorectal cancer",
         samplesize.outcome = 407746,
         ncase.outcome = 4161,
         ncontrol.outcome = 403585)


readr::write_csv(outcome_raw, file.path(result_dir, "CRC_outcome_for_selected_IVs.csv"))
message("Outcome SNPs extracted: ", nrow(outcome_raw))

# ---- 格式化、等位基因协调与 MR ----
names(iv_raw)

exposure_dat <- iv_raw |>
  select(id.exposure, exposure, SNP, beta.exposure, se.exposure,
         effect_allele.exposure, other_allele.exposure, eaf.exposure,
         pval.exposure, samplesize.exposure)

outcome_dat <- outcome_raw 


harmonised <- harmonise_data(exposure_dat, outcome_dat) |>
  filter(mr_keep)

if (!nrow(harmonised)) stop("IV 与 outcome 没有可协调的 SNP 交集。")
readr::write_csv(harmonised, file.path(result_dir, "harmonised_selected_IVs.csv"))

mr_results <- harmonised |>
  group_by(id.exposure, exposure, id.outcome, outcome) |>
  group_split(.keep = TRUE) |>
  purrr::map_dfr(function(dat) {
    tryCatch(
      mr(dat, method_list = c("mr_ivw", "mr_egger_regression",
                              "mr_weighted_median", "mr_weighted_mode",
                              "mr_simple_mode")),
      error = function(e) tibble()
    )
  }) |>
  left_join(harmonised |>
              count(id.exposure, exposure, id.outcome, outcome, name = "nsnp"),
            by = c("id.exposure", "exposure", "id.outcome", "outcome")) |>
  group_by(method) |>
  mutate(pval_fdr = p.adjust(pval, method = "fdr")) |>
  ungroup()
readr::write_csv(mr_results, file.path(result_dir, "MR_all_methods.csv"))

# ---- IVW 主结果与 DEG 交集 ----
ivw_results <- mr_results |>
  filter(method == "Inverse variance weighted") |>
  mutate(OR = exp(b), OR_low = exp(b - 1.96 * se),
         OR_high = exp(b + 1.96 * se),
         MR_significant = pval < 0.05)
readr::write_csv(ivw_results, file.path(result_dir, "MR_IVW_results.csv"))
deg_file <- "D:/Bioinfor/colorectal/results/GSE196006_GSE146009_GSE110224/DEG_Tumor_vs_Normal_significant.csv"
deg <- readr::read_csv(deg_file, show_col_types = FALSE) |>
  transmute(GeneSymbol = toupper(str_trim(GeneSymbol)),
            DEG_logFC = logFC, DEG_FDR = adj.P.Val, DEG_Category = Category) |>
  filter(!is.na(GeneSymbol), GeneSymbol != "") |>
  distinct(GeneSymbol, .keep_all = TRUE)
mr_deg_intersection <- ivw_results |>
  filter(MR_significant) |>
  mutate(GeneSymbol = toupper(str_trim(exposure))) |>
  inner_join(deg, by = "GeneSymbol") |>
  arrange(pval_fdr)
readr::write_csv(mr_deg_intersection,
                 file.path(result_dir, "MR_IVW_FDR05_DEG_intersection.csv"))
######-------------------------韦恩图绘制--------------------------------------
library(VennDiagram)
library(dplyr)
library(stringr)
# ---------- 1. 读取数据 ----------
mr_dir     <- "D:/Bioinfor/colorectal/MR"
deg_file   <- "D:/Bioinfor/colorectal/results/GSE196006_GSE146009_GSE110224/DEG_Tumor_vs_Normal_significant.csv"
result_dir <- file.path('D:/Bioinfor/colorectal/MR/pQTL')

mr_res <- read.csv(file.path(result_dir, "MR_IVW_results.csv"))
mr_sig <- mr_res %>%
  filter(MR_significant == TRUE) %>%
  pull(exposure) %>% toupper() %>% unique()

deg <- read.csv(deg_file) %>%
  mutate(GeneSymbol = toupper(trimws(GeneSymbol))) %>%
  filter(!is.na(GeneSymbol), GeneSymbol != "") %>%
  distinct(GeneSymbol, .keep_all = TRUE)

deg_up   <- deg %>% filter(Category == "Up")   %>% pull(GeneSymbol)
deg_down <- deg %>% filter(Category == "Down") %>% pull(GeneSymbol)

# ---------- 2. 公共参数（重点修改标签位置） ----------
venn_common <- list(
  scaled = FALSE,
  ext.percent = 0.70,
  ext.text = FALSE,
  
  fill   = c("#4e79a7", "#b07aa1"),
  alpha  = 0.55,
  col    = c("#4e79a7", "#b07aa1"),
  lwd    = 2.5,
  lty    = "solid",
  
  cex          = 1.3,
  fontface     = "bold",
  cat.cex      = 0.70,   # 标签进一步缩小
  cat.fontface = "bold",
  cat.col      = c("#4e79a7", "#b07aa1"),
  # =========关键参数=========
  cat.pos = c(-20, 20), # 两个分类标签角度，左圈-20°，右圈+20°，向内收，不贴画布边缘
  cat.dist = 0.08,      # 标签距离圆圈的距离，向内拉近
  margin = 0.12,        # 大幅增加四周留白！！
  print.mode = "raw",
  main.cex = 0.9,
  main.fontface = "bold"
)

# ---------- 3. 图A ----------
do.call(venn.diagram, c(
  list(
    x = list(DEG_up = deg_up, MR_significant = mr_sig),
    filename = file.path(result_dir, "Venn_A_DEGup_MRsig.png"),
    height = 12, width = 18, units = "cm", resolution = 300,
    compression = "lzw",
    main = "(A) DEG_up vs MR_significant"
  ),
  venn_common
))

# ---------- 4. 图B ----------
do.call(venn.diagram, c(
  list(
    x = list(DEG_down = deg_down, MR_significant = mr_sig),
    filename = file.path(result_dir, "Venn_B_DEGdown_MRsig.png"),
    height = 12, width = 18, units = "cm", resolution = 300,
    compression = "lzw",
    main = "(B) DEG_down vs MR_significant"
  ),
  venn_common
))

# ---------- 5. 交集基因输出 ----------
inter_up   <- intersect(deg_up, mr_sig)
inter_down <- intersect(deg_down, mr_sig)
cat("上调交集 (", length(inter_up), "):", inter_up, "\n")
cat("下调交集 (", length(inter_down), "):", inter_down, "\n")
#####-------绘制森林图-----------------------###################################
library(tidyverse)
library(patchwork)

mr_dat <- tibble(
  exposure = c("GREM1","HHIP","LAMC2","ITIH5","NTRK3","ANGPTL1","UGDH","KLB" ),
  nsnp     = c(7,2,4,6,2,2,2,4),
  method   = rep("IVW", 8),
  pval     = c(6.638436e-06, 4.769232e-05, 8.704422e-04, 2.605977e-03, 1.390666e-02, 4.471245e-02, 4.123032e-02, 3.671420e-02),
  or       = c(1.3067461, 1.7850550, 0.7587666, 0.8091456, 0.5044520, 1.1918740, 1.1953791, 1.1518175),
  or_lower = c(1.1631594, 1.3500924, 0.6449527, 0.7049400, 0.2924220, 1.0041489, 1.0071252, 1.0087617),
  or_upper = c(1.4680579, 2.3601505, 0.8926652, 0.9287552, 0.8702214,1.4146942, 1.4188218, 1.3151605)
) %>%
  mutate(
    `OR(95% CI)` = paste0(sprintf("%.3f", or), " (",
                          sprintf("%.3f", or_lower), " to ",
                          sprintf("%.3f", or_upper), ")"),
    pval_text = sprintf("%.3f", pval),
    col = "#E63946",
    y = nrow(.):1
  )

# 固定每行高度，避免森林图和文本列错位
row_height <- 0.8
y_min <- min(mr_dat$y) - 0.5
y_max <- max(mr_dat$y) + 0.5

# -------- 左侧表格 --------
p_table_left <- mr_dat %>%
  select(y, exposure, nsnp, method, pval_text) %>%
  mutate(
    nsnp = as.character(nsnp),
    pval_text = as.character(pval_text)
  ) %>%
  pivot_longer(-y, names_to = "colname", values_to = "val") %>%
  mutate(colname = factor(colname,
                          levels = c("exposure", "nsnp", "method", "pval_text"),
                          labels = c("exposure", "nsnp", "method", "pval"))) %>%
  ggplot(aes(x = colname, y = y)) +
  geom_rect(data = ~ filter(.x, y %% 2 == 0),
            aes(xmin = -Inf, xmax = Inf, ymin = y - 0.5, ymax = y + 0.5),
            fill = "#f0f0f0", inherit.aes = FALSE) +
  geom_text(aes(label = val), size = 3.2, hjust = 0.5) +
  scale_x_discrete(limits = c("exposure", "nsnp", "method", "pval")) +
  scale_y_continuous(limits = c(y_min, y_max), expand = c(0, 0)) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 9),
    axis.title = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  )

# -------- 中间森林图：修复每行显示不全 --------
p_forest <- ggplot(mr_dat, aes(y = y)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = or_lower, xmax = or_upper),
                 height = 0.12, linewidth = 0.6) +
  geom_point(aes(x = or, color = col), size = 2.5) +
  scale_color_identity() +
  scale_x_continuous(trans = "log",
                     breaks = c(0.2, 0.6, 1.0, 1.4, 1.8, 2.1, 2.5),
                     limits = c(0.2, 2.6)) +
  scale_y_continuous(limits = c(y_min, y_max), expand = c(0, 0)) +
  labs(x = "Odds Ratio") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 9),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 10),
    plot.margin = margin(0, 0, 0, 0)
  )

# -------- 右侧 OR(95% CI) 文本列 --------
p_table_right <- mr_dat %>%
  select(y, `OR(95% CI)`) %>%
  pivot_longer(-y, names_to = "colname", values_to = "val") %>%
  ggplot(aes(x = colname, y = y)) +
  geom_rect(data = ~ filter(.x, y %% 2 == 0),
            aes(xmin = -Inf, xmax = Inf, ymin = y - 0.5, ymax = y + 0.5),
            fill = "#f0f0f0", inherit.aes = FALSE) +
  geom_text(aes(label = val), size = 3.2, hjust = 0.5) +
  scale_x_discrete(limits = "OR(95% CI)") +
  scale_y_continuous(limits = c(y_min, y_max), expand = c(0, 0)) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 9),
    axis.title = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  )

# -------- 拼接 --------
p_final <- wrap_plots(
  p_table_left,
  p_forest,
  p_table_right,
  ncol = 3,
  widths = c(4, 2.5, 3.5)
)

print(p_final)

ggsave("MR_forest_Fig5.png", p_final,
       width = 18, height = 7, dpi = 300)

##############################################################################

summary_tbl <- tibble(
  metric = c("IV proteins", "IV SNPs", "Outcome SNPs extracted",
             "Harmonised SNP rows", "MR proteins", "IVW FDR<0.05",
             "MR-DEG intersection"),
  value = c(n_distinct(iv_raw$exposure), n_distinct(iv_raw$SNP),
            nrow(outcome_raw), nrow(harmonised), n_distinct(mr_results$exposure),
            n_distinct(ivw_results$exposure[ivw_results$MR_significant]),
            nrow(mr_deg_intersection))
)
readr::write_csv(summary_tbl, file.path(result_dir, "MR_run_summary.csv"))
message("MR completed. Results written to: ", result_dir)
# ---- 计数 ----
cat("IVW FDR<0.05 蛋白数:", n_distinct(mr_results$exposure[mr_results$method=="Inverse variance weighted" & mr_results$pval_fdr < 0.05]), "\n")
cat("IVW raw p<0.05 蛋白数:", n_distinct(ivw_results$exposure[ivw_results$MR_significant]), "\n")
cat("MR-DEG 交集基因数:", nrow(mr_deg_intersection), "\n")

mr_deg_intersection |>
  select(exposure, pval, pval_fdr, OR, OR_low, OR_high, DEG_logFC, DEG_FDR, DEG_Category)

summary_tbl <- tibble(
  metric = c("IV proteins", "IV SNPs", "Outcome SNPs extracted",
             "Harmonised SNP rows", "MR proteins", "IVW raw p<0.05",
             "MR-DEG intersection"),
  value = c(n_distinct(iv_raw$exposure), n_distinct(iv_raw$SNP),
            nrow(outcome_raw), nrow(harmonised), n_distinct(mr_results$exposure),
            n_distinct(ivw_results$exposure[ivw_results$MR_significant]),
            nrow(mr_deg_intersection))
)
# ==============================================================================
# pQTL -> 结直肠癌 MR 敏感性分析
#
# 目的：按每个蛋白保留的 SNP 数量，安全地运行主 MR、异质性检验、
#       Egger 截距、留一法和 MR-PRESSO，并统一导出表格与图片。
# ==============================================================================

# 1. 参数 ----------------------------------------------------------------------
target_proteins <- c(
  "GREM1", "HHIP", "LAMC2", "ITIH5",
  "NTRK3", "ANGPTL1", "UGDH", "KLB"
)
output_dir <- "D:/Bioinfor/colorectal/MR/pQTL/pQTL_heterogeneity_review"
set.seed(20240918)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 2. 依赖包 --------------------------------------------------------------------
required_packages <- c(
  "dplyr", "purrr", "readr", "tidyr", "ggplot2",
  "writexl", "TwoSampleMR", "MRPRESSO", "ragg"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "缺少以下 R 包：", paste(missing_packages, collapse = ", "),
    "。请安装后重新运行。"
  )
}

# 3. 辅助函数 ------------------------------------------------------------------
empty_result <- function(columns) {
  result <- as.data.frame(
    setNames(replicate(length(columns), logical(0), simplify = FALSE), columns),
    stringsAsFactors = FALSE
  )
  result
}

first_existing_value <- function(object, candidate_names, default = NA_real_) {
  if (is.null(object)) {
    return(default)
  }
  
  matching_names <- intersect(candidate_names, names(object))
  if (length(matching_names) == 0) {
    return(default)
  }
  
  value <- object[[matching_names[[1]]]]
  if (length(value) == 0) default else value[[1]]
}

save_ggplot <- function(plot_object, file_stem, width, height) {
  png_file <- file.path(output_dir, paste0(file_stem, ".png"))
  pdf_file <- file.path(output_dir, paste0(file_stem, ".pdf"))
  
  ggplot2::ggsave(
    filename = png_file,
    plot = plot_object,
    device = ragg::agg_png,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    background = "white"
  )
  
  ggplot2::ggsave(
    filename = pdf_file,
    plot = plot_object,
    device = grDevices::cairo_pdf,
    width = width,
    height = height,
    units = "in",
    bg = "white"
  )
}

# 4. 读取、检查和清理数据 -------------------------------------------------------
if (!file.exists(input_file)) {
  stop("找不到输入文件：", input_file)
}

input_file <- "D:/Bioinfor/colorectal/MR/pQTL/harmonised_selected_IVs.csv"

harmonised_raw <- readr::read_csv(
  input_file,
  guess_max = Inf,
  show_col_types = FALSE,
  progress = FALSE
)

required_columns <- c(
  "SNP", "exposure", "outcome", "id.exposure", "id.outcome",
  "beta.exposure", "se.exposure", "beta.outcome", "se.outcome", "mr_keep"
)

missing_columns <- setdiff(required_columns, names(harmonised_raw))
if (length(missing_columns) > 0) {
  stop("输入文件缺少必要列：", paste(missing_columns, collapse = ", "))
}

harmonised_data <- harmonised_raw |>
  dplyr::mutate(
    exposure = toupper(trimws(as.character(exposure))),
    id.exposure = as.character(id.exposure),
    id.outcome = as.character(id.outcome),
    outcome = as.character(outcome),
    mr_keep = as.logical(mr_keep),
    dplyr::across(
      dplyr::all_of(c(
        "beta.exposure", "se.exposure", "beta.outcome", "se.outcome",
        "pval.exposure", "pval.outcome", "eaf.exposure", "eaf.outcome"
      )),
      as.numeric
    )
  ) |>
  dplyr::filter(exposure %in% target_proteins, mr_keep %in% TRUE) |>
  dplyr::distinct(exposure, SNP, .keep_all = TRUE) |>
  dplyr::mutate(
    exposure = factor(exposure, levels = target_proteins)
  ) |>
  dplyr::arrange(exposure, SNP) |>
  dplyr::mutate(exposure = as.character(exposure))

critical_numeric_columns <- c(
  "beta.exposure", "se.exposure", "beta.outcome", "se.outcome"
)

invalid_numeric_rows <- harmonised_data |>
  dplyr::filter(
    dplyr::if_any(
      dplyr::all_of(critical_numeric_columns),
      ~ is.na(.x) | !is.finite(.x)
    )
  )

if (nrow(invalid_numeric_rows) > 0) {
  stop("目标数据中存在缺失或非有限的 beta/se，无法安全运行 MR。")
}

missing_proteins <- setdiff(target_proteins, unique(harmonised_data$exposure))
if (length(missing_proteins) > 0) {
  stop("以下目标蛋白没有可用工具变量：", paste(missing_proteins, collapse = ", "))
}

protein_data_list <- split(
  harmonised_data,
  factor(harmonised_data$exposure, levels = target_proteins),
  drop = TRUE
)

# 5. 建立方法可用性表 -----------------------------------------------------------
analysis_overview <- harmonised_data |>
  dplyr::count(exposure, name = "n_snp") |>
  dplyr::mutate(
    exposure = factor(exposure, levels = target_proteins)
  ) |>
  dplyr::arrange(exposure) |>
  dplyr::mutate(
    exposure = as.character(exposure),
    primary_mr = dplyr::case_when(
      n_snp == 1 ~ "Wald ratio",
      n_snp == 2 ~ "IVW",
      n_snp >= 3 ~ "IVW + MR-Egger + weighted median"
    ),
    heterogeneity = dplyr::case_when(
      n_snp == 1 ~ "不可估计",
      n_snp == 2 ~ "仅 IVW Q",
      n_snp >= 3 ~ "IVW Q + MR-Egger Q"
    ),
    egger_intercept_availability = ifelse(n_snp >= 3, "可估计", "不可估计（需要 >= 3 SNP）"),
    leave_one_out = ifelse(n_snp >= 3, "可估计", "不可估计（需要 >= 3 SNP）"),
    mr_presso = ifelse(n_snp >= 4, "运行", "不运行（需要 >= 4 SNP）")
  )

# 6. 主 MR 分析 -----------------------------------------------------------------
run_mr <- function(protein_data) {
  number_of_snps <- nrow(protein_data)
  
  method_list <- if (number_of_snps == 1) {
    "mr_wald_ratio"
  } else if (number_of_snps == 2) {
    "mr_ivw"
  } else {
    c("mr_ivw", "mr_egger_regression", "mr_weighted_median")
  }
  
  TwoSampleMR::mr(protein_data, method_list = method_list)
}

mr_results <- purrr::map_dfr(protein_data_list, run_mr) |>
  dplyr::mutate(
    exposure = factor(exposure, levels = target_proteins)
  ) |>
  dplyr::arrange(exposure, method) |>
  dplyr::mutate(exposure = as.character(exposure))

# 7. Cochran Q 异质性检验 -------------------------------------------------------
run_heterogeneity <- function(protein_data) {
  number_of_snps <- nrow(protein_data)
  if (number_of_snps < 2) {
    return(NULL)
  }
  
  method_list <- if (number_of_snps == 2) {
    "mr_ivw"
  } else {
    c("mr_ivw", "mr_egger_regression")
  }
  
  TwoSampleMR::mr_heterogeneity(
    protein_data,
    method_list = method_list
  )
}

heterogeneity_results <- purrr::map_dfr(
  protein_data_list,
  run_heterogeneity
) |>
  dplyr::mutate(
    exposure = factor(exposure, levels = target_proteins)
  ) |>
  dplyr::arrange(exposure, method) |>
  dplyr::mutate(exposure = as.character(exposure))

# 8. MR-Egger 截距检验 ----------------------------------------------------------
run_egger_intercept <- function(protein_data) {
  if (nrow(protein_data) < 3) {
    return(NULL)
  }
  TwoSampleMR::mr_pleiotropy_test(protein_data)
}

egger_results <- purrr::map_dfr(
  protein_data_list,
  run_egger_intercept
) |>
  dplyr::mutate(
    exposure = factor(exposure, levels = target_proteins)
  ) |>
  dplyr::arrange(exposure) |>
  dplyr::mutate(exposure = as.character(exposure))

# 9. IVW 留一法 -----------------------------------------------------------------
run_leave_one_out <- function(protein_data) {
  if (nrow(protein_data) < 3) {
    return(NULL)
  }
  
  TwoSampleMR::mr_leaveoneout(
    protein_data,
    method = TwoSampleMR::mr_ivw
  )
}

leave_one_out_results <- purrr::map_dfr(
  protein_data_list,
  run_leave_one_out
) |>
  dplyr::mutate(
    exposure = factor(exposure, levels = target_proteins)
  ) |>
  dplyr::arrange(exposure, SNP) |>
  dplyr::mutate(exposure = as.character(exposure))

# 10. MR-PRESSO -----------------------------------------------------------------
run_mr_presso <- function(protein_data) {
  current_protein <- unique(protein_data$exposure)[[1]]
  number_of_snps <- nrow(protein_data)
  
  if (number_of_snps < 4) {
    return(tibble::tibble(
      exposure = current_protein,
      n_snp = number_of_snps,
      global_rssobs = NA_real_,
      global_pvalue = NA_real_,
      outlier_count = NA_integer_,
      status = "not estimable: n_snp < 4"
    ))
  }
  
  presso_data <- protein_data |>
    dplyr::transmute(
      BetaExposure = beta.exposure,
      BetaOutcome = beta.outcome,
      SdExposure = se.exposure,
      SdOutcome = se.outcome
    ) |>
    as.data.frame()
  
  tryCatch(
    {
      presso_result <- MRPRESSO::mr_presso(
        BetaOutcome = "BetaOutcome",
        BetaExposure = "BetaExposure",
        SdOutcome = "SdOutcome",
        SdExposure = "SdExposure",
        OUTLIERtest = TRUE,
        DISTORTIONtest = TRUE,
        data = presso_data,
        NbDistribution = 5000,
        SignifThreshold = 0.05,
        seed = 20240918
      )
      
      presso_tests <- presso_result[["MR-PRESSO results"]]
      global_test <- presso_tests[["Global Test"]]
      outlier_test <- presso_tests[["Outlier Test"]]
      
      outlier_p_column <- intersect(
        c("Pvalue", "p.value", "P.value"),
        names(outlier_test)
      )
      
      outlier_count <- if (length(outlier_p_column) == 0) {
        global_p <- suppressWarnings(as.numeric(first_existing_value(global_test, c("Pvalue", "p.value", "P.value"))))
        if (is.finite(global_p) && global_p >= 0.05) 0L else NA_integer_
      } else {
        sum(outlier_test[[outlier_p_column[[1]]]] < 0.05, na.rm = TRUE)
      }
      
      tibble::tibble(
        exposure = current_protein,
        n_snp = number_of_snps,
        global_rssobs = as.numeric(first_existing_value(
          global_test,
          c("RSSobs", "RSS.obs")
        )),
        global_pvalue = as.numeric(first_existing_value(
          global_test,
          c("Pvalue", "p.value", "P.value")
        )),
        outlier_count = as.integer(outlier_count),
        status = "ok"
      )
    },
    error = function(error_condition) {
      tibble::tibble(
        exposure = current_protein,
        n_snp = number_of_snps,
        global_rssobs = NA_real_,
        global_pvalue = NA_real_,
        outlier_count = NA_integer_,
        status = paste0("error: ", conditionMessage(error_condition))
      )
    }
  )
}

mr_presso_results <- purrr::map_dfr(
  protein_data_list,
  run_mr_presso
) |>
  dplyr::mutate(
    exposure = factor(exposure, levels = target_proteins)
  ) |>
  dplyr::arrange(exposure) |>
  dplyr::mutate(exposure = as.character(exposure))

# 11. 汇总每个蛋白的主要结果 ---------------------------------------------------
primary_mr_results <- mr_results |>
  dplyr::filter(method %in% c("Wald ratio", "Inverse variance weighted")) |>
  dplyr::transmute(
    exposure,
    primary_method = method,
    primary_beta = b,
    primary_se = se,
    primary_ci_low = b - 1.96 * se,
    primary_ci_high = b + 1.96 * se,
    primary_pvalue = pval
  )

ivw_heterogeneity <- heterogeneity_results |>
  dplyr::filter(method == "Inverse variance weighted") |>
  dplyr::transmute(
    exposure,
    ivw_Q = Q,
    ivw_Q_df = Q_df,
    ivw_Q_pvalue = Q_pval
  )

egger_summary <- egger_results |>
  dplyr::transmute(
    exposure,
    egger_intercept = egger_intercept,
    egger_intercept_se = se,
    egger_intercept_pvalue = pval
  )

summary_results <- analysis_overview |>
  dplyr::left_join(primary_mr_results, by = "exposure") |>
  dplyr::left_join(ivw_heterogeneity, by = "exposure") |>
  dplyr::left_join(egger_summary, by = "exposure") |>
  dplyr::left_join(
    mr_presso_results |>
      dplyr::select(
        exposure,
        presso_global_pvalue = global_pvalue,
        presso_outlier_count = outlier_count,
        presso_status = status
      ),
    by = "exposure"
  )

# 12. 绘图 ----------------------------------------------------------------------
forest_plot_data <- summary_results |>
  dplyr::mutate(
    exposure = factor(exposure, levels = rev(target_proteins)),
    significance = ifelse(primary_pvalue < 0.05, "P < 0.05", "P >= 0.05")
  )

primary_forest_plot <- ggplot2::ggplot(
  forest_plot_data,
  ggplot2::aes(
    x = primary_beta,
    y = exposure,
    color = significance
  )
) +
  ggplot2::geom_vline(
    xintercept = 0,
    linewidth = 0.45,
    linetype = 2,
    color = "grey45"
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = primary_ci_low, xmax = primary_ci_high),
    orientation = "y",
    width = 0.20,
    linewidth = 0.65
  ) +
  ggplot2::geom_point(size = 2.7) +
  ggplot2::scale_color_manual(
    values = c("P < 0.05" = "#C0392B", "P >= 0.05" = "#365F91"),
    name = NULL
  ) +
  ggplot2::labs(
    title = "Primary Mendelian randomization estimates",
    subtitle = "IVW estimates with 95% confidence intervals",
    x = "Effect estimate (beta, 95% CI)",
    y = NULL
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 12),
    plot.subtitle = ggplot2::element_text(color = "grey35"),
    legend.position = "top",
    axis.text.y = ggplot2::element_text(face = "bold")
  )

#save_ggplot(primary_forest_plot, "primary_MR_forest", 6.8, 4.8)
ggsave(
  filename = file.path(result_dir, "primary_MR_forest.png"),
  plot     = primary_forest_plot,
  width    = 6.8,
  height   = 4.8,
  dpi      = 300
)

leave_one_out_plot_data <- leave_one_out_results |>
  dplyr::filter(SNP != "All") |>
  dplyr::mutate(
    exposure = factor(exposure, levels = target_proteins),
    SNP = factor(SNP, levels = rev(unique(SNP)))
  )

leave_one_out_plot <- ggplot2::ggplot(
  leave_one_out_plot_data,
  ggplot2::aes(x = b, y = SNP)
) +
  ggplot2::geom_vline(
    xintercept = 0,
    linewidth = 0.4,
    linetype = 2,
    color = "grey50"
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se),
    orientation = "y",
    width = 0.16,
    linewidth = 0.45,
    color = "#527AA3"
  ) +
  ggplot2::geom_point(size = 1.8, color = "#B03A2E") +
  ggplot2::facet_wrap(~ exposure, scales = "free", ncol = 2) +
  ggplot2::labs(
    title = "IVW leave-one-out analysis",
    x = "IVW effect after omitting one SNP",
    y = "Omitted SNP"
  ) +
  ggplot2::theme_bw(base_size = 9.5) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "#EAF0F6", color = "grey60"),
    strip.text = ggplot2::element_text(face = "bold"),
    plot.title = ggplot2::element_text(face = "bold", size = 12)
  )

#save_ggplot(leave_one_out_plot, "leave_one_out", 9.0, 7.0)
ggsave(
  filename = file.path(result_dir, "leave_one_out.png"),
  plot     = leave_one_out_plot,
  width    = 9.0,
  height   = 7.0,
  dpi      = 300
)
# 13. 导出 Excel 与 CSV ---------------------------------------------------------
excel_sheets <- list(
  Summary = summary_results,
  Analysis_overview = analysis_overview,
  MR_results = mr_results,
  Heterogeneity = heterogeneity_results,
  Egger_intercept = egger_results,
  MR_PRESSO = mr_presso_results,
  Leave_one_out = leave_one_out_results
)

writexl::write_xlsx(
  excel_sheets,
  path = file.path(output_dir, "pQTL_MR_sensitivity_results.xlsx")
)

readr::write_csv(
  summary_results,
  file.path(output_dir, "pQTL_MR_sensitivity_summary.csv"),
  na = ""
)

# 14. 最终完整性检查 ------------------------------------------------------------
stopifnot(
  nrow(harmonised_data) == sum(analysis_overview$n_snp),
  nrow(analysis_overview) == length(target_proteins),
  setequal(analysis_overview$exposure, target_proteins),
  nrow(primary_mr_results) == length(target_proteins),
  all(analysis_overview$n_snp >= 2),
  all(summary_results$primary_method == "Inverse variance weighted")
)

message("分析完成。输出目录：", output_dir)
message("纳入蛋白：", nrow(analysis_overview), "；纳入 SNP：", nrow(harmonised_data))
message("MR 结果：", nrow(mr_results), "；异质性结果：", nrow(heterogeneity_results))
message("Egger 截距：", nrow(egger_results), "；留一法：", nrow(leave_one_out_results))


