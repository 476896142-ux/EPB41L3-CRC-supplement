suppressPackageStartupMessages({
  library(arrow); library(data.table); library(dplyr); library(tidyr); library(stringr)
  library(MungeSumstats); library(TwoSampleMR);library(MendelR);library(QTLMR)
})
set.seed(123)
project_dir <- Sys.getenv("CRC_MR_DIR", "D:/Bioinfor/colorectal/MR/eQTL")
output_dir <- file.path(project_dir, "Mr_colorectal_result")
sig_file <- file.path(project_dir, "Colon_Sigmoid.v10.eQTLs.signif_pairs.parquet")
tra_file <- file.path(project_dir, "Colon_Transverse.v10.eQTLs.signif_pairs.parquet")
#ld_ref <- file.path(project_dir, "g1000_eur", "g1000_eur")
#plink_path <- file.path(project_dir, "LD", "plink_win64", "plink.exe")
ld_ref <- "D:/Bioinfor/colorectal/MR/g1000_eur/g1000_eur"  
plink_path <- "D:/Bioinfor/colorectal/MR/LD/plink_win64/plink.exe"
p_thresh <- 5e-6; 
r2_thresh <- 0.1; 
kb_window <- 10000; 
f_thresh <- 10
message("[1/6] 检查输入路径和 PLINK reference")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
required_files <- c(sig_file, tra_file, plink_path, paste0(ld_ref, c(".bed", ".bim", ".fam")))
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) stop("缺少输入文件：\n", paste(missing_files, collapse = "\n"))

message("[2/6] 读取并合并两个组织")
required_cols <- c("gene_id", "variant_id", "pval_nominal", "af", "slope", "slope_se","ma_count")

eqtl_data_sig <- arrow::read_parquet(sig_file, col_select = required_cols) |> mutate(tissue = "Colon_Sigmoid") |> as.data.table()

eqtl_data_tra <- arrow::read_parquet(tra_file, col_select = required_cols) |> mutate(tissue = "Colon_Transverse") |> as.data.table()

missing_cols <- setdiff(required_cols, union(names(eqtl_data_sig), names(eqtl_data_tra)))
if (length(missing_cols)) stop("eQTL 缺少必要列：", paste(missing_cols, collapse = ", "))
eqtl_all <- rbindlist(list(eqtl_data_sig, eqtl_data_tra), use.names = TRUE, fill = TRUE)
eqtl_all[, n_tissue := uniqueN(tissue), by = .(gene_id, variant_id)]
eqtl_dedup <- eqtl_all[order(pval_nominal), .SD[1L], by = .(gene_id, variant_id)]

message("[3/6] 整理字段并转换 rsID")

eqtl_data <- eqtl_dedup %>%
  separate(variant_id, into = c("CHR", "BP", "A1", "A2", "Build"), sep = "_",
           remove = FALSE, extra = "merge", fill = "right") %>%
  mutate(CHR = str_remove(as.character(CHR), "^chr"), BP = as.integer(BP),
         A1 = toupper(A1), A2 = toupper(A2), EAF = as.numeric(af),
         BETA = as.numeric(slope), SE = as.numeric(slope_se), P = as.numeric(pval_nominal),
         ma_count = as.numeric(ma_count),
         maf = pmin(EAF, 1 - EAF),
         samplesize = ma_count / (2 * maf))%>%
  filter(!is.na(CHR), !is.na(BP), !is.na(A1), !is.na(A2), is.finite(BETA),
         is.finite(SE), SE > 0, is.finite(P), P > 0, P <= 1,
         is.finite(ma_count), ma_count >= 0, maf > 0, is.finite(samplesize)
  )
View(eqtl_data)
make_key <- function(chr, bp, a1, a2) paste0(tolower(sub("^chr", "", chr)), ":", bp, ":",
                                             toupper(a1), ":", toupper(a2))
eqtl_data <- mutate(eqtl_data, join_key = make_key(CHR, BP, A1, A2))

formatted <- MungeSumstats::format_sumstats(
  as.data.frame(select(eqtl_data, all_of(c("CHR", "BP", "A1", "A2", "BETA", "P","samplesize")))),
  ref_genome = "GRCh38",
  nThread = max(1L, parallel::detectCores() - 1L), dbSNP = 155, return_data = TRUE,
  on_ref_genome = TRUE) |>
  mutate(join_key = make_key(CHR, BP, A1, A2)) |>
  left_join(select(eqtl_data, join_key, gene_id, tissue, n_tissue, EAF, SE), by = "join_key") |>
  filter(!is.na(gene_id))

if (!nrow(formatted)) stop("rsID 转换后无法回接 gene_id，请检查 variant_id 格式")

exposure_dat <- TwoSampleMR::format_data(
  as.data.frame(formatted),
  type = "exposure", snp_col = "SNP", beta_col = "BETA", se_col = "SE",
  effect_allele_col = "A2", other_allele_col = "A1", pval_col = "P", eaf_col = "EAF",samplesize_col = "N",
  phenotype_col = "gene_id") |>
  mutate(id.exposure = paste0("Colorectal eQTL: ", exposure), F_stat = (beta.exposure / se.exposure)^2) |>
  filter(pval.exposure <= p_thresh, is.finite(F_stat), F_stat >= f_thresh, !is.na(SNP), SNP != "") |>
  distinct(id.exposure, SNP, .keep_all = TRUE)

if (!nrow(exposure_dat)) stop("显著性和 F 统计量筛选后没有工具变量")

write.table(exposure_dat, file.path(output_dir, "Colorectal_exposure_preclump.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

clump_res <- ieugwasr::ld_clump(
  dat = exposure_dat |>
    rename(rsid = SNP, pval = pval.exposure),    # ← 改这两列名
  clump_kb  = kb_window,
  clump_r2  = r2_thresh,
  clump_p   = p_thresh,
  bfile     = ld_ref,
  plink_bin = plink_path
)
View(clump_res)
# clump 完改回 TwoSampleMR 标准列名
exp_dat_clump <- clump_res |>
  rename(SNP = rsid, pval.exposure = pval)
message("Clumping 完成: ", n_distinct(exp_dat_clump$exposure), " 个基因, ",
        nrow(exp_dat_clump), " 个 SNP")
# ========== 3. 弱工具筛选（F > 10） ==========
exp_dat_final <- exp_dat_clump %>%
  filter(F_stat > f_thresh)
# ========== 4. 结果汇总与导出 ==========
exp_dat_final <- exp_dat_final %>%
  left_join(
    formatted %>% select(SNP, CHR) %>% distinct(SNP, .keep_all = TRUE),
    by = "SNP"
  ) %>%
  rename(chr.exposure = CHR)
#readr::write_csv(exp_dat_final, "IVs_selected_colorectal.csv")
readr::write_csv(exp_dat_final,file.path(output_dir, "IVs_selected_colorectal.csv"))
View(exp_dat_final)
message("Clumping 完成: ", n_distinct(exp_dat_clump$exposure), " 个基因, ",
        nrow(exp_dat_clump), " 个 SNP")
exp_dat<- readr::read_csv(file.path(output_dir, "IVs_selected_colorectal.csv"), show_col_types = FALSE)
outcome_raw <- extract_outcome_data(
  snps     = exp_dat$SNP,        # 字符向量，API 路径按 rsID 查，不受 build 影响
  outcomes = "ebi-a-GCST90013862",
  proxies  = F
) |>
  mutate(id.outcome = "ebi-a-GCST90013862",
         outcome = "Colorectal cancer",
         samplesize.outcome = 407746,
         ncase.outcome = 4161,
         ncontrol.outcome = 403585)
View(outcome_raw)

harmonised_data <- harmonise_data(exp_dat, outcome_raw) |>
  filter(mr_keep)

# 增加file参数
setwd("D:/Bioinfor/colorectal/MR/eQTL")
save(eqtl_all, eqtl_dedup, eqtl_data, exp_dat, exp_dat_clump, exp_dat_final, outcome_raw, harmonised_data, file = "eqtl_result.RData")
###########--------------MR分析--------------------------############################
#MR分析 IVW为主，MR‑Egger,weighted median,simple mode,weighted mode
mr_res <- mr(harmonised_data, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_simple_mode","mr_weighted_mode"))
head(mr_res)
write.table(mr_res, file.path(output_dir, "mr_res.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

#save(mr_res,file.path(output_dir,"mr_res.RData"))
save(mr_res, file = file.path(output_dir, "mr_res.RData"))
load("mr_res.RData")
names(mr_res)
# ===== 加载依赖 =====
library(dplyr)
library(readr)
library(AnnotationDbi)
library(org.Hs.eg.db)

# ===== 1) 读取 DEG =====
deg_file   <- "D:/Bioinfor/colorectal/results/GSE196006_GSE146009_GSE110224/DEG_Tumor_vs_Normal_significant.csv"
deg <- readr::read_csv(deg_file, show_col_types = FALSE) |>
  transmute(GeneSymbol = toupper(str_trim(GeneSymbol)),
            DEG_logFC = logFC, DEG_FDR = adj.P.Val, DEG_Category = Category) |>
  filter(!is.na(GeneSymbol), GeneSymbol != "") |>
  distinct(GeneSymbol, .keep_all = TRUE)

# ===== 2) 构建 gene_map =====
# 先给 mr_res 加一列去版本号的暴露 ID，后续 join 统一用它
mr_res$exposure_clean <- sub("\\..*$", "", mr_res$exposure)

gene_ids <- unique(mr_res$exposure_clean)
gene_map <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys    = gene_ids,
                                  keytype = "ENSEMBL",
                                  columns = c("SYMBOL", "GENENAME")) |>
  filter(!is.na(SYMBOL) & !is.na(GENENAME)) |> 
  distinct(ENSEMBL, .keep_all = TRUE)

# ===== 3) MR 结果：不做 FDR，直接用 pval < 0.05 筛选 =====
mr_sig <- mr_res %>%
  filter(pval < 0.05)

#mr_res <- mr_res %>% mutate(fdr_p = p.adjust(pval, method = "BH")) # BH FDR校正

# FDR<0.05筛选，只有一个1个SNP，舍弃
#mr_sig_fdr <- mr_res %>% filter(fdr_p < 0.05)

# ===== 4) 注释：用 exposure_clean 与 gene_map 匹配 =====
mr_sig_annot <- mr_sig %>%
  left_join(gene_map, by = c("exposure_clean" = "ENSEMBL")) %>%
  mutate(GeneSymbol = toupper(coalesce(SYMBOL, exposure_clean))) %>%
  filter(!is.na(GeneSymbol), GeneSymbol != "") %>%
  distinct(GeneSymbol, method, exposure_clean, .keep_all = TRUE)

head(mr_sig_annot)   # 此时 SYMBOL / GENENAME / GeneSymbol 应显示基因符号
# ===== 5) 与 DEG 取交集（保留双边信息） =====
names(mr_sig_annot)
ivw_eqtl_results <- mr_sig_annot |>
  filter(method == "Inverse variance weighted") |>
  mutate(OR = exp(b), OR_low = exp(b - 1.96 * se),
         OR_high = exp(b + 1.96 * se),
         MR_significant = pval < 0.05)

eqtl_deg_intersection <- ivw_eqtl_results |>
  filter(MR_significant) |>
  inner_join(deg, by = "GeneSymbol")

readr::write_csv(eqtl_deg_intersection, file.path(output_dir, "eqtl_deg_intersection.csv"))
##########---------韦恩图-------------------------#############################
library(VennDiagram)
library(dplyr)
# ---------- 1. 读取数据 ----------

mr_res <- read.csv(file.path(output_dir, "eqtl_deg_intersection.csv"))
mr_sig <- mr_res %>%
  filter(MR_significant == TRUE) %>%
  pull(GeneSymbol) %>% toupper() %>% unique()

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
    filename = file.path(output_dir, "Venn_A_DEGup_MRsig.png"),
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
    filename = file.path(output_dir, "Venn_B_DEGdown_MRsig.png"),
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
print(mr_res$OR_high)
mr_dat <- tibble(
  exposure = c("CSE1L","ZNF239","KDELR3","LAMC2","FPR1","EPHB1","WDR72","EPB41L3","MYO15B","GPER1","FABP2","MS4A8"),
  nsnp     = c(2,4,2,2,3,3,2,4,2,3,2,2),
  method   = rep("IVW", 12),
  pval     = c(0.0023951966, 0.0431908240, 0.0248061900, 0.0478810707, 0.0187892830, 0.0082401633, 0.0192900218, 0.0339193943,0.0408254440, 0.0315525149, 0.0003294808, 0.0126803729),
  or       = c(0.7820784, 1.0524239, 1.0677817, 0.9162047, 0.8844206, 0.9463075, 0.7090681, 0.8244231, 1.1566779, 1.0761603,0.9039688, 1.1639122),
  or_lower = c(0.6673296, 1.0015645, 1.0083452, 0.8401155, 0.7982963, 0.9083475, 0.5316388, 0.6897095, 1.0060887, 1.0065086,0.8555025, 1.0329655),
  or_upper = c(0.9165584, 1.1058659, 1.1307217, 0.9991854, 0.9798365, 0.9858539, 0.9457128, 0.9854489, 1.3298070, 1.1506320,0.9551810, 1.3114588)
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
  dplyr::select(y, exposure, nsnp, method, pval_text) %>%
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
  geom_errorbar(aes(xmin = or_lower, xmax = or_upper),
                height = 0.12, linewidth = 0.6, orientation = "y") +
  geom_point(aes(x = or, color = col), size = 2.5) +
  scale_color_identity() +
  scale_x_continuous(trans = "log",
                     breaks = c(0.2, 0.6, 1.0, 1.4, 1.8,2.5),
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
  dplyr::select(y, `OR(95% CI)`) %>%
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

ggsave(filename = file.path(output_dir,"MR_forest_Fig5.png"),
       plot=p_final,
       width = 18, height = 7, dpi = 300)
# eQTL -> CRC：主 MR、异质性、多效性及敏感性分析 ----
# 输入：原始 eqtl_result.RData + 已有 eqtl_deg_intersection.csv 的显著基因名单。
# 输出：独立目录中的 Excel、完整分析对象及森林图/留一法图/漏斗图。
# 保留上游筛选口径，不重新选择工具变量，不清空用户工作区，不修改原文件。
rm(list = ls())
setwd("D:/Bioinfor/colorectal/MR/eQTL")
# 1. 参数与依赖 ----
project_dir <- "D:/Bioinfor/colorectal/MR/eQTL"
input_rdata <- file.path(project_dir, "eQTL/eqtl_result.RData")
selection_csv <- file.path(project_dir, "Mr_colorectal_results/eqtl_deg_intersection.csv")

outcome_id <- "ebi-a-GCST90013862"
presso_simulations <- 5000L
random_seed <- 20240918L
set.seed(random_seed)

packages <- c("dplyr", "purrr", "readr", "tibble", "TwoSampleMR", "MRPRESSO",
              "ggplot2", "ggrepel", "ggsci", "patchwork", "ragg", "writexl")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing packages: ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
source("./eqtl_plot_helpers.R", encoding = "UTF-8")

# 2. 独立环境加载，避免 load() 覆盖分析变量 ----
stopifnot(file.exists(input_rdata), file.exists(selection_csv))
input_environment <- new.env(parent = emptyenv())
load(input_rdata, envir = input_environment)
if (!exists("harmonised_data", envir = input_environment, inherits = FALSE)) {
  stop("RData does not contain harmonised_data")
}
harmonised_raw <- input_environment$harmonised_data
rm(input_environment)
selection <- readr::read_csv(selection_csv, show_col_types = FALSE, guess_max = Inf)
stopifnot(all(c("exposure", "GeneSymbol", "MR_significant") %in% names(selection)))
selection <- selection |>
  mutate(GeneSymbol = toupper(trimws(GeneSymbol)),
         ensg_clean = sub("\\.[0-9]+$", "", trimws(exposure)))
target_genes <- selection |>
  filter(MR_significant %in% TRUE, !is.na(GeneSymbol), nzchar(GeneSymbol)) |>
  distinct(GeneSymbol) |>
  pull(GeneSymbol)
if (!length(target_genes)) stop("No selected genes")

# 只保留映射所需两列；避免 exposure.x/exposure.y 和一对多合并膨胀。
gene_map <- selection |>
  filter(GeneSymbol %in% target_genes) |>
  distinct(ensg_clean, GeneSymbol)
if (anyDuplicated(gene_map$ensg_clean) || anyDuplicated(gene_map$GeneSymbol)) {
  stop("Ambiguous ENSG-to-symbol mapping; resolve before pooling instruments")
}
required <- c("id.exposure", "id.outcome", "SNP", "mr_keep", "beta.exposure",
              "se.exposure", "beta.outcome", "se.outcome")
if (!all(required %in% names(harmonised_raw))) stop("Missing harmonised columns")
harmonised_selected <- harmonised_raw |>
  mutate(ensg_clean = sub("\\.[0-9]+$", "", trimws(sub("^.*:", "", id.exposure)))) |>
  select(-any_of("GeneSymbol")) |>
  inner_join(gene_map, by = "ensg_clean", relationship = "many-to-one") |>
  filter(id.outcome == outcome_id, mr_keep %in% TRUE)

# 冲突的重复 SNP 不能直接保留首行；完全一致的分析记录才去重。
analysis_columns <- c(required, "GeneSymbol", "effect_allele.exposure",
                      "other_allele.exposure", "effect_allele.outcome", "other_allele.outcome")
duplicate_check <- harmonised_selected |>
  select(any_of(analysis_columns)) |>
  distinct() |>
  count(GeneSymbol, SNP) |>
  filter(n > 1)
if (nrow(duplicate_check)) stop("Conflicting duplicate gene/SNP records")
harmonised <- harmonised_selected |>
  distinct(GeneSymbol, SNP, .keep_all = TRUE) |>
  mutate(original_id_exposure = id.exposure,
         id.exposure = GeneSymbol, exposure = GeneSymbol,
         outcome = "Colorectal cancer") |>
  arrange(GeneSymbol, SNP)
missing_genes <- setdiff(target_genes, harmonised$GeneSymbol)
if (length(missing_genes)) stop("No instruments for: ", paste(missing_genes, collapse = ", "))
numeric_columns <- c("beta.exposure", "se.exposure", "beta.outcome", "se.outcome")
stopifnot(all(vapply(harmonised[numeric_columns], is.numeric, logical(1))))
invalid_rows <- harmonised |>
  filter(if_any(all_of(numeric_columns), ~ !is.finite(.x)),
         .preserve = FALSE)
if (nrow(invalid_rows) || any(harmonised$se.exposure <= 0 | harmonised$se.outcome <= 0) ||
    any(harmonised$beta.exposure == 0) || anyNA(harmonised$SNP) || any(!nzchar(harmonised$SNP))) {
  stop("Invalid beta, standard error, or SNP identifier")
}
overview <- harmonised |>
  count(GeneSymbol, name = "n_snp") |>
  mutate(primary_method = if_else(n_snp == 1, "Wald ratio", "Inverse variance weighted"),
         heterogeneity_eligible = n_snp >= 2,
         egger_eligible = n_snp >= 3, loo_eligible = n_snp >= 3,
         presso_eligible = n_snp >= 4)

# 3. 每项方法分别捕获错误，保留已完成的其他方法及警告 ----
run_checked <- function(gene, method, eligible, compute) {
  warnings <- character()
  if (!eligible) return(list(value = NULL, status = tibble(
    GeneSymbol = gene, analysis = method, status = "not_estimable",
    detail = "Insufficient SNP count", warnings = "")))
  result <- tryCatch(withCallingHandlers(compute(), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = identity)
  failed <- inherits(result, "error")
  list(value = if (failed) NULL else result,
       status = tibble(GeneSymbol = gene, analysis = method,
                       status = if (failed) "error" else "ok",
                       detail = if (failed) conditionMessage(result) else "",
                       warnings = paste(unique(warnings), collapse = " | ")))
}
run_gene <- function(gene_data) {
  gene <- gene_data$GeneSymbol[[1]]
  n_snp <- nrow(gene_data)
  methods <- if (n_snp == 1) "mr_wald_ratio" else if (n_snp == 2) "mr_ivw" else
    c("mr_ivw", "mr_egger_regression", "mr_weighted_median")
  het_methods <- if (n_snp >= 3) c("mr_ivw", "mr_egger_regression") else "mr_ivw"
  list(
    mr = run_checked(gene, "MR", TRUE, function() TwoSampleMR::mr(gene_data, method_list = methods)),
    het = run_checked(gene, "Heterogeneity", n_snp >= 2,
                      function() TwoSampleMR::mr_heterogeneity(gene_data, method_list = het_methods)),
    egger = run_checked(gene, "Egger_intercept", n_snp >= 3,
                        function() TwoSampleMR::mr_pleiotropy_test(gene_data)),
    loo = run_checked(gene, "Leave_one_out", n_snp >= 3,
                      function() TwoSampleMR::mr_leaveoneout(gene_data, method = TwoSampleMR::mr_ivw)),
    presso = run_checked(gene, "MR_PRESSO", n_snp >= 4, function() MRPRESSO::mr_presso(
      BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
      SdOutcome = "se.outcome", SdExposure = "se.exposure",
      data = as.data.frame(gene_data), OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
      NbDistribution = presso_simulations, SignifThreshold = 0.05, seed = random_seed))
  )
}
analysis <- purrr::map(split(harmonised, harmonised$GeneSymbol), run_gene)
collect_results <- function(method, columns) {
  empty <- as_tibble(setNames(lapply(columns, function(column) {
    if (column %in% c("exposure", "method", "SNP")) character() else numeric()
  }), columns))
  purrr::map_dfr(analysis, ~ .x[[method]]$value) |>
    bind_rows(empty)
}
mr <- collect_results("mr", c("exposure", "method", "b", "se", "pval"))
het <- collect_results("het", c("exposure", "method", "Q", "Q_df", "Q_pval"))
egger <- collect_results("egger", c("exposure", "egger_intercept", "se", "pval"))
loo <- collect_results("loo", c("exposure", "SNP", "b", "se", "p"))
status <- purrr::map_dfr(analysis, ~ purrr::map_dfr(.x, "status"))

# 4. MR-PRESSO 正确读取带空格的 Global Test；保留 <0.0002 等原始文本 ----
presso <- purrr::imap_dfr(analysis, function(gene_result, gene) {
  tests <- gene_result$presso$value[["MR-PRESSO results"]]
  global_p <- tests[["Global Test"]][["Pvalue"]]
  if (is.null(global_p)) global_p <- NA_character_
  tibble(GeneSymbol = gene, presso_global_p = as.character(global_p),
         presso_status = gene_result$presso$status$status)
})
primary <- mr |>
  filter(method %in% c("Inverse variance weighted", "Wald ratio")) |>
  transmute(GeneSymbol = exposure, beta = b, se, p_value = pval,
            OR = exp(b), OR_low = exp(b - 1.96 * se), OR_high = exp(b + 1.96 * se))
summary_table <- overview |>
  left_join(primary, by = "GeneSymbol", relationship = "one-to-one") |>
  left_join(het |> filter(method == "Inverse variance weighted") |>
              transmute(GeneSymbol = exposure, Q, Q_df, Q_p = Q_pval), by = "GeneSymbol") |>
  left_join(egger |> transmute(GeneSymbol = exposure, egger_intercept,
                               egger_se = se, egger_p = pval), by = "GeneSymbol") |>
  left_join(presso, by = "GeneSymbol")

# 保留原始精度；不在已筛选显著子集上重新计算全扫描 FDR。
paper_table <- summary_table |>
  transmute(Gene = GeneSymbol, SNPs = n_snp, Method = primary_method,
            OR_95CI = sprintf("%.3f (%.3f-%.3f)", OR, OR_low, OR_high),
            P = p_value, Q, Q_df, Q_P = Q_p,
            Egger_intercept = egger_intercept, Egger_P = egger_p,
            PRESSO_global_P = presso_global_p, PRESSO_status = presso_status)
audit <- tibble(stage = c("Raw harmonised rows", "Selected kept rows", "After deduplication", "Genes"),
                n = c(nrow(harmonised_raw), nrow(harmonised_selected), nrow(harmonised), nrow(overview)))
versions <- tibble(package = c("R", packages),
                   version = c(as.character(getRversion()), vapply(packages,
                                                                   function(pkg) as.character(utils::packageVersion(pkg)), character(1))))
writexl::write_xlsx(list(Summary = summary_table, Paper_table = paper_table,
                         Availability = overview, MR = mr, Heterogeneity = het, Egger = egger,
                         Leave_one_out = loo, MR_PRESSO = presso, Analysis_status = status,
                         Harmonised = harmonised, Selection = selection, Data_audit = audit, Versions = versions),
                    file.path(output_dir, "eQTL_MR_sensitivity_results.xlsx"))
# 嵌套分析对象保留完整 PRESSO 主效应、离群点、distortion 信息，便于复核。
saveRDS(analysis, file.path(output_dir, "eQTL_analysis_objects.rds"))
readr::write_excel_csv(summary_table, file.path(output_dir, "eQTL_MR_summary.csv"), na = "")

# 5. 验证与绘图 ----
if (any(status$status == "error")) {
  stop("Some analyses failed; completed results saved. See Analysis_status sheet.")
}
stopifnot(nrow(primary) == nrow(overview), all(is.finite(primary$beta)),
          all(is.finite(primary$OR_low)), all(is.finite(primary$OR_high)),
          all(primary$p_value >= 0 & primary$p_value <= 1),
          nrow(het) == sum(overview$n_snp >= 2) + sum(overview$n_snp >= 3),
          nrow(egger) == sum(overview$n_snp >= 3),
          nrow(loo) == sum(overview$n_snp[overview$n_snp >= 3] + 1))
make_eqtl_plots(summary_table, mr, loo, harmonised, output_dir)
message("Completed: ", nrow(overview), " genes; ", nrow(harmonised), " instrument rows; ",
        nrow(mr), " MR; ", nrow(het), " heterogeneity; ", nrow(egger), " Egger; ", nrow(loo), " LOO.")
message("Outputs: ", output_dir)
