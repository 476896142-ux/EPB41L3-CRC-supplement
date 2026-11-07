# 全基因 GO/KEGG GSEA ----
# 输入：全量 limma 表；不按 P 值、FDR 或 |logFC| 预筛选。
# 中文注释；原始脚本和显著 DEG 结果不覆盖。
# clusterProfiler 4.20.0 的 gseGO/gseKEGG 未转发 eps/seed 等参数，
# 因此直接调用同一 enrichit 后端，确保 eps=0 和种子实际生效。
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(GO.db)
  library(enrichit)
})
input_file <- "D:/Bioinfor/colorectal/results/GSE196006_GSE146009_GSE110224/DEG_Tumor_vs_Normal_limma.csv"
output_dir <- "E:/codex/2026-10-31/dui/outputs/GSEA_all_genes"
annotation_dir <- file.path(output_dir, "annotations")
dir.create(annotation_dir, recursive = TRUE, showWarnings = FALSE)
random_seed <- 20261031L
set.seed(random_seed)

# 1. 整理排序向量 ----
expression_data <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("GeneSymbol", "logFC") %in% names(expression_data)))
input_rows <- nrow(expression_data)
valid_rows <- is.finite(expression_data$logFC) & !is.na(expression_data$GeneSymbol) &
  nzchar(trimws(expression_data$GeneSymbol))
write.csv(expression_data[!valid_rows, ], file.path(output_dir, "excluded_invalid_rows.csv"), row.names = FALSE)
expression_data <- expression_data[valid_rows, ]
expression_data$GeneSymbol <- trimws(expression_data$GeneSymbol)
# 多行同名基因保留绝对 logFC 最大的一行；同分按名称排序。
expression_data <- expression_data[order(-abs(expression_data$logFC), expression_data$GeneSymbol), ]
duplicate_symbols <- expression_data[duplicated(expression_data$GeneSymbol), ]
write.csv(duplicate_symbols, file.path(output_dir, "excluded_duplicate_symbols.csv"), row.names = FALSE)
expression_data <- expression_data[!duplicated(expression_data$GeneSymbol), ]
gene_mapping <- unique(AnnotationDbi::select(org.Hs.eg.db,
  keys = expression_data$GeneSymbol, keytype = "SYMBOL", columns = "ENTREZID"))
unmapped_symbols <- unique(gene_mapping$SYMBOL[is.na(gene_mapping$ENTREZID)])
write.csv(data.frame(GeneSymbol = unmapped_symbols), file.path(output_dir, "unmapped_symbols.csv"), row.names = FALSE)
gene_mapping <- gene_mapping[!is.na(gene_mapping$ENTREZID), ]
mapping_counts <- table(gene_mapping$SYMBOL)
ambiguous_symbols <- names(mapping_counts)[mapping_counts > 1L]
write.csv(gene_mapping[gene_mapping$SYMBOL %in% ambiguous_symbols, ],
  file.path(output_dir, "excluded_ambiguous_symbols.csv"), row.names = FALSE)
# 一对多映射予以剔除，避免把一个表达值重复计入多个基因。
gene_mapping <- gene_mapping[!gene_mapping$SYMBOL %in% ambiguous_symbols, ]
rank_data <- merge(expression_data[, c("GeneSymbol", "logFC")], gene_mapping,
  by.x = "GeneSymbol", by.y = "SYMBOL", sort = FALSE)
rank_data <- rank_data[order(-abs(rank_data$logFC), rank_data$GeneSymbol, rank_data$ENTREZID), ]
duplicate_entrez <- rank_data[duplicated(rank_data$ENTREZID), ]
write.csv(duplicate_entrez, file.path(output_dir, "excluded_duplicate_entrez.csv"), row.names = FALSE)
rank_data <- rank_data[!duplicated(rank_data$ENTREZID), ]
rank_data <- rank_data[order(-rank_data$logFC, rank_data$ENTREZID), ]
gene_list <- setNames(rank_data$logFC, rank_data$ENTREZID)
stopifnot(length(gene_list) > 1000, all(is.finite(gene_list)),
  !anyDuplicated(names(gene_list)), all(diff(gene_list) <= 0))
write.csv(rank_data, file.path(output_dir, "all_genes_ranked.csv"), row.names = FALSE)
mapping_summary <- data.frame(input_rows = input_rows, invalid_rows = sum(!valid_rows),
  duplicate_symbols = nrow(duplicate_symbols), unmapped_symbols = length(unmapped_symbols),
  ambiguous_symbols = length(ambiguous_symbols), duplicate_entrez_rows = nrow(duplicate_entrez),
  ranked_genes = length(gene_list), positive_logFC = sum(gene_list > 0),
  negative_logFC = sum(gene_list < 0), zero_logFC = sum(gene_list == 0), tied_logFC = sum(duplicated(gene_list)))
write.csv(mapping_summary, file.path(output_dir, "mapping_summary.csv"), row.names = FALSE)
message("排序完成：", length(gene_list), " 个基因；开始准备数据库。")

# 2. 冻结 GO 和 KEGG 注释，使后续重跑不受在线数据库更新影响 ----
go_snapshot <- file.path(annotation_dir, "GO_ALL_gson.rds")
if (!file.exists(go_snapshot)) {
  saveRDS(clusterProfiler::gson_GO(org.Hs.eg.db, keytype = "ENTREZID", ont = "ALL"), go_snapshot)
}
go_annotation <- readRDS(go_snapshot)
link_file <- file.path(annotation_dir, "kegg_pathway_genes.tsv")
name_file <- file.path(annotation_dir, "kegg_pathway_names.tsv")
if (!file.exists(link_file)) download.file("https://rest.kegg.jp/link/hsa/pathway", link_file, mode = "wb")
if (!file.exists(name_file)) download.file("https://rest.kegg.jp/list/pathway/hsa", name_file, mode = "wb")
kegg_links <- read.delim(link_file, header = FALSE, col.names = c("gsid", "gene"))
kegg_links$gsid <- sub("^path:", "", kegg_links$gsid)
kegg_links$gene <- sub("^hsa:", "", kegg_links$gene)
kegg_names <- read.delim(name_file, header = FALSE, col.names = c("gsid", "name"), quote = "")
kegg_names$gsid <- sub("^path:", "", kegg_names$gsid)
kegg_names$name <- sub(" - Homo sapiens \\(human\\)$", "", kegg_names$name)
kegg_annotation <- gson::gson(gsid2gene = unique(kegg_links), gsid2name = kegg_names,
  species = "Homo sapiens", gsname = "KEGG", keytype = "ENTREZID",
  accessed_date = as.character(as.Date(file.info(link_file)$mtime)))
saveRDS(kegg_annotation, file.path(annotation_dir, "KEGG_gson.rds"))
write.csv(AnnotationDbi::metadata(org.Hs.eg.db), file.path(annotation_dir, "org_Hs_eg_db_metadata.csv"), row.names = FALSE)

# 3. 全部检验后分别对 GO（BP/CC/MF 合并）与 KEGG 做 BH 校正 ----
run_gsea <- function(annotation, collection) {
  message("运行 ", collection, " 全基因 GSEA……")
  set.seed(random_seed)
  result <- enrichit::gsea_gson(geneList = gene_list, gson = annotation,
    exponent = 1, minGSSize = 10, maxGSSize = 500,
    pvalueCutoff = NULL, pAdjustMethod = "BH", method = "multilevel",
    eps = 0, seed = random_seed, nPermSimple = 10000, sampleSize = 101,
    scoreType = "std", verbose = TRUE)
  if (is.null(result) || !nrow(result@result)) stop(collection, " 未产生结果。")
  result@setType <- collection
  result@params$method <- "multilevel"
  result@params$eps <- 0
  result@params$seed <- random_seed
  result@params$nPermSimple <- 10000L
  result@params$sampleSize <- 101L
  result@params$scoreType <- "std"
  result@params$pvalueCutoff <- 1
  table <- result@result
  if (collection == "GO") {
    ontology <- AnnotationDbi::Ontology(GO.db::GOTERM)
    table$ONTOLOGY <- unname(ontology[table$ID])
    stopifnot(!anyNA(table$ONTOLOGY))
  }
  table$Count <- vapply(strsplit(table$core_enrichment, "/", fixed = TRUE), function(x) sum(nzchar(x)), integer(1))
  table$GeneRatio <- table$Count / table$setSize
  table <- table[order(table$p.adjust, -abs(table$NES), table$ID), ]
  result@result <- table
  write.csv(table, file.path(output_dir, paste0(collection, "_GSEA_all_tested.csv")), row.names = FALSE)
  saveRDS(result, file.path(output_dir, paste0(collection, "_GSEA.rds")))
  if (any(!is.finite(table$pvalue) | !is.finite(table$NES))) {
    stop(collection, " 存在未成功估计的 P 值或 NES；已保存诊断结果，停止下游绘图。")
  }
  stopifnot(all(table$pvalue >= 0 & table$pvalue <= 1),
    isTRUE(all.equal(table$p.adjust, p.adjust(table$pvalue, "BH"), tolerance = 1e-12)),
    all(table$setSize >= 10 & table$setSize <= 500), all(table$Count <= table$setSize))
  significant <- table[table$p.adjust < 0.05, ]
  write.csv(significant, file.path(output_dir, paste0(collection, "_GSEA_FDR05.csv")), row.names = FALSE)
  message(collection, "：检验 ", nrow(table), " 个条目，FDR < 0.05：", nrow(significant))
  result
}
gsea_kegg <- run_gsea(kegg_annotation, "KEGG")
gsea_go <- run_gsea(go_annotation, "GO")

# 4. 最小可重复性记录，不输出持续运行日志 ----
versions <- data.frame(package = c("R", "clusterProfiler", "enrichit", "org.Hs.eg.db", "GO.db", "enrichplot", "ggplot2"),
  version = c(as.character(getRversion()), vapply(c("clusterProfiler", "enrichit", "org.Hs.eg.db", "GO.db", "enrichplot", "ggplot2"),
    function(p) as.character(packageVersion(p)), character(1))))
write.csv(versions, file.path(output_dir, "software_versions.csv"), row.names = FALSE)
provenance_files <- c(input_file, go_snapshot, link_file, name_file)
write.csv(data.frame(file = provenance_files, md5 = unname(tools::md5sum(provenance_files))),
  file.path(output_dir, "input_annotation_checksums.csv"), row.names = FALSE)
source(file.path(output_dir, "GSEA_all_genes_plots.R"), encoding = "UTF-8")
message("分析与绘图完成：", output_dir)
