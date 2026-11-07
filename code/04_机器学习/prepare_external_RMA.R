# 外部数据独立 RMA；不按表达或差异阈值删掉探针，不读取结局参与归一化。
suppressPackageStartupMessages({library(affy);library(hgu133plus2cdf)})
project <- 'D:/Bioinfor/colorectal/GEO_ML_hub10'
out <- 'E:/codex/2026-10-31/dui/outputs/GEO_ML_hub10'
files <- sort(list.files(file.path(project,'raw/GSE20916_CEL'),pattern='CEL.gz$',full.names=TRUE))
stopifnot(length(files)==60)
raw <- ReadAffy(filenames=files)
normalized <- exprs(rma(raw))
colnames(normalized) <- sub('.CEL.gz','',basename(files),fixed=TRUE)
saveRDS(normalized,file.path(out,'GSE20916_RMA_probe_expression.rds'))
ann <- read.csv(file.path(project,'raw/hgu133plus2_unique_probe_map.csv'))
genes <- c('NTRK3','EPHB1','KLB','GREM1','HHIP','ANGPTL1','LAMC2','WDR72','GPER1','EPB41L3')
gene_expr <- vapply(genes,function(g) {
 probes <- intersect(ann$PROBEID[ann$SYMBOL==g],rownames(normalized));stopifnot(length(probes)>0)
 apply(normalized[probes,,drop=FALSE],2,median)
},numeric(ncol(normalized)))
rownames(gene_expr) <- colnames(normalized)
write.csv(gene_expr,file.path(out,'GSE20916_expression_RMA.csv'))
cat('外部 RMA 完成：',nrow(gene_expr),'样本，',ncol(gene_expr),'基因。\n')
