# 从既有单细胞原始计数构建患者均衡参考；保留作者细胞标签。
suppressPackageStartupMessages({library(Seurat);library(Matrix)})
set.seed(20261102)
root <- "D:/Bioinfor/colorectal/GSE225857_spatial"
source <- "D:/Bioinfor/colorectal/single_raw/results"
sc <- readRDS(file.path(source,"seurat_GSE132465.rds"))
meta <- read.csv(file.path(source,"annotation_GSE132465.csv"),row.names=1)
key <- function(x)gsub("-",".",x,fixed=TRUE)
idx <- match(key(colnames(sc)),key(rownames(meta)))
stopifnot(!anyNA(idx));meta<-meta[idx,];rownames(meta)<-colnames(sc)
# 仅以原发肿瘤细胞为参考；过少的罕见亚群并入作者的大类。
meta$type <- meta$Cell_type
fine <- meta$Cell_type %in% c("Stromal cells","Myeloids")
meta$type[fine] <- paste(meta$Cell_type[fine],meta$Cell_subtype[fine],sep=": ")
meta <- meta[meta$Class=="Tumor",]
tab <- table(meta$type); small <- names(tab)[tab<50]
meta$type[meta$type %in% small] <- meta$Cell_type[meta$type %in% small]
# 肿瘤肥大细胞仅 3 个，不足以估计参考表达谱。
meta <- meta[meta$type!="Mast cells",]
groups <- split(rownames(meta),interaction(meta$Patient,meta$type,drop=TRUE))
cells <- unlist(lapply(groups,function(x)if(length(x)>100)sample(x,100)else x),use.names=FALSE)
counts <- LayerData(sc,layer="counts")[,cells,drop=FALSE]
meta<-meta[cells,];rm(sc);gc()
saveRDS(list(counts=counts,meta=meta),file.path(root,"results","single_cell_reference.rds"))
write.csv(as.data.frame(table(meta$Patient,meta$type)),file.path(root,"results","reference_cell_counts.csv"),row.names=FALSE)
message("Reference ready: ",ncol(counts)," cells, ",length(unique(meta$type))," types")
