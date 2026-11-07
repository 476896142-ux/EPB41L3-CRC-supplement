# 读取真实空间矩阵和坐标；所有三基因统计使用原始检出值，不插补。
suppressPackageStartupMessages({library(Matrix);library(Seurat);library(jsonlite);library(data.table)})
root <- "D:/Bioinfor/colorectal/GSE225857_spatial"
genes <- c("EPB41L3","KLB","ANGPTL1")
qc <- list(); detection <- list()
for(i in 1:4) {
  sample <- paste0("C",i); prefix <- paste0("GSM",7058755+i,"_",sample)
  file <- function(suffix) file.path(root,"raw",paste0(prefix,suffix))
  message("Reading ", sample)
  features <- read.delim(gzfile(file(".features.tsv.gz")),header=FALSE)
  barcodes <- read.delim(gzfile(file(".barcodes.tsv.gz")),header=FALSE)[,1]
  counts <- as(readMM(gzfile(file(".matrix.mtx.gz"))),"CsparseMatrix")
  stopifnot(nrow(counts)==nrow(features),ncol(counts)==length(barcodes))
  # 同名基因求和，保留真实基因名；固定条码顺序与坐标对应。
  symbols <- as.character(features[,2]); uniq <- unique(symbols)
  counts <- sparseMatrix(i=match(symbols,uniq),j=seq_along(symbols),x=1)%*%counts
  dimnames(counts) <- list(uniq,barcodes)
  pos <- read.csv(gzfile(file("_tissue_positions_list.csv.gz")),header=FALSE)
  names(pos) <- c("barcode","in_tissue","array_row","array_col","pixel_row","pixel_col")
  pos <- pos[match(barcodes,pos$barcode),];stopifnot(!anyNA(pos$barcode))
  pos$nCount <- Matrix::colSums(counts);pos$nFeature <- Matrix::colSums(counts>0)
  # FFPE 为探针测定，不套用单细胞线粒体百分比阈值。
  pos$keep <- pos$in_tissue==1 & pos$nCount>=500 & pos$nFeature>=200
  write.csv(pos,file.path(root,"results",paste0(sample,"_spot_QC.csv")),row.names=FALSE)
  qc[[sample]] <- data.frame(sample,available=ncol(counts),in_tissue=sum(pos$in_tissue==1),retained=sum(pos$keep),median_UMI=median(pos$nCount[pos$keep]),median_genes=median(pos$nFeature[pos$keep]))
  counts <- counts[,pos$keep,drop=FALSE];pos <- pos[pos$keep,];rownames(pos)<-pos$barcode
  obj <- CreateSeuratObject(counts,meta.data=pos,project=sample)
  obj <- NormalizeData(obj,verbose=FALSE)
  for(g in genes) {
    present <- g %in% rownames(counts)
    detection[[paste(sample,g)]] <- data.frame(sample,gene=g,measured=present,n_spots=ncol(counts),positive=if(present)sum(counts[g,]>0)else NA,pct_positive=if(present)100*mean(counts[g,]>0)else NA,mean_log=if(present)mean(LayerData(obj,layer="data")[g,])else NA)
  }
  sf <- fromJSON(paste(readLines(gzfile(file("_scalefactors_json.json.gz")),warn=FALSE),collapse=""))
  saveRDS(list(object=obj,scale=sf),file.path(root,"results",paste0(sample,"_spatial.rds")))
  rm(obj,counts);gc()
}
write.csv(rbindlist(qc),file.path(root,"results","QC_summary.csv"),row.names=FALSE)
write.csv(rbindlist(detection),file.path(root,"results","three_gene_detection.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(root,"results","sessionInfo.txt"))
