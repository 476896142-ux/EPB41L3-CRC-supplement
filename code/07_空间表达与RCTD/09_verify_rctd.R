# 实际拟合结果核验：坐标、条码、权重、收敛和靶基因排除。
suppressPackageStartupMessages({library(spacexr);library(SummarizedExperiment);library(data.table)})
root<-"D:/Bioinfor/colorectal/GSE225857_spatial"
summary<-list()
for(s in paste0("C",1:4)) {
  z<-readRDS(file.path(root,"results",paste0(s,"_RCTD.rds")))
  raw<-t(as.matrix(assay(z,"weights")))
  w<-as.matrix(read.csv(file.path(root,"results",paste0(s,"_celltype_weights.csv")),row.names=1,check.names=FALSE))
  cv<-read.csv(file.path(root,"results",paste0(s,"_RCTD_convergence.csv")))
  md<-metadata(z)
  genes<-md$internal_vars$gene_list_reg
  coord<-SpatialExperiment::spatialCoords(z)
  stopifnot(identical(rownames(w),rownames(raw)),identical(rownames(w),cv$barcode),!anyNA(w),!anyNA(coord),all(rowSums(w)>0),min(w)>=0,max(abs(rowSums(w)-1))<1e-8,!any(c("EPB41L3","KLB","ANGPTL1")%in%genes),ncol(w)==17)
  stopifnot(max(abs(w-pmax(raw,0)/rowSums(pmax(raw,0))))<1e-12)
  summary[[s]]<-data.frame(sample=s,spots=nrow(w),types=ncol(w),converged=sum(cv$converged),sigma=md$internal_vars$sigma,regression_genes=length(genes),raw_sum_min=min(rowSums(raw)),raw_sum_max=max(rowSums(raw)),normalized_sum_error=max(abs(rowSums(w)-1)),targets_excluded=TRUE)
}
summary<-rbindlist(summary);stopifnot(sum(summary$spots)==12430)
write.csv(summary,file.path(root,"rctd_results","validation_checks.csv"),row.names=FALSE)
print(summary)
