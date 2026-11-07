# 分阶段执行 spacexr full 模式，保留检查点及每点收敛信息。
# 计算函数与 runRctd/full 相同；分批仅用于恢复计算，不改变拟合参数。
suppressPackageStartupMessages({library(spacexr);library(SummarizedExperiment);library(Matrix)})
set.seed(20261102)
root<-"D:/Bioinfor/colorectal/GSE225857_spatial"
bulk_only<-"bulk"%in%commandArgs(trailingOnly=TRUE)
if(!bulk_only)stopifnot(spacexr:::verify_Q_cache())
samples<-intersect(paste0("C",1:4),commandArgs(trailingOnly=TRUE))
if(!length(samples))samples<-paste0("C",1:4)
for(s in samples) {
  # 每张切片固定独立种子，使并行分派和断点恢复不改变校准抽样。
  set.seed(20261102+match(s,paste0("C",1:4)))
  checkpoint<-file.path(root,"results",paste0(s,"_rctd_checkpoints"));dir.create(checkpoint,showWarnings=FALSE)
  input<-readRDS(file.path(root,"results",paste0(s,"_RCTD_input.rds")))
  stopifnot(!any(c("EPB41L3","KLB","ANGPTL1")%in%rownames(input$spatial_experiment)))
  config<-spacexr:::createRctdConfig(input,max_cores=1,max_multi_types=4,confidence_threshold=5,doublet_threshold=20)
  bulkfile<-file.path(checkpoint,"bulk.rds");sigmafile<-file.path(checkpoint,"calibrated.rds")
  if(file.exists(bulkfile))config<-readRDS(bulkfile)else {
    message(s," fitting platform effect ",Sys.time());config<-spacexr:::fitBulk(config);saveRDS(config,bulkfile)
  }
  if(bulk_only)next
  if(file.exists(sigmafile))config<-readRDS(sigmafile)else {
    message(s," selecting sigma ",Sys.time());config<-spacexr:::chooseSigmaC(config);saveRDS(config,sigmafile)
  }
  iv<-spacexr:::internal_vars(config);info<-spacexr:::cell_type_info(config)$renorm
  spacexr:::set_likelihood_vars(iv$Q_mat,iv$X_vals)
  sr<-spacexr:::spatialRNA(config);gl<-iv$gene_list_reg
  beads<-t(as.matrix(spacexr:::counts(sr)[gl,,drop=FALSE]));umi<-spacexr:::nUMI(sr)
  chunks<-split(seq_len(nrow(beads)),ceiling(seq_len(nrow(beads))/200))
  fits<-vector("list",length(chunks))
  for(b in seq_along(chunks)) {
    cachefile<-file.path(checkpoint,sprintf("batch_%03d.rds",b));idx<-chunks[[b]]
    if(file.exists(cachefile))fits[[b]]<-readRDS(cachefile)else {
      fits[[b]]<-spacexr:::decompose_batch(umi[idx],info[[1]],beads[idx,,drop=FALSE],gl,constrain=FALSE,max_cores=1,MIN.CHANGE=spacexr:::config(config)$MIN_CHANGE_REG)
      saveRDS(fits[[b]],cachefile)
    }
    message(s," batch ",b,"/",length(chunks)," ",Sys.time())
  }
  fits<-unlist(fits,recursive=FALSE)
  result<-spacexr:::create_spe_full(config,fits)
  raw<-t(as.matrix(assay(result,"weights")))
  # 官方 full 模式为非负但不强制总和为1；保留原始值后按点归一化。
  stopifnot(!anyNA(raw),all(rowSums(raw)>0),min(raw)>-1e-8)
  weights<-pmax(raw,0)/rowSums(pmax(raw,0))
  stopifnot(max(abs(rowSums(weights)-1))<1e-8)
  write.csv(raw,file.path(root,"results",paste0(s,"_celltype_weights_raw.csv")))
  write.csv(weights,file.path(root,"results",paste0(s,"_celltype_weights.csv")))
  write.csv(data.frame(barcode=rownames(raw),converged=vapply(fits,function(x)x$converged,logical(1)),raw_weight_sum=rowSums(raw)),file.path(root,"results",paste0(s,"_RCTD_convergence.csv")),row.names=FALSE)
  saveRDS(result,file.path(root,"results",paste0(s,"_RCTD.rds")))
  message("Completed ",s,": ",nrow(weights)," spots; convergence = ",mean(vapply(fits,function(x)x$converged,logical(1))))
  rm(input,config,iv,info,sr,beads,fits,result,raw,weights);gc()
}
