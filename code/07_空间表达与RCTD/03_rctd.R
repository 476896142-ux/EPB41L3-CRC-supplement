# RCTD 全模式适用于一个 Visium 点包含多种细胞的情况。
suppressPackageStartupMessages({library(Seurat);library(spacexr);library(SpatialExperiment);library(SummarizedExperiment);library(Matrix)})
set.seed(20261102)
root <- "D:/Bioinfor/colorectal/GSE225857_spatial"
ref <- readRDS(file.path(root,"results","single_cell_reference.rds"))
targets <- c("EPB41L3","KLB","ANGPTL1")
args <- commandArgs(trailingOnly=TRUE)
prepare_only <- "prepare" %in% args
samples <- intersect(paste0("C",1:4),args)
if(!length(samples))samples <- paste0("C",1:4)
for(sample in samples) {
  output <- file.path(root,"results",paste0(sample,"_RCTD.rds"))
  if(file.exists(output))next
  message("RCTD ",sample," ",Sys.time())
  z <- readRDS(file.path(root,"results",paste0(sample,"_spatial.rds")))
  scounts <- LayerData(z$object,layer="counts")
  common <- setdiff(intersect(rownames(scounts),rownames(ref$counts)),targets)
  # 排除三靶基因，防止表达定位与参考特征循环论证。
  reference <- SummarizedExperiment(assays=list(counts=ref$counts[common,]),colData=S4Vectors::DataFrame(cell_type=factor(ref$meta$type),nUMI=Matrix::colSums(ref$counts)))
  coords <- as.matrix(z$object[[]][,c("pixel_col","pixel_row")])
  spatial <- SpatialExperiment(assays=list(counts=scounts[common,]),colData=S4Vectors::DataFrame(nUMI=Matrix::colSums(scounts)),spatialCoords=coords)
  input <- createRctd(spatial,reference,UMI_min=500,ref_n_cells_min=25)
  saveRDS(input,file.path(root,"results",paste0(sample,"_RCTD_input.rds")))
  message("Prepared ",sample,"; model spots = ",ncol(input$spatial_experiment))
  if(prepare_only)next
  # 官方实现使用 fork 并行，Windows 下使用单核保证一致性。
  result <- runRctd(input,rctd_mode="full",max_cores=1)
  saveRDS(result,output)
  weights <- t(as.matrix(assay(result,"weights")))
  write.csv(weights,file.path(root,"results",paste0(sample,"_celltype_weights_raw.csv")))
  stopifnot(!anyNA(weights),all(rowSums(weights)>0),min(weights)>-1e-8)
  weights <- pmax(weights,0)/rowSums(pmax(weights,0))
  write.csv(weights,file.path(root,"results",paste0(sample,"_celltype_weights.csv")))
  message("Completed ",sample," ",Sys.time())
  rm(z,scounts,spatial,input,result,weights);gc()
}
