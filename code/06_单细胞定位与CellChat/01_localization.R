# 三基因定位：分别分析两队列；保留原作者注释，所有质控合格细胞用于汇总。
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(patchwork)
set.seed(20261102)
root <- "D:/Bioinfor/colorectal/single_raw"
out <- file.path(root, "three_gene_analysis")
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
genes <- c("EPB41L3", "KLB", "ANGPTL1")
save_plot <- function(p, name, w, h) {
  ggsave(file.path(out, "figures", paste0(name, ".pdf")), p, width=w, height=h, device=cairo_pdf)
  ggsave(file.path(out, "figures", paste0(name, ".png")), p, width=w, height=h, dpi=300)
}
all_summary <- list()
qc_summary <- list()
for (id in c("GSE132465", "GSE132257")) {
  message("Processing ", id)
  sc <- readRDS(file.path(root, "results", paste0("seurat_", id, ".rds")))
  # 按原始条码匹配，修复旧对象中横杠被转为点的问题；不依赖文件行顺序。
  meta <- read.csv(file.path(root,"results",paste0("annotation_",id,".csv")), row.names=1, check.names=FALSE)
  key <- function(x) gsub("-", ".", x, fixed=TRUE)
  stopifnot(!anyDuplicated(key(rownames(meta))))
  idx <- match(key(colnames(sc)), key(rownames(meta)))
  stopifnot(!anyNA(idx), all(genes %in% rownames(sc)))
  meta <- meta[idx,,drop=FALSE]; rownames(meta) <- colnames(sc)
  sc <- AddMetaData(sc, meta)
  sc$dataset <- id
  sc$percent.mt <- PercentageFeatureSet(sc, pattern="^MT-")
  before <- sc[[]]
  before$keep <- before$nFeature_RNA > 200 & before$nFeature_RNA < 6000 & before$percent.mt < 20 & before$Cell_type != "Unspecified"
  qc_summary[[id]] <- before %>% group_by(dataset,Patient,Sample,Class) %>% summarise(before=n(), after=sum(keep),.groups="drop")
  sc <- subset(sc, cells=rownames(before)[before$keep])
  sc <- NormalizeData(sc, verbose=FALSE)
  sc$celltype <- sc$Cell_type
  # 主队列保留髓系与基质亚群，避免把所有基质细胞直接标成 CAF。
  sc$communication_type <- sc$Cell_type
  if (id == "GSE132465") {
    fine <- sc$Cell_subtype
    use <- sc$Cell_type %in% c("Stromal cells","Myeloids")
    sc$communication_type[use] <- paste(sc$Cell_type[use], fine[use], sep=": ")
  }
  dat <- GetAssayData(sc, layer="data")[genes,,drop=FALSE]
  long <- bind_rows(lapply(genes,function(g) {
    data.frame(sc[[]],gene=g,expression=as.numeric(dat[g,]),check.names=FALSE)
  }))
  patient <- long %>% group_by(dataset,Patient,Class,Cell_type,gene) %>%
    summarise(n_cells=n(),mean_expression=mean(expression),pct_positive=100*mean(expression>0),.groups="drop")
  write.csv(patient,file.path(out,paste0(id,"_patient_expression.csv")),row.names=FALSE)
  bycell <- long %>% group_by(dataset,Class,Cell_type,gene) %>%
    summarise(n_cells=n(),patients=n_distinct(Patient),mean_expression=mean(expression),pct_positive=100*mean(expression>0),.groups="drop")
  all_summary[[id]] <- bycell
  if(id == "GSE132465") {
    subtype <- long %>% group_by(Class,Cell_type,Cell_subtype,gene) %>% summarise(n_cells=n(),patients=n_distinct(Patient),mean_expression=mean(expression),pct_positive=100*mean(expression>0),.groups="drop")
    write.csv(subtype,file.path(out,"GSE132465_subtype_expression.csv"),row.names=FALSE)
  }
  saveRDS(sc,file.path(out,paste0(id,"_qc_normalized.rds")),compress=FALSE)
  p <- ggplot(bycell,aes(gene,Cell_type,size=pct_positive,colour=mean_expression)) +
    geom_point() + facet_wrap(~Class) + scale_colour_viridis_c(option="C") +
    scale_size_area(max_size=9,limits=c(0,100)) + theme_bw(base_size=11) +
    labs(x=NULL,y=NULL,size="Detected (%)",colour="Mean log\nexpression",title=id)
  save_plot(p,paste0(id,"_expression_dotplot"),9,4)
  p_patient <- ggplot(patient %>% filter(n_cells>=20),aes(Cell_type,mean_expression,colour=Class)) +
    geom_boxplot(outlier.shape=NA,position=position_dodge(width=.7)) +
    geom_point(position=position_jitterdodge(jitter.width=.12,dodge.width=.7),size=1.2,alpha=.7) +
    facet_wrap(~gene,scales="free_y",ncol=1) + scale_colour_manual(values=c(Normal="#389B98",Tumor="#CE6657")) +
    theme_bw(base_size=11) + theme(axis.text.x=element_text(angle=25,hjust=1)) +
    labs(x=NULL,y="Mean log expression per donor",title=paste(id,"donor-level expression"))
  save_plot(p_patient,paste0(id,"_donor_expression"),10,9)
  # UMAP 仅用于定位展示；按样本及细胞类型抽样，每层最多 250 个，统计保留全量。
  groups <- split(colnames(sc),interaction(sc$Sample,sc$Cell_type,drop=TRUE))
  chosen <- unlist(lapply(groups,function(x) if(length(x)>250) sample(x,250) else x),use.names=FALSE)
  vis <- subset(sc,cells=chosen)
  rm(sc,long,dat);gc()
  vis <- FindVariableFeatures(vis,nfeatures=2000,verbose=FALSE)
  vis <- ScaleData(vis,features=VariableFeatures(vis),verbose=FALSE)
  vis <- RunPCA(vis,npcs=30,verbose=FALSE,seed.use=20261102)
  vis <- RunUMAP(vis,dims=1:30,seed.use=20261102,verbose=FALSE)
  p_umap <- DimPlot(vis,group.by="Cell_type",label=TRUE,repel=TRUE,raster=FALSE,pt.size=.12) + NoLegend() + ggtitle("Cell types")
  p_features <- FeaturePlot(vis,features=genes,order=TRUE,cols=c("#E7E9EC","#AF3D65"),raster=FALSE,pt.size=.12,combine=FALSE)
  p_all <- wrap_plots(c(list(p_umap),p_features),ncol=2) + plot_annotation(title=id,tag_levels="A")
  save_plot(p_all,paste0(id,"_localization"),12,10)
  saveRDS(vis,file.path(out,paste0(id,"_visualization.rds")))
  rm(vis);gc()
}
write.csv(bind_rows(qc_summary),file.path(out,"QC_by_sample.csv"),row.names=FALSE)
write.csv(bind_rows(all_summary),file.path(out,"three_gene_expression.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
