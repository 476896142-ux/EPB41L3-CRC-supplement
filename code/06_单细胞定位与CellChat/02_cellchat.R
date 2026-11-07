# 在主队列内分别推断肿瘤、正常通讯；不将技术重复队列混入主网络。
library(Seurat)
library(CellChat)
library(dplyr)
library(ggplot2)
library(patchwork)
future::plan("sequential")
set.seed(20261102)
out <- "D:/Bioinfor/colorectal/single_raw/three_gene_analysis"
dir.create(file.path(out,"cellchat"),showWarnings=FALSE)
sc <- readRDS(file.path(out,"GSE132465_qc_normalized.rds"))
# 使用已发布的亚群注释；SPP1+ 与 Pro-inflammatory 保留原名，不直接等同于 M2/M1。
meta <- sc[[]]
meta$group <- meta$Cell_type
stromal <- meta$Cell_type == "Stromal cells"
meta$group[stromal] <- meta$Cell_subtype[stromal]
myeloid <- meta$Cell_type == "Myeloids"
meta$group[myeloid] <- paste0("Myeloid ",meta$Cell_subtype[myeloid])
expr <- GetAssayData(sc,layer="data")
rm(sc);gc()
palette <- setNames(grDevices::hcl.colors(length(unique(meta$group)),"Dark 3"),sort(unique(meta$group)))
plots <- list(); summaries <- list(); max_strength <- 0
for (condition in c("Normal","Tumor")) {
  path <- file.path(out,"cellchat",paste0(condition,"_cellchat.rds"))
  if (!file.exists(path)) {
    m <- meta[meta$Class == condition,,drop=FALSE]
    # 按患者×细胞群分层抽样控制内存，避免高细胞数患者占绝大多数。
    cell_file <- file.path(out,"cellchat",paste0(condition,"_cells.csv"))
    strata <- split(rownames(m),interaction(m$Patient,m$group,drop=TRUE))
    # 固定已保存的细胞列表，使从中断处恢复也能复现同一网络。
    keep <- if(file.exists(cell_file)) rownames(read.csv(cell_file,row.names=1)) else
      unlist(lapply(strata,function(x) if(length(x)>200) sample(x,200) else x),use.names=FALSE)
    m <- m[keep,,drop=FALSE]
    counts <- table(m$group)
    m <- m[m$group %in% names(counts)[counts>=10],,drop=FALSE]
    m$group <- droplevels(factor(m$group))
    write.csv(m,file.path(out,"cellchat",paste0(condition,"_cells.csv")))
    cc <- createCellChat(object=expr[,rownames(m)],meta=m,group.by="group")
    cc@DB <- CellChatDB.human
    cc <- subsetData(cc)
    cc <- identifyOverExpressedGenes(cc)
    cc <- identifyOverExpressedInteractions(cc)
    # raw.use=TRUE 指未投影的标准化数据，不是 UMI 原始计数。
    cc <- computeCommunProb(cc,type="triMean",raw.use=TRUE,population.size=FALSE,nboot=100,seed.use=20261102)
    cc <- filterCommunication(cc,min.cells=10)
    cc <- computeCommunProbPathway(cc)
    cc <- aggregateNet(cc)
    saveRDS(cc,path)
  } else cc <- readRDS(path)
  pairs <- subsetCommunication(cc)
  write.csv(pairs,file.path(out,"cellchat",paste0(condition,"_LR_pairs.csv")),row.names=FALSE)
  summaries[[condition]] <- data.frame(condition=condition,cells=length(cc@idents),groups=nlevels(cc@idents),LR_edges=nrow(pairs),pathways=length(cc@netP$pathways))
  mat <- as.data.frame(as.table(cc@net$weight));names(mat)<-c("Sender","Receiver","Strength")
  max_strength <- max(max_strength,mat$Strength)
  plots[[condition]] <- ggplot(mat,aes(Receiver,Sender,fill=Strength)) + geom_tile(colour="white",linewidth=.15) +
    scale_fill_viridis_c(option="C") + coord_equal() + theme_minimal(base_size=9) +
    theme(axis.text.x=element_text(angle=60,hjust=1),panel.grid=element_blank()) + labs(title=condition,x=NULL,y=NULL)
  pathway <- data.frame(pathway=cc@netP$pathways,strength=apply(cc@netP$prob,3,sum))
  write.csv(pathway,file.path(out,"cellchat",paste0(condition,"_pathway_strength.csv")),row.names=FALSE)
  rm(cc);gc()
}
write.csv(bind_rows(summaries),file.path(out,"cellchat","summary.csv"),row.names=FALSE)
plots <- lapply(plots,function(p)p+scale_fill_viridis_c(option="C",limits=c(0,max_strength)))
p <- wrap_plots(plots,ncol=2) + plot_annotation(title="Predicted cell-cell communication",tag_levels="A")
ggsave(file.path(out,"figures","CellChat_networks.pdf"),p,width=17,height=8,device=cairo_pdf)
ggsave(file.path(out,"figures","CellChat_networks.png"),p,width=17,height=8,dpi=300)
# 数据库覆盖核对：包括受体复合物和辅助因子，不把没有边误写成没有功能。
db <- CellChatDB.human
coverage <- bind_rows(lapply(c("EPB41L3","KLB","ANGPTL1"),function(g) {
  bind_rows(lapply(names(db),function(part) {
    z <- db[[part]]
    if(!is.data.frame(z)) return(NULL)
    hit <- apply(z,1,function(x) any(as.character(x)==g,na.rm=TRUE))
    if(!any(hit)) return(NULL)
    data.frame(gene=g,table=part,row=rownames(z)[hit])
  }))
}))
write.csv(coverage,file.path(out,"cellchat","candidate_database_coverage.csv"),row.names=FALSE)
