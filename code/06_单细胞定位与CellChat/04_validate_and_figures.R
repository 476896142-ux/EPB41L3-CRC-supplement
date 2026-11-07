# 输出可用于论文的组合图、原始计数核验及标志基因核对。
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(patchwork)
out <- "D:/Bioinfor/colorectal/single_raw/three_gene_analysis"
genes <- c("EPB41L3","KLB","ANGPTL1")
colors <- c("B cells"="#C58760","Epithelial cells"="#5F8FA8","Mast cells"="#B87393",
            "Myeloids"="#508D79","Stromal cells"="#8270A4","T cells"="#B5A24C")
savefig <- function(p,name,w,h) {
  ggsave(file.path(out,"figures",paste0(name,".pdf")),p,width=w,height=h,device=cairo_pdf)
  ggsave(file.path(out,"figures",paste0(name,".png")),p,width=w,height=h,dpi=300)
}
checks<-list()
for(id in c("GSE132465","GSE132257")) {
  sc<-readRDS(file.path(out,paste0(id,"_qc_normalized.rds")))
  raw<-read.csv(file.path(out,paste0(id,"_raw_three_genes.csv")),check.names=FALSE)
  key<-function(x)gsub("-",".",x,fixed=TRUE)
  idx<-match(key(colnames(sc)),key(raw$cell));stopifnot(!anyNA(idx))
  actual<-as.matrix(GetAssayData(sc,layer="counts")[genes,,drop=FALSE])
  expected<-t(as.matrix(raw[idx,genes]))
  stopifnot(all(actual==expected))
  checks[[id]]<-data.frame(dataset=id,cells=ncol(sc),target_UMI_exact_match=TRUE,
    mt_genes=sum(grepl("^MT-",rownames(sc))),max_mt=max(sc$percent.mt),
    min_features=min(sc$nFeature_RNA),max_features=max(sc$nFeature_RNA))
  markers<-intersect(c("EPCAM","KRT8","COL1A1","DCN","PECAM1","VWF","LST1","TYROBP","CD3D","MS4A1","JCHAIN","TPSAB1"),rownames(sc))
  p<-DotPlot(sc,features=markers,group.by="Cell_type",dot.scale=7)+RotatedAxis()+
    scale_colour_viridis_c(option="C")+theme_bw(base_size=10)+labs(x=NULL,y=NULL,title=paste(id,"annotation markers"))
  savefig(p,paste0(id,"_annotation_markers"),10,4)
  if(id=="GSE132257") {
    x<-as.matrix(GetAssayData(sc,layer="data")[genes,,drop=FALSE])
    z<-bind_rows(lapply(genes,function(g)data.frame(sc[[]],gene=g,value=x[g,]))) %>%
      group_by(Patient,Status,Class,Cell_type,gene) %>%
      summarise(n_cells=n(),mean_expression=mean(value),pct_positive=100*mean(value>0),.groups="drop")
    write.csv(z,file.path(out,"GSE132257_processing_strata.csv"),row.names=FALSE)
  }
  rm(sc);gc()
  vis<-readRDS(file.path(out,paste0(id,"_visualization.rds")))
  a<-DimPlot(vis,group.by="Cell_type",cols=colors,label=TRUE,repel=TRUE,raster=FALSE,pt.size=.1)+
    NoLegend()+labs(title="Cellular landscape",x="UMAP 1",y="UMAP 2")
  features<-FeaturePlot(vis,features=genes,order=TRUE,cols=c("#E6E8EB","#AE3F6A"),raster=FALSE,pt.size=.12,combine=FALSE)
  features<-lapply(features,function(p)p+labs(x="UMAP 1",y="UMAP 2"))
  tab<-read.csv(file.path(out,"three_gene_expression.csv")) %>% filter(dataset==id)
  e<-ggplot(tab,aes(factor(gene,levels=genes),Cell_type,size=pct_positive,colour=mean_expression))+
    geom_point()+facet_wrap(~Class)+scale_colour_viridis_c(option="C")+scale_size_area(max_size=8)+
    theme_bw(base_size=11)+labs(x=NULL,y=NULL,size="Detected (%)",colour="Mean log\nexpression",title="Expression across cell types")
  p<-wrap_plots(c(list(a),features,list(e)),design="AB\nCD\nEE",heights=c(1,1,.65))+
    plot_annotation(title=paste(id,"| Three-gene cellular localization"),tag_levels="A")
  savefig(p,paste0(id,"_Figure_single_cell"),13,13)
  rm(vis);gc()
}
write.csv(bind_rows(checks),file.path(out,"cache_validation.csv"),row.names=FALSE)
