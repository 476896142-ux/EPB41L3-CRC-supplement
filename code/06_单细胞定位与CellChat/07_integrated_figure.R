# 整合单细胞定位与 CellChat；不改变细胞、表达值、通讯概率或筛选条件。
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
out <- "D:/Bioinfor/colorectal/single_raw/three_gene_analysis"
# 复用上一版通讯绘图对象，保持环图与气泡图的数据一致。
cc_env <- new.env()
sys.source(file.path(out,"06_cellchat_Huang.R"),envir=cc_env)
vis <- readRDS(file.path(out,"GSE132465_visualization.rds"))
genes <- c("EPB41L3","KLB","ANGPTL1")
xy <- as.data.frame(Embeddings(vis,"umap"));names(xy)<-c("x","y")
xy$celltype <- vis$Cell_type
expression <- as.matrix(GetAssayData(vis,layer="data")[genes,,drop=FALSE])
upper <- ceiling(max(expression)*10)/10
palette <- c("B cells"="#C58760","Epithelial cells"="#5F8FA8","Mast cells"="#B87393",
             "Myeloids"="#508D79","Stromal cells"="#8270A4","T cells"="#B5A24C")
base_theme <- theme_classic(base_size=10)+theme(plot.title=element_text(face="bold",hjust=.5),
  plot.margin=margin(4,5,4,5),legend.title=element_text(size=8),legend.text=element_text(size=8))
xrange <- range(xy$x);yrange<-range(xy$y)
umap <- function(p) p+coord_fixed(xlim=xrange,ylim=yrange)+
  scale_x_continuous(expand=expansion(mult=.025))+scale_y_continuous(expand=expansion(mult=.025))+
  labs(x="UMAP 1",y="UMAP 2")+base_theme
centres<-xy %>% group_by(celltype) %>% summarise(x=median(x),y=median(y),.groups="drop")
a<-umap(ggplot(xy,aes(x,y))+geom_point(aes(colour=celltype),size=.18)+
  scale_colour_manual(values=palette,guide="none")+
  ggrepel::geom_text_repel(data=centres,aes(label=celltype),size=2.7,seed=20261102,max.overlaps=Inf)+
  labs(title="Cellular landscape"))
features<-lapply(genes,function(g) {
  z<-xy;z$value<-as.numeric(expression[g,]);positive<-z[z$value>0,];positive<-positive[order(positive$value),]
  umap(ggplot()+geom_point(data=z[z$value==0,],aes(x,y),colour="#AEBFCB",size=.21)+
    geom_point(data=positive,aes(x,y,colour=value),size=.28)+
    scale_colour_gradientn(colours=c("#F3D18A","#E99C43","#CA553B","#84284F"),limits=c(0,upper),
      name="Expression (>0)",guide=guide_colourbar(barwidth=unit(22,"mm"),barheight=unit(2,"mm"),title.position="top"))+
    labs(title=g))+theme(legend.position="bottom")
})
tab<-read.csv(file.path(out,"three_gene_expression.csv")) %>% filter(dataset=="GSE132465")
e<-ggplot(tab,aes(factor(gene,levels=genes),Cell_type,size=pct_positive,colour=mean_expression))+
  geom_point()+facet_wrap(~Class)+scale_colour_viridis_c(option="C")+scale_size_area(max_size=7)+
  theme_bw(base_size=10)+theme(plot.margin=margin(5,5,5,5),legend.text=element_text(size=8),legend.title=element_text(size=9))+
  labs(x=NULL,y=NULL,size="Detected (%)",colour="Mean log\nexpression",title="Expression across cell types")
tag<-function(p,label)p+labs(tag=label)+theme(plot.tag=element_text(face="bold",size=14))
single<-wrap_plots(list(tag(a,"A"),tag(features[[1]],"B"),tag(features[[2]],"C"),tag(features[[3]],"D"),tag(e,"E")),
  design="AB\nCD\nEE",heights=c(1,1,.65))
f<-cc_env$a+coord_equal(xlim=c(-1.4,1.4),ylim=c(-1.2,1.2),clip="off")+theme(plot.margin=margin(2,5,2,5))
g<-cc_env$b+coord_equal(xlim=c(-1.4,1.4),ylim=c(-1.2,1.2),clip="off")+theme(plot.margin=margin(2,5,2,5))
h<-cc_env$c+theme(plot.margin=margin(3,5,3,5),axis.text.x=element_text(angle=45,hjust=1,size=9))
combined<-wrap_plots(list(tag(a,"A"),tag(features[[1]],"B"),tag(features[[2]],"C"),tag(features[[3]],"D"),
  tag(e,"E"),tag(f,"F"),tag(g,"G"),tag(h,"H")),design="ABCD\nEEEE\nFFGG\nHHHH",heights=c(3.3,1.8,5.4,5.9))
savefigure<-function(p,name,width,height) {
  ggsave(file.path(out,"figures",paste0(name,".pdf")),p,width=width,height=height,device=cairo_pdf,bg="white")
  ggsave(file.path(out,"figures",paste0(name,".png")),p,width=width,height=height,dpi=300,bg="white")
}
savefigure(single,"GSE132465_localization_recolored",12,11)
savefigure(combined,"Figure_single_cell_CellChat_integrated",16,16.4)
write.csv(data.frame(gene=genes,positive_cells=rowSums(expression>0),visualized_cells=ncol(expression),
  color_min=0,color_max=upper,zero_colour="#AEBFCB"),file.path(out,"integrated_figure_display_settings.csv"),row.names=FALSE)
