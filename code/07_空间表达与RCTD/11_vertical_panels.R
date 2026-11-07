# 为上下组图重新标注底部面板；仅使用既有结果表，不重新分析。
suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
root<-"D:/Bioinfor/colorectal/GSE225857_spatial"
out<-file.path(root,"combined_figures");dir.create(out,showWarnings=FALSE)
genes<-c("EPB41L3","KLB","ANGPTL1")
d<-read.csv(file.path(root,"results","three_gene_detection.csv"))
m<-read.csv(file.path(root,"results","Moran_spatial_autocorrelation.csv"))
d$gene<-factor(d$gene,levels=genes);m$gene<-factor(m$gene,levels=genes)
p1<-ggplot(d,aes(sample,gene,size=pct_positive,fill=pct_positive))+geom_point(shape=21,colour="#52606D")+scale_size_area(max_size=14)+scale_fill_gradient(low="#F5D18C",high="#84284F")+theme_bw(base_size=11)+labs(x=NULL,y=NULL,size="Detected (%)",fill="Detected (%)",title="B  Spatial detection")
p2<-ggplot(m,aes(sample,gene,fill=Moran_I))+geom_tile(colour="white")+geom_text(aes(label=paste0(sprintf("%.3f",Moran_I),ifelse(FDR<.05," *",""))),size=3.5)+scale_fill_gradient2(low="#456DA5",mid="white",high="#C65D50",midpoint=0)+theme_minimal(base_size=11)+theme(panel.grid=element_blank())+labs(x=NULL,y=NULL,fill="Moran's I",title="C  Spatial autocorrelation",caption="* BH-adjusted permutation P < 0.05")
ggsave(file.path(out,"S_bottom_panels.pdf"),p1+p2,width=12,height=3.5,device=cairo_pdf)
# 在独立环境读取原有热图代码，仅执行绘图对象构建，避免覆盖旧图。
env<-new.env();code<-readLines(file.path(root,"scripts","10_combine_R2_R3.R"),warn=FALSE)
code<-code[!grepl("^ggsave|^stopifnot",code)]
eval(parse(text=code),envir=env)
p_r<-(env$p_a+labs(title="B  Gene-cell associations"))+(env$p_b+labs(title="C  Marker-expression concordance"))+plot_layout(widths=c(1.5,1))
ggsave(file.path(root,"combined_figures","R_bottom_panels.pdf"),p_r,width=12,height=6,device=cairo_pdf)
