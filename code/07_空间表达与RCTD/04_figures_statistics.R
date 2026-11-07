# 组织图像上展示真实表达和 Seurat 匹配分数；切片分别汇总，不将空间点当作患者。
suppressPackageStartupMessages({library(Seurat);library(Matrix);library(ggplot2);library(patchwork);library(png);library(RANN);library(data.table)})
set.seed(20261102)
root <- "D:/Bioinfor/colorectal/GSE225857_spatial"
genes <- c("EPB41L3","KLB","ANGPTL1")
savefig <- function(p,name,w,h) {
  ggsave(file.path(root,"figures",paste0(name,".pdf")),p,width=w,height=h,device=cairo_pdf)
  ggsave(file.path(root,"figures",paste0(name,".png")),p,width=w,height=h,dpi=300,bg="white")
}
objects <- lapply(paste0("C",1:4),function(s)readRDS(file.path(root,"results",paste0(s,"_spatial.rds"))))
names(objects)<-paste0("C",1:4)
limits <- sapply(genes,function(g)max(sapply(objects,function(z)max(LayerData(z$object,layer="data")[g,]))))
safe_cor <- function(x,y)if(sd(x)>0 && sd(y)>0)cor(x,y,method="spearman")else NA_real_
baseplot <- function(s,z,opacity=1) {
  f<-file.path(root,"raw",paste0("GSM",7058755+as.integer(sub("C","",s)),"_",s,"_tissue_lowres_image.png.gz"))
  con<-gzfile(f,"rb");img<-readPNG(readBin(con,"raw",n=20000000));close(con)
  img<-img*opacity+(1-opacity)
  # Visium FFPE 标准捕获点直径为55微米，利用作者缩放参数换算比例尺。
  bar_length<-500/55*z$scale$spot_diameter_fullres*z$scale$tissue_lowres_scalef
  ggplot()+annotation_raster(img,xmin=0,xmax=dim(img)[2],ymin=-dim(img)[1],ymax=0)+
    annotate("rect",xmin=10,xmax=30+bar_length,ymin=-dim(img)[1]+8,ymax=-dim(img)[1]+49,fill="white",alpha=.9)+
    annotate("segment",x=20,xend=20+bar_length,y=-dim(img)[1]+20,yend=-dim(img)[1]+20,linewidth=.65)+
    annotate("text",x=20+bar_length/2,y=-dim(img)[1]+37,label="500 um",size=2.7)+
    coord_fixed(xlim=c(0,dim(img)[2]),ylim=c(-dim(img)[1],0),expand=FALSE)+theme_void(base_size=11)+
    theme(plot.title=element_text(face="bold",size=12),legend.key.height=grid::unit(.25,"cm"),legend.key.width=grid::unit(.8,"cm"),legend.title=element_text(size=9),legend.position="bottom",plot.margin=margin(2,2,2,2))
}
plots<-list();moran<-list();associations<-list();composition<-list();maps<-list();niches<-list()
for(s in names(objects)) {
  z<-objects[[s]];md<-z$object[[]];dat<-LayerData(z$object,layer="data")
  md$x<-md$pixel_col*z$scale$tissue_lowres_scalef;md$y<--md$pixel_row*z$scale$tissue_lowres_scalef
  plots[[paste(s,"HE")]]<-baseplot(s,z)+labs(title=paste(s,"H&E"))
  # 距离阈值保留真正相邻点，避免跨越组织空白连接。
  kn<-nn2(as.matrix(md[,c("array_row","array_col")]),k=7)
  ii<-rep(seq_len(nrow(md)),each=6);jj<-as.vector(t(kn$nn.idx[,-1]));dd<-as.vector(t(kn$nn.dists[,-1]))
  keep<-dd<=2.1;W<-sparseMatrix(i=ii[keep],j=jj[keep],x=1,dims=c(nrow(md),nrow(md)))
  W<-((W+t(W))>0)*1;diag(W)<-0
  weight_sum<-sum(W)
  for(g in genes) {
    md$value<-as.numeric(dat[g,]); pos<-md[md$value>0,]
    plots[[paste(s,g)]]<-baseplot(s,z,.3)+geom_point(data=md,aes(x,y),colour="#CBD0D5",size=.65,alpha=.65)+
      geom_point(data=pos,aes(x,y,colour=value),size=.9)+
      scale_colour_gradientn(colours=c("#F5D18C","#E88A46","#BC4359","#6D225C"),limits=c(0,limits[g]))+
      labs(title=paste(s,g),colour="Log expression")
    x<-md$value-mean(md$value)
    stat<-function(v)if(sum(v*v)>0)nrow(md)/weight_sum*sum(v*(W%*%v))/sum(v*v)else NA_real_
    obs<-stat(x);null<-replicate(999,stat(sample(x)))
    moran[[paste(s,g)]]<-data.frame(sample=s,gene=g,Moran_I=obs,p_permutation=if(is.finite(obs))(1+sum(null>=obs))/(1000)else NA)
  }
  wf<-file.path(root,"results",paste0(s,"_celltype_scores.csv"))
  if(!file.exists(wf))next
  weights<-as.matrix(read.csv(wf,row.names=1,check.names=FALSE));weights<-weights[match(rownames(md),rownames(weights)),,drop=FALSE]
  stopifnot(!anyNA(weights))
  for(g in genes)for(ct in colnames(weights)) {
    associations[[paste(s,g,ct)]]<-data.frame(sample=s,gene=g,cell_type=ct,rho=safe_cor(as.numeric(dat[g,]),weights[,ct]))
  }
  composition[[s]]<-data.frame(sample=s,cell_type=colnames(weights),mean_mapping_score=colMeans(weights))
  for(label in c("Myeloids","Stromal cells","Epithelial cells")) {
    md$value<-rowSums(weights[,startsWith(colnames(weights),label),drop=FALSE])
    maps[[paste(s,label)]]<-baseplot(s,z,.3)+geom_point(data=md,aes(x,y,colour=value),size=.9)+
      scale_colour_gradient(low="#EEF0F3",high=c(Myeloids="#C65D50",`Stromal cells`="#258F8A",`Epithelial cells`="#6B5D9C")[[label]],limits=c(0,1))+
      labs(title=paste(s,label),colour="Mapping score")
  }
  # 描述性邻域相关，不将其写为配体-受体因果证据。
  stromal<-rowSums(weights[,startsWith(colnames(weights),"Stromal cells"),drop=FALSE])
  myeloid<-rowSums(weights[,startsWith(colnames(weights),"Myeloids"),drop=FALSE])
  neigh<-as.numeric(W%*%myeloid)/pmax(Matrix::rowSums(W),1)
  niches[[s]]<-data.frame(sample=s,same_spot_rho=safe_cor(stromal,myeloid),neighbor_rho=safe_cor(stromal,neigh),n_spots=nrow(md))
}
savefig(wrap_plots(plots,ncol=4),"Figure_spatial_expression",16,15)
mt<-rbindlist(moran);mt$FDR<-p.adjust(mt$p_permutation,"BH")
write.csv(mt,file.path(root,"results","Moran_spatial_autocorrelation.csv"),row.names=FALSE)
if(length(associations)) {
  a<-rbindlist(associations);write.csv(a,file.path(root,"results","gene_celltype_correlations.csv"),row.names=FALSE)
  write.csv(rbindlist(composition),file.path(root,"results","celltype_mapping_summary.csv"),row.names=FALSE)
  write.csv(rbindlist(niches),file.path(root,"results","stromal_myeloid_neighborhood.csv"),row.names=FALSE)
  p<-ggplot(a,aes(gene,cell_type,fill=rho))+geom_tile(colour="white")+facet_wrap(~sample,nrow=1)+
    scale_fill_gradient2(low="#456DA5",mid="white",high="#C65D50",limits=c(-1,1),na.value="#DDDDDD")+theme_minimal(base_size=10)+
    theme(panel.grid=element_blank(),axis.text.x=element_text(angle=45,hjust=1))+labs(x=NULL,y=NULL,fill="Spearman rho")
  savefig(p,"Figure_celltype_associations",14,7)
  savefig(wrap_plots(maps,ncol=3),"Figure_spatial_celltypes",13,15)
}
