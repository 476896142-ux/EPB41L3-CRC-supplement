# 汇总 RCTD 相对权重，并核对其与细胞标志表达的一致性。
suppressPackageStartupMessages({library(Seurat);library(Matrix);library(ggplot2);library(patchwork);library(png);library(data.table);library(RANN)})
root<-"D:/Bioinfor/colorectal/GSE225857_spatial"
out<-file.path(root,"rctd_results");dir.create(out,showWarnings=FALSE)
genes<-c("EPB41L3","KLB","ANGPTL1")
groups<-c(Epithelial="Epithelial cells",Stromal="Stromal cells",Myeloid="Myeloids",Smooth_muscle="Stromal cells: Smooth muscle cells",SPP1_myeloid="Myeloids: SPP1+")
safe_cor<-function(x,y)if(sd(x)>0&&sd(y)>0)cor(x,y,method="spearman")else NA_real_
savefig<-function(p,n,w,h){ggsave(file.path(out,paste0(n,".pdf")),p,width=w,height=h,device=cairo_pdf);ggsave(file.path(out,paste0(n,".png")),p,width=w,height=h,dpi=300,bg="white")}
plots<-list();quality<-list();corr<-list();proportions<-list();marker_check<-list();neighbor<-list()
for(s in paste0("C",1:4)) {
  z<-readRDS(file.path(root,"results",paste0(s,"_spatial.rds")));md<-z$object[[]];dat<-LayerData(z$object,layer="data")
  w<-as.matrix(read.csv(file.path(root,"results",paste0(s,"_celltype_weights.csv")),row.names=1,check.names=FALSE));w<-w[rownames(md),,drop=FALSE]
  stopifnot(!anyNA(w),all(abs(rowSums(w)-1)<1e-8),min(w)>=0)
  conv<-read.csv(file.path(root,"results",paste0(s,"_RCTD_convergence.csv")));conv<-conv[match(rownames(md),conv$barcode),]
  quality[[s]]<-data.frame(sample=s,spots=nrow(w),converged=sum(conv$converged),convergence_rate=mean(conv$converged),min_weight=min(w),max_weight=max(w),max_sum_error=max(abs(rowSums(w)-1)))
  proportions[[s]]<-data.frame(sample=s,cell_type=colnames(w),mean_weight=colMeans(w))
  broad<-sapply(groups,function(prefix)rowSums(w[,startsWith(colnames(w),prefix),drop=FALSE]))
  write.csv(broad,file.path(out,paste0(s,"_broad_weights.csv")))
  marker<-as.matrix(read.csv(file.path(root,"results",paste0(s,"_marker_mean_scores.csv")),row.names=1,check.names=FALSE));marker<-marker[rownames(md),,drop=FALSE]
  for(ct in colnames(marker))marker_check[[paste(s,ct)]]<-data.frame(sample=s,cell_type=ct,RCTD_marker_rho=safe_cor(broad[,ct],marker[,ct]),mean_weight=mean(broad[,ct]))
  for(g in genes)for(ct in names(groups)) {
    depth<-rank(md$nCount);x<-rank(as.numeric(dat[g,]));y<-rank(broad[,ct])
    corr[[paste(s,g,ct)]]<-data.frame(sample=s,gene=g,cell_type=ct,rho=safe_cor(as.numeric(dat[g,]),broad[,ct]),depth_adjusted_rank_rho=if(sd(y)>0)cor(resid(lm(x~depth)),resid(lm(y~depth)))else NA_real_)
  }
  # 基质与周围髓系权重的描述性关系；不将点作为独立患者进行P值检验。
  kn<-nn2(as.matrix(md[,c("array_row","array_col")]),k=7);ii<-rep(seq_len(nrow(md)),each=6);jj<-as.vector(t(kn$nn.idx[,-1]));dd<-as.vector(t(kn$nn.dists[,-1]));keep<-dd<=2.1
  W<-sparseMatrix(i=ii[keep],j=jj[keep],x=1,dims=c(nrow(md),nrow(md)));W<-((W+t(W))>0)*1;diag(W)<-0
  adj<-as.numeric(W%*%broad[,"Myeloid"])/pmax(Matrix::rowSums(W),1)
  neighbor[[s]]<-data.frame(sample=s,same_spot_rho=safe_cor(broad[,"Stromal"],broad[,"Myeloid"]),neighbor_rho=safe_cor(broad[,"Stromal"],adj),mean_myeloid=mean(broad[,"Myeloid"]))
  md$x<-md$pixel_col*z$scale$tissue_lowres_scalef;md$y<--md$pixel_row*z$scale$tissue_lowres_scalef
  f<-file.path(root,"raw",paste0("GSM",7058755+as.integer(sub("C","",s)),"_",s,"_tissue_lowres_image.png.gz"));con<-gzfile(f,"rb");img<-readPNG(readBin(con,"raw",n=20000000));close(con);img<-img*.3+.7
  bar<-500/55*z$scale$spot_diameter_fullres*z$scale$tissue_lowres_scalef
  colors<-c(Myeloid="#C65D50",Stromal="#258F8A",Epithelial="#6B5D9C")
  for(ct in names(colors)) {
    md$value<-broad[,ct]
    plots[[paste(s,ct)]]<-ggplot()+annotation_raster(img,xmin=0,xmax=dim(img)[2],ymin=-dim(img)[1],ymax=0)+
      geom_point(data=md,aes(x,y,colour=value),size=.85)+
      annotate("rect",xmin=10,xmax=30+bar,ymin=-dim(img)[1]+8,ymax=-dim(img)[1]+49,fill="white",alpha=.9)+
      annotate("segment",x=20,xend=20+bar,y=-dim(img)[1]+20,yend=-dim(img)[1]+20,linewidth=.65)+annotate("text",x=20+bar/2,y=-dim(img)[1]+37,label="500 um",size=2.7)+
      coord_fixed(xlim=c(0,dim(img)[2]),ylim=c(-dim(img)[1],0),expand=FALSE)+scale_colour_gradient(low="#F0F1F4",high=colors[[ct]],limits=c(0,1))+
      theme_void(base_size=11)+theme(plot.title=element_text(face="bold"),legend.position="bottom",legend.key.height=grid::unit(.25,"cm"),legend.key.width=grid::unit(.75,"cm"),legend.title=element_text(size=9))+labs(title=paste(s,ct),colour="Relative weight")
  }
}
q<-rbindlist(quality);a<-rbindlist(corr);p<-rbindlist(proportions);mc<-rbindlist(marker_check)
write.csv(q,file.path(out,"fit_quality.csv"),row.names=FALSE);write.csv(a,file.path(out,"gene_RCTD_correlations.csv"),row.names=FALSE);write.csv(p,file.path(out,"celltype_weights_summary.csv"),row.names=FALSE)
write.csv(rbindlist(neighbor),file.path(out,"stromal_myeloid_neighborhood.csv"),row.names=FALSE)
write.csv(mc,file.path(out,"RCTD_marker_concordance.csv"),row.names=FALSE)
savefig(wrap_plots(plots,ncol=3),"Figure_RCTD_spatial",13,15)
a$gene<-factor(a$gene,levels=genes)
hp<-ggplot(a,aes(gene,cell_type,fill=depth_adjusted_rank_rho))+geom_tile(colour="white")+geom_text(aes(label=sprintf("%.2f",depth_adjusted_rank_rho)),size=3)+facet_wrap(~sample,nrow=1)+scale_y_discrete(labels=c(Epithelial="Epithelial",Stromal="Stromal",Myeloid="Myeloid",Smooth_muscle="Smooth muscle",SPP1_myeloid="SPP1+ myeloid"))+scale_fill_gradient2(low="#456DA5",mid="white",high="#C65D50",limits=c(-1,1),na.value="#DDDDDD")+theme_minimal(base_size=11)+theme(panel.grid=element_blank(),axis.text.x=element_text(angle=45,hjust=1))+labs(x=NULL,y=NULL,fill="Depth-adjusted\nrank correlation")
savefig(hp,"Figure_RCTD_gene_associations",12,4)
cp<-ggplot(mc,aes(sample,cell_type,fill=RCTD_marker_rho))+geom_tile(colour="white")+geom_text(aes(label=sprintf("%.2f",RCTD_marker_rho)),size=4)+scale_y_discrete(labels=c(Epithelial="Epithelial",Stromal="Stromal",Myeloid="Myeloid",Smooth_muscle="Smooth muscle"))+scale_fill_gradient(low="#F4F0EB",high="#258F8A",limits=c(0,1))+theme_minimal(base_size=12)+theme(panel.grid=element_blank())+labs(x=NULL,y=NULL,fill="Spearman rho",title="RCTD marker concordance")
savefig(cp,"Figure_RCTD_marker_concordance",6,3.8)
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
