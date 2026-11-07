# 参考Huang Fig.8的三联布局，全部数据来自本研究；不复用原文图像。
suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
root <- 'D:/Bioinfor/colorectal/TCGA_CIBERSORT_final'
dest <- file.path(root,'figures/Huang_style');dir.create(dest,recursive=TRUE,showWarnings=FALSE)
cb <- readRDS(file.path(root,'results/cibersort.rds'))
genes <- c('EPB41L3','KLB','ANGPTL1')
stars <- function(p) ifelse(is.na(p),'NA',ifelse(p<.001,'***',ifelse(p<.01,'**',ifelse(p<.05,'*','ns'))))
short <- function(s) {
 s<-gsub('T cells CD4 memory','CD4 memory',s,fixed=TRUE)
 s<-gsub('T cells regulatory (Tregs)','Tregs',s,fixed=TRUE)
 s<-gsub('T cells follicular helper','Tfh',s,fixed=TRUE)
 s<-gsub('T cells gamma delta','T gamma delta',s,fixed=TRUE)
 s<-gsub('T cells CD4 naive','CD4 naive',s,fixed=TRUE)
 s<-gsub('T cells CD8','CD8 T cells',s,fixed=TRUE)
 s<-gsub('Dendritic cells','DC',s,fixed=TRUE)
 s
}
theme_set(theme_classic(base_size=9,base_family='Arial'))
draw <- function(scores,m,co,dif,label,prefix,note) {
 cells<-rownames(scores);n<-length(cells);tum<-m$Type=='01'
 stopifnot(identical(colnames(scores),m$Sample),all(cells%in%dif$Cell))
 long<-do.call(rbind,lapply(cells,function(cell)data.frame(Cell=cell,Group=m$Group,Score=as.numeric(scores[cell,]))))
 long$Cell<-factor(long$Cell,levels=cells)
 sig<-data.frame(Cell=factor(cells,levels=cells),Label=stars(dif$adj.P.Val[match(cells,dif$Cell)]))
 sig$Y<-max(scores)+.04*diff(range(scores))
 a<-ggplot(long,aes(Cell,Score,fill=Group))+geom_boxplot(position=position_dodge(.78),width=.68,outlier.size=.35,linewidth=.25)+
  geom_text(data=sig,aes(Cell,Y,label=Label),inherit.aes=FALSE,size=2.6)+
  scale_fill_manual(values=c('Adjacent normal'='#247C74','Tumour'='#B96551'),name=NULL)+
  scale_x_discrete(labels=short)+scale_y_continuous(expand=expansion(mult=c(.03,.12)))+
  labs(title=paste('A',label,'tissue comparison'),subtitle=sprintf('Adjacent normal n=%d; tumour n=%d',sum(!tum),sum(tum)),x=NULL,
       y='Relative fraction')+
  theme(legend.position='top',axis.text.x=element_text(angle=55,hjust=1,size=7),plot.title=element_text(face='bold',size=11),plot.subtitle=element_text(size=9))
 # 下三角：免疫细胞/签名之间的Spearman相关。连线：三基因与各指标的相关。
 cm<-cor(t(scores[,tum,drop=FALSE]),method='spearman')
 grid<-expand.grid(i=seq_len(n),j=seq_len(n));grid$Rho<-cm[cbind(grid$i,grid$j)]
 grid<-grid[grid$i>=grid$j,];grid$x<-grid$j;grid$y<--grid$i
 nodes<-data.frame(Gene=genes,x=n+c(4,5.5,4),y=-n*c(.15,.48,.82))
 edges<-co;edges$x<-match(edges$Cell,cells);edges$y<--edges$x
 edges$xend<-nodes$x[match(edges$Gene,genes)];edges$yend<-nodes$y[match(edges$Gene,genes)]
 edges$Kind<-ifelse(edges$FDR<.05,ifelse(edges$Rho>0,'Positive','Negative'),'Not significant')
 b<-ggplot(grid,aes(x,y))+geom_tile(aes(fill=Rho),width=.96,height=.96)+
  scale_fill_gradient2(low='#347B93',mid='white',high='#BC6652',limits=c(-1,1),name='Cell-cell rho')+
  geom_curve(data=edges[edges$Kind=='Not significant',],aes(x=x,y=y,xend=xend,yend=yend),inherit.aes=FALSE,colour='grey80',linewidth=.15,alpha=.45,curvature=.08)+
  geom_curve(data=edges[edges$Kind!='Not significant',],aes(x=x,y=y,xend=xend,yend=yend,colour=Kind,linewidth=abs(Rho)),inherit.aes=FALSE,curvature=.08)+
  scale_colour_manual(values=c(Positive='#B96551',Negative='#247C74'),name='Gene association\nFDR < 0.05')+
  scale_linewidth_continuous(range=c(.3,1.1),limits=c(0,1),breaks=c(.2,.4,.6),name='Absolute rho')+
  geom_point(data=data.frame(x=seq_len(n),y=-seq_len(n)),aes(x,y),inherit.aes=FALSE,size=.7,colour='#53616B')+
  geom_point(data=nodes,aes(x,y),inherit.aes=FALSE,size=2,colour='#53616B')+
  geom_text(data=nodes,aes(x+.4,y,label=Gene),inherit.aes=FALSE,hjust=0,size=2.5,fontface='bold')+
  scale_x_continuous(breaks=seq_len(n),labels=short(cells),limits=c(.4,n+9),expand=c(0,0))+
  scale_y_continuous(breaks=-seq_len(n),labels=short(cells),limits=c(-n-.6,0),expand=c(0,0))+
  coord_fixed(clip='off')+labs(title='B  Immune correlations and gene associations',subtitle=paste('Tumours only: n =',sum(tum)),x=NULL,y=NULL)+
  theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text=element_text(size=5.5),axis.text.x=element_text(angle=55,hjust=1),
        plot.title=element_text(face='bold',size=10),plot.subtitle=element_text(size=8),legend.position='right',legend.title=element_text(size=7),legend.text=element_text(size=6),legend.key.height=grid::unit(.25,'cm'))
 co$Cell<-factor(co$Cell,levels=cells);co$Gene<-factor(co$Gene,levels=rev(genes))
 co$Label<-paste0(sprintf('%.2f',co$Rho),'\n',ifelse(co$FDR<.05,stars(co$FDR),''))
 c<-ggplot(co,aes(Cell,Gene,fill=Rho))+geom_tile(colour='white',linewidth=.3)+geom_text(aes(label=Label),size=2.5,lineheight=.85)+
  scale_fill_gradient2(low='#347B93',mid='white',high='#BC6652',limits=c(-1,1),name='Spearman rho')+
  scale_x_discrete(labels=short)+labs(title='C  Three-gene associations with immune profiles',subtitle='Numbers: Spearman rho; asterisks: BH-adjusted significance',x=NULL,y=NULL)+
  theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(angle=55,hjust=1,size=8),axis.text.y=element_text(face='bold'),plot.title=element_text(face='bold',size=11),plot.subtitle=element_text(size=8))
 p<-((a|b)/c)+plot_layout(heights=c(1.55,1))+plot_annotation(caption=note,theme=theme(plot.caption=element_text(size=8,hjust=0)))
 for(ext in c('pdf','png','tiff')) {
  path<-file.path(dest,paste0('Figure_',prefix,'_Huang_style.',ext))
  if(ext=='pdf')ggsave(path,p,width=15,height=9.5,device=cairo_pdf)else if(ext=='tiff')ggsave(path,p,width=15,height=9.5,dpi=600,compression='lzw',bg='white')else ggsave(path,p,width=15,height=9.5,dpi=450,bg='white')
 }
 write.csv(data.frame(Cell=cells,cm,check.names=FALSE),file.path(dest,paste0(prefix,'_cell_correlations.csv')),row.names=FALSE)
 write.csv(edges,file.path(dest,paste0(prefix,'_network_edges.csv')),row.names=FALSE)
}
keep<-cb$qc$Pass
draw(t(cb$raw[keep,1:22]),droplevels(cb$metadata[keep,]),
 read.csv(file.path(root,'results/CIBERSORT_QC_correlations.csv')),
 read.csv(file.path(root,'results/CIBERSORT_QC_differences.csv')),'CIBERSORT','CIBERSORT',
 'CIBERSORT fit P < 0.05. A: site-adjusted model, patient-cluster robust SE, BH across 22 tests. B/C: BH across 66 gene-cell tests.\nGrey links: FDR >= 0.05. * FDR < 0.05; ** < 0.01; *** < 0.001. Relative fractions do not measure absolute cell abundance.')
