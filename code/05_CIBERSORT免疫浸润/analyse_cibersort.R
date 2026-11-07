# CIBERSORT组间差异、肿瘤内三基因相关及全样本敏感性分析。
suppressPackageStartupMessages({library(limma);library(ggplot2);library(patchwork)})
root <- 'D:/Bioinfor/colorectal/TCGA_CIBERSORT_final';out<-file.path(root,'results')
fig<-file.path(root,'figures');dir.create(fig,showWarnings=FALSE)
cib<-readRDS(file.path(out,'cibersort.rds'))
m<-cib$metadata
wt<-function(x,n)write.csv(x,file.path(out,paste0(n,'.csv')),row.names=FALSE)
cs<-t(cib$raw[,1:22]);genes<-c('EPB41L3','KLB','ANGPTL1')
analyse<-function(scores,keep,label) {
 d<-droplevels(m[keep,]);s<-scores[,keep,drop=FALSE];v<-apply(s,1,sd)>0
 # QC后仅剩1对样本，无法稳定估计共同组内相关；统一使用患者聚类稳健标准误。
 dif<-do.call(rbind,lapply(rownames(s)[v],function(cell){
  dd<-d;dd$Score<-s[cell,];fit<-lm(Score~Group+Site,dd)
  cov<-sandwich::vcovCL(fit,cluster=dd$Patient,type='HC1')
  beta<-coef(fit)['GroupTumour'];se<-sqrt(cov['GroupTumour','GroupTumour']);df<-length(unique(dd$Patient))-1
  data.frame(Cell=cell,Difference=unname(beta),CI_low=beta-qt(.975,df)*se,CI_high=beta+qt(.975,df)*se,P.Value=2*pt(-abs(beta/se),df))
 }));dif$adj.P.Val<-p.adjust(dif$P.Value,'BH',n=nrow(scores))
 wt(dif,paste0(label,'_differences'))
 cc<-list();tum<-keep & m$Type=='01'
 for(g in genes)for(cell in rownames(scores)) {
  if(sd(scores[cell,tum])==0) {rho<-NA;p<-NA} else {z<-cor.test(cib$genes[g,tum],scores[cell,tum],method='spearman',exact=FALSE);rho<-unname(z$estimate);p<-z$p.value}
  cc[[paste(g,cell)]]<-data.frame(Gene=g,Cell=cell,N=sum(tum),Rho=rho,P=p)
 }
 co<-do.call(rbind,cc);co$FDR<-p.adjust(co$P,'BH',n=length(genes)*nrow(scores));wt(co,paste0(label,'_correlations'))
 pairs<-intersect(d$Patient[d$Type=='01'],d$Patient[d$Type=='11'])
 ti<-match(paste(pairs,'01'),paste(d$Patient,d$Type));ni<-match(paste(pairs,'11'),paste(d$Patient,d$Type))
 pp<-do.call(rbind,lapply(rownames(s),function(cell){
  delta<-s[cell,ti]-s[cell,ni];z<-delta[delta!=0];n<-length(z);r<-rank(abs(z));w<-sum(r[z>0]);ties<-table(abs(z))
  variance<-n*(n+1)*(2*n+1)/24-sum(ties^3-ties)/48;a<-w-n*(n+1)/4
  p<-if(length(pairs)<5) NA_real_ else if(variance>0)2*pnorm(-abs((a-sign(a)*.5)/sqrt(variance)))else 1
  data.frame(Cell=cell,N_pairs=length(pairs),Median_difference=median(delta),P=p)
 }));pp$FDR<-p.adjust(pp$P,'BH');wt(pp,paste0(label,'_paired'))
 wt(data.frame(N=nrow(d),Tumours=sum(tum),Normal=sum(keep & m$Type=='11'),Pairs=length(pairs),Constant_scores=sum(!v)),paste0(label,'_counts'))
 list(dif=dif,cor=co,data=d,scores=s)
}
res<-analyse(cs,cib$qc$Pass,'CIBERSORT_QC')
all<-analyse(cs,rep(TRUE,nrow(m)),'CIBERSORT_all')
theme_set(theme_classic(base_size=9,base_family='Arial'))
long<-do.call(rbind,lapply(rownames(res$scores),function(cell)data.frame(res$data,Cell=cell,Score=as.numeric(res$scores[cell,]))))
labs<-setNames(paste0(res$dif$Cell,'\nFDR ',format.pval(res$dif$adj.P.Val,digits=2,eps=.001)),res$dif$Cell)
p1<-ggplot(long,aes(Group,Score,fill=Group))+geom_boxplot(outlier.size=.3,linewidth=.3)+
 facet_wrap(~Cell,scales='free_y',ncol=4,labeller=as_labeller(labs))+
 scale_fill_manual(values=c('Adjacent normal'='#247C74','Tumour'='#B96551'))+
 labs(title='A  CIBERSORT relative immune fractions',subtitle=paste('Fit P < 0.05:',sum(cib$qc$Pass),'samples; site-adjusted, patient-cluster robust SE'),x=NULL,y='Relative fraction')+
 theme(legend.position='none',axis.text.x=element_text(angle=30,hjust=1),strip.background=element_blank(),strip.text=element_text(size=7))
co<-res$cor;co$Gene<-factor(co$Gene,levels=genes);co$Cell<-factor(co$Cell,levels=rev(rownames(cs)))
co$Label<-ifelse(is.na(co$Rho),'NA',sprintf('%.2f%s',co$Rho,ifelse(co$FDR<.001,'***',ifelse(co$FDR<.01,'**',ifelse(co$FDR<.05,'*','')))))
p2<-ggplot(co,aes(Gene,Cell,fill=Rho))+geom_tile(colour='white')+geom_text(aes(label=Label),size=2.5)+
 scale_fill_gradient2(low='#347B93',mid='white',high='#BC6652',limits=c(-1,1),name='rho',na.value='grey90')+
 labs(title='B  Gene-fraction associations within tumours',subtitle='BH correction across 66 tests',x=NULL,y=NULL)+theme(axis.line=element_blank(),axis.ticks=element_blank())
ggsave(file.path(fig,'Figure_CIBERSORT.pdf'),p1/p2+plot_layout(heights=c(1.5,1)),width=11,height=17,device=cairo_pdf)
ggsave(file.path(fig,'Figure_CIBERSORT.png'),p1/p2+plot_layout(heights=c(1.5,1)),width=11,height=17,dpi=350,bg='white')
