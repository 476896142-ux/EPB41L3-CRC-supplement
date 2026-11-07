# 非线性结果同时展示 OS 与 PFI，不按显著性选择终点。
try(Sys.setlocale('LC_CTYPE','Chinese (Simplified)_China.utf8'),silent=TRUE)
suppressPackageStartupMessages({library(data.table);library(survival);library(ggplot2)})
root<-'D:/Bioinfor/colorectal/EPB41L3_bulk_validation'
d<-as.data.frame(fread('D:/Bioinfor/colorectal/TCGA_three_gene_prognosis/results/matched_clinical_survival.csv'))
d$Stage<-factor(d$Stage,levels=c('I','II','III','IV'));d$Sex<-factor(d$Sex,levels=c('FEMALE','MALE'));d$Site<-factor(d$Site,levels=c('COAD','READ'))
curves<-list();phs<-list();events<-list()
for(e in c('OS','PFI')) {
 z<-d;z$Time<-z[[paste0(e,'_time')]];z$Event<-z[[paste0(e,'_event')]]
 z<-z[is.finite(z$Time)&z$Time>0&z$Event%in%c(0,1)&complete.cases(z[,c('Age','Sex','Stage','Site','EPB41L3')]),]
 fit<-coxph(Surv(Time,Event)~splines::ns(EPB41L3,df=3)+Age+Sex+Site+Stage,data=z,x=TRUE)
 grid<-z[rep(1,200),];grid$EPB41L3<-seq(quantile(z$EPB41L3,.01),quantile(z$EPB41L3,.99),length.out=200)
 reference<-grid;reference$EPB41L3<-median(z$EPB41L3)
 contrast<-model.matrix(fit,grid)-model.matrix(fit,reference)
 lp<-as.numeric(contrast%*%coef(fit));se<-sqrt(rowSums((contrast%*%vcov(fit))*contrast))
 curves[[e]]<-data.table(Endpoint=e,Expression_SD=grid$EPB41L3,HR=exp(lp),Low=exp(lp-1.96*se),High=exp(lp+1.96*se),Reference_SD=median(z$EPB41L3))
 pt<-cox.zph(fit)$table;phs[[e]]<-data.table(Endpoint=e,Term=rownames(pt),P=pt[,'p'])
 z$Bin<-cut(z$EPB41L3,quantile(z$EPB41L3,seq(0,1,.2)),include.lowest=TRUE)
 events[[e]]<-as.data.table(z)[,.(N=.N,Events=sum(Event),Min=min(EPB41L3),Max=max(EPB41L3)),by=Bin][,Endpoint:=e]
}
cc<-rbindlist(curves);fwrite(cc,file.path(root,'results/nonlinear_curves.csv'))
fwrite(rbindlist(phs),file.path(root,'results/nonlinear_PH.csv'));fwrite(rbindlist(events),file.path(root,'results/expression_quintile_event_audit.csv'))
tests<-fread(file.path(root,'results/nonlinearity.csv'))
old_overall<-fread('D:/Bioinfor/colorectal/TCGA_three_gene_prognosis/results/spline_overall_associations.csv')[Gene=='EPB41L3']
old_nonlin<-fread('D:/Bioinfor/colorectal/TCGA_three_gene_prognosis/results/gene_nonlinearity_tests.csv')[Predictor=='EPB41L3']
tests[,Original_family_overall_FDR:=old_overall$FDR[match(Endpoint,old_overall$Endpoint)]]
tests[,Original_family_nonlinearity_FDR:=old_nonlin$FDR[match(Endpoint,old_nonlin$Endpoint)]]
stopifnot(all(abs(tests$Overall_P-old_overall$P[match(tests$Endpoint,old_overall$Endpoint)])<1e-8))
tests[,Overall_FDR:=p.adjust(Overall_P,'BH')];tests[,Nonlinearity_FDR:=p.adjust(Nonlinearity_P,'BH')]
tests[,All_four_test_FDR_overall:=p.adjust(c(Overall_P,Nonlinearity_P),'BH')[seq_len(.N)]]
tests[,All_four_test_FDR_nonlinear:=p.adjust(c(Overall_P,Nonlinearity_P),'BH')[.N+seq_len(.N)]]
fwrite(tests,file.path(root,'results/nonlinearity_adjusted.csv'))
labels<-tests[,.(Endpoint,Label=sprintf('Overall P = %.3g\nNonlinearity P = %.3g',Overall_P,Nonlinearity_P))]
p<-ggplot(cc,aes(Expression_SD,HR,color=Endpoint,fill=Endpoint))+geom_hline(yintercept=1,linetype=2,color='grey55')+
 geom_ribbon(aes(ymin=Low,ymax=High),alpha=.18,color=NA)+geom_line(linewidth=.9)+
 geom_text(data=labels,aes(x=-Inf,y=Inf,label=Label),inherit.aes=FALSE,hjust=-.05,vjust=1.2,size=3.5)+
 facet_wrap(~Endpoint)+scale_y_log10()+scale_color_manual(values=c(OS='#387CA3',PFI='#CF8744'))+
 scale_fill_manual(values=c(OS='#387CA3',PFI='#CF8744'))+
 labs(x='EPB41L3 expression (SD units)',y='Adjusted hazard ratio vs median (95% CI)',title='Exploratory nonlinear associations',subtitle='Natural cubic spline, 3 df; age, sex, site and stage adjusted')+
 theme_classic(base_size=11)+theme(legend.position='none',strip.background=element_blank(),strip.text=element_text(face='bold'))
ggsave(file.path(root,'Figure_prognosis_nonlinear.pdf'),p,width=9,height=4.5)
ggsave(file.path(root,'Figure_prognosis_nonlinear.png'),p,width=9,height=4.5,dpi=300)

# 独立检查：新线性模型须精确复现旧结果，防止缺失分期被作为有效类别。
old<-fread(file.path(root,'results/original_prognosis_family_preserved.csv'))
new<-fread(file.path(root,'results/prognosis_models.csv'))
check<-merge(new,old[Model%in%c('Unadjusted','Adjusted','Stage_stratified')],by=c('Endpoint','Model'))
stopifnot(all(check$N.x==check$N.y),all(abs(check$HR.x-check$HR.y)<1e-8),all(abs(check$P.x-check$P.y)<1e-8))
inp<-readRDS(file.path(root,'results/analysis_input.rds'))
markers<-read.delim(file.path(root,'reference/genes.txt'),check.names=FALSE)
for(cohort in names(inp$expression)) {
 x<-inp$expression[[cohort]];s<-inp$scores[[cohort]]
 for(pop in c('Monocytic lineage','Fibroblasts','Endothelial cells')) {
  genes<-unique(markers[markers[['Cell population']]==pop,'HUGO symbols']);genes<-setdiff(intersect(genes,rownames(x)),'EPB41L3')
  stopifnot(max(abs(colMeans(x[genes,,drop=FALSE])-s[[pop]]))<1e-10)
 }
 stopifnot(all(abs(s$EPB41L3-x['EPB41L3',])<1e-10))
}
writeLines('PASS: scores independently recomputed; original linear prognosis reproduced; sample identifiers checked by main script.',file.path(root,'results/verification.txt'))
print(tests);print(rbindlist(phs));print(rbindlist(events))
