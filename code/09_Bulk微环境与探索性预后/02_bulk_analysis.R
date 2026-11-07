# 队列内计算微环境评分，再比较关联方向；不使用 ssGSEA。
try(Sys.setlocale('LC_CTYPE','Chinese (Simplified)_China.utf8'),silent=TRUE)
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(survival)})
root <- 'D:/Bioinfor/colorectal/EPB41L3_bulk_validation'
geo <- 'E:/codex/2026-10-31/dui/outputs/GEO_ML_hub10'
raw <- 'D:/Bioinfor/colorectal/GEO_ML_hub10/raw'
dir.create(file.path(root,'results'),showWarnings=FALSE)
write_out <- function(x,name) fwrite(x,file.path(root,'results',paste0(name,'.csv')))
source(file.path(root,'reference/MCPcounter.R'))
markers <- read.delim(file.path(root,'reference/genes.txt'),check.names=FALSE)
markers <- markers[markers[['HUGO symbols']]!='EPB41L3',]
needed <- unique(c('EPB41L3',markers[['HUGO symbols']]))
primary <- c('Monocytic lineage','Fibroblasts','Endothelial cells')
set.seed(20261106)

# TCGA：仅保留原发肿瘤，沿用已核查的样本和临床表。
tcga <- fread('D:/Bioinfor/colorectal/TCGA_\u673a\u5668\u5b66\u4e60\u7279\u5f81\u7b5b\u9009/TCGA_RAW/HiSeqV2')
rn <- tcga[[1]]; tx <- as.matrix(tcga[rn %in% needed,-1,with=FALSE]);rownames(tx)<-rn[rn %in% needed]
tm <- fread('D:/Bioinfor/colorectal/TCGA_CIBERSORT_final/results/sample_metadata.csv')
tm <- tm[Group=='Tumour'];stopifnot(!anyDuplicated(tm$Patient))
tx <- tx[,tm$Sample,drop=FALSE]

# GSE44076：完整平台表中提取标记，多探针已取中位数。
a <- fread(file.path(root,'GSE44076_marker_expression.csv'))
ax <- t(as.matrix(a[,-1,with=FALSE]));colnames(ax)<-a$Sample
am <- fread(file.path(geo,'GSE44076_metadata.csv'))[Class=='Tumor']
stopifnot(!anyDuplicated(am$Patient));ax<-ax[,am$Sample,drop=FALSE]

# GSE20916：沿用非显微切割样本的全平台 RMA 结果。
b <- readRDS(file.path(geo,'GSE20916_RMA_probe_expression.rds'))
ann <- fread(file.path(raw,'hgu133plus2_unique_probe_map.csv'))
ann <- ann[SYMBOL %in% needed & PROBEID %in% rownames(b)]
bx <- t(vapply(split(ann$PROBEID,ann$SYMBOL),function(p) apply(b[p,,drop=FALSE],2,median),numeric(ncol(b))))
colnames(bx)<-colnames(b)
bm <- fread(file.path(geo,'GSE20916_sample_inclusion.csv'))
bm <- bm[tolower(as.character(Included))=='true' & tissue=='adenocarcinoma']
bm$Sample<-bm$geo_accession;stopifnot(nrow(bm)>10);bx<-bx[,bm$Sample,drop=FALSE]
datasets <- list(TCGA=tx,GSE44076=ax,GSE20916=bx)
metadata <- list(TCGA=tm,GSE44076=am,GSE20916=bm)
coverage<-list();scores<-list();results<-list();sens<-list()
for(cohort in names(datasets)) {
 x<-datasets[[cohort]];stopifnot('EPB41L3'%in%rownames(x),all(is.finite(x)),!anyDuplicated(colnames(x)))
 sc<-MCPcounter.estimate(x,featuresType='HUGO_symbols',genes=markers)
 d<-data.table(Sample=colnames(x),Cohort=cohort,EPB41L3=as.numeric(x['EPB41L3',]),t(sc))
 scores[[cohort]]<-d
 for(pop in unique(markers[['Cell population']])) {
  gs<-unique(markers[markers[['Cell population']]==pop,'HUGO symbols'])
  coverage[[length(coverage)+1]]<-data.table(Cohort=cohort,Population=pop,Available=sum(gs%in%rownames(x)),Total=length(gs),Missing=paste(setdiff(gs,rownames(x)),collapse=';'))
 }
 for(pop in rownames(sc)) {
  xx<-d$EPB41L3;yy<-d[[pop]];n<-length(xx)
  ct<-cor.test(xx,yy,method='spearman',exact=FALSE)
  boot<-replicate(2000,{ix<-sample.int(n,n,replace=TRUE);cor(xx[ix],yy[ix],method='spearman')})
  ci<-quantile(boot,c(.025,.975),na.rm=TRUE)
  results[[length(results)+1]]<-data.table(Cohort=cohort,Population=pop,N=n,Rho=unname(ct$estimate),Low=ci[1],High=ci[2],P=ct$p.value,Family=ifelse(pop%in%primary,'Primary','Exploratory'))
 }
 # 单核谱系的偏秩相关：控制两类基质评分，不等同于纯度校正。
 z<-data.frame(x=rank(d$EPB41L3),y=rank(d[['Monocytic lineage']]),fib=rank(d$Fibroblasts),endo=rank(d[['Endothelial cells']]))
 r<-cor(resid(lm(x~fib+endo,z)),resid(lm(y~fib+endo,z)))
 df<-nrow(z)-4
 sens[[cohort]]<-data.table(Cohort=cohort,N=nrow(z),Partial_rho=r,P=2*pt(-abs(r*sqrt(df/(1-r*r))),df))
}
assoc<-rbindlist(results);assoc[,FDR:=p.adjust(P,'BH'),by=Family]
ss<-rbindlist(sens);ss[,FDR:=p.adjust(P,'BH')]
all_scores<-rbindlist(scores,fill=TRUE)
write_out(assoc,'bulk_correlations');write_out(rbindlist(coverage),'marker_coverage');write_out(ss,'stromal_adjusted_monocytic');write_out(all_scores,'sample_scores')
write_out(rbindlist(lapply(names(metadata),function(k) data.table(Cohort=k,N=nrow(metadata[[k]]),Patient_identifiers=ifelse(k=='GSE20916','Not available; sample-level','Unique tumor patients')))),'cohort_audit')
saveRDS(list(expression=datasets,metadata=metadata,scores=scores),file.path(root,'results/analysis_input.rds'))

# 临床校正后的单核谱系关联：只在有完整分期资料的队列实施。
clinical<-fread('D:/Bioinfor/colorectal/TCGA_three_gene_prognosis/results/matched_clinical_survival.csv')
clinical[,Stage:=factor(Stage,levels=c('I','II','III','IV'))]
clinical[,Sex:=factor(Sex,levels=c('FEMALE','MALE'))]
clinical[,Site:=factor(Site,levels=c('COAD','READ'))]
cm<-merge(scores$TCGA,clinical[,.(Sample,Age,Sex,Site,Stage)],by='Sample')
ca<-copy(scores$GSE44076);mi<-match(ca$Sample,am$Sample)
ca[,`:=`(Age=as.numeric(am$age[mi]),Sex=am$gender[mi],Site=am$location[mi],Stage=as.character(am$stage[mi]))]
clin_res<-list()
for(k in c('TCGA','GSE44076')) {
 d<-as.data.frame(if(k=='TCGA') cm else ca)
 d<-d[complete.cases(d[,c('Age','Sex','Site','Stage')]),]
 d$x<-rank(d$EPB41L3);d$y<-rank(d[['Monocytic lineage']]);d$Age<-rank(d$Age)
 for(v in c('Sex','Site','Stage')) d[[v]]<-factor(d[[v]])
 vars<-c('Age',c('Sex','Site','Stage')[vapply(d[c('Sex','Site','Stage')],nlevels,integer(1))>1])
 f<-reformulate(vars,'x');g<-reformulate(vars,'y')
 rx<-resid(lm(f,d));ry<-resid(lm(g,d));r<-cor(rx,ry);df<-nrow(d)-qr(model.matrix(f,d))$rank-1
 clin_res[[k]]<-data.table(Cohort=k,N=nrow(d),Partial_rho=r,P=2*pt(-abs(r*sqrt(df/(1-r*r))),df),Covariates=paste(vars,collapse='+'))
}
cr<-rbindlist(clin_res);cr[,FDR:=p.adjust(P,'BH')];write_out(cr,'clinical_adjusted_monocytic')

# 固定连续 EPB41L3 模型；沿用原始 380 人标准化表达，不寻找切点。
surv_rows<-list();ph_rows<-list();spline_rows<-list()
surv_data<-merge(clinical,scores$TCGA[,.(Sample,Monocytic=`Monocytic lineage`,Fibroblasts,Endothelial=`Endothelial cells`)],by='Sample')
for(endpoint in c('OS','PFI')) {
 d<-as.data.frame(surv_data);d$Time<-d[[paste0(endpoint,'_time')]];d$Event<-d[[paste0(endpoint,'_event')]]
 d<-d[is.finite(d$Time)&d$Time>0&d$Event%in%c(0,1),]
 for(v in c('Sex','Site','Stage')) d[[v]]<-factor(d[[v]])
 formulas<-list(Unadjusted=Surv(Time,Event)~EPB41L3,
  Adjusted=Surv(Time,Event)~EPB41L3+Age+Sex+Site+Stage,
  Stage_stratified=Surv(Time,Event)~EPB41L3+Age+Sex+Site+strata(Stage),
  Microenvironment_adjusted=Surv(Time,Event)~EPB41L3+Age+Sex+Site+Stage+Monocytic+Fibroblasts+Endothelial)
 for(model in names(formulas)) {
  fit<-coxph(formulas[[model]],d,x=TRUE);s<-summary(fit);ci<-s$conf.int['EPB41L3',]
  surv_rows[[length(surv_rows)+1]]<-data.table(Endpoint=endpoint,Model=model,N=fit$n,Events=fit$nevent,HR=ci[1],Low=ci[3],High=ci[4],P=s$coefficients['EPB41L3','Pr(>|z|)'])
  ph<-cox.zph(fit)$table
  ph_rows[[length(ph_rows)+1]]<-data.table(Endpoint=endpoint,Model=model,Term=rownames(ph),P=ph[,'p'])
 }
 dd<-d[complete.cases(d[,c('EPB41L3','Age','Sex','Site','Stage')]),]
 linear<-coxph(formulas$Adjusted,dd)
 spline<-coxph(Surv(Time,Event)~splines::ns(EPB41L3,df=3)+Age+Sex+Site+Stage,dd)
 null<-coxph(Surv(Time,Event)~Age+Sex+Site+Stage,dd)
 spline_rows[[endpoint]]<-data.table(Endpoint=endpoint,N=nrow(dd),Overall_P=anova(null,spline,test='Chisq')[2,'Pr(>|Chi|)'],Nonlinearity_P=anova(linear,spline,test='Chisq')[2,'Pr(>|Chi|)'])
}
sr<-rbindlist(surv_rows);sr[,FDR:=p.adjust(P,'BH'),by=Model]
write_out(sr,'prognosis_models');write_out(rbindlist(ph_rows),'proportional_hazards');write_out(rbindlist(spline_rows),'nonlinearity')
old<-fread('D:/Bioinfor/colorectal/TCGA_three_gene_prognosis/results/cox_target_associations.csv')[Predictor=='EPB41L3']
write_out(old,'original_prognosis_family_preserved')

# 出图：队列用统一配色；相关系数和 HR 均附置信区间。
cols<-c(TCGA='#387CA3',GSE44076='#CF8744',GSE20916='#459C8C')
theme_set(theme_classic(base_size=11)+theme(plot.title=element_text(face='bold'),legend.position='bottom'))
main<-assoc[Family=='Primary'];main[,Population:=factor(Population,levels=rev(primary))]
p1<-ggplot(main,aes(Rho,Population,color=Cohort))+geom_vline(xintercept=0,linetype=2,color='grey65')+
 geom_errorbar(aes(xmin=Low,xmax=High),orientation='y',width=.18,position=position_dodge(.6))+
 geom_point(position=position_dodge(.6),size=2.8)+scale_color_manual(values=cols)+coord_cartesian(xlim=c(-1,1))+
 labs(title='EPB41L3 and microenvironment scores',x='Spearman correlation (95% bootstrap CI)',y=NULL)
long<-rbindlist(lapply(primary,function(pop) data.table(all_scores[,.(Sample,Cohort,EPB41L3)],Population=pop,Score=all_scores[[pop]])))
long[,Population:=factor(Population,levels=primary)]
long[,Cohort:=factor(Cohort,levels=c('TCGA','GSE44076','GSE20916'))]
# 为避免跨平台绝对量比较，散点坐标用各队列内百分位秩。
long[,`:=`(Expression_rank=rank(EPB41L3)/.N,Score_rank=rank(Score)/.N),by=.(Cohort,Population)]
p2<-ggplot(long,aes(Expression_rank,Score_rank,color=Cohort))+geom_point(alpha=.4,size=.8)+
 geom_smooth(method='lm',formula=y~x,se=FALSE,linewidth=.65)+facet_grid(Population~Cohort)+scale_color_manual(values=cols)+
 labs(x='EPB41L3 expression: within-cohort rank',y='MCP-counter score: within-cohort rank')+theme(legend.position='none',strip.background=element_blank(),strip.text=element_text(size=9))
combined<-(p1/p2)+plot_layout(heights=c(1,2))+plot_annotation(tag_levels='A')
ggsave(file.path(root,'Figure_bulk_microenvironment.pdf'),combined,width=10,height=10)
ggsave(file.path(root,'Figure_bulk_microenvironment.png'),combined,width=10,height=10,dpi=300)
sr[,Model:=factor(Model,levels=rev(names(formulas)))]
p3<-ggplot(sr,aes(HR,Model,color=Endpoint))+geom_vline(xintercept=1,linetype=2,color='grey55')+
 geom_errorbar(aes(xmin=Low,xmax=High),orientation='y',width=.2,position=position_dodge(.5))+
 geom_point(position=position_dodge(.5),size=2.5)+scale_x_log10()+scale_color_manual(values=c(OS='#387CA3',PFI='#CF8744'))+
 labs(title='Exploratory prognosis analysis',x='Hazard ratio per 1 SD of EPB41L3 (95% CI)',y=NULL)
ggsave(file.path(root,'Figure_prognosis_supplement.pdf'),p3,width=9,height=4.5)
ggsave(file.path(root,'Figure_prognosis_supplement.png'),p3,width=9,height=4.5,dpi=300)
writeLines(capture.output(sessionInfo()),file.path(root,'results/sessionInfo.txt'))
print(main);print(ss);print(cr);print(sr);print(rbindlist(spline_rows))
