# 基于已锁定结果绘图：不重新选基因，不改变概率、阈值或原始统计量。
suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(tidyr)
 library(patchwork);library(glmnet);library(pROC);library(ggbeeswarm)})
base <- 'E:/codex/2026-10-31/dui/outputs/GEO_ML_hub10'
selection <- file.path(base,'nested_selection'); out <- file.path(selection,'publication')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
read_table <- function(folder,name) read.csv(file.path(folder,paste0(name,'.csv')))
write_source <- function(table,name) write.csv(table,file.path(out,paste0('Source_',name,'.csv')),row.names=FALSE)
theme_set(theme_classic(base_size=8,base_family='Arial')+
 theme(plot.title=element_text(size=9,face='bold'),plot.subtitle=element_text(size=7),
       legend.title=element_text(size=7),legend.text=element_text(size=6.5),
       legend.position='bottom',plot.margin=margin(5,5,5,5)))
save_figure <- function(plot,name,width,height) {
 ggsave(file.path(out,paste0(name,'.pdf')),plot,width=width,height=height,units='mm',device=cairo_pdf)
 ggsave(file.path(out,paste0(name,'.png')),plot,width=width,height=height,units='mm',dpi=600,bg='white')
 ggsave(file.path(out,paste0(name,'.tiff')),plot,width=width,height=height,units='mm',dpi=600,compression='lzw',bg='white')
}
locked <- readRDS(file.path(selection,'locked_development_models.rds'))
choice <- locked$choices; candidates <- read_table(selection,'final_candidate_summary')
panel <- strsplit(choice$SelectedGenes[choice$Strategy=='Panel3to5'],';',fixed=TRUE)[[1]]
stopifnot(setequal(panel,c('EPB41L3','KLB','ANGPTL1')))
lambda <- sort(unique(candidates$Lambda[!is.na(candidates$Lambda)]),decreasing=TRUE)
coefficients <- as.matrix(coef(locked$path$lasso,s=lambda))[-1,,drop=FALSE]
paths <- data.frame(Gene=rep(rownames(coefficients),times=length(lambda)),
 Lambda=rep(lambda,each=nrow(coefficients)),Coefficient=as.vector(coefficients))
paths$Highlight <- ifelse(paths$Gene %in% panel,paths$Gene,'Other candidates')
gene_colors <- c(EPB41L3='#247C74',KLB='#47769E',ANGPTL1='#9B719D','Other candidates'='#B5BBC1')
selected_lambda <- choice$Lambda[choice$Strategy=='LASSO_1SE']
write_source(paths,'A_LASSO_paths')
pA <- ggplot(paths,aes(log10(Lambda),Coefficient,group=Gene,color=Highlight))+
 geom_line(linewidth=.45)+geom_vline(xintercept=log10(selected_lambda),linetype=2,color='grey30')+
 scale_color_manual(values=gene_colors,name=NULL)+guides(color=guide_legend(nrow=2,byrow=TRUE))+
 labs(title='LASSO coefficient paths',subtitle='Dashed line: selected lambda; 1 retained gene',
      x=expression(log[10](lambda)),y='Standardized coefficient')
lasso <- candidates[candidates$Method=='LASSO',];write_source(lasso,'B_LASSO_CV')
pB <- ggplot(lasso,aes(log10(Lambda),MeanAUC))+
 geom_errorbar(aes(ymin=pmax(0,MeanAUC-SE_AUC),ymax=pmin(1,MeanAUC+SE_AUC)),width=.03,color='grey65')+
 geom_line(linewidth=.4,color='#47769E')+geom_point(size=1,color='#47769E')+
 geom_vline(xintercept=log10(selected_lambda),linetype=2,color='#B7644C')+
 labs(title='LASSO cross-validation',subtitle='Unconstrained 1-SE selection: EPB41L3',
      x=expression(log[10](lambda)),y='Inner-CV AUC (mean +/- SE)')
svm <- candidates[candidates$Method=='SVM_RFE',];write_source(svm,'C_SVM_RFE_CV')
pC <- ggplot(svm,aes(K,MeanAUC,color=factor(C),linetype=factor(C)))+
 annotate('rect',xmin=2.5,xmax=5.5,ymin=-Inf,ymax=Inf,fill='#E2EFEC',alpha=.65)+
 geom_errorbar(aes(ymin=MeanAUC-SE_AUC,ymax=pmin(1,MeanAUC+SE_AUC)),width=.12,alpha=.45)+
 geom_line(linewidth=.45)+geom_point(size=1)+
 geom_point(data=svm[svm$ID==choice$ID[choice$Strategy=='Panel3to5'],],
            shape=21,fill='white',color='black',size=2.4)+
 scale_color_manual(values=c('#247C74','#47769E','#B7644C'),name='C')+
 scale_linetype_discrete(name='C')+scale_x_continuous(breaks=1:10)+
 labs(title='SVM-RFE feature selection',subtitle='Target range: 3-5 genes; selected: 3 (C = 0.1)',
      x='Number of retained genes',y='Inner-CV AUC (mean +/- SE)')
frequency <- read_table(selection,'gene_selection_stability') %>% filter(Strategy=='Panel3to5') %>%
 mutate(Final_panel=ifelse(Gene %in% panel,'Selected','Other')) %>% arrange(Frequency,Gene)
frequency$Gene <- factor(frequency$Gene,levels=frequency$Gene)
write_source(frequency,'D_selection_frequency')
pD <- ggplot(frequency,aes(Frequency,Gene,fill=Final_panel))+
 geom_col(width=.65)+geom_text(aes(label=paste0(round(Frequency*100),'%')),hjust=-.15,size=2.2)+
 scale_fill_manual(values=c(Selected='#247C74',Other='#C3C9CE'),name=NULL)+
 scale_x_continuous(limits=c(0,1.16),breaks=c(0,.5,1),labels=c('0','50','100'))+
 labs(title='Gene-selection stability',subtitle='3-5 gene strategy; 25 outer-training fits',
      x='Selection frequency (%)',y=NULL)
external <- read_table(selection,'external_reused_cohort_predictions') %>%
 filter(Strategy %in% c('Adaptive_1SE','Panel3to5'))
roc_data <- bind_rows(lapply(split(external,external$Strategy),function(data) {
 curve <- pROC::roc(data$Observed,data$Probability,levels=c('Normal','Tumor'),direction='<',quiet=TRUE)
 data.frame(Strategy=data$Strategy[1],FPR=1-curve$specificities,TPR=curve$sensitivities,AUC=as.numeric(auc(curve)))
}))
write_source(roc_data,'E_external_ROC')
roc_labels <- c(Adaptive_1SE='1 gene (AUC 0.972)',Panel3to5='3 genes (AUC 0.999)')
pE <- ggplot(roc_data,aes(FPR,TPR,color=Strategy,linetype=Strategy))+
 geom_abline(slope=1,intercept=0,color='grey75',linetype=2)+geom_path(linewidth=.6)+
 scale_color_manual(values=c(Adaptive_1SE='#47769E',Panel3to5='#247C74'),labels=roc_labels,name=NULL)+
 scale_linetype_manual(values=c(Adaptive_1SE='dashed',Panel3to5='solid'),labels=roc_labels,name=NULL)+
 coord_equal(xlim=c(0,1),ylim=c(0,1))+guides(color=guide_legend(nrow=2),linetype=guide_legend(nrow=2))+
 labs(title='External discrimination',subtitle='GSE20916 (previously evaluated cohort)',
      x='1 - Specificity',y='Sensitivity')
panel_pred <- external[external$Strategy=='Panel3to5',]
confusion <- as.data.frame(table(Observed=factor(panel_pred$Observed,levels=c('Normal','Tumor')),
 Predicted=factor(ifelse(panel_pred$Probability>=.5,'Tumor','Normal'),levels=c('Normal','Tumor'))))
confusion <- confusion %>% group_by(Observed) %>% mutate(RowPercent=100*Freq/sum(Freq)) %>% ungroup()
stopifnot(sum(confusion$Freq)==60,confusion$Freq[confusion$Observed=='Tumor' & confusion$Predicted=='Tumor']==9)
write_source(confusion,'F_confusion_matrix')
pF <- ggplot(confusion,aes(Predicted,Observed,fill=Freq))+
 geom_tile(color='white',linewidth=1)+
 geom_text(aes(label=sprintf('%d\n(%.0f%%)',Freq,RowPercent),color=Freq>15),size=3.1)+
 scale_fill_gradient(low='#F0F4F5',high='#285F72',guide='none')+
 scale_color_manual(values=c('FALSE'='#222222','TRUE'='white'),guide='none')+
 scale_y_discrete(limits=c('Tumor','Normal'))+coord_equal()+
 labs(title='Three-gene classification',subtitle='Threshold 0.5: sensitivity 25%; specificity 100%',
      x='Predicted class',y='Observed class')+theme(axis.line=element_blank(),axis.ticks=element_blank())
main <- wrap_plots(pA,pB,pC,pD,pE,pF,ncol=2,widths=c(1,1),heights=c(1,1,1))+
 plot_annotation(tag_levels='A',theme=theme(plot.tag=element_text(size=11,face='bold')))
save_figure(main,'Figure_ML_three_gene_panel',190,245)

# 补充图保留原始四模型及十基因 RF SHAP，避免冒充最终三基因模型解释。
model_colors <- c(glmnet='#247C74',rf='#47769E',svmRadial='#B7644C',nnet='#9B719D')
model_names <- c(glmnet='GLMNET',rf='RF',svmRadial='Radial SVM',nnet='NNET')
make_roc <- function(data,title,subtitle) {
 curves <- bind_rows(lapply(split(data,data$Model),function(data) {
  curve <- pROC::roc(data$Observed,data$Probability,levels=c('Normal','Tumor'),direction='<',quiet=TRUE)
  data.frame(Model=data$Model[1],FPR=1-curve$specificities,TPR=curve$sensitivities,AUC=as.numeric(auc(curve)))
 }))
 entries <- unique(curves[,c('Model','AUC')]); labels <- setNames(sprintf('%s: %.3f',model_names[entries$Model],entries$AUC),entries$Model)
 ggplot(curves,aes(FPR,TPR,color=Model,linetype=Model))+geom_abline(slope=1,intercept=0,linetype=2,color='grey75')+
 geom_path(linewidth=.6)+scale_color_manual(values=model_colors,labels=labels,name=NULL)+
 scale_linetype_manual(values=c(glmnet='solid',rf='dashed',svmRadial='dotdash',nnet='dotted'),labels=labels,name=NULL)+
 coord_equal(xlim=c(0,1),ylim=c(0,1))+guides(color=guide_legend(nrow=2),linetype=guide_legend(nrow=2))+
 labs(title=title,subtitle=subtitle,x='1 - Specificity',y='Sensitivity')
}
holdout <- read_table(base,'holdout_predictions'); old_external <- read_table(base,'external_GSE20916_predictions')
sA <- make_roc(holdout,'Four-model internal ROC','10 genes; patient-held-out samples (n = 50)')
sB <- make_roc(old_external,'Four-model external ROC','10 genes; GSE20916 (n = 60)')
importance <- read_table(base,'ML_gene_ranking')
sC <- ggplot(importance,aes(MeanImportance,reorder(Gene,MeanImportance)))+geom_col(fill='#47769E',width=.65)+
 labs(title='Permutation importance',subtitle='Mean across four 10-gene models',x='Mean decrease in held-out AUC',y=NULL)
shap <- read_table(base,'RF_SHAP_long'); shap_rank <- read_table(base,'RF_SHAP_ranking')
shap$Gene <- factor(shap$Gene,levels=rev(shap_rank$Gene))
sD <- ggplot(shap,aes(SHAP,Gene,color=ExpressionScaled))+geom_vline(xintercept=0,linetype=2,color='grey70')+
 geom_quasirandom(size=.6,width=.25,alpha=.7,orientation='y')+
 scale_color_gradient(low='#47769E',high='#B7644C',name='Relative expression')+
 labs(title='Random forest SHAP',subtitle='Original 10-gene RF; internal held-out samples',x='SHAP value (tumour probability)',y=NULL)
supplement <- wrap_plots(sA,sB,sC,sD,ncol=2,widths=c(1,1),heights=c(1,1))+plot_annotation(tag_levels='A',
 theme=theme(plot.tag=element_text(size=11,face='bold')))
save_figure(supplement,'Figure_S1_four_models_SHAP',190,190)
write_source(holdout,'S1_internal_predictions');write_source(old_external,'S1_external_predictions')
write_source(importance,'S1_importance');write_source(shap,'S1_SHAP')
cat('主图六联图、补充图四联图及全部作图数据已导出。\n')
