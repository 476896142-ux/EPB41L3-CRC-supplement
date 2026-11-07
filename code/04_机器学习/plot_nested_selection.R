# 结果可视化与实数据核查；不从外部结果反向选择基因。
suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(patchwork)})
project <- 'D:/Bioinfor/colorectal/GEO_ML_hub10'
input <- 'E:/codex/2026-10-31/dui/outputs/GEO_ML_hub10'
out <- file.path(input,'nested_selection')
source(file.path(project,'selection_helpers.R'))
read_result <- function(name) read.csv(file.path(out,paste0(name,'.csv')))
save_figure <- function(plot,name,width,height) {
  ggsave(file.path(out,paste0(name,'.png')),plot,width=width,height=height,dpi=300,bg='white')
  ggsave(file.path(out,paste0(name,'.pdf')),plot,width=width,height=height,device=cairo_pdf)
}
stability <- read_result('gene_selection_stability')
choices <- read_result('locked_final_selection')
candidates <- read_result('final_candidate_summary')
performance <- read_result('nested_performance_summary')
external <- read_result('external_reused_cohort_metrics')
outer_choices <- read_result('outer_selected_models')
predictions <- read_result('nested_oof_predictions')
metadata <- read.csv(file.path(input,'GSE44076_metadata.csv'))
expression <- as.matrix(read.csv(file.path(input,'GSE44076_expression.csv'),row.names=1))
folds <- readRDS(file.path(out,'outer_folds.rds'))
correlations <- list()

# 验证全部 25 个外层及 250 个内层划分的患者隔离和选择规则。
for(index in seq_along(folds)) {
  fold <- folds[[index]]; name <- names(folds)[index]
  stopifnot(length(intersect(metadata$Patient[fold$train],metadata$Patient[fold$test]))==0)
  correlation <- as.data.frame(as.table(cor(expression[fold$train,,drop=FALSE],method='spearman')))
  names(correlation) <- c('Gene1','Gene2','Spearman_rho')
  correlations[[index]] <- cbind(OuterFold=name,correlation)
  inner <- readRDS(file.path(out,paste0('inner_',name,'.rds')))
  for(split in inner$folds) stopifnot(length(intersect(metadata$Patient[fold$train][split$train],
                                                     metadata$Patient[fold$train][split$test]))==0)
  for(strategy in unique(outer_choices$Strategy)) {
    expected <- choose_one_se(inner$summary,strategy)$ID
    observed <- outer_choices$ID[outer_choices$OuterFold==name & outer_choices$Strategy==strategy]
    stopifnot(identical(expected,observed))
  }
}
write.csv(bind_rows(correlations),file.path(out,'outer_training_gene_correlations.csv'),row.names=FALSE)
stopifnot(nrow(predictions)==246*5*5,all(table(predictions$Sample,predictions$Strategy)==5),
          all(performance$AUC_mean>=0 & performance$AUC_mean<=1))
write.csv(data.frame(Check=c('outer_patient_separation','inner_patient_separation','inner_only_selection_rule',
                            'one_oof_prediction_per_sample_per_repeat','ten_initial_candidates'),
                     Passed=TRUE,Count=c(25,250,125,246*5*5,10)),file.path(out,'verification_checks.csv'),row.names=FALSE)

# A：稳定性；频率来自重复外层训练，不能解释为统计显著性或因果概率。
order_genes <- stability %>% filter(Strategy=='Panel3to5') %>% arrange(Frequency,Gene) %>% pull(Gene)
display <- stability %>% filter(Strategy %in% c('LASSO_1SE','SVMRFE_1SE','Panel3to5'))
display$Gene <- factor(display$Gene,levels=order_genes)
display$Strategy <- factor(display$Strategy,levels=c('LASSO_1SE','SVMRFE_1SE','Panel3to5'),
                           labels=c('LASSO\n1-SE','SVM-RFE\n1-SE','SVM-RFE\n3-5 genes'))
p_stability <- ggplot(display,aes(Strategy,Gene,fill=Frequency))+
 geom_tile(color='white',linewidth=.5)+geom_text(aes(label=paste0(round(Frequency*100),'%')),size=3.3)+
 scale_fill_gradient(low='#F4F6F7',high='#57A497',limits=c(0,1),breaks=c(0,.5,1),name='Selection frequency')+
 labs(title='Selection stability',subtitle='25 outer-training fits',x=NULL,y=NULL)+theme_minimal(base_size=11)+
 theme(panel.grid=element_blank(),legend.position='bottom')

# B：所有候选 C 值的内层曲线；阴影表示用户提出的 3–5 基因范围。
svm <- candidates[candidates$Method=='SVM_RFE',]
p_size <- ggplot(svm,aes(K,MeanAUC,color=factor(C),group=C))+
 annotate('rect',xmin=2.5,xmax=5.5,ymin=-Inf,ymax=Inf,fill='#DDECE8',alpha=.4)+
 geom_errorbar(aes(ymin=pmax(0,MeanAUC-SE_AUC),ymax=pmin(1,MeanAUC+SE_AUC)),width=.12,alpha=.5)+
 geom_line(linewidth=.7)+geom_point(size=1.7)+
 geom_point(data=svm[svm$ID==choices$ID[choices$Strategy=='Panel3to5'],],
            shape=21,fill='white',color='black',size=3)+scale_x_continuous(breaks=1:10)+
 scale_color_manual(values=c('#368A83','#557CA4','#CC7964'),name='SVM C')+
 labs(title='SVM-RFE size comparison',subtitle='Development inner CV; mean AUC +/- fold SE',
      x='Number of retained genes',y='Inner-CV AUC')+theme_bw(base_size=11)+theme(legend.position='bottom')

# C：外层验证评价筛选流程；各折的实际基因可以不同。
labels <- c(LASSO_1SE='LASSO 1-SE',SVMRFE_1SE='SVM-RFE 1-SE',Adaptive_1SE='Adaptive 1-SE',
            Panel3to5='SVM-RFE 3-5',SVM10='SVM 10 genes')
performance$Display <- factor(labels[performance$Strategy],levels=unname(labels))
p_outer <- ggplot(performance,aes(Display,AUC_mean))+
 geom_errorbar(aes(ymin=AUC_mean-AUC_SD,ymax=AUC_mean+AUC_SD),width=.15,color='#557CA4')+
 geom_point(size=2.5,color='#557CA4')+coord_flip()+
 labs(title='Nested evaluation of selection workflows',subtitle='Mean +/- SD across five repeats; not a confidence interval',
      x=NULL,y='Outer out-of-fold AUC')+theme_bw(base_size=11)
figure <- (p_stability | p_size) / p_outer + plot_layout(heights=c(1.4,1))+
 plot_annotation(tag_levels='A',title='Reduction of 10 candidate hub genes')
save_figure(figure,'Fig1_Nested_feature_selection',12,10)

# LASSO 参数与基因数一起展示；不会把岭回归写成 LASSO。
lasso <- candidates[candidates$Method=='LASSO',]
p_lasso <- ggplot(lasso,aes(log10(Lambda),MeanAUC))+
 geom_errorbar(aes(ymin=pmax(0,MeanAUC-SE_AUC),ymax=pmin(1,MeanAUC+SE_AUC)),width=.03,color='grey60')+
 geom_line(color='#368A83')+geom_point(aes(color=MeanGenes),size=2)+
 geom_vline(xintercept=log10(choices$Lambda[choices$Strategy=='LASSO_1SE']),linetype=2,color='#CC7964')+
 scale_color_gradient(low='#D5E8E4',high='#26796E',name='Mean nonzero\ngenes')+
 labs(title='LASSO development cross-validation',x='log10(lambda)',y='Inner-CV AUC')+theme_bw()
save_figure(p_lasso,'Fig2_LASSO_selection',8,5)

# 固定组合的外部概率：保留 0.5 阈值失效的真实现象。
ext_pred <- read_result('external_reused_cohort_predictions') %>%
 filter(Strategy %in% c('Adaptive_1SE','Panel3to5','SVM10'))
ext_pred$Strategy <- factor(ext_pred$Strategy,levels=c('Adaptive_1SE','Panel3to5','SVM10'),
                            labels=c('Locked 1-gene model','Locked 3-gene panel','Locked 10-gene model'))
ext_pred$Observed <- factor(ext_pred$Observed,levels=c('Normal','Tumor'))
p_external <- ggplot(ext_pred,aes(Observed,Probability,color=Strategy))+
 geom_hline(yintercept=.5,linetype=2,color='grey40')+geom_boxplot(outlier.shape=NA,width=.5)+
 geom_point(position=position_jitter(width=.12,height=0,seed=5),size=1.3,alpha=.7)+
 facet_wrap(~Strategy,nrow=1)+scale_y_continuous(limits=c(0,1))+
 scale_color_manual(values=c('#557CA4','#368A83','#CC7964'),guide='none')+
 labs(title='GSE20916: reused external cohort',subtitle='Models locked before external evaluation; threshold = 0.5',
      x=NULL,y='Predicted tumour probability')+theme_bw()
save_figure(p_external,'Fig3_Reused_external_probabilities',10,4)
cat('25 外层、250 内层患者隔离及全部选择规则检查通过；图件完成。\n')
