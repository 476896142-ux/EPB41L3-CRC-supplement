# 十基因精简：所有排序、标准化及调参均放在训练折内。
project <- 'D:/Bioinfor/colorectal/GEO_ML_hub10'
input <- 'E:/codex/2026-10-31/dui/outputs/GEO_ML_hub10'
out <- file.path(input,'nested_selection')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
source(file.path(project,'selection_helpers.R'))
write_result <- function(value,name) write.csv(value,file.path(out,paste0(name,'.csv')),row.names=FALSE)
genes <- c('NTRK3','EPHB1','KLB','GREM1','HHIP','ANGPTL1','LAMC2','WDR72','GPER1','EPB41L3')
expression <- as.matrix(read.csv(file.path(input,'GSE44076_expression.csv'),row.names=1)[,genes])
metadata <- read.csv(file.path(input,'GSE44076_metadata.csv'))
outcome <- factor(metadata$Class,levels=c('Normal','Tumor')); patient <- metadata$Patient
stopifnot(identical(rownames(expression),metadata$Sample),nrow(expression)==246,
          length(unique(patient))==148,all(is.finite(expression)))
strategies <- c('LASSO_1SE','SVMRFE_1SE','Adaptive_1SE','Panel3to5','SVM10')
outer_folds <- group_folds(patient,outcome,repeats=5,seed=20261102)
saveRDS(outer_folds,file.path(out,'outer_folds.rds'))
outer_predictions <- list(); selected_rows <- list(); inner_summaries <- list(); membership <- list()

# 1. 外层只用于评估；每次外层训练中重新进行内层选择。
for(index in seq_along(outer_folds)) {
  fold <- outer_folds[[index]]; fold_name <- names(outer_folds)[index]
  cat('外层',index,'/',length(outer_folds),fold_name,'\n')
  training <- expression[fold$train,,drop=FALSE]; y_training <- outcome[fold$train]
  tuning <- inner_cv(training,y_training,patient[fold$train],repeats=2,seed=1000+index*100)
  saveRDS(tuning,file.path(out,paste0('inner_',fold_name,'.rds')))
  inner_summaries[[index]] <- cbind(OuterFold=fold_name,tuning$summary)
  path <- fit_paths(training,y_training,seed=5000+index)
  for(strategy in strategies) {
    choice <- choose_one_se(tuning$summary,strategy)
    fitted <- path$metadata[path$metadata$ID==choice$ID,]
    stopifnot(fitted$NGenes>0)
    probability <- predict_path(path,expression[fold$test,,drop=FALSE],choice$ID)
    selected_rows[[length(selected_rows)+1]] <- cbind(OuterFold=fold_name,Repeat=fold$repeat_id,
      choice,ActualNGenes=fitted$NGenes,SelectedGenes=fitted$Genes)
    outer_predictions[[length(outer_predictions)+1]] <- data.frame(OuterFold=fold_name,Repeat=fold$repeat_id,
      Strategy=strategy,Sample=rownames(expression)[fold$test],Patient=patient[fold$test],
      Observed=outcome[fold$test],Probability=probability)
    membership[[length(membership)+1]] <- data.frame(OuterFold=fold_name,Strategy=strategy,Gene=genes,
      Selected=genes %in% strsplit(fitted$Genes,';',fixed=TRUE)[[1]])
  }
}
predictions <- bind_rows(outer_predictions); selections <- bind_rows(selected_rows)
write_result(predictions,'nested_oof_predictions'); write_result(selections,'outer_selected_models')
write_result(bind_rows(inner_summaries),'all_inner_candidate_summaries')
write_result(bind_rows(membership),'outer_gene_membership')
stability <- bind_rows(membership) %>% group_by(Strategy,Gene) %>%
  summarise(Selected_count=sum(Selected),Outer_fits=n(),Frequency=mean(Selected),.groups='drop')
write_result(stability,'gene_selection_stability')
# 每次重复中，每份样本恰好有一次折外预测；不能把五次重复当作五倍样本量。
repeat_metrics <- predictions %>% group_by(Strategy,Repeat) %>% group_modify(~{
  stopifnot(nrow(.x)==246,!anyDuplicated(.x$Sample))
  metrics(.x$Observed,.x$Probability)
}) %>% ungroup()
write_result(repeat_metrics,'nested_per_repeat_metrics')
performance <- repeat_metrics %>% group_by(Strategy) %>%
  summarise(across(c(AUC,Sensitivity,Specificity,Brier),list(mean=mean,SD=sd)),.groups='drop')
write_result(performance,'nested_performance_summary')

# 2. 用全开发队列执行相同选择规则，保存锁定模型后才读取外部表达和标签。
cat('外层完成；开始全开发队列的最终内层交叉验证。\n')
final_tuning <- inner_cv(expression,outcome,patient,repeats=5,seed=9000)
saveRDS(final_tuning,file.path(out,'final_inner_cv.rds'))
write_result(final_tuning$summary,'final_candidate_summary')
final_path <- fit_paths(expression,outcome,seed=10000)
final_choices <- bind_rows(lapply(strategies,function(strategy) choose_one_se(final_tuning$summary,strategy)))
positions <- match(final_choices$ID,final_path$metadata$ID)
final_choices$ActualNGenes <- final_path$metadata$NGenes[positions]
final_choices$SelectedGenes <- final_path$metadata$Genes[positions]
stopifnot(all(final_choices$ActualNGenes>0),
 final_choices$ActualNGenes[final_choices$Strategy=='Panel3to5'] %in% 3:5)
write_result(final_choices,'locked_final_selection')
saveRDS(list(path=final_path,choices=final_choices,locked_at=Sys.time(),genes=genes),
        file.path(out,'locked_development_models.rds'))
final_membership <- bind_rows(lapply(strategies,function(strategy) data.frame(Strategy=strategy,Gene=genes,
 Selected=genes %in% strsplit(final_choices$SelectedGenes[final_choices$Strategy==strategy],';',fixed=TRUE)[[1]])))
write_result(final_membership,'final_gene_membership')

# 3. GSE20916 仅作复用队列验证，不据其结果重新筛选或修改阈值。
external <- as.matrix(read.csv(file.path(input,'GSE20916_expression_RMA.csv'),row.names=1)[,genes])
external_meta <- read.csv(file.path(project,'raw/GSE20916_metadata.csv'))
external_meta <- external_meta[match(rownames(external),external_meta$geo_accession),]
external_outcome <- factor(ifelse(external_meta$tissue=='adenocarcinoma','Tumor','Normal'),levels=levels(outcome))
stopifnot(nrow(external)==60,all(is.finite(external)),!anyNA(external_meta$geo_accession),
 all(external_meta$tissue %in% c('adenocarcinoma','normal colon')),
 !any(grepl('(micro)',external_meta$title,fixed=TRUE)))
external_predictions <- list(); external_metrics <- list()
for(strategy in strategies) {
  id <- final_choices$ID[final_choices$Strategy==strategy]
  probability <- predict_path(final_path,external,id)
  external_predictions[[strategy]] <- data.frame(Strategy=strategy,Sample=rownames(external),
    Observed=external_outcome,Probability=probability)
  external_metrics[[strategy]] <- cbind(Strategy=strategy,metrics(external_outcome,probability))
}
write_result(bind_rows(external_predictions),'external_reused_cohort_predictions')
write_result(bind_rows(external_metrics),'external_reused_cohort_metrics')
packages <- c('glmnet','e1071','dplyr','ggplot2')
write_result(data.frame(Package=c('R',packages),Version=c(as.character(getRversion()),
 vapply(packages,function(package) as.character(packageVersion(package)),character(1)))),'software_versions')
cat('全部完成；锁定组合：\n');print(final_choices[,c('Strategy','Method','ActualNGenes','SelectedGenes')])
print(performance);print(bind_rows(external_metrics))
