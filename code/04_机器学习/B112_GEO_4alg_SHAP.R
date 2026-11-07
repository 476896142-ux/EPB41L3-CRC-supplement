# 四算法评估 10 个 PPI 交集基因；保留原脚本，本文件独立运行。
# 训练集：患者分组 10 折 × 5 次交叉验证；留出集：仅用于评估及解释。
suppressPackageStartupMessages({
  library(caret); library(dplyr); library(tidyr); library(ggplot2)
  library(pROC); library(fastshap); library(ggbeeswarm)
})
project <- 'D:/Bioinfor/colorectal/GEO_ML_hub10'
out <- 'E:/codex/2026-10-31/dui/outputs/GEO_ML_hub10'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
write_csv <- function(x,name) write.csv(x,file.path(out,paste0(name,'.csv')),row.names=FALSE)
save_plot <- function(p,name,w=8,h=6) {
  ggsave(file.path(out,paste0(name,'.png')),p,width=w,height=h,dpi=300,bg='white')
  ggsave(file.path(out,paste0(name,'.pdf')),p,width=w,height=h,device=cairo_pdf)
}

# 开发队列仅用 GSE44076；外部队列 GSE20916 不参与以下划分和调参。
genes <- c('NTRK3','EPHB1','KLB','GREM1','HHIP','ANGPTL1','LAMC2','WDR72','GPER1','EPB41L3')
meta <- read.csv(file.path(out,'GSE44076_metadata.csv'),check.names=FALSE)
expr <- read.csv(file.path(out,'GSE44076_expression.csv'),row.names=1,check.names=FALSE)
stopifnot(identical(meta$Sample,rownames(expr)),all(genes %in% names(expr)))
X <- expr[,genes,drop=FALSE]; y <- factor(meta$Class,levels=c('Tumor','Normal'));patient <- meta$Patient
stopifnot(nrow(X)==246,all(is.finite(as.matrix(X))),all(vapply(X,sd,numeric(1))>0))
patient_class <- tapply(as.character(y),patient,function(z) paste(sort(unique(z)),collapse='+'))
set.seed(20261101)
train_patients <- unlist(lapply(split(names(patient_class),patient_class),function(z) sample(z,floor(.8*length(z)))))
train_i <- which(patient %in% train_patients); test_i <- which(!patient %in% train_patients)
xtr <- X[train_i,,drop=FALSE]; xte <- X[test_i,,drop=FALSE]; ytr <- y[train_i]; yte <- y[test_i]
stopifnot(length(intersect(patient[train_i],patient[test_i]))==0,all(table(yte)>0))
write_csv(data.frame(Sample=rownames(X),Patient=patient,Class=y,Set=ifelse(seq_along(y)%in%train_i,'Train','Holdout')),'sample_split')
cat('训练/留出样本：',length(train_i),'/',length(test_i),'\n'); print(table(Set=ifelse(seq_along(y)%in%train_i,'Train','Holdout'),Class=y))
# 2. 四模型共享相同的患者分组折；标准化仅在每个训练折内拟合。
index <- list(); index_out <- list(); set.seed(123)
for(r in 1:5) {
 fold <- setNames(integer(length(train_patients)),train_patients)
 strata <- split(train_patients,patient_class[train_patients])
 for(z in strata) fold[z] <- sample(rep(1:10,length.out=length(z)))
 for(k in 1:10) {
  nm <- paste0('Fold',k,'.Rep',r)
  validation <- which(fold[patient[train_i]]==k)
  index_out[[nm]] <- validation; index[[nm]] <- setdiff(seq_along(ytr),validation)
  stopifnot(all(table(ytr[validation])>0),length(intersect(patient[train_i][validation],patient[train_i][index[[nm]]]))==0)
 }
}
ctrl <- trainControl(method='repeatedcv',number=10,repeats=5,index=index,indexOut=index_out,
 classProbs=TRUE,summaryFunction=twoClassSummary,savePredictions='final',allowParallel=FALSE)
saveRDS(list(train=train_i,test=test_i,index=index,indexOut=index_out),file.path(out,'resampling_indices.rds'))
grids <- list(glmnet=expand.grid(alpha=c(0,.5,1),lambda=10^seq(-4,-1,length.out=10)),
 rf=expand.grid(mtry=c(2,4,6)),svmRadial=expand.grid(sigma=c(.01,.05,.1),C=c(.5,1,2)),
 nnet=expand.grid(size=c(3,5,8),decay=c(.1,.5,1)))
fits <- list(); performance <- list(); predictions <- list()
auc <- function(obs,prob) as.numeric(pROC::auc(pROC::roc(obs,prob,levels=c('Normal','Tumor'),direction='<',quiet=TRUE)))
# 3. 必须四模型全部成功；不再静默跳过失败模型。
for(m in names(grids)) {
 cat('开始训练：',m,'\n'); set.seed(123)
 extra <- switch(m,rf=list(ntree=500),nnet=list(trace=FALSE,maxit=500,MaxNWts=2000),list())
 fit <- do.call(caret::train,c(list(x=xtr,y=ytr,method=m,trControl=ctrl,metric='ROC',
  tuneGrid=grids[[m]],preProcess=c('center','scale')),extra))
 fits[[m]] <- fit; saveRDS(fit,file.path(out,paste0('fit_',m,'.rds')))
 prob <- predict(fit,xte,type='prob')[,'Tumor']
 cls <- factor(ifelse(prob>=.5,'Tumor','Normal'),levels=levels(y))
 performance[[m]] <- data.frame(Model=m,CV_AUC_mean=mean(fit$resample$ROC),CV_AUC_SD=sd(fit$resample$ROC),
  Holdout_AUC=auc(yte,prob),Sensitivity=mean(cls[yte=='Tumor']=='Tumor'),Specificity=mean(cls[yte=='Normal']=='Normal'))
 predictions[[m]] <- data.frame(Sample=rownames(xte),Observed=yte,Probability=prob,Model=m)
 write_csv(fit$results,paste0('tuning_',m)); write_csv(fit$bestTune,paste0('best_parameters_',m))
 cat(m,'完成，留出集 AUC=',round(auc(yte,prob),4),'\n')
}
perf <- bind_rows(performance); pred <- bind_rows(predictions)
write_csv(perf,'model_performance'); write_csv(pred,'holdout_predictions')
# 4. 使用同一留出集置换重要性（AUC 降幅），避免把 SVM 的单基因 ROC 当模型重要性。
# 解释结果不参与调参；配对患者共享同一数据划分。保留负值，表示没有正向贡献证据。
set.seed(456); importance <- list()
for(m in names(fits)) for(gene in genes) {
 delta <- replicate(20,{
  xp <- xte; xp[[gene]] <- sample(xp[[gene]])
  perf$Holdout_AUC[perf$Model==m]-auc(yte,predict(fits[[m]],xp,type='prob')[,'Tumor'])
 })
 importance[[paste(m,gene)]] <- data.frame(Model=m,Gene=gene,Importance=mean(delta),PermutationSD=sd(delta))
}
imp <- bind_rows(importance)
ranked <- imp %>% group_by(Gene) %>% summarise(MeanImportance=mean(Importance),N_models=n_distinct(Model),.groups='drop') %>% arrange(desc(MeanImportance))
write_csv(imp,'permutation_importance_by_model'); write_csv(ranked,'ML_gene_ranking')
# PPI 已用于候选筛选，等权融合只作探索性补充，不按指定基因排名调整权重。
ppi <- read.csv('E:/codex/2026-10-31/dui/outputs/Hub_Degree_MCC_MNC/consensus_hubs.csv')
combined <- ranked %>% mutate(Degree=ppi$Degree[match(Gene,ppi$node_name)]) %>%
 mutate(ML_z=as.numeric(scale(MeanImportance)),PPI_z=as.numeric(scale(Degree)),Combined=ML_z+PPI_z) %>% arrange(desc(Combined))
stopifnot(!anyNA(combined),nrow(combined)==10)
write_csv(combined,'ML_PPI_equal_weight_exploratory'); saveRDS(combined,file.path(out,'combined_gene_score.rds'))
# 5. 留出集 ROC、统一定义的重要性和探索性融合排名。
curves <- lapply(names(fits),function(m) {
 z <- pred[pred$Model==m,]; rr <- pROC::roc(z$Observed,z$Probability,levels=c('Normal','Tumor'),direction='<',quiet=TRUE)
 data.frame(FPR=1-rr$specificities,TPR=rr$sensitivities,Model=m)
})
cols <- c(glmnet='#368A83',rf='#557CA4',svmRadial='#CC7964',nnet='#9278AA')
labels <- setNames(sprintf('%s (AUC %.3f)',perf$Model,perf$Holdout_AUC),perf$Model)
p <- ggplot(bind_rows(curves),aes(FPR,TPR,color=Model))+geom_path(linewidth=.9)+
 geom_abline(slope=1,intercept=0,linetype=2,color='grey70')+coord_equal()+
 scale_color_manual(values=cols,labels=labels)+labs(title='Patient-held-out ROC curves',x='1 - Specificity',y='Sensitivity')+theme_bw()
save_plot(p,'Fig1_ROC_4models')
p <- ggplot(ranked,aes(MeanImportance,reorder(Gene,MeanImportance)))+geom_col(fill='#557CA4',width=.7)+
 labs(title='Mean held-out permutation importance',x='Mean decrease in AUC (four models)',y=NULL)+theme_bw()
save_plot(p,'Fig2_ML_importance')
p <- ggplot(combined,aes(Combined,reorder(Gene,Combined)))+geom_col(fill='#CC7964',width=.7)+
 labs(title='Exploratory equal-weight ML + PPI ranking',x='ML z-score + Degree z-score',y=NULL)+theme_bw()
save_plot(p,'Fig3_ML_PPI_exploratory')
# 6. RF 的 Monte Carlo SHAP，不是 TreeSHAP；训练集作为背景，解释留出样本。
wrapper <- function(object,newdata) predict(object,newdata,type='prob')[,'Tumor']
set.seed(789)
shap <- fastshap::explain(fits$rf,X=xtr,newdata=xte,pred_wrapper=wrapper,nsim=100,adjust=TRUE)
saveRDS(shap,file.path(out,'RF_SHAP.rds'))
shap_error <- max(abs(rowSums(shap)+attr(shap,'baseline')-wrapper(fits$rf,xte)))
stopifnot(shap_error<1e-6)
shap_df <- as.data.frame(shap) %>% mutate(Sample=rownames(xte)) %>% pivot_longer(-Sample,names_to='Gene',values_to='SHAP')
expression_df <- xte %>% mutate(Sample=rownames(xte)) %>% pivot_longer(-Sample,names_to='Gene',values_to='Expression')
shap_df <- left_join(shap_df,expression_df,by=c('Sample','Gene')) %>% group_by(Gene) %>%
 mutate(ExpressionScaled=if(max(Expression)>min(Expression)) (Expression-min(Expression))/(max(Expression)-min(Expression)) else .5) %>% ungroup()
shap_rank <- shap_df %>% group_by(Gene) %>% summarise(MeanAbsSHAP=mean(abs(SHAP)),.groups='drop') %>% arrange(MeanAbsSHAP)
shap_df$Gene <- factor(shap_df$Gene,levels=shap_rank$Gene)
write_csv(shap_df,'RF_SHAP_long'); write_csv(shap_rank %>% arrange(desc(MeanAbsSHAP)),'RF_SHAP_ranking')
p <- ggplot(shap_df,aes(SHAP,Gene,color=ExpressionScaled))+geom_vline(xintercept=0,color='grey75',linetype=2)+
 geom_quasirandom(size=1.5,width=.3,alpha=.8)+scale_color_gradient(low='#457BA0',high='#CC6F68',name='Relative\nexpression')+
 labs(title='Random forest SHAP: held-out samples',x='SHAP value for tumour probability',y=NULL)+theme_bw()
save_plot(p,'Fig4_RF_SHAP',9,7)
write_csv(data.frame(Package=c('R',names(grids),'caret','fastshap'),Version=c(as.character(getRversion()),
 vapply(c('glmnet','randomForest','kernlab','nnet','caret','fastshap'),function(z)as.character(packageVersion(z)),character(1)))),'software_versions')
saveRDS(list(fits=fits,X=X,y=y,genes=genes,cols=cols),file.path(out,'development_state.rds'))
cat('开发阶段完成；SHAP 加和误差=',shap_error,'\n'); print(perf)

# 7. 固定最佳参数后在全部 GSE44076 重拟合；仅此后读取 GSE20916 结局。
# 外部表达由原始 CEL 独立 RMA 得到，未与开发队列合并校正，也不利用外部标签调阈值。
if(!file.exists(file.path(out,'GSE20916_expression_RMA.csv'))) stop('开发阶段已完成；外部原始数据 RMA 尚未完成，运行后可从第 7 节恢复。')
validation <- read.csv(file.path(out,'GSE20916_expression_RMA.csv'),row.names=1,check.names=FALSE)
val_meta <- read.csv(file.path(project,'raw/GSE20916_metadata.csv'),check.names=FALSE)
val_meta <- val_meta[match(rownames(validation),val_meta$geo_accession),]
xval <- validation[,genes,drop=FALSE]
yval <- factor(ifelse(val_meta$tissue=='adenocarcinoma','Tumor','Normal'),levels=levels(y))
stopifnot(nrow(xval)==60,!anyNA(val_meta$geo_accession),all(is.finite(as.matrix(xval))),
 all(val_meta$tissue %in% c('adenocarcinoma','normal colon')),!any(grepl('(micro)',val_meta$title,fixed=TRUE)))
external_results <- list(); external_pred <- list()
for(m in names(fits)) {
 set.seed(123)
 extra <- switch(m,rf=list(ntree=500),nnet=list(trace=FALSE,maxit=500,MaxNWts=2000),list())
 final <- do.call(caret::train,c(list(x=X,y=y,method=m,trControl=trainControl(method='none',classProbs=TRUE),
 tuneGrid=fits[[m]]$bestTune,preProcess=c('center','scale')),extra))
 saveRDS(final,file.path(out,paste0('final_all_development_',m,'.rds')))
 pv <- predict(final,xval,type='prob')[,'Tumor']; cls <- ifelse(pv>=.5,'Tumor','Normal')
 rr <- pROC::roc(yval,pv,levels=c('Normal','Tumor'),direction='<',quiet=TRUE); ci <- pROC::ci.auc(rr)
 external_results[[m]] <- data.frame(Model=m,External_AUC=as.numeric(pROC::auc(rr)),CI_low=ci[1],CI_high=ci[3],
 Sensitivity=mean(cls[yval=='Tumor']=='Tumor'),Specificity=mean(cls[yval=='Normal']=='Normal'),
  Brier=mean((pv-as.numeric(yval=='Tumor'))^2),Probability_min=min(pv),Probability_max=max(pv),
  CI_note=if(as.numeric(pROC::auc(rr))==1) 'Degenerate DeLong CI at perfect AUC; do not interpret as certainty' else 'Sample-level DeLong CI; patient pairing unknown')
 external_pred[[m]] <- data.frame(Sample=rownames(xval),Observed=yval,Probability=pv,Model=m)
}
write_csv(bind_rows(external_results),'external_GSE20916_performance');write_csv(bind_rows(external_pred),'external_GSE20916_predictions')
curves <- lapply(external_pred,function(z) {
 rr <- pROC::roc(z$Observed,z$Probability,levels=c('Normal','Tumor'),direction='<',quiet=TRUE)
 data.frame(FPR=1-rr$specificities,TPR=rr$sensitivities,Model=z$Model[1])
})
ep <- bind_rows(external_results)
p <- ggplot(bind_rows(curves),aes(FPR,TPR,color=Model))+geom_path(linewidth=.9)+
 geom_abline(slope=1,intercept=0,linetype=2,color='grey70')+coord_equal()+
 scale_color_manual(values=cols,labels=setNames(sprintf('%s (AUC %.3f)',ep$Model,ep$External_AUC),ep$Model))+
 labs(title='External validation: GSE20916',x='1 - Specificity',y='Sensitivity')+theme_bw()
save_plot(p,'Fig5_External_ROC')
cat('GSE44076 开发 + GSE20916 外部验证全部完成。\n');print(ep)
