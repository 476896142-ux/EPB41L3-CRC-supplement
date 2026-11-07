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

state <- readRDS(file.path(out,'development_state.rds')); list2env(state,envir=environment())
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
