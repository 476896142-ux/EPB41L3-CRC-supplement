# 核查真实输入、患者划分和外部概率尺度；不改变模型或阈值。
suppressPackageStartupMessages(library(ggplot2))
out <- 'E:/codex/2026-10-31/dui/outputs/GEO_ML_hub10'
development <- read.csv(file.path(out,'GSE44076_expression.csv'),row.names=1)
external <- read.csv(file.path(out,'GSE20916_expression_RMA.csv'),row.names=1)
split <- read.csv(file.path(out,'sample_split.csv'))
state <- readRDS(file.path(out,'development_state.rds'))
stopifnot(nrow(development)==246,nrow(external)==60,ncol(external)==10,
 identical(names(development),names(external)),all(is.finite(as.matrix(external))),
 !anyDuplicated(rownames(external)),length(intersect(rownames(development),rownames(external)))==0,
 length(intersect(split$Patient[split$Set=='Train'],split$Patient[split$Set=='Holdout']))==0)
for(model in names(state$fits)) stopifnot(identical(setdiff(names(state$fits[[model]]$trainingData),'.outcome'),names(development)))
audit <- data.frame(Gene=names(development),Development_median=sapply(development,median),
 External_median=sapply(external,median),Development_SD=sapply(development,sd),External_SD=sapply(external,sd))
audit$External_median_in_development_z <- (audit$External_median-sapply(development,mean))/audit$Development_SD
write.csv(audit,file.path(out,'cross_platform_expression_audit.csv'),row.names=FALSE)
pred <- read.csv(file.path(out,'external_GSE20916_predictions.csv'))
stopifnot(nrow(pred)==240,all(table(pred$Model)==60),all(pred$Probability>=0 & pred$Probability<=1))
pred$Observed <- factor(pred$Observed,levels=c('Normal','Tumor'))
colors <- c(glmnet='#368A83',rf='#557CA4',svmRadial='#CC7964',nnet='#9278AA')
p <- ggplot(pred,aes(Observed,Probability,color=Model))+
 geom_hline(yintercept=.5,linetype=2,color='grey40')+
 geom_boxplot(width=.45,outlier.shape=NA,alpha=.2)+
 geom_point(position=position_jitter(width=.12,height=0,seed=42),size=1.2,alpha=.7)+
 facet_wrap(~Model,nrow=1)+scale_color_manual(values=colors,guide='none')+
 scale_y_continuous(limits=c(0,1))+
 labs(title='GSE20916: predicted tumour probability',subtitle='Dashed line: fixed classification threshold of 0.5',
      x=NULL,y='Predicted probability')+theme_bw(base_size=11)
ggsave(file.path(out,'Fig6_External_probabilities.png'),p,width=10,height=4,dpi=300,bg='white')
ggsave(file.path(out,'Fig6_External_probabilities.pdf'),p,width=10,height=4,device=cairo_pdf)
cat('输入维度、10 基因、患者分组和外部预测检查通过。\n');print(audit)
