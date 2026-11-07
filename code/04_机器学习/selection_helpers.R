# 患者分组、训练折标准化、LASSO/SVM-RFE 和一标准误差选模。
suppressPackageStartupMessages({library(glmnet);library(e1071);library(dplyr)})
lambda_grid <- 10^seq(0,-4,length.out=30)
cost_grid <- c(.1,1,10)

group_folds <- function(patient,outcome,repeats,seed) {
  set.seed(seed)
  strata <- tapply(as.character(outcome),patient,function(v) paste(sort(unique(v)),collapse='+'))
  folds <- list()
  for(repetition in seq_len(repeats)) {
    allocation <- setNames(integer(length(strata)),names(strata))
    for(ids in split(names(strata),strata)) allocation[ids] <- sample(rep(1:5,length.out=length(ids)))
    for(fold in 1:5) {
      test <- which(allocation[patient]==fold); train <- setdiff(seq_along(patient),test)
      stopifnot(length(intersect(patient[train],patient[test]))==0,all(table(outcome[test])>0))
      folds[[paste0('R',repetition,'F',fold)]] <- list(train=train,test=test,repeat_id=repetition)
    }
  }
  folds
}

auc_fixed <- function(outcome,probability) {
  positive <- outcome=='Tumor'; n1 <- sum(positive); n0 <- sum(!positive)
  (sum(rank(probability,ties.method='average')[positive])-n1*(n1+1)/2)/(n1*n0)
}
metrics <- function(outcome,probability) {
  stopifnot(all(is.finite(probability)),all(probability>=0 & probability<=1))
  data.frame(AUC=auc_fixed(outcome,probability),
    Sensitivity=mean(probability[outcome=='Tumor']>=.5),
    Specificity=mean(probability[outcome=='Normal']<.5),
    Brier=mean((probability-as.numeric(outcome=='Tumor'))^2))
}
standardize <- function(expression,center,scale) sweep(sweep(expression,2,center,'-'),2,scale,'/')

fit_paths <- function(expression,outcome,seed) {
  set.seed(seed)
  center <- colMeans(expression); scale <- apply(expression,2,sd)
  stopifnot(all(scale>0))
  standardized <- standardize(expression,center,scale)
  lasso <- glmnet(standardized,as.numeric(outcome=='Tumor'),family='binomial',alpha=1,
                  lambda=lambda_grid,standardize=FALSE,control=list(thresh=1e-9,maxit=100000))
  coefficients <- as.matrix(coef(lasso,s=lambda_grid))[-1,,drop=FALSE]
  gene_sets <- lapply(seq_along(lambda_grid),function(i) rownames(coefficients)[abs(coefficients[,i])>1e-8])
  metadata <- data.frame(ID=paste0('L',seq_along(lambda_grid)),Method='LASSO',Lambda=lambda_grid,
                         C=NA_real_,K=NA_integer_,NGenes=lengths(gene_sets),
                         Genes=vapply(gene_sets,paste,collapse=';',FUN.VALUE=character(1)))
  svm_models <- list()
  for(cost in cost_grid) {
    retained <- colnames(standardized)
    for(size in 10:1) {
      id <- paste0('S',cost,'K',size)
      fit <- e1071::svm(x=standardized[,retained,drop=FALSE],y=outcome,kernel='linear',
                        cost=cost,scale=FALSE,probability=TRUE)
      svm_models[[id]] <- fit
      metadata <- rbind(metadata,data.frame(ID=id,Method='SVM_RFE',Lambda=NA_real_,C=cost,K=size,
                                            NGenes=size,Genes=paste(retained,collapse=';')))
      # 权重在训练折标准化空间计算；并列时按基因名确定删除次序。
      if(size>1) {
        weight <- drop(t(fit$coefs)%*%fit$SV)
        retained <- retained[-order(weight^2,retained)[1]]
      }
    }
  }
  list(center=center,scale=scale,lasso=lasso,svm=svm_models,metadata=metadata)
}

predict_path <- function(path,expression,id) {
  standardized <- standardize(expression,path$center,path$scale)
  row <- path$metadata[path$metadata$ID==id,]
  if(row$Method=='LASSO') return(drop(predict(path$lasso,newx=standardized,s=row$Lambda,type='response')))
  genes <- strsplit(row$Genes,';',fixed=TRUE)[[1]]
  prediction <- predict(path$svm[[id]],standardized[,genes,drop=FALSE],probability=TRUE)
  unname(attr(prediction,'probabilities')[,'Tumor'])
}

inner_cv <- function(expression,outcome,patient,repeats,seed) {
  folds <- group_folds(patient,outcome,repeats,seed); records <- list()
  for(index in seq_along(folds)) {
    fold <- folds[[index]]
    path <- fit_paths(expression[fold$train,,drop=FALSE],outcome[fold$train],seed+index)
    for(id in path$metadata$ID) {
      probability <- predict_path(path,expression[fold$test,,drop=FALSE],id)
      records[[length(records)+1]] <- cbind(Fold=names(folds)[index],path$metadata[path$metadata$ID==id,],
                                          metrics(outcome[fold$test],probability))
    }
  }
  records <- bind_rows(records)
  summary <- records %>% group_by(ID,Method,Lambda,C,K) %>%
    summarise(MeanAUC=mean(AUC),SE_AUC=sd(AUC)/sqrt(n()),MeanGenes=mean(NGenes),
              MeanBrier=mean(Brier),MeanSensitivity=mean(Sensitivity),MeanSpecificity=mean(Specificity),.groups='drop')
  list(summary=summary,records=records,folds=folds)
}

choose_one_se <- function(summary,strategy) {
  candidates <- summary[summary$MeanGenes>0,]
  if(strategy=='LASSO_1SE') candidates <- candidates[candidates$Method=='LASSO',]
  if(strategy=='SVMRFE_1SE') candidates <- candidates[candidates$Method=='SVM_RFE',]
  if(strategy=='Panel3to5') candidates <- candidates[candidates$Method=='SVM_RFE' & candidates$K %in% 3:5,]
  if(strategy=='SVM10') candidates <- candidates[candidates$Method=='SVM_RFE' & candidates$K==10,]
  stopifnot(nrow(candidates)>0)
  best <- candidates[order(-candidates$MeanAUC,candidates$MeanBrier,candidates$ID)[1],]
  eligible <- candidates[candidates$MeanAUC>=best$MeanAUC-best$SE_AUC-1e-12,]
  regularization <- ifelse(eligible$Method=='LASSO',-eligible$Lambda,eligible$C)
  selected <- eligible[order(eligible$MeanGenes,eligible$MeanBrier,-eligible$MeanAUC,regularization,eligible$ID)[1],]
  selected$Strategy <- strategy; selected$BestAUC <- best$MeanAUC; selected$OneSE_cutoff <- best$MeanAUC-best$SE_AUC
  selected
}
