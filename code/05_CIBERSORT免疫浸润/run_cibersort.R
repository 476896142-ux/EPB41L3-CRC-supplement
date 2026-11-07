# TCGA原发肿瘤与癌旁组织的CIBERSORT分析；直接读取原始表达矩阵和临床注释。
suppressPackageStartupMessages(library(CIBERSORT))
root <- 'D:/Bioinfor/colorectal/TCGA_CIBERSORT_final'
dir.create(file.path(root,'results'),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(root,'reference'),showWarnings=FALSE)
wt <- function(x,n) write.csv(x,file.path(root,'results',paste0(n,'.csv')),row.names=FALSE)
input <- 'D:/Bioinfor/colorectal/TCGA_\u673a\u5668\u5b66\u4e60\u7279\u5f81\u7b5b\u9009/TCGA_RAW/HiSeqV2'
x <- as.matrix(read.delim(input,row.names=1,check.names=FALSE))
clin <- read.delim(file.path(dirname(input),'TCGA.COADREAD.sampleMap_COADREAD_clinicalMatrix'),check.names=FALSE)
m <- data.frame(Sample=colnames(x),Patient=substr(colnames(x),1,12),Type=substr(colnames(x),14,15))
m <- m[m$Type %in% c('01','11'),]
m$Group <- factor(ifelse(m$Type=='01','Tumour','Adjacent normal'),levels=c('Adjacent normal','Tumour'))
ix <- match(m$Sample,clin$sampleID)
m$Site <- factor(ifelse(clin$`_primary_disease`[ix]=='colon adenocarcinoma','COAD',ifelse(clin$`_primary_disease`[ix]=='rectum adenocarcinoma','READ',NA)))
stopifnot(!anyNA(m$Site),!anyDuplicated(paste(m$Patient,m$Type)))
x <- x[,m$Sample,drop=FALSE]
genes <- x[c('EPB41L3','KLB','ANGPTL1'),,drop=FALSE]
sig <- as.matrix(read.delim(file.path(root,'reference/LM22.txt'),row.names=1,check.names=FALSE))
stopifnot(ncol(sig)==22,!anyDuplicated(rownames(sig)),all(is.finite(sig)),all(sig>=0))
wt(data.frame(Gene=rownames(sig),Present=rownames(sig) %in% rownames(x)),'LM22_coverage')
# 输入为log2(RSEM normalized counts+1)，还原后再传入，避免自动反变换。
linear <- pmax(2^x-1,0);stopifnot(max(linear)>50,all(is.finite(linear)))
wt(data.frame(Sample=m$Sample,Patient=m$Patient,Group=m$Group,Site=m$Site),'sample_metadata')
# 仅在局部环境替换并行调度；保留已安装版本的SVR、置换及P值实现。
# 避免Windows中每次拟合重新配置multicore；不会改动安装的软件包。
env <- new.env(parent=asNamespace('CIBERSORT'))
for(n in c('cibersort','CoreAlg','doPerm')) {
 f <- get(n,asNamespace('CIBERSORT'));environment(f)<-env;assign(n,f,env)
 writeLines(deparse(get(n,asNamespace('CIBERSORT'))),file.path(root,'reference',paste0(n,'_installed.txt')))
}
env$enableParallel <- function(...) invisible(NULL)
env$future_map <- function(.x,.f,...) lapply(.x,.f,...)
# 已有结果通过样本、基因表达和输入指纹核对后复用；需要重算时设force=TRUE。
force <- '--force' %in% commandArgs(TRUE)
cache <- file.path(root,'results/cibersort.rds')
fingerprint <- tools::md5sum(c(input,file.path(root,'reference/LM22.txt')))
fpfile <- file.path(root,'reference/input_fingerprint.rds')
if(file.exists(cache) && file.exists(fpfile) && !force) {
 stopifnot(identical(unname(fingerprint),unname(readRDS(fpfile))))
 previous <- readRDS(cache)
 stopifnot(identical(m$Sample,previous$metadata$Sample),identical(as.character(m$Site),as.character(previous$metadata$Site)),identical(as.character(m$Group),as.character(previous$metadata$Group)),isTRUE(all.equal(genes,previous$genes)))
 ans <- previous$raw
 cat('Verified raw inputs and reused completed CIBERSORT deconvolution.\n')
} else {
 set.seed(20261102)
 ans <- env$cibersort(sig_matrix=sig,mixture_file=linear,perm=1000,QN=FALSE)
 saveRDS(fingerprint,fpfile)
}
stopifnot(identical(rownames(ans),m$Sample),all(is.finite(ans)),all(ans[,1:22]>=0),max(abs(rowSums(ans[,1:22])-1))<1e-8)
wt(data.frame(Sample=rownames(ans),ans,check.names=FALSE),'CIBERSORT_all_samples')
qc <- data.frame(m,Fit_P=ans[,'P-value'],Correlation=ans[,'Correlation'],RMSE=ans[,'RMSE'],Pass=ans[,'P-value']<.05)
wt(qc,'fit_QC')
saveRDS(list(metadata=m,raw=ans,genes=genes,qc=qc),file.path(root,'results/cibersort.rds'))
wt(data.frame(Package=c('R','CIBERSORT','e1071'),Version=c(as.character(getRversion()),as.character(packageVersion('CIBERSORT')),as.character(packageVersion('e1071')))),'versions')
print(table(qc$Group,qc$Pass));cat('LM22 markers:',sum(rownames(sig)%in%rownames(x)),'/',nrow(sig),'\n')
