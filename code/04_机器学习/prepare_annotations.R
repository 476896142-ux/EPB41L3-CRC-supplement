# 固定候选基因，保留仅映射到单一 SYMBOL 的探针；不使用分组或差异结果筛探针。
suppressPackageStartupMessages({library(AnnotationDbi);library(hgu219.db);library(hgu133plus2.db)})
raw_dir <- 'D:/Bioinfor/colorectal/GEO_ML_hub10/raw'
for(platform in c('hgu219','hgu133plus2')) {
  database <- get(paste0(platform,'.db'))
  mapping <- unique(AnnotationDbi::select(database,keys=keys(database,keytype='PROBEID'),
                                        keytype='PROBEID',columns='SYMBOL'))
  mapping <- mapping[!is.na(mapping$SYMBOL) & nzchar(mapping$SYMBOL),]
  counts <- table(mapping$PROBEID)
  mapping <- mapping[mapping$PROBEID %in% names(counts[counts==1]),]
  write.csv(mapping,file.path(raw_dir,paste0(platform,'_unique_probe_map.csv')),row.names=FALSE)
}
