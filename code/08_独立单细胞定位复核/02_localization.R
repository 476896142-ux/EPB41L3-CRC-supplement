# GSE146771：原作者注释下的EPB41L3定位；以患者等权汇总，描述性复核。
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
out <- 'D:/Bioinfor/colorectal/single_raw/GSE146771_localization'
d <- as.data.table(read.csv(gzfile(file.path(out,'cells_fixed_markers.csv.gz')),check.names=FALSE))
genes <- c('EPB41L3','SPP1','LST1','TYROBP','FCER1G','CD68','C1QA','C1QB','C1QC','APOE','TREM2','LPL','NLRP3','IL1B','FCN1','S100A8','S100A9','EPCAM','COL1A1','CD3D','MS4A1')
stopifnot(nrow(d)==43817,!anyDuplicated(d$CellName))
write <- function(x,n) fwrite(x,file.path(out,paste0(n,'.csv')))
m <- d[Global_Cluster=='Myeloid cell']
counts <- m[,.(cells=.N),by=.(Sample,Tissue,Sub_Cluster)]
write(counts,'patient_subtype_counts')
long <- melt(m,id.vars=c('CellName','Sample','Tissue','Sub_Cluster'),measure.vars=genes,variable.name='gene',value.name='expression')
patient <- long[,.(cells=.N,detected=sum(expression>0),detection=mean(expression>0),mean_expression=mean(expression)),by=.(Sample,Tissue,Sub_Cluster,gene)]
write(patient,'patient_marker_summary')
summary <- rbindlist(lapply(c(1,10,20,30),function(cut) patient[Tissue=='T'&cells>=cut,.(min_cells=cut,patients=.N,total_cells=sum(cells),mean_expression=mean(mean_expression),detection=mean(detection),min_detection=min(detection),max_detection=max(detection)),by=.(Sub_Cluster,gene)]))
write(summary,'patient_equal_localization')
epb <- patient[gene=='EPB41L3'&Sub_Cluster=='hM13_TAM-SPP1']
write(epb,'SPP1_patient_EPB41L3')
pooled <- long[Tissue=='T'&gene%in%c('EPB41L3','SPP1'),.(cells=.N,detected=sum(expression>0),detection=mean(expression>0),mean_expression=mean(expression)),by=.(Sub_Cluster,gene)]
write(pooled,'pooled_descriptive_only')
# 原作者的巨噬细胞标签；不将Monolike-FCN1自动归并进其他巨噬细胞。
macro <- c('hM08_Macro-NLRP3','hM09_Macro-PLTP','hM10_Macro-IL1B','hM12_TAM-C1QC','hM13_TAM-SPP1')
cmp <- d[Tissue=='T'&Sub_Cluster%in%macro,.(cells=.N,mean_expression=mean(EPB41L3),detection=mean(EPB41L3>0)),by=.(Sample,group=ifelse(Sub_Cluster=='hM13_TAM-SPP1','SPP1 TAM','Other macrophages'))]
good <- cmp[cells>=20,.(groups=uniqueN(group)),by=Sample][groups==2]$Sample
cmp <- cmp[Sample%in%good];write(cmp,'matched_patient_macrophage_comparison')
w <- dcast(cmp,Sample~group,value.var=c('mean_expression','detection','cells'))
if(nrow(w)) {
 w[,expression_difference:=get('mean_expression_SPP1 TAM')-get('mean_expression_Other macrophages')]
 w[,detection_difference:=get('detection_SPP1 TAM')-get('detection_Other macrophages')]
}
write(w,'matched_patient_differences')
theme_set(theme_classic(base_size=11)+theme(plot.title=element_text(face='bold',size=12),axis.text=element_text(colour='#273B45'),legend.title=element_text(size=10)))
label <- function(x)sub('^hM[0-9]+_','',x)
shown <- c('EPB41L3','SPP1','LST1','CD68','C1QC','APOE','TREM2','NLRP3','IL1B','FCN1')
a <- summary[min_cells==20&gene%in%shown]
a[,gene:=factor(gene,levels=shown)]
a[,subtype:=factor(label(Sub_Cluster),levels=rev(label(sort(unique(a$Sub_Cluster)))))]
p1 <- ggplot(a,aes(gene,subtype,size=detection*100,colour=mean_expression))+geom_point()+scale_size_area(max_size=6,limits=c(0,100))+scale_colour_gradient(low='#DBE9EC',high='#17697B')+labs(title='A  Myeloid expression profiles',subtitle='Tumor; patient-equal means; >=20 cells per subtype',x=NULL,y=NULL,size='Detected (%)',colour='Published\nexpression')+theme(axis.text.x=element_text(angle=55,hjust=1))
b <- patient[Tissue=='T'&gene=='EPB41L3']
b[,subtype:=factor(label(Sub_Cluster),levels=rev(label(sort(unique(b$Sub_Cluster)))))]
b[,coverage:=factor(ifelse(cells>=20,'>=20 cells','<20 cells'),levels=c('>=20 cells','<20 cells'))]
p2 <- ggplot(b,aes(Sample,subtype))+geom_point(aes(size=100*detection,colour=coverage),alpha=.85)+scale_colour_manual(values=c('>=20 cells'='#C97646','<20 cells'='#BDC6CA'))+scale_size_area(max_size=6,limits=c(0,100))+labs(title='B  EPB41L3 detection by patient',subtitle='Each point represents one patient and subtype',x=NULL,y=NULL,colour='Coverage',size='Detected (%)')+theme(axis.text.x=element_text(angle=55,hjust=1))
p3 <- ggplot(cmp,aes(group,100*detection,group=Sample))+geom_line(colour='#B0BDC3')+geom_point(aes(colour=group),size=3)+geom_text(aes(label=ifelse(group=='SPP1 TAM',Sample,'')),nudge_x=.13,size=3,hjust=0)+scale_colour_manual(values=c('Other macrophages'='#277C8D','SPP1 TAM'='#C97646'),guide='none')+scale_x_discrete(expand=expansion(add=c(.4,1)))+labs(title='C  Within-patient macrophage comparison',subtitle=paste(length(good),'patients with >=20 cells in both groups'),x=NULL,y='EPB41L3 detected (%)')
e <- epb[Tissue=='T']
p4 <- ggplot(e,aes(Sample,cells,fill=cells>=20))+geom_col(width=.7)+geom_hline(yintercept=20,linetype=2,colour='#677C83')+geom_text(aes(label=cells),vjust=-.4,size=3.4)+scale_fill_manual(values=c('FALSE'='#C8D1D5','TRUE'='#C97646'),guide='none')+scale_y_continuous(expand=expansion(mult=c(0,.12)))+labs(title='D  SPP1 TAM coverage',subtitle='Original author labels; tumor samples',x=NULL,y='Cells per patient')+theme(axis.text.x=element_text(angle=55,hjust=1))
p <- (p1|p2)/(p3|p4)+plot_annotation(title='Independent localization of EPB41L3 in GSE146771',subtitle='10x cohort | Published annotations | Descriptive patient-level validation',theme=theme(plot.title=element_text(size=17,face='bold'),plot.subtitle=element_text(size=11)))
ggsave(file.path(out,'Figure_GSE146771_localization.pdf'),p,width=15,height=11,device=cairo_pdf)
ggsave(file.path(out,'Figure_GSE146771_localization.png'),p,width=15,height=11,dpi=220,bg='white')
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
print(epb);print(summary[min_cells==20&gene=='EPB41L3']);print(w)

