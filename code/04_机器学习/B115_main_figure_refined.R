# 直接读取已核验绘图源数据；不重新筛选、不重新训练。
suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
args <- commandArgs(trailingOnly=TRUE)
src <- if(length(args)) args[1] else 'E:/codex/2026-10-31/dui/outputs/GEO_ML_hub10/nested_selection/publication'
out <- file.path(src,'refined_v2'); dir.create(out,recursive=TRUE,showWarnings=FALSE)
read_source <- function(name) read.csv(file.path(src,paste0('Source_',name,'.csv')))
teal <- '#247C74'; blue <- '#426A94'; purple <- '#956A96'; orange <- '#C17C4E'
gene_cols <- c(EPB41L3=teal,KLB=blue,ANGPTL1=purple,'Other candidates'='#BFC6CD')
theme_set(theme_classic(base_size=8,base_family='Arial')+
 theme(axis.line=element_line(linewidth=.3,color='#59636E'),axis.ticks=element_line(linewidth=.25),
 axis.text=element_text(color='#303A43',size=7),axis.title=element_text(size=8),
 plot.title=element_text(size=9,face='bold',color='#202D38',margin=margin(b=4)),
 plot.subtitle=element_text(size=7,color='#64717B',margin=margin(b=6)),
 legend.position='bottom',legend.text=element_text(size=6.7),legend.title=element_text(size=7),
 legend.key.height=unit(3,'mm'),legend.key.width=unit(5,'mm'),legend.spacing.y=unit(0,'mm'),
 legend.margin=margin(0,0,0,0),plot.margin=margin(7,8,7,7),
 plot.tag=element_text(size=12,face='bold',color='#202D38')))

paths <- read_source('A_LASSO_paths'); cv <- read_source('B_LASSO_CV')
# 已锁定的 LASSO 惩罚值，仅用于标记。
selected_lambda <- 0.204335971785694
A <- ggplot(paths,aes(log10(Lambda),Coefficient,group=Gene,color=Highlight))+
 geom_hline(yintercept=0,color='#E8EBEE',linewidth=.3)+geom_line(linewidth=.55)+
 geom_vline(xintercept=log10(selected_lambda),linetype=2,linewidth=.35,color='#67717B')+
 scale_color_manual(values=gene_cols,name=NULL)+guides(color=guide_legend(nrow=2,byrow=TRUE))+
 labs(title='LASSO coefficient paths',subtitle='10 candidate genes',x=expression(log[10](lambda)),y='Coefficient')
B <- ggplot(cv,aes(log10(Lambda),MeanAUC))+
 geom_errorbar(aes(ymin=pmax(0,MeanAUC-SE_AUC),ymax=pmin(1,MeanAUC+SE_AUC)),width=.035,color='#BBC5CF',linewidth=.3)+
 geom_line(color=blue,linewidth=.5)+geom_point(color=blue,size=1)+
 geom_vline(xintercept=log10(selected_lambda),linetype=2,color=orange,linewidth=.45)+
 annotate('text',x=-1.1,y=.73,label='1-SE selection\nEPB41L3',hjust=1,size=2.5,color='#56626C')+
 scale_y_continuous(limits=c(.48,1.01),breaks=seq(.5,1,.1))+
 labs(title='LASSO cross-validation',subtitle='Mean AUC and standard error',x=expression(log[10](lambda)),y='Inner-CV AUC')

svm <- read_source('C_SVM_RFE_CV')
C <- ggplot(svm,aes(K,MeanAUC,color=factor(C),linetype=factor(C)))+
 annotate('rect',xmin=2.5,xmax=5.5,ymin=-Inf,ymax=Inf,fill=teal,alpha=.06)+
 geom_errorbar(aes(ymin=MeanAUC-SE_AUC,ymax=pmin(1,MeanAUC+SE_AUC)),width=.10,alpha=.38,linewidth=.3)+
 geom_line(linewidth=.5)+geom_point(aes(shape=factor(C)),size=1.25)+
 geom_point(data=svm[svm$C==.1 & svm$K==3,],shape=21,size=2.7,stroke=.6,fill='white',color='#202D38')+
 scale_color_manual(values=c(teal,blue,orange),name='C')+scale_linetype_discrete(name='C')+scale_shape_discrete(name='C')+
 scale_x_continuous(breaks=1:10)+
 labs(title='SVM-RFE selection',subtitle='3-5 gene search range; selected: 3',x='Number of genes',y='Inner-CV AUC')

freq <- read_source('D_selection_frequency')
freq$Gene <- factor(freq$Gene,levels=freq$Gene)
freq$Colour <- ifelse(freq$Final_panel=='Selected',as.character(freq$Gene),'Other candidates')
freq$Panel <- ifelse(freq$Final_panel=='Selected','Final panel','Other')
D <- ggplot(freq,aes(Frequency*100,Gene))+
 geom_segment(aes(x=0,xend=Frequency*100,yend=Gene),color='#DCE2E6',linewidth=.6)+
 geom_point(aes(color=Colour,shape=Panel),size=2.4,stroke=.7)+
 geom_text(aes(label=paste0(round(Frequency*100),'%')),nudge_x=4,hjust=0,size=2.35,color='#39454E')+
 scale_color_manual(values=gene_cols,guide='none')+scale_shape_manual(values=c('Final panel'=16,Other=1),name=NULL)+
 scale_x_continuous(limits=c(0,115),breaks=c(0,25,50,75,100),expand=expansion(mult=c(.025,0)))+
 labs(title='Selection stability',subtitle='25 outer-training fits',x='Selection frequency (%)',y=NULL)+
 theme(axis.ticks.y=element_blank(),axis.line.y=element_blank())

roc <- read_source('E_external_ROC')
labels <- c(Adaptive_1SE='EPB41L3  |  AUC 0.9722',Panel3to5='Three genes  |  AUC 0.9988')
E <- ggplot(roc,aes(FPR,TPR,color=Strategy,linetype=Strategy))+
 geom_abline(slope=1,intercept=0,linetype=2,color='#CDD3D8',linewidth=.35)+geom_path(linewidth=.65)+
 scale_color_manual(values=c(Adaptive_1SE=blue,Panel3to5=teal),labels=labels,name=NULL)+
 scale_linetype_manual(values=c(Adaptive_1SE='dashed',Panel3to5='solid'),labels=labels,name=NULL)+
 scale_x_continuous(breaks=c(0,.25,.5,.75,1))+scale_y_continuous(breaks=c(0,.25,.5,.75,1))+
 coord_fixed(xlim=c(-.01,1.01),ylim=c(-.01,1.01),expand=FALSE,clip='off')+
 guides(color=guide_legend(nrow=2,byrow=TRUE),linetype=guide_legend(nrow=2,byrow=TRUE))+
 labs(title='External discrimination',subtitle='GSE20916 | 36 tumour / 24 normal',x='1 - Specificity',y='Sensitivity')

conf <- read_source('F_confusion_matrix')
# 色深按真实类别内百分比映射；每格同时显示原始计数。
F <- ggplot(conf,aes(Predicted,Observed,fill=RowPercent))+
 geom_tile(color='white',linewidth=1.5)+
 geom_text(aes(label=sprintf('%d\n%.0f%%',Freq,RowPercent),color=RowPercent>=60),size=3.8,lineheight=1.25)+
 scale_fill_gradient(low='#F0F4F6',high=teal,limits=c(0,100),guide='none')+
 scale_color_manual(values=c('FALSE'='#303D47','TRUE'='white'),guide='none')+
 scale_x_discrete(limits=c('Normal','Tumor'),labels=c('Normal','Tumour'),expand=c(0,0))+
 scale_y_discrete(limits=c('Tumor','Normal'),labels=c('Tumour','Normal'),expand=c(0,0))+
 coord_fixed()+labs(title='External classification',subtitle='Fixed threshold = 0.5',x='Predicted class',y='Observed class')+
 theme(axis.line=element_blank(),axis.ticks=element_blank())

# 等宽双列；紧凑六联图，保留原 A–F 顺序和正文引用。
figure <- wrap_plots(A,B,C,D,E,F,ncol=2,heights=c(1,1,1.08))+plot_annotation(tag_levels='A')
for(ext in c('pdf','png','tiff')) {
 file <- file.path(out,paste0('Figure_ML_refined.',ext))
 if(ext=='pdf') ggsave(file,figure,width=190,height=215,units='mm',device=cairo_pdf,bg='white')
 else if(ext=='tiff') ggsave(file,figure,width=190,height=215,units='mm',dpi=600,compression='lzw',bg='white')
 else ggsave(file,figure,width=190,height=215,units='mm',dpi=600,bg='white')
}
cat('Refined figure saved to:',out,'\n')

