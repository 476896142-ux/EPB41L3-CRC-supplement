# 从保存的结果表重绘；不重新计算富集，不改变网络边或布局。
library(ggplot2)
library(cowplot)
library(stringr)
args <- commandArgs(trailingOnly=TRUE)
out <- if(length(args)) args[1] else 'E:/codex/2026-10-31/dui/outputs/Figure_GSEA_Hub_4panel'
d <- read.csv(file.path(out,'source_data/lollipop_selected_terms.csv'), check.names=FALSE)
n <- read.csv(file.path(out,'source_data/cytoscape_nodes_double_ring.csv'))
e <- read.csv(file.path(out,'source_data/cytoscape_edges_double_ring.csv'))
stopifnot(nrow(d)==20, nrow(n)==36, nrow(e)==100, sum(n$Is_CEG=='CEG')==16,
          all(d$p.adjust<0.05), all(abs(d$GeneRatio-d$Count/d$setSize)<1e-10),
          all(e$source %in% n$GeneSymbol),all(e$target %in% n$GeneSymbol))
# 统一配色：类别色用于 A；B 的颜色仅表示节点身份。
category <- c('GO BP'='#4B9562','GO CC'='#368DA9','GO MF'='#9A6BA5','KEGG'='#536B9B')
d$Group <- factor(d$Group, levels=names(category))
d$Label <- str_wrap(d$Description, width=39)
a <- ggplot(d,aes(y=y))+
 geom_segment(aes(x=0,xend=GeneRatio,yend=y,color=Group),linewidth=.65)+
 geom_point(aes(x=GeneRatio,size=Count,fill=neg_log10_FDR),shape=21,color='#647078',stroke=.35)+
 geom_text(aes(x=-.025,label=Label,color=Group),hjust=1,size=3.1,lineheight=.95)+
 scale_color_manual(values=category,name=NULL)+
 scale_fill_gradientn(colors=c('#338FB4','#C0C5BB','#EBA139'),name=expression(-log[10](FDR)))+
 scale_size_area(max_size=7.2,name='Leading-edge genes',breaks=c(50,150,250))+
 scale_x_continuous(limits=c(0,.72),breaks=seq(0,.6,.2),expand=expansion(mult=c(0,.02)),oob=scales::oob_keep)+
 scale_y_continuous(limits=c(min(d$y)-.7,max(d$y)+.7),expand=c(0,0))+
 coord_cartesian(clip='off')+
 labs(x='Leading-edge gene ratio',y=NULL)+
 guides(color=guide_legend(order=1,nrow=1,override.aes=list(linewidth=1.5)),
        fill=guide_colorbar(order=2,barwidth=grid::unit(28,'mm'),barheight=grid::unit(3,'mm')),
        size=guide_legend(order=3,nrow=1))+
 theme_classic(base_size=10,base_family='sans')+
 theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
       axis.text.x=element_text(color='#37434B'),axis.title.x=element_text(margin=margin(t=8)),
       legend.position='bottom',legend.box='vertical',legend.margin=margin(0,0,0,0),
       legend.spacing.y=grid::unit(1,'mm'),legend.key.height=grid::unit(4,'mm'),
       plot.margin=margin(12,16,6,212))

d <- read.csv(file.path(out,'source_data/rankings_and_membership.csv'))
methods <- c('Degree','MCC','MNC')
stopifnot(sum(d$Intersection)==10, sum(d$Degree_selected)==11, sum(d$MCC_selected)==10, sum(d$MNC_selected)==16)
# 使用相同坐标展示完整网络背景，确保三个面板可直接比较。
e$x <- n$X[match(e$source,n$GeneSymbol)]; e$y <- n$Y[match(e$source,n$GeneSymbol)]
e$xend <- n$X[match(e$target,n$GeneSymbol)]; e$yend <- n$Y[match(e$target,n$GeneSymbol)]
outer <- subset(n,Is_CEG!='CEG'); inner <- subset(n,Is_CEG=='CEG')
outer$ang <- atan2(outer$Y,outer$X)
outer$lx <- 3.02*cos(outer$ang); outer$ly <- 3.02*sin(outer$ang)
plots <- list()
for(m in methods) {
 inn <- inner
 inn$score <- d[[m]][match(inn$GeneSymbol,d$node_name)]
 inn$selected <- d[[paste0(m,'_selected')]][match(inn$GeneSymbol,d$node_name)]
 cutoff <- min(inn$score[inn$selected])
 p <- ggplot()+
 geom_segment(data=e,aes(x=x,y=y,xend=xend,yend=yend),color='#9AA8AF',alpha=.33,linewidth=.33)+
 geom_segment(data=outer,aes(x=X,y=Y,xend=lx,yend=ly),color='#ACB7BE',linewidth=.25)+
 geom_point(data=outer,aes(X,Y),shape=21,size=4.2,fill='#DCE5EB',color='#97A7B2',stroke=.4)+
 geom_point(data=inn[!inn$selected,],aes(X,Y),shape=21,size=12.5,fill='#EBEEF0',color='#AAB2B8',stroke=.6)+
 geom_point(data=inn[inn$selected,],aes(X,Y,fill=score),shape=21,size=12.5,color='#9C6047',stroke=.65)+
 geom_text(data=inn,aes(X,Y,label=GeneSymbol),color='#25333C',fontface='bold',size=2.75)+
 geom_text(data=outer,aes(lx,ly,label=GeneSymbol,hjust=ifelse(lx>0,0,1)),size=2.5,color='#63717A')+
 scale_fill_gradient(low='#FCE5A0',high='#DE684B',name=paste(m,'score'))+
 annotate('text',x=0,y=.10,label=paste(sum(inn$selected),'CEGs selected'),size=3.2,fontface='bold',color='#42515B')+
 annotate('text',x=0,y=-.18,label=paste('Cut-off:',cutoff),size=2.8,color='#63717A')+
 coord_equal(xlim=c(-3.65,4.05),ylim=c(-3.4,3.45),clip='off')+
 labs(title=NULL,subtitle=NULL)+
 theme_void(base_family='sans')+
 theme(plot.title=element_text(face='bold',size=15),plot.subtitle=element_text(size=9,color='#63717A'),
 legend.position='bottom',legend.title=element_text(size=9),legend.text=element_text(size=8),plot.margin=margin(12,12,12,12))
 plots[[m]] <- p
}

# 统一四面板标题与字号；网络使用统一暖色渐变，保留各自评分色标。
wrap_panel <- function(p,letter,title) {
 heading <- ggdraw()+draw_label(letter,x=.015,hjust=0,fontface='bold',size=18)+
 draw_label(title,x=.065,hjust=0,fontface='bold',size=15,color='#2E3B43')
 plot_grid(heading,p,ncol=1,rel_heights=c(.065,.935))
}
a <- a+theme(plot.margin=margin(6,18,8,212))
pa <- wrap_panel(a,'A','GO and KEGG enrichment')
pb <- wrap_panel(plots$Degree,'B','Degree')
pc <- wrap_panel(plots$MCC,'C','MCC')
pd <- wrap_panel(plots$MNC,'D','MNC')
body <- plot_grid(pa,pb,pc,pd,ncol=2,rel_heights=c(1,1))
fig <- plot_grid(body,ggdraw()+draw_label('Consensus CEGs (n = 10): NTRK3, EPHB1, KLB, GREM1, HHIP, ANGPTL1, LAMC2, WDR72, GPER1, EPB41L3',size=11,color='#42515B'),ncol=1,rel_heights=c(1,.025))
ggsave(file.path(out,'GSEA_Degree_MCC_MNC.png'),fig,width=18.5,height=16,dpi=300,device=ragg::agg_png,bg='white')
ggsave(file.path(out,'GSEA_Degree_MCC_MNC.pdf'),fig,width=18.5,height=16,device=cairo_pdf,bg='white')
ggsave(file.path(out,'GSEA_Degree_MCC_MNC.svg'),fig,width=18.5,height=16,device=svglite::svglite,bg='white')
cat('Four-panel figure exported; 20 terms and 36 nodes/100 edges preserved.\n')
