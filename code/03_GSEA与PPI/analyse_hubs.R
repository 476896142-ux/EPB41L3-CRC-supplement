library(igraph)
out <- 'E:/codex/2026-10-31/dui/outputs/Hub_Degree_MCC_MNC'
# 原 CSV 每条数据记录末尾有一个空字段；删除末尾分隔符后显式读取，防止首列错作行名。
lines <- readLines(file.path(out,'source_data/cytoHubba_original.csv'))
d <- read.csv(text=paste(sub(',[[:space:]]*$','',lines),collapse='\n'),check.names=FALSE)
n <- read.csv(file.path(out,'source_data/cytoscape_nodes_double_ring.csv'))
e <- read.csv(file.path(out,'source_data/cytoscape_edges_double_ring.csv'))
stopifnot(!anyDuplicated(d$node_name),setequal(d$node_name,n$GeneSymbol))
g <- graph_from_data_frame(e[,c('source','target')],directed=FALSE,vertices=n$GeneSymbol)
stopifnot(is_simple(g),vcount(g)==36,ecount(g)==100)
cl <- max_cliques(g)
calc <- data.frame(node_name=V(g)$name,Degree=as.numeric(degree(g)),MCC=0,MNC=0)
for(i in seq_len(nrow(calc))) {
 name <- calc$node_name[i]
 calc$MCC[i] <- sum(vapply(cl,function(cq) if(name %in% as_ids(cq)) factorial(length(cq)-1) else 0,numeric(1)))
 nb <- neighbors(g,name)
 calc$MNC[i] <- if(length(nb)) max(components(induced_subgraph(g,nb))$csize) else 0
}
audit <- merge(d[,c('node_name','Degree','MCC','MNC')],calc,by='node_name',suffixes=c('_CSV','_recomputed'))
for(m in c('Degree','MCC','MNC')) audit[[paste0(m,'matches')]] <- audit[[paste0(m,'_CSV')]]==audit[[paste0(m,'_recomputed')]]
write.csv(audit,file.path(out,'score_verification.csv'),row.names=FALSE)
print(colSums(audit[,grepl('matches',names(audit))]))
stopifnot(all(as.matrix(audit[,grepl('matches',names(audit))])))
write.csv(d,file.path(out,'cytoHubba_cleaned.csv'),row.names=FALSE)
eligible <- d$node_name %in% n$GeneSymbol[n$Is_CEG=='CEG']
for(m in c('Degree','MCC','MNC')) {
 d[[paste0(m,'_rank')]] <- NA_real_
 d[[paste0(m,'_rank')]][eligible] <- rank(-d[[m]][eligible],ties.method='min')
 cutoff <- sort(d[[m]][eligible],decreasing=TRUE)[10]
 d[[paste0(m,'_selected')]] <- eligible & d[[m]]>=cutoff
 cat(m,'cutoff',cutoff,'n',sum(d[[paste0(m,'_selected')]]),'\n')
}
d$Is_CEG <- n$Is_CEG[match(d$node_name,n$GeneSymbol)]
d$Intersection <- d$Degree_selected & d$MCC_selected & d$MNC_selected
write.csv(d,file.path(out,'rankings_and_membership.csv'),row.names=FALSE)
print(d[d$Intersection,c('node_name','Degree','MCC','MNC','Is_CEG')])

# 仅 CEG 参与排名，分数来自完整 36 节点网络；各方法保留第 10 位的全部并列。
library(ggplot2)
library(cowplot)
methods <- c('Degree','MCC','MNC')
hubs <- d[d$Intersection,]
hubs <- hubs[order(-hubs$MCC,hubs$node_name),]
write.csv(hubs,file.path(out,'consensus_hubs.csv'),row.names=FALSE)
for(m in methods) {
 selected <- d[d[[paste0(m,'_selected')]],]
 selected <- selected[order(-selected[[m]],selected$node_name),]
 write.csv(selected,file.path(out,paste0(m,'_CEG_Top10_with_ties.csv')),row.names=FALSE)
}
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
 geom_point(data=inn[!inn$selected,],aes(X,Y),shape=21,size=10.5,fill='#EBEEF0',color='#AAB2B8',stroke=.6)+
 geom_point(data=inn[inn$selected,],aes(X,Y,fill=score),shape=21,size=10.5,color='#9C6047',stroke=.65)+
 geom_text(data=inn,aes(X,Y,label=GeneSymbol),color='#25333C',fontface='bold',size=2.35)+
 geom_text(data=outer,aes(lx,ly,label=GeneSymbol,hjust=ifelse(lx>0,0,1)),size=2.1,color='#63717A')+
 scale_fill_gradient(low='#FCE5A0',high='#DE684B',name=paste(m,'score'))+
 annotate('text',x=0,y=.10,label=paste(sum(inn$selected),'CEGs selected'),size=3.2,fontface='bold',color='#42515B')+
 annotate('text',x=0,y=-.18,label=paste('Cut-off:',cutoff),size=2.8,color='#63717A')+
 coord_equal(xlim=c(-3.65,4.05),ylim=c(-3.4,3.45),clip='off')+
 labs(title=m,subtitle='CEG Top 10, including boundary ties')+
 theme_void(base_family='sans')+
 theme(plot.title=element_text(face='bold',size=15),plot.subtitle=element_text(size=9,color='#63717A'),
 legend.position='bottom',legend.title=element_text(size=9),legend.text=element_text(size=8),plot.margin=margin(12,12,12,12))
 plots[[m]] <- p
 ggsave(file.path(out,paste0(m,'_PPI.png')),p,width=7.4,height=7.2,dpi=300,device=ragg::agg_png,bg='white')
 ggsave(file.path(out,paste0(m,'_PPI.pdf')),p,width=7.4,height=7.2,device=cairo_pdf,bg='white')
}
combined <- plot_grid(plotlist=plots,nrow=1,labels=c('A','B','C'),label_size=17)
footer <- ggdraw()+draw_label(paste0('Consensus CEGs (n = ',nrow(hubs),'): ',paste(hubs$node_name,collapse=', ')),size=12,fontface='bold')
combined <- plot_grid(combined,footer,ncol=1,rel_heights=c(1,.08))
ggsave(file.path(out,'Degree_MCC_MNC_PPI.png'),combined,width=22.2,height=7.8,dpi=300,device=ragg::agg_png,bg='white')
ggsave(file.path(out,'Degree_MCC_MNC_PPI.pdf'),combined,width=22.2,height=7.8,device=cairo_pdf,bg='white')
ggsave(file.path(out,'Degree_MCC_MNC_PPI.svg'),combined,width=22.2,height=7.8,device=svglite::svglite,bg='white')
# 交集矩阵同时保留所有 CEG，展示并列导致的入选差异。
ceg <- d[d$Is_CEG=='CEG',]
ceg <- ceg[order(-ceg$MCC,ceg$node_name),]
mat <- do.call(rbind,lapply(methods,function(m)data.frame(Gene=ceg$node_name,Method=m,
 Selected=ceg[[paste0(m,'_selected')]],Rank=ceg[[paste0(m,'_rank')]],Score=ceg[[m]])))
mat$Gene <- factor(mat$Gene,levels=rev(ceg$node_name));mat$Method <- factor(mat$Method,levels=methods)
mat$Text <- paste0('Rank ',mat$Rank,'\n',mat$Score)
p <- ggplot(mat,aes(Method,Gene))+geom_tile(aes(fill=Selected),color='white',linewidth=1)+
 geom_text(aes(label=Text),size=3,lineheight=.95)+
 scale_fill_manual(values=c('FALSE'='#EFF2F4','TRUE'='#B6D5CC'),name='Selected')+
 labs(x=NULL,y=NULL,title='CEG ranking and three-method intersection',subtitle='Rank within 16 CEGs; scores computed on the full network')+
 theme_minimal(base_size=11)+theme(panel.grid=element_blank(),plot.title=element_text(face='bold'))
ggsave(file.path(out,'CEG_selection_matrix.png'),p,width=7,height=9,dpi=300,device=ragg::agg_png,bg='white')
ggsave(file.path(out,'CEG_selection_matrix.pdf'),p,width=7,height=9,device=cairo_pdf,bg='white')
