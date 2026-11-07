# 参考 Huang Fig.9G-I：互作数量环图、强度环图、配体受体气泡图。
# 直接使用已完成的肿瘤 CellChat；不重新估计概率或调整显著性阈值。
library(CellChat)
library(ggplot2)
library(dplyr)
library(patchwork)
out <- "D:/Bioinfor/colorectal/single_raw/three_gene_analysis"
cc <- readRDS(file.path(out,"cellchat","Tumor_cellchat.rds"))
groups <- rownames(cc@net$count)
short <- function(x) {
  x <- sub("Myeloid ","",x,fixed=TRUE)
  x <- sub("Epithelial cells","Epithelial",x,fixed=TRUE)
  x <- sub("Smooth muscle cells","Smooth muscle",x,fixed=TRUE)
  x <- sub("Enteric glial cells","Enteric glia",x,fixed=TRUE)
  x <- sub("Pro-inflammatory","Inflammatory myeloid",x,fixed=TRUE)
  x <- sub("SPP1+","SPP1+ myeloid",x,fixed=TRUE)
  x
}
angle <- pi/2-seq(0,2*pi,length.out=length(groups)+1)[- (length(groups)+1)]
nodes <- data.frame(group=groups,x=cos(angle),y=sin(angle),label=short(groups),n=as.numeric(table(cc@idents)[groups]))
colors <- setNames(grDevices::hcl.colors(length(groups),"Dark 3"),groups)
circle <- function(mat,title,legend_title) {
  edges <- as.data.frame(as.table(mat)); names(edges)<-c("source","target","value")
  edges <- edges[edges$value>0,,drop=FALSE]
  i<-match(edges$source,nodes$group); j<-match(edges$target,nodes$group)
  edges$x<-nodes$x[i];edges$y<-nodes$y[i];edges$xend<-nodes$x[j];edges$yend<-nodes$y[j]
  loops<-edges[edges$source==edges$target,,drop=FALSE]
  loopdata<-bind_rows(lapply(seq_len(nrow(loops)),function(k) {
    a<-seq(0,2*pi,length.out=50)
    data.frame(source=loops$source[k],value=loops$value[k],x=1.07*loops$x[k]+.07*cos(a),y=1.07*loops$y[k]+.07*sin(a))
  }))
  p<-ggplot()+geom_curve(data=edges[edges$source!=edges$target,],aes(x,y,xend=xend,yend=yend,colour=source,linewidth=value),
    curvature=.14,alpha=.32,arrow=grid::arrow(length=grid::unit(1.1,"mm"),type="closed"))+
    geom_path(data=loopdata,aes(x,y,group=source,colour=source,linewidth=value),alpha=.45)+
    geom_point(data=nodes,aes(x,y,size=n,fill=group),shape=21,colour="white",stroke=.5)+
    geom_text(data=nodes,aes(x=1.16*x,y=1.16*y,label=label),size=2.7)+
    scale_colour_manual(values=colors,guide="none")+scale_fill_manual(values=colors,guide="none")+
    scale_linewidth_continuous(range=c(.12,1.5),name=legend_title)+
    scale_size_area(max_size=9,name="Sampled cells")+
    coord_equal(xlim=c(-1.5,1.5),ylim=c(-1.3,1.3),clip="off")+
    guides(linewidth=guide_legend(order=1,nrow=1),size=guide_legend(order=2,nrow=1,override.aes=list(fill="#81939F",colour="white")))+
    theme_void(base_size=10)+theme(legend.position="bottom",legend.box="vertical",legend.text=element_text(size=8),legend.title=element_text(size=9),plot.title=element_text(face="bold",hjust=.5),plot.margin=margin(8,15,8,15))+
    labs(title=title)
  p
}
a<-circle(cc@net$count,"Number of interactions","Count")
b<-circle(cc@net$weight,"Interaction strength","Strength")
lr<-read.csv(file.path(out,"cellchat","Tumor_LR_pairs.csv"))
senders<-c("Myofibroblasts","Smooth muscle cells","Myeloid SPP1+")
receivers<-c("Myeloid SPP1+","Myeloid Pro-inflammatory","T cells","Epithelial cells")
focus<-lr %>% filter(source %in% senders,target %in% receivers)
top<-focus %>% group_by(interaction_name,interaction_name_2) %>% summarise(total=sum(prob),.groups="drop") %>%
  arrange(desc(total),interaction_name) %>% slice_head(n=20)
focus<-focus %>% filter(interaction_name %in% top$interaction_name)
focus$pair<-paste(short(focus$source),short(focus$target),sep=" -> ")
pair_order<-unlist(lapply(senders,function(s)paste(short(s),short(receivers),sep=" -> ")))
focus$pair<-factor(focus$pair,levels=pair_order)
focus$LR<-factor(focus$interaction_name_2,levels=rev(top$interaction_name_2))
focus$p_group<-factor(ifelse(focus$pval<.01,"P < 0.01","0.01 <= P < 0.05"),levels=c("0.01 <= P < 0.05","P < 0.01"))
c<-ggplot(focus,aes(pair,LR,colour=prob,size=p_group)) + geom_point(alpha=.95)+
  scale_colour_gradientn(colours=c("#567BA1","#8FB9AD","#E7D8A3","#C85D60"),name="Communication\nprobability")+
  scale_size_manual(values=c(2,3.6),drop=FALSE,name="Permutation P")+
  scale_x_discrete(drop=FALSE)+theme_bw(base_size=10)+
  theme(axis.text.x=element_text(angle=48,hjust=1),panel.grid.minor=element_blank(),panel.grid.major=element_line(colour="#EFEFEF"),plot.title=element_text(face="bold"))+
  labs(x=NULL,y=NULL,title="Selected stromal and myeloid ligand-receptor interactions")
fig <- (a|b)/c + plot_layout(heights=c(1,1.2)) + plot_annotation(tag_levels="A")
for(ext in c("pdf","png")) {
  target<-file.path(out,"figures",paste0("Figure_CellChat_Huang_style.",ext))
  if(ext=="pdf")ggsave(target,fig,width=14,height=13,device=cairo_pdf,bg="white") else
    ggsave(target,fig,width=14,height=13,dpi=300,bg="white")
}
write.csv(focus,file.path(out,"cellchat","Huang_style_displayed_LR.csv"),row.names=FALSE)
write.csv(nodes,file.path(out,"cellchat","Huang_style_node_key.csv"),row.names=FALSE)
