# 通讯结果图：展示完整网络概况及按关联强度选出的代表链，不强制纳入候选基因。
library(dplyr)
library(ggplot2)
library(patchwork)
out<-"D:/Bioinfor/colorectal/single_raw/three_gene_analysis"
savefig<-function(p,name,w,h) {
  ggsave(file.path(out,"figures",paste0(name,".pdf")),p,width=w,height=h,device=cairo_pdf,bg="white")
  ggsave(file.path(out,"figures",paste0(name,".png")),p,width=w,height=h,dpi=300,bg="white")
}
s<-read.csv(file.path(out,"scmlnet","summary.csv"))
s$direction<-paste(s$sender,s$receiver,sep=" -> ")
p1<-ggplot(s,aes(reorder(direction,paths),paths,fill=receiver))+geom_col(width=.65)+coord_flip()+
  scale_fill_manual(values=c(Stromal="#8270A4",Myeloid="#508D79",Epithelial="#5F8FA8"))+
  theme_classic(base_size=11)+labs(x=NULL,y="Complete predicted L-R-TF-target paths",fill="Receiver",title="Multilayer network overview")
savefig(p1,"scMLnet_direction_summary",9,4.5)
paths<-read.csv(file.path(out,"scmlnet","all_paths.csv"))
if(nrow(paths)) {
  tag<-"Stromal_to_Myeloid"
  z<-paths %>% filter(sender=="Stromal",receiver=="Myeloid")
  if(nrow(z)) {
    rt<-read.csv(file.path(out,"scmlnet",paste0(tag,"_RecTF_correlation.csv")))
    tt<-read.csv(file.path(out,"scmlnet",paste0(tag,"_TFTarget_correlation.csv")))
    z$score<-pmin(rt$tau[match(paste(z$Receptor,z$TF,sep="_"),rt$key)],tt$tau[match(paste(z$TF,z$Target,sep="_"),tt$key)])
    z<-z %>% arrange(desc(score),Ligand,Receptor,TF,Target) %>% slice_head(n=12)
    write.csv(z,file.path(out,"scmlnet","displayed_representative_paths.csv"),row.names=FALSE)
    layers<-c("Ligand","Receptor","TF","Target")
    nodes<-bind_rows(lapply(seq_along(layers),function(i) {
      v<-sort(unique(z[[layers[i]]]))
      data.frame(gene=v,layer=layers[i],x=i,y=if(length(v)==1).5 else seq(0,1,length.out=length(v)),key=paste(i,v))
    }))
    edges<-bind_rows(lapply(1:3,function(i) {
      pairs<-unique(z[,layers[c(i,i+1)]])
      a<-match(paste(i,pairs[[1]]),nodes$key);b<-match(paste(i+1,pairs[[2]]),nodes$key)
      data.frame(x=nodes$x[a],y=nodes$y[a],xend=nodes$x[b],yend=nodes$y[b])
    }))
    p2<-ggplot()+geom_segment(data=edges,aes(x=x,y=y,xend=xend,yend=yend),colour="#BEC5CD",linewidth=.6,alpha=.7)+
      geom_point(data=nodes,aes(x,y,colour=layer),size=3)+
      geom_label(data=nodes,aes(x,y,label=gene),size=3,linewidth=0,nudge_y=.035,fill="white")+
      scale_colour_manual(values=c(Ligand="#8270A4",Receptor="#508D79",TF="#C58760",Target="#5F8FA8"))+
      scale_x_continuous(breaks=1:4,labels=c("Stromal\nligand","Myeloid\nreceptor","Myeloid\nTF","Myeloid\ntarget"),expand=expansion(mult=.12))+
      scale_y_continuous(expand=expansion(mult=.12))+theme_void(base_size=11)+
      theme(axis.text.x=element_text(),legend.position="none")+
      labs(title="Representative stromal-to-myeloid paths",subtitle="12 paths ranked by the weaker of the two Kendall correlations")
    savefig(p2,"scMLnet_representative_network",11,7)
  }
}
# 各队列 CellChat 通路总强度：仅描述，不作为患者层面的组间显著性检验。
pathways<-bind_rows(lapply(c("Normal","Tumor"),function(k) {
  z<-read.csv(file.path(out,"cellchat",paste0(k,"_pathway_strength.csv")));z$condition<-k;z
}))
top<-pathways %>% group_by(pathway) %>% summarise(total=sum(strength),.groups="drop") %>% arrange(desc(total)) %>% slice_head(n=15)
p3<-ggplot(filter(pathways,pathway %in% top$pathway),aes(strength,factor(pathway,levels=rev(top$pathway)),colour=condition))+
  geom_point(size=3)+scale_colour_manual(values=c(Normal="#389B98",Tumor="#CE6657"))+
  theme_classic(base_size=11)+labs(x="Summed predicted communication strength",y=NULL,colour=NULL,title="Major signaling pathways")
savefig(p3,"CellChat_top_pathways",9,6)
lr<-read.csv(file.path(out,"cellchat","Tumor_LR_pairs.csv"))
is_myeloid<-function(x)grepl("^Myeloid",x)
is_stromal<-function(x)grepl("Stromal|fibroblast|muscle|Pericyte|ECs|glial",x)
focus<-lr[(is_myeloid(lr$source)&is_stromal(lr$target))|(is_stromal(lr$source)&is_myeloid(lr$target)),]
focus<-focus %>% arrange(desc(prob),source,target,interaction_name) %>% slice_head(n=25)
write.csv(focus,file.path(out,"cellchat","displayed_stromal_myeloid_edges.csv"),row.names=FALSE)
if(nrow(focus)) {
  focus$label<-paste(focus$source,"->",focus$target,"|",focus$interaction_name)
  p4<-ggplot(focus,aes(prob,reorder(label,prob),colour=pathway_name))+geom_point(size=3)+
    theme_classic(base_size=10)+labs(x="Predicted communication probability",y=NULL,colour="Pathway",title="Stromal-myeloid communication in tumors")
  savefig(p4,"CellChat_stromal_myeloid",13,9)
}
