rm(list = ls())
setwd("D:/Bioinfor/colorectal/PPI")
#1.加载差异分析结果，并且根据标识提取得到需要的差异基因列表
eqtl <- read.csv("D:/Bioinfor/colorectal/MR/eQTL/Mr_colorectal_result/eqtl_deg_intersection.csv", row.names = NULL)
decode<-read.csv("D:/Bioinfor/colorectal/MR/pQTL/MR_IVW_FDR05_DEG_intersection.csv",row.names = NULL)
names(eqtl)
names(decode)
eg_all<-read.csv("D:/Bioinfor/colorectal/results/GSE196006_GSE146009_GSE110224/DEG_Tumor_vs_Normal_significant.csv", row.names = NULL)
names(eg_all)
output_dir <- "D:/Bioinfor/colorectal/PPI"
# ===== 合并：decode(8基因) ，eqtl (12基因)

# ===== 提取两个基因列表 =====
g_decode <- unique(decode$GeneSymbol)
g_eqtl   <- unique(eqtl$GeneSymbol)

# ===== 合并（并集，去重） =====
all_genes <- union(g_decode, g_eqtl) # 全部基因

library(dplyr)

# 以 all_genes 为基准，左连接 eg_all 的 DEG 信息
gene_info <- left_join(data.frame(GeneSymbol = all_genes), 
                       eg_all[,c("GeneSymbol","Category")], 
                       by = "GeneSymbol")

# 检查
head(gene_info)
table(gene_info$Category, useNA = "ifany")   # UP / DOWN 分布（NA 表示不在 eg_all）
nrow(gene_info)                              # 应 = length(all_genes)

# 保存
write.csv(gene_info, file.path(output_dir,"gene_info.csv" ))

gene_up= gene_info[gene_info$Category=='Up','GeneSymbol']
gene_down=gene_info[gene_info$Category=='Down','GeneSymbol']
gene_diff = c(gene_up,gene_down)

#2.制作string的输入数据
write.table(gene_diff,
            file="diffgene.txt",
            row.names = F,
            col.names = F,
            quote = F)
#############################-------sting网页Version: 12.5版；0.15；20————#############################
library(igraph)
library(dplyr)

# ===== 1. 读取边表（不要加 comment.char = "#"） =====
edges <- read.delim("D:/Bioinfor/colorectal/PPI/string_interactions_short_39_nodes.tsv",
                    stringsAsFactors = FALSE)
names(edges)[1:2] <- c("node1", "node2")      # 第一列原名是 #node1

# ===== 2. 构建网络 =====
net <- graph_from_data_frame(edges[, c("node1", "node2")], directed = FALSE)
cat("节点数:", vcount(net), " 边数:", ecount(net), "\n")

# ===== 3. 中心性计算 =====
deg <- degree(net)                             # 先定义 deg
hub <- data.frame(
  GeneSymbol  = V(net)$name,
  Degree      = as.numeric(deg),
  Betweenness = betweenness(net, weights = NA),
  Closeness   = closeness(net, weights = NA)
) %>% arrange(desc(Degree))

# ===== 4. 标记 CEGs =====
cegs <- c("EPHB1","ITIH5","GPER1","LAMC2","MYO15B","MS4A8","CSE1L","HHIP",
          "ANGPTL1","KLB","FABP2","WDR72","EPB41L3","NTRK3","FPR1","GREM1")
hub$Is_CEG <- ifelse(hub$GeneSymbol %in% cegs, "CEG", "Neighbor")

head(hub, 36)
write.csv(hub, file.path(output_dir,'hub_genes_36nodes.csv'))

# ===== 5. 画优化后的网络图 =====
V(net)$size  <- sqrt(deg) * 5 + 4              # 节点大小（平方根压缩，不会过大）
V(net)$color <- ifelse(V(net)$name %in% cegs, "#EA6668", "#C0C0C0")
V(net)$frame.color <- "white"
V(net)$label <- ifelse(deg >= 5, V(net)$name, NA)   # 只标 degree≥5 的 hub
V(net)$label.cex <- 0.75
V(net)$label.color <- "#1A1B1C"

l <- layout_with_kk(net)                        # Kamada-Kawai 均匀布局

png("D:/Bioinfor/colorectal/PPI/ppi_network_optimized.png",
    width = 2000, height = 1600, res = 200)
par(mar = c(1, 1, 2, 1))
plot(net, layout = l,
     edge.width = 1.2, edge.color = "#999999", edge.curved = 0.1,
     vertex.label.dist = 0.6, vertex.label.degree = -pi/2,
     main = "CEGs PPI network  (red = CEG, gray = neighbor)")
dev.off()
#####

