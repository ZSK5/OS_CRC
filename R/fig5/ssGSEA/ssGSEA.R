######################ssGSEA.R######################==
rm(list=ls())   #清空环境变量
options(stringsAsFactors = F)

#1. R包及相关文件的准备
##1.1 R包的安装及读取
if(!require("data.table")) install.packages("data.table",update = F,ask = F)
if(!require("GSVA")) BiocManager::install("GSVA",update = F,ask = F)
library(data.table)
library(GSVA)
#1.2 准备细胞marker
cellMarker <- data.table::fread("cellMarker.csv",data.table = F)
colnames(cellMarker)[2] <- "celltype"

type <- split(cellMarker,cellMarker$celltype)

cellMarker <- lapply(type, function(x){
  dd = x$Metagene
  unique(dd)
})
#将list中每个celltype中的基因进行合并
save(cellMarker,file = "cellMarker_ssGSEA.Rdata")
##1.3 表达量矩阵的准备
###行是基因，列是样本
expr <- data.table::fread("tpms_mrna.txt",data.table = F)   #读取表达文件
rownames(expr) <- expr[,1]   #将第一列作为行名
expr <- expr[,-1]   #去除第一列
expr <- as.matrix(expr)   #将expr转换为矩阵格式

#2. 使用ssGSEA量化免疫浸润
gsva_data <- gsva(expr,cellMarker, method = "ssgsea")

#3. 简单作图可视化展示
library(pheatmap)
pheatmap(gsva_data,
         cluster_rows = TRUE,   #行聚类
         cluster_cols = TRUE,   #列聚类，可以看出样本之间的区分度
         annotation_legend=TRUE,   # 显示注释
         show_rownames = T,   # 显示行名
         show_colnames = F,   # 显示行名
         scale = "row",   #以行来标准化
         color =colorRampPalette(c("blue", "white","red"))(100),   #调色
         #filename = "heatmap_F.pdf",   #是否保存
         fontsize = 10)
