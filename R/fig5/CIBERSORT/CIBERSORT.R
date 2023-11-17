######################CIBERSORT.R######################==
rm(list=ls())   #清空环境变量
options(stringsAsFactors = F)

#1. 准备工作
##1.1 相关依赖性R包安装
#install.packages('e1071')
#install.packages('parallel')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("preprocessCore", version = "3.8")

#1.2 函数的下载及准备
#在线分析：https://cibersort.stanford.edu/
#函数下载：https://content.cruk.cam.ac.uk/fmlab/sivakumar2016/Cibersort.R
source("source.R")   #注释文件

#2. CIBERSORT计算
sig_matrix <- "LM22.txt"   #注释文件名
mixture_file = 'cibersort_input.txt'   #表达数据文件名，需要修改
res_cibersort <- CIBERSORT(sig_matrix, mixture_file, perm=100, QN=TRUE)
save(res_cibersort,file = "res_cibersort.Rdata")   #保存中间文件

#3. CIBERSORT分析结果的可视化展示
# 读取CIBERSORT结果
rm(list=ls())
load("res_cibersort.Rdata")
res_cibersort <- res_cibersort[,1:22]   #取前22列为细胞丰度数据
ciber.res <- res_cibersort[,colSums(res_cibersort) > 0]   #去除丰度全为0的细胞
##3.1 barplot图
mycol <- ggplot2::alpha(rainbow(ncol(ciber.res)), 0.7) #创建彩虹色板（带70%透明度）
par(bty="o", mgp = c(2.5,0.3,0), mar = c(2.1,4.1,2.1,10.1),tcl=-.25,las = 1,xpd = F)
barplot(as.matrix(t(ciber.res)),
        border = NA, # 柱子无边框写
        names.arg = rep("",nrow(ciber.res)), # 无横坐标样本名
        yaxt = "n", # 先不绘制y轴
        ylab = "Relative percentage", # 修改y轴名称
        col = mycol) # 采用彩虹色板
axis(side = 2, at = c(0,0.2,0.4,0.6,0.8,1), # 补齐y轴添加百分号
     labels = c("0%","20%","40%","60%","80%","100%"))
legend(par("usr")[2]-20, # 这里-20要根据实际出图的图例位置情况调整
       par("usr")[4], 
       legend = colnames(ciber.res), 
       xpd = T,
       fill = mycol,
       cex = 0.8, 
       border = NA, 
       y.intersp = 1,
       x.intersp = 0.2,
       bty = "n")
dev.off()   #关闭画板

##3.2 相关性热图
M <- round(cor(ciber.res),2) # 计算相关性矩阵并保留两位小数
#绘制方法一：pheatmap
pheatmap::pheatmap(M)
#绘制方法一：corrplot
library(corrplot)
corrplot.mixed(M,
               lower.col = "black", # 左下方字体颜色为黑色
               upper = "color", # 右上方以颜色方块法绘制
               tl.pos = "lt", # 标签出现在左侧和顶部
               number.cex = 0.8) # 左下方字号为0.8
dev.off()   #关闭画板

##3.3 PCA图
library(FactoMineR)#画主成分分析图需要加载这两个包
library(factoextra)
#输入数据：exp和group_list
dat=as.data.frame(ciber.res)
#加载表达谱
expr <- read.table("tpms_mrna.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
target.gene <- "BTK"   #选取目标基因
target.gene.expr <- as.numeric(log2(expr[target.gene,] + 1))
names(target.gene.expr) <- colnames(expr)
identical(names(target.gene.expr), rownames(ciber.res))
#[1] TRUE
group_list=ifelse(target.gene.expr < median(target.gene.expr), "Low", "High")
group_list = factor(group_list,
                    levels = c("Low","High"))
table(group_list)

dat.pca <- PCA(dat, graph = FALSE)
pca_plot <- fviz_pca_ind(dat.pca,
                         geom.ind = "point", #只展现点
                         col.ind = group_list, #分组的颜色
                         palette = c("#00AFBB", "#E7B800"),
                         addEllipses = TRUE, # Concentration ellipses
                         legend.title = "Groups"
)
pca_plot
dev.off()   #关闭画板

##3.4 免疫细胞相关性分析
# install.packages("ggpubr")
library(ggpubr)
##循环生成相关性散点图
p.cutoff <- 0.05 # 相关性检验p值的阈值
outTab <- NULL
for (i in colnames(ciber.res)) {
  message(paste0("--analysis of ",i," done..."))
  # 构建数据库，包括基因表达以及对应细胞的丰度
  dat <- data.frame(gene.expr = as.numeric(target.gene.expr),
                    cell.score = as.numeric(ciber.res[names(target.gene.expr),i]),
                    stringsAsFactors = F)
  
  cor.res <- cor.test(dat$gene.expr,dat$cell.score, method = "pearson") # 相关性分析
  outTab <- rbind.data.frame(outTab,
                             data.frame(gene = target.gene,
                                        cell = i,
                                        rho = cor.res$estimate,
                                        pval = cor.res$p.value,
                                        stringsAsFactors = F),
                             stringsAsFactors = F)
  
  if(cor.res$p.value < p.cutoff) { # 若检验p值小于阈值则出图
    sp <- ggscatter(dat, 
                    x = "gene.expr", # x为基因表达
                    y = "cell.score", # y为细胞丰度
                    xlab = paste0(target.gene," expression level"),
                    ylab = i,
                    size = 0.8,
                    color = "blue", # 散点颜色
                    add = "reg.line",  # 添加回归线
                    add.params = list(color = "red", fill = "grey40"), # 设置回归线以及回归置信区间的颜色
                    conf.int = TRUE) + 
      stat_cor(method = "pearson", # 相关性分析方法
               label.x = min(dat$gene.expr), # 结果标注的x位置
               label.y = max(dat$cell.score)) # 结果标注的y位置
    ggsave(file = paste0("correlation scatter plot between expression of ", target.gene," and ", i, ".pdf"), width = 4, height = 4)
  }
}
write.table(outTab,paste0("correlation results between ",target.gene," and cibersort cells.txt"),sep = "\t",row.names = F,col.names = T,quote = F)

##3.5 小提琴图
rm(list=ls())
###3.5.1 加载R包
#install.packages("vioplot")
library(vioplot)

###3.5.2 加载表达谱
expr <- read.table("tpms_mrna.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
target.gene <- "BTK"   #选取目标基因

###3.5.3 取出目标基因表达
target.gene.expr <- as.numeric(log2(expr[target.gene,] + 1))
names(target.gene.expr) <- colnames(expr)

#取出高低表达组
hsam <- names(target.gene.expr[target.gene.expr > median(target.gene.expr)]) # 高于中位数的样本
lsam <- names(target.gene.expr[target.gene.expr <= median(target.gene.expr)]) # 低于中位数的样本

n.hsam <- length(hsam) # 高表达样本数目
n.lsam <- length(lsam) # 低表达样本数目

###3.5.4 加载CIBERSOciber免疫丰度数据
load("res_cibersort.Rdata")
res_cibersort <- res_cibersort[,1:22]   #取前22列为细胞丰度数据
ciber.res <- res_cibersort[,colSums(res_cibersort) > 0]   #去除丰度全为0的细胞
ciber.res <- ciber.res[c(lsam,hsam),] # 根据高低样本排序表达谱

###3.5.5 绘制小提琴图
par(bty = "o", mgp = c(2,.5,0), mar = c(9,3,2,2), las = 1, xpd = F)
x <- c(1:ncol(ciber.res))
y <- c(1:ncol(ciber.res))
plot(x,y,
     xlim=c(0,ncol(ciber.res)*3-3),ylim=c(min(ciber.res),max(ciber.res)+0.02),
     main = "",
     xlab = "", ylab = "Fraction",
     pch = 21,
     col = "white",
     xaxt = "n")
legend("topright", legend = c("Low", "High"), fill = c("blue","red"), bty = "n", x.intersp = 0.2, y.intersp = 0.8) # 绘制图例

# 循环每一个丰度不为0的细胞绘制小提琴图
for(i in 1:ncol(ciber.res)){
  lsam.expr <- ciber.res[1:n.lsam,i]
  hsam.expr <- ciber.res[(n.lsam+1):(n.lsam+n.hsam),i]
  
  vioplot(lsam.expr,at = 3*(i-1),lty = 1,add = T,col = "blue")
  vioplot(hsam.expr,at = 3*(i-1)+1,lty = 1,add = T,col = "red")
  
  wt <- wilcox.test(lsam.expr,hsam.expr) # 非参数检验
  p <- round(wt$p.value,3) 
  mx <- max(c(lsam.expr,hsam.expr))
  
  lines(c(x=3*(i-1) + 0.2,x = 3*(i-1) + 0.8),c(mx,mx))
  text(x = 3*(i-1) + 0.5,y = mx + 0.02,labels = ifelse(p<0.001,paste0("P<0.001"),paste0("P=",p)),cex = 0.8)
  text(seq(1,ncol(ciber.res)*3-2,3), -0.03, labels=colnames(ciber.res),cex = 1, srt = 45, pos=2, xpd = T)
}
dev.off()
