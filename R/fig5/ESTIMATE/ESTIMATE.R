######################ESTIMATE包评估肿瘤纯度.R######################==
rm(list=ls())   #清空环境变量
options(stringsAsFactors = F)

#1. 准备工作
##1.1 R包的安装与读取
#library(utils)
#rforge <- "http://r-forge.r-project.org"
#install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)
##1.2 读取表达文件
inputFile="tpms_mrna.txt"
expr <- read.table(inputFile,sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#对数转化并保存到当前目录下
expr <- log2(expr + 1)   #对数转换
write.table(expr, "tpms_mrna_log2transformed.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

#2. 估计各类免疫得分
filterCommonGenes(input.f = "tpms_mrna_log2transformed.txt",   #输入文件名
                  output.f = "tpms_mrna_log2transformed.gct",   #输出文件名
                  id = "GeneSymbol")   #行名为gene symbol
estimateScore("tpms_mrna_log2transformed.gct",   #刚才的输出文件名
              "tpm_mrna_estimate_score.txt",   #新的输出文件名
              platform="affymetrix")   #默认平台

#3. 输出每个样品的打分
est <- read.table("tpm_mrna_estimate_score.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
est <- est[,-1]   #移除第一列
colnames(est) <- est[1,]   #设置列名
est <- as.data.frame(t(est[-1,]))
rownames(est) <- colnames(expr)
write.table(est, file = "tpm_mrna_estimate_score.txt",sep = "\t",row.names = T,col.names = NA,quote = F) # 保存并覆盖得分
