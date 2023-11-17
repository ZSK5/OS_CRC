######################42_corrplot######################==
rm(list=ls())   #清空环境变量
options(stringsAsFactors = F)

#1. R包和数据准备
##1.1 读取R包
#install.packages("corrplot")
library(corrplot)
##1.2 读取数据
load("res_cibersort.Rdata")
res_cibersort <- res_cibersort[,1:22]   #取前22列为细胞丰度数据
ciber.res <- res_cibersort[,colSums(res_cibersort) > 0]   #去除丰度全为0的细胞

#2. 计算相关系数
cor_ciber <- cor(ciber.res, method = 'pearson')
cor_ciber <- round(cor(ciber.res),2)

#3. 简单绘制
corrplot(cor_ciber)

#4.调整参数
##4.1 可视化方法
corrplot(cor_ciber, method = "square")   #方形
corrplot(cor_ciber, method = "ellipse")   #椭圆形
corrplot(cor_ciber, method = "number")   #数字
corrplot(cor_ciber, method = "pie")   #饼图
corrplot(cor_ciber, method = "shade")   #阴影
corrplot(cor_ciber, method = "color")   #颜色

##4.2 矩阵显示类型
corrplot(cor_ciber, type = "lower")   #下三角矩阵
corrplot(cor_ciber, type = "upper")   #上三角矩阵

##4.3 绘制组合图形
corrplot(cor_ciber, method = "square", 
         type = "lower") # 下三角矩阵
corrplot(cor_ciber, method = "pie", 
         type = "upper") # 上三角矩阵

##4.4 相关矩阵排序
corrplot(cor_ciber, order = "AOE") # 特征向量角序
corrplot(cor_ciber, order = "FPC") # 第一主成分顺序
corrplot(cor_ciber, order = "hclust") # 按层次聚类
corrplot(cor_ciber, order = "alphabet") # 按字母顺序

##4.5 设置矩阵颜色
col = colorRampPalette(c('navy', 'white', 'yellow'))(40)
corrplot(cor_ciber, method = "color", col = col) # 矩阵颜色
corrplot(cor_ciber, method = "circle", bg = "grey") # 背景颜色

##4.6 设置文本标签属性
corrplot(cor_ciber, tl.pos = "n") # 不显示文本标签
corrplot(cor_ciber, tl.pos = "lt") # 在左边和顶部显示
corrplot(cor_ciber, tl.cex = 0.5) # 设置文本标签的缩放倍数

##4.7 设置显著水平
res1 <- cor.mtest(ciber.res, conf.level = .95)
corrplot(cor_ciber, p.mat = res1$p, 
         sig.level = .05)   #设置p值>0.05的不显示
corrplot(cor_ciber, p.mat = res1$p, 
         insig = "blank")   #设置p值>0.05的相关系数为空白
corrplot(cor_ciber, p.mat = res1$p, 
         insig = "p-value")
