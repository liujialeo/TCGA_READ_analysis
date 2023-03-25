################################################
################################################
### 作者：liujialeo
### 更新时间：2023-02-03
### Github: https://github.com/liujialeo

rm(list = ls())
### 本节任务: 转录组数据变成清洁数据
###################################################################
### 一个基因，或者多个基因的差异作图
### 需要清洁数据
load("output/READ_exprSet_vst.Rdata")
load("output/READ_allDiff.Rdata")
load("output/READ_metadata.Rdata")
###删除正常样本
exprSet_vst_t <- as.data.frame(t(exprSet_vst))
test <- exprSet_vst_t[1:10,1:10]
exprSet_vst_t <- cbind(rownames(exprSet_vst_t),exprSet_vst_t)
colnames(exprSet_vst_t)[1] <- "TCGA_id"
library(dplyr)
exprSet_vst_t <- inner_join(metadata,exprSet_vst_t)
exprSet_vst_t <- exprSet_vst_t%>%
  filter(sample =="cancer")
exprSet_vst_t <- exprSet_vst_t[,-2]
rownames(exprSet_vst_t) <- exprSet_vst_t[,1]
exprSet_vst_t <- exprSet_vst_t[,-1]
exprSet_vst <- as.data.frame(t(exprSet_vst_t))

exprSet <- exprSet_vst

### 保存数据，画热图
save(exprSet,file = "output/exprSet_heatmap.Rdata")

### 行列转置
exprSet <- t(exprSet)
exprSet <- as.data.frame(exprSet)
test <- exprSet[,1:10]

### 7.添加分组
colnames(metadata) <- c("sample","group")
exprSet <- cbind(group= metadata$group,exprSet)
test <- exprSet[,1:10]
save(exprSet,file = "output/exprSet_tidy.Rdata")

###########################################################

rm(list = ls())
library(RColorBrewer)
display.brewer.all()
###获取自定义颜色
colset <- brewer.pal(10,"Paired") ##brewer.pal选择颜色的最小种类是3
load(file = "output/exprSet_tidy.Rdata")

## steal plot
my_comparisons <- list(
  c("cancer", "normal")
)
library(ggpubr)
ggboxplot(
  exprSet, x = "group", y = "OR51E1",
  color = "group", palette = c("#00AFBB", "#E7B800"),
  add = "jitter"
)+
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

## 改写成函数
diffplot <- function(gene){
  my_comparisons <- list(
    c("treat", "con")
  )
  library(ggpubr)
  ggboxplot(
    exprSet, x = "group", y = gene,
    color = "group", palette = c("#00AFBB", "#E7B800"),
    add = "jitter"
  )+
    stat_compare_means(comparisons = my_comparisons, method = "t.test")
}

### AGR3,ESR1,SLC4A10, ALPP,VSIR,PLA2G2F
diffplot("NFKBIE")
diffplot("ALPP")

## 多个基因作图查看
## 先把基因提取出来
genelist <- c("OR51E1","FFAR2","FFAR3")
## 再提取表达量，使用名称选取行
data <- exprSet[,c("group",genelist)]
## 用pivot_longer调整数据，数据变长，增加的是行
library(tidyr)
data <- data %>% 
  pivot_longer(cols=-1,
               names_to= "gene",
               values_to = "expression")
## 多基因作图
## 作图
ggplot(data = data,aes(x=group,y=expression,fill=group))+
  geom_boxplot()+
  geom_jitter(size= 1)+
  theme_bw()+
  #facet_grid(.~gene)+
  facet_wrap(~gene,nrow = 1)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test")+
  scale_fill_manual(values= colset[c(1,2)] )

ggsave(filename = "./output/multi_genes_expression.pdf", width = 15, height = 5 )
