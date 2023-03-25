################################################
################################################
### 作者：liujialeo
### 更新时间：2023-02-03
### Github: https://github.com/liujialeo

rm(list = ls())
dir.create("output/OR51E1_associated")
### 加载表达数据，这是vst标准化后的数据
load(file = "output/exprSet_heatmap.Rdata")
test <- exprSet[1:10,1:10]
##1.设定容器
correlation <- data.frame()
##2.准备数据
data <- as.data.frame(t(exprSet)) 
test <- data[1:10,1:10]
##3.获取基因列表
genelist <- colnames(data)
##4.指定基因
gene <- "OR51E1"
genedata <- as.numeric(data[,gene])
##5.开始for循环
for(i in 1:length(genelist)){
  ## 1.指示
  print(i)
  ## 2.计算
  dd = cor.test(genedata,as.numeric(data[,i]),method="spearman")
  ## 3.填充
  correlation[i,1] = gene
  correlation[i,2] = genelist[i]
  correlation[i,3] = dd$estimate
  correlation[i,4] = dd$p.value
}

colnames(correlation) <- c("gene1","gene2","cor","p.value")
save(correlation, file = "./output/OR51E1_association.Rdata")
## 排序
library(dplyr)

genelist = correlation$cor
names(genelist) = correlation$gene2
genelist = sort(genelist,decreasing = T)
head(genelist)
##############################
## 基因集 
library(msigdbr)
h_df <- msigdbr(species = "Homo sapiens") %>% 
  filter(gs_cat == "H") %>% 
  dplyr::select(gs_name,gene_symbol)

#################################
## enrichr
## 选取一定数目上调基因
upgene = head(names(genelist),1000)
## 选取一定数目下调基因
dngene = tail(names(genelist),1000)
## 输入文件
gcSample = list(up = upgene,down = dngene,all = c(upgene,dngene))

## 万能代码！！
library(clusterProfiler)
xx <- compareCluster(gcSample, fun="enricher",TERM2GENE = h_df)


## 作图
library(ggplot2)
dotplot(xx)
dotplot(xx)+ggplot2::scale_color_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"))



ggsave("output/OR51E1_associated/GO_enrichment.pdf")
#################################
## GSEA
y <- GSEA(genelist,TERM2GENE = h_df,pvalueCutoff = 0.05)
yd <- as.data.frame(y)
library(ggplot2)
dotplot(y,showCategory=15,split=".sign")+facet_grid(~.sign)
### 自定义画图

ggplot(y, showCategory = 38, aes(NES, forcats::fct_reorder(Description, NES))) + 
  #geom_segment(aes(xend=0, yend = Description)) + ##原点上连接直线
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_colour_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3")) +
  #scale_colour_gradientn(colours= c("#f7ca64", "#46bac2","#7e62a3")) +
  #scale_colour_gradientn(colours=c("#7fc97f","#984ea3")) +
  #scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  #scale_color_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE))+
  scale_size_continuous(range=c(2, 10)) +
  theme_bw() + 
  xlab("Normalized Enrichment Score") +
  ylab(NULL)

ggsave(filename = "./output/OR51E1_associated/GSEAplot.pdf",width = 8,height = 9)

library(enrichplot)
library(export)
pathway.id = "HALLMARK_PI3K_AKT_MTOR_SIGNALING"
gseaplot2(y,color = "green",geneSetID = pathway.id,pvalue_table = T)
ggsave(filename = "./output/OR51E1_associated/HALLMARK_PI3K_AKT_MTOR_SIGNALING.pdf")
graph2pdf(file="./output/OR51E1_associated/HALLMARK_PI3K_AKT_MTOR_SIGNALING.pdf")

####在图中标记出特定的基因
pathway.id = "HALLMARK_PI3K_AKT_MTOR_SIGNALING"
pg2 <- gseaplot2(y, color = "green",geneSetID = pathway.id,pvalue_table = T)
g1 <- c("RALB")
pg2[[1]] = pg2[[1]]  + geom_gsea_gene(g1, geom=ggrepel::geom_text_repel)
pg2
ggsave(filename = "./output/OR51E1_associated/HALLMARK_PI3K_AKT_MTOR_SIGNALING_labeled.pdf")

#####基因之间两两相关性
exprSet <- as.data.frame(t(exprSet))
test <- exprSet[1:10,1:10]
library(ggstatsplot)
library(cowplot)

p1 <- ggscatterstats(data = exprSet,
               y = OR51E1,
               x = ATF1,
               centrality.para = "mean",
               margins = "both",
               xfill = "#CC79A7",
               yfill = "#009E73",
               marginal.type = "histogram",
               title = "Relationship between OR51E1 and ATF1")

p2 <- ggscatterstats(data = exprSet,
                     y = OR51E1,
                     x = NGF,
                     centrality.para = "mean",
                     margins = "both",
                     xfill = "#CC79A7",
                     yfill = "#009E73",
                     marginal.type = "histogram",
                     title = "Relationship between OR51E1 and NGF")

p3 <- ggscatterstats(data = exprSet,
                     y = OR51E1,
                     x = RALB,
                     centrality.para = "mean",
                     margins = "both",
                     xfill = "#CC79A7",
                     yfill = "#009E73",
                     marginal.type = "histogram",
                     title = "Relationship between OR51E1 and RALB")

p4 <- ggscatterstats(data = exprSet,
                     y = OR51E1,
                     x = MYD88,
                     centrality.para = "mean",
                     margins = "both",
                     xfill = "#CC79A7",
                     yfill = "#009E73",
                     marginal.type = "histogram",
                     title = "Relationship between OR51E1 and MYD88")

plot_grid(p1,p2,p3,p4, nrow = 2, labels = LETTERS[1:4])
ggsave(filename = "./output/OR51E1_associated/OR51E1_correlation_with_4genes.pdf")
## R 包和网站 GTBA
## 如何用?
## 单基因批量相关性分析的妙用
## https://mp.weixin.qq.com/s/TfE2koPhSkFxTWpb7TlGKA
## 单基因批量相关性分析的GSEA
## https://mp.weixin.qq.com/s/sZJPW8OWaLNBiXXrs7UYFw