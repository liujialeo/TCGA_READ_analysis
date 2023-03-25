################################################
################################################
### 作者：liujialeo
### 更新时间：2023-02-03
### Github: https://github.com/liujialeo


rm(list = ls())
load(file = "output/READ_RNASEQ_exprdf.Rdata")
class(expr_df)

exprSet <- as.data.frame(expr_df)
### 查看分组
### https://dwz.cn/WVgQUqfw
### 样本名称
TCGA_id <- colnames(exprSet)[-1]
table(substring(TCGA_id,14,15))
### 我们发现了7个转移的样本，本次分析，我们关注的是癌症和癌旁，先把转移的样本去掉
### 原发和转移的对比作为家庭作业

# TCGA_id <- TCGA_id[substring(TCGA_id,14,15)!="06"]

exprSet <- cbind(exprSet$gene_id,exprSet[,TCGA_id])
TCGA_id <- colnames(exprSet)[-1]
table(substring(TCGA_id,14,15))
colnames(exprSet)[1] <- "symbol"
### 创建metadata
sample <- ifelse(substring(TCGA_id,14,15)=="01","cancer","normal")
sample <- factor(sample,levels = c("normal","cancer"),ordered = F)
metadata <- data.frame(TCGA_id,sample) 
save(metadata,file = "output/READ_metadata.Rdata")
##去重
library(dplyr)
library(tibble)
exprSet <- exprSet %>% 
  ## rowMeans求出行的平均数(这边的.代表上面传入的数据)
  ## .[,-1]表示去掉出入数据的第一列，然后求行的平均值
  mutate(rowMean =rowMeans(.[,-1])) %>% 
  ## 把表达量的平均值按从大到小排序
  arrange(desc(rowMean)) %>% 
  ## 去重，symbol留下第一个
  distinct(symbol,.keep_all = T) %>% 
  ## 反向选择去除rowMean这一列
  dplyr::select(-rowMean)



library(DESeq2)
### 第一列有名称，所以tidy=TRUE
dds <-DESeqDataSetFromMatrix(countData=exprSet, 
                             colData=metadata, 
                             design=~sample,
                             tidy=TRUE)
nrow(dds)
### 如果一个基因在所有样本中的counts数小于等于1，我们就把他删掉
dds <- dds[rowSums(counts(dds))>1,]
nrow(dds)

### 数据标准化用于看聚类
### https://dwz.cn/xJTuI4aO
### 很耗时间
if(T){
  vsd <- vst(dds, blind = FALSE)
}
save(vsd,file = "output/READ_vsd.Rdata")
#load(file = "output/READ_vsd.Rdata")
### PCAf
plotPCA(vsd, "sample")
ggsave(filename = "output/PCA_READ.pdf")
### 保存数据用于热图
exprSet_vst <- as.data.frame(assay(vsd))
test <- exprSet_vst[1:10,1:10]
save(exprSet_vst,file = "output/READ_exprSet_vst.Rdata")

### 最困难的一步来了.燃烧你的小电脑。
### 这里还使用了并行化处理来解决速度的问题
### Deseq2 更新后速度大幅度提升
### https://dwz.cn/bb1jDs12

library(BiocParallel)
dds <- DESeq(dds,parallel = T)
save(dds,file="output/dds_very_long.Rdata")

load(file="output/dds_very_long.Rdata")

###################################################################
### 提取标准化后的数据，注意，我们data warngling 用的就是这个数据
normalized_counts <- as.data.frame(counts(dds, normalized=TRUE))
save(normalized_counts,file = "output/READ_normalized_counts.Rdata")
### 单独作图简易展示
plotCounts(dds, gene = "OR51E1", intgroup=c("sample"))
### 还可以把数据返回
plotdata <- plotCounts(dds, gene = "OR51E1", intgroup=c("sample","paire_info"),returnData = T)
library(ggplot2)
ggplot(plotdata,aes(x=sample,y=count,col=sample))+
  geom_jitter()+
  theme_bw()

#######################################################################
### logFC矫正,注意顺序哈，?号
contrast <- c("sample","cancer","normal")
dd1 <- results(dds, contrast=contrast, alpha = 0.05)
plotMA(dd1, ylim=c(-5,5))
### 发现样本需要做logFC校正
###dd2 <- lfcShrink(dds, contrast=contrast, res=dd1)
###plotMA(dd2, ylim=c(-5,5))

summary(dd1, alpha = 0.05)
library(dplyr)
library(tibble)
library(tidyr)
### 导出差异分析的结果
res <- dd1 %>% 
  data.frame() %>% 
  rownames_to_column("gene_id") 

### 基因注释
## 准备注释文件
### https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files
### 下载后先解压gencode.v22.annotation.gtf.gz
### gtf1 <- rtracklayer::import('gencode.v22.annotation.gtf')
### gtf_df <- as.data.frame(gtf1)
### save(gtf_df,file = "gtf_df.Rdata")

# load(file = "gtf_df.Rdata")
# ## 注释差异基因
# allDiff <- gtf_df %>% 
#   ##筛选gene
#   dplyr::filter(type=="gene") %>%
#   ## 筛选基因名称和类型
#   dplyr::select(c(gene_name,gene_id,gene_type)) %>% 
#   ## 合并
#   dplyr::inner_join(res,by ="gene_id")
allDiff <- res
save(allDiff,file = "output/READ_allDiff.Rdata")
### 现在有了差异基因列表，热图，火山图，GO，KEGG，GSEA，pathview都可以做了
### 自行完成