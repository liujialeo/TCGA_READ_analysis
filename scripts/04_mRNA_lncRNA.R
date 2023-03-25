################################################
################################################
### 作者：liujialeo
### 更新时间：2023-02-03
### Github: https://github.com/liujialeo

### 分别提取mRNA和lncRNA
rm(list = ls())
load(file = "BRCA_normalized_counts.Rdata")
expr_df <- cbind(gene_id= rownames(normalized_counts),normalized_counts,stringsAsFactors=F)
## 准备注释文件
### https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files
### 下载后先解压gencode.v22.annotation.gtf.gz
### gtf1 <- rtracklayer::import('gencode.v22.annotation.gtf')
### gtf_df <- as.data.frame(gtf1)
### save(gtf_df,file = "gtf_df.Rdata")

load(file = "gtf_df.Rdata")
test <- gtf_df[1:100,]

library(dplyr)
library(tidyr)
## 提取mRNA
mRNA_exprSet <- gtf_df %>% 
  #筛选gene,和编码指标
  dplyr::filter(type=="gene",gene_type=="protein_coding") %>%
  dplyr::select(c(gene_name,gene_id,gene_type)) %>% 
  dplyr::inner_join(expr_df,by ="gene_id") %>% 
  tidyr::unite(gene_id,gene_name,gene_id,gene_type,sep = "|")
test <- mRNA_exprSet[1:10,1:10]
save(mRNA_exprSet,file = "BRCA_mRNA_exprSet.Rdata")

### lncRNA
### 如何定义非编码RNA呢？
### https://www.gencodegenes.org/pages/biotypes.html
ncRNA <- c("3prime_overlapping_ncrna","antisense","lincRNA",
           "macro_lncRNA","non_coding",
           "processed_transcript","sense_intronic",
           "sense_overlapping")

LncRNA_exprSet <- gtf_df %>% 
  dplyr::filter(type=="gene",gene_type %in% ncRNA) %>% 
  dplyr::select(c(gene_name,gene_id,gene_type)) %>% 
  dplyr::inner_join(expr_df,by ="gene_id") %>% 
  tidyr::unite(gene_id,gene_name,gene_id,gene_type,sep = "|")

test <- LncRNA_exprSet[1:10,1:10]
save(LncRNA_exprSet,file = "BRCA_lncRNA_exprSet.Rdata")

### 可以集体合并，同时提取编码和非编码的信息
mytype<- c("protein_coding","3prime_overlapping_ncrna","antisense","lincRNA",
           "macro_lncRNA","non_coding",
           "processed_transcript","sense_intronic",
           "sense_overlapping")
all_exprSet <- gtf_df %>% 
  dplyr::filter(type=="gene",gene_type %in% mytype) %>% 
  dplyr::select(c(gene_name,gene_id,gene_type)) %>% 
  dplyr::inner_join(expr_df,by ="gene_id") %>% 
  tidyr::unite(gene_id,gene_name,gene_id,gene_type,sep = "|")
save(all_exprSet,file = "BRCA_all_exprSet.Rdata")

#########################################################################
### 那现在能干什么呢？
### 1.调整数据格式到清洁数据对接ggplot2
### 三大步，基因注释，行列转换，分组信息
### 2.单基因的批量相关性分析
### 参考a:https://dwz.cn/ebZiHkEK
### 参考b:https://dwz.cn/nQqbZQaH
### 3.单基因GSEA分析，神器
### 可以注释任意基因
### 核心点，在于用单个基因的表达量把样本分成高表达和低表达
### 之后就变成了两组RNA-seq分析，接着进行GSEA分析