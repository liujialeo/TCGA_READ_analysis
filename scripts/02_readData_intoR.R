################################################
################################################
### 作者：liujialeo
### 更新时间：2023-02-03
### Github: https://github.com/liujialeo


rm(list = ls())
# 02_读入数据
### 原始文件在rawdata
### 首先把所有数据读入在一个文件夹中
### 创建目的文件夹，可以右击创建，也可以dir.create函数创建
dir.create("data_in_one")
dir("rawdata/")
### 使用for循环来批量做，回顾for循环的要点
for (dirname in dir("rawdata/")){  
  ## 要查看的单个文件夹的绝对路径
  mydir <- paste0(getwd(),"/rawdata/",dirname)
  ##找到对应文件夹中的文件并提取名称，pattern表示模式，可以是正则表达式
  file <- list.files(mydir,pattern = "*.counts")
  ## 当前文件的绝对路径是
  myfile <- paste0(mydir,"/",file)
  #复制这个文件到目的文件夹
  file.copy(myfile,"data_in_one")  
}

### 在终端也可以简单实现这个效果：cp rawdata/*/*.gz data_in_one

## 文件和TCGA ID的对应关系
metadata <- jsonlite::fromJSON("metadata.cart.2023-01-14.json")
library(dplyr)
metadata_id <- metadata %>% 
  dplyr::select(c(file_name,associated_entities)) 

## 提取对应的文件
naid_df <- data.frame()
for (i in 1:nrow(metadata)){
  naid_df[i,1] <- metadata_id$file_name[i]
  naid_df[i,2] <- metadata_id$associated_entities[i][[1]]$entity_submitter_id
}

colnames(naid_df) <- c("filename","TCGA_id")

################################################################################
## 读入数据
### 1.存储文件读入的顺序
files <- dir("data_in_one")

## 构建函数
myfread <- function(files){
  data.table::fread(paste0("data_in_one/",files))[-(1:4),4]
}

## 测试文件以及时间，大概2秒一个
system.time(test <- myfread(files[1]))

###############################################
## lapply 的使用方法
## 1.for 循环正常速度是1222*2 大概40分钟，所以不适合
## 2.此处使用的是lapply+ function的形式，目的是为了提速
## 3.提速的实现方式是，lapply的并行化处理，使用的是future.apply这个包

library(future.apply)
plan(multisession)
### 用system.time返回时间
### 我的台式机12个线程大概需要440秒
system.time(f <- future_lapply(files,myfread))

## 列表变成数据框，回想晾衣杆的结构
expr_df <- as.data.frame(do.call(cbind,f))

## 为了把文件名称和TCGAid对应起来
rownames(naid_df) <- naid_df[,1]
naid_df <- naid_df[files,]

## 命名
colnames(expr_df) <- naid_df$TCGA_id

test <- expr_df[1:100,1:4]
## 加上一列
gene_id <- data.table::fread(paste0("data_in_one/",files[1]))[-(1:4),2]

expr_df <- cbind(gene_id=gene_id,expr_df)
test <- expr_df[1:100,1:4]

dir.create("output")
save(expr_df,file = "output/READ_RNASEQ_exprdf.Rdata")
