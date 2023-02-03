################################################
################################################
### 作者：liujialeo
### 更新时间：2023-02-03
### Github: https://github.com/liujialeo

### 1.从GDC下载counts数据
## https://portal.gdc.cancer.gov/repository
dir.create("rawdata")
manifest <- "gdc_manifest_20230114_023951.txt"
rawdata <- "rawdata"
command <- sprintf("./gdc-client download -m %s -d %s",
                   manifest,
                   rawdata)
system(command = command)

### 备选方案
### 2.下载数据，使用TCGA GDC以及https://xena.ucsc.edu/
### 注意事项，下载的数据是log之后的数据，所以不能直接使用Deseq2来做差异分析。