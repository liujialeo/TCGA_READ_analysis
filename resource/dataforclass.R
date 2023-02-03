load(file = "resource/GSE153250_ESR1_74186928.Rds")

metadata = kosadata$metadata
exprSet = kosadata$exprSet
metadata = metadata[,c("GSM","group")]

write.csv(exprSet,file = "exprSet_counts.csv")

exprSet <- data.table::fread("exprSet_counts.csv",data.table = F)
colnames(exprSet)[1] <- "geneid"
write.csv(exprSet,file = "exprSet_counts.csv",row.names = F)
exprSet <- data.table::fread("exprSet_counts.csv",data.table = F)


### 读取metadata
metadata <- data.table::fread("data/metadata.txt",data.table = F,header = F)
 
rownames(metadata) <- metadata[,1]
metadata <- metadata[colnames(exprSet)[-1],]
metadata$group <- ifelse(grepl("ESR1",metadata$V2),"treat","control")


