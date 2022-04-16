rm(list = ls())

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")

library(edgeR)

setwd("I:\\Batch_1")
x1 <-read.table('counts_anno.xls',header = T)
rownames(x1) <- x1[,1]
x1 <- x1[,-1]
colnames <- colnames(x1)

DEG_Statistics <-as.data.frame(c(1:2))
rownames(DEG_Statistics) <- c("Down","Up" )

for(i in 4:79){
  y <- as.matrix(x1[,c(1,2,3,i)])
  group <- c(rep(1,3),2)
  y <- DGEList(counts=y, group = group)
  keep <- rowSums(cpm(y)>1) >= 1
  y <- y[keep, , keep.lib.sizes=FALSE]
  y_bcv <- calcNormFactors(y)
  bcv <- 0.01
  et <- exactTest(y_bcv, dispersion = bcv ^ 2)
  results <- cbind(y$counts,et$table)
  results=results[order(-results[,5]),]
  write.csv(results,file=paste(colnames[i],"ALL_genes.csv",sep=""))
  print(i)
  i <- i+1
}
