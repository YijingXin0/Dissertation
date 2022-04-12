rm(list=ls())
gc()

setwd('K:\\compound')
getwd()

expreset <- read.table('compound.txt',header = T)
expreset_agg <- aggregate(x = expreset[,2:7],by = list(expreset$gene_id),FUN = mean)
rownames(expreset_agg) <- expreset_agg$Group.1
expreset_agg <- expreset_agg[,-1]
expreset_agg_round <- round(expreset_agg)
condition <- factor(rep(c('control', 'treat'), each = 3))
coldata <- data.frame(row.names=colnames(expreset_agg_round), condition)

library(DESeq2)
dds <- DESeqDataSetFromMatrix(expreset_agg_round, DataFrame(condition), design= ~ condition )
head(dds)
dds <- DESeq(dds)
res <- results(dds)
summary(res)

diff_gene_deseq2 <-subset(res,pvalue< 0.05 & (log2FoldChange > 0.584962500721156 | log2FoldChange < -1))
a <- data.frame(diff_gene_deseq2, stringsAsFactors = FALSE, check.names = FALSE)
write.table(a, 'compound.DEGs.txt', col.names = NA, sep = '\t', quote = FALSE)
