# This is an optional step
rm(list=ls())
options(stringsAsFactors = F)
load("DEG_result.Rdata")

# compare three methods difference
d1 <- data.frame(gene=rownames(Deseq_DEG),
                 deseq=Deseq_DEG$log2FoldChange)
d2 <- data.frame(gene=rownames(edgeR_DEG),
                 edgeR=edgeR_DEG$logFC)
d3 <- data.frame(gene=rownames(limma_voom_DEG),
                 limma=limma_voom_DEG$F)

tmp <- merge(d2, d3, by='gene')
plot(tmp[,2:3])
