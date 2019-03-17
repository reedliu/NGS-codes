### ---------------
###
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-3-3
### Email: jieandze1314@gmail.com
### Blog: https://reedliu.github.io/
### CAAS/AGIS/SDAU 
###
### ---------------

rm(list = ls())
options(stringsAsFactors = F)
# for salmon
count <- read.csv("Salmon.isoform.counts.matrix", sep = "\t", row.names = 1)
# for RSEM
count <- read.csv("RSEM.gene.counts.matrix", sep = "\t", row.names = 1)

count <- apply(count, 1:2, round)
count <- as.data.frame(count)
countData <- count[apply(count, 1, sum) > 0 ,]
(colData <- data.frame(row.names =colnames(countData), condition=factor(rep(c("SS","LFC2","XXCAD","RS","DU"),each=3))))
group_list <- rep(c("SS","LFC2","XXCAD","RS","DU"),each=3)

##################
# Deseq2 
##################
suppressMessages(library(DESeq2))

dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData = colData,
  design = ~ condition)
dds <- DESeq(dds)

# for salmon
if(F){
  # for DU vs SS
  res <- results(dds, contrast = c("condition","DU", "SS"))
  resOrdered <- res[order(res$padj),] 
  salmon_DESeq2_DU_SS <- as.data.frame(resOrdered)
  View(salmon_DESeq2_DU_SS)
  salmon_DESeq2_DU_SS <- na.omit(salmon_DESeq2_DU_SS)
  # for LFC2 vs SS
  res <- results(dds, contrast = c("condition","LFC2", "SS"))
  resOrdered <- res[order(res$padj),] 
  salmon_DESeq2_LFC2_SS <- as.data.frame(resOrdered)
  View(salmon_DESeq2_LFC2_SS)
  salmon_DESeq2_LFC2_SS <- na.omit(salmon_DESeq2_LFC2_SS)
  # for RS vs SS
  res <- results(dds, contrast = c("condition","RS", "SS"))
  resOrdered <- res[order(res$padj),] 
  salmon_DESeq2_RS_SS <- as.data.frame(resOrdered)
  View(salmon_DESeq2_RS_SS)
  salmon_DESeq2_RS_SS <- na.omit(salmon_DESeq2_RS_SS)
  # for XXCAD vs SS
  res <- results(dds, contrast = c("condition","XXCAD", "SS"))
  resOrdered <- res[order(res$padj),] 
  salmon_DESeq2_XXCAD_SS <- as.data.frame(resOrdered)
  View(salmon_DESeq2_XXCAD_SS)
  salmon_DESeq2_XXCAD_SS <- na.omit(salmon_DESeq2_XXCAD_SS)
  
}
save("salmon_DESeq2_DU_SS","salmon_DESeq2_LFC2_SS",
    "salmon_DESeq2_RS_SS","salmon_DESeq2_XXCAD_SS", 
    file = "DESeq2-salmon.Rdata")

# for RSEM
if(T){
  # for DU vs SS
  res <- results(dds, contrast = c("condition","DU", "SS"))
  resOrdered <- res[order(res$padj),] 
  RSEM_DESeq2_DU_SS <- as.data.frame(resOrdered)
  View(RSEM_DESeq2_DU_SS)
  RSEM_DESeq2_DU_SS <- na.omit(RSEM_DESeq2_DU_SS)
  # for LFC2 vs SS
  res <- results(dds, contrast = c("condition","LFC2", "SS"))
  resOrdered <- res[order(res$padj),] 
  RSEM_DESeq2_LFC2_SS <- as.data.frame(resOrdered)
  RSEM_DESeq2_LFC2_SS <- na.omit(RSEM_DESeq2_LFC2_SS)
  # for RS vs SS
  res <- results(dds, contrast = c("condition","RS", "SS"))
  resOrdered <- res[order(res$padj),] 
  RSEM_DESeq2_RS_SS <- as.data.frame(resOrdered)
  RSEM_DESeq2_RS_SS <- na.omit(RSEM_DESeq2_RS_SS)
  # for XXCAD vs SS
  res <- results(dds, contrast = c("condition","XXCAD", "SS"))
  resOrdered <- res[order(res$padj),] 
  RSEM_DESeq2_XXCAD_SS <- as.data.frame(resOrdered)
  RSEM_DESeq2_XXCAD_SS <- na.omit(RSEM_DESeq2_XXCAD_SS)
  
}
save("RSEM_DESeq2_DU_SS","RSEM_DESeq2_LFC2_SS",
     "RSEM_DESeq2_RS_SS","RSEM_DESeq2_XXCAD_SS", 
     file = "DESeq2-RSEM.Rdata")

##################
# Plots
##################

if (F){
  # Heatmap
  countData = countData 
  top30_gene=head(rownames(DEG),30)
  top30_matrix=countData[top30_gene,]
  top30_matrix=t(scale(t(top30_matrix)))
  
  nrDEG=DEG
  library(pheatmap)
  choose_gene=head(rownames(nrDEG),50) ## 50 maybe better
  choose_matrix=countData[choose_gene,]
  choose_matrix=t(scale(t(choose_matrix)))
  pheatmap(choose_matrix,filename = 'DEG_top50_heatmap.png')
  
  # more complex one
  if(T){
    library(ComplexHeatmap)
    #install.packages("dendextend")
    library(dendextend)
    library(dendsort)
    e_mean = tail(sort(apply(countData,1,mean)),30)
    #top30_cluster_cols = hclust(dist(t(top30_matrix)))
    sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...))) 
    top30_cluster_cols = sort_hclust(hclust(dist(t(top30_matrix))))
    top30_cluster_rows <- sort_hclust(hclust(dist(top30_matrix)))
    base_mean = e_mean
    ef_lable = HeatmapAnnotation(text = anno_text(colnames(countData), rot = 45,offset = unit(23, "mm"), gp = gpar(fontsize	= 8.5)))
    p <- Heatmap(top30_matrix, name = "expresssion", 
                 column_title = "Samples", 
                 column_title_gp = gpar(fontsize = 14, fontface = "bold"),
                 row_title = "Genes",
                 row_title_gp = gpar(fontsize = 14, fontface = "bold"),
                 row_names_gp = gpar(fontsize = 8.5),
                 row_names_side = "left",
                 cluster_rows = color_branches(top30_cluster_rows, k = 4),
                 cluster_columns = color_branches(top30_cluster_cols, k = 3),
                 show_column_names = F,
                 bottom_annotation = ef_lable, bottom_annotation_height = unit(0.5, "cm"))+
      Heatmap(base_mean, name = "base_mean", show_row_names = FALSE, width = unit(5, "mm"))
    print(p)
  }
  
  # Vocalno
  DEG = RSEM_DESeq2_DU_SS
  DEG = RSEM_DESeq2_LFC2_SS
  DEG = RSEM_DESeq2_RS_SS
  DEG = RSEM_DESeq2_XXCAD_SS
  if(T){
    library(ggplot2)
    log2FoldChange_cutoff <- with(DEG,mean(abs(log2FoldChange)) + 2*sd(abs(log2FoldChange)) )
    
    DEG$change = as.factor(
      ifelse(DEG$pvalue < 0.05 & abs(DEG$log2FoldChange) > log2FoldChange_cutoff,
             ifelse(DEG$log2FoldChange > log2FoldChange_cutoff ,'UP','DOWN'),'NOT')
    )
    this_tile <- paste0('RSEM DESeq2 for XXCAD vs SS',
                        '\nCutoff for log2FoldChange is ',round(log2FoldChange_cutoff,3),
                        '\nThe number of up gene is ',nrow(DEG[DEG$change =='UP',]) ,
                        '\nThe number of down gene is ',nrow(DEG[DEG$change =='DOWN',])
    )
    g= ggplot(data=DEG, 
              aes(x=log2FoldChange, y=-log10(pvalue), 
                  color=change)) +
      geom_point(alpha=0.4, size=1.75) +
      theme_set(theme_set(theme_bw(base_size=20)))+
      xlab("log2 fold change") + ylab("-log10 p-value") +
      ggtitle( this_tile ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
      scale_colour_manual(values = c('blue','black','red')) 
    ggsave(g, filename = "RSEM_deseq2_XXCAD_SS.png", width = 10, height = 7 )
    
  }
}




