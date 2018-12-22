######################################
# featureCounts result => row (genes), column (samples)
######################################

count <- read.csv("hisat2_counts.txt", sep = "\t")
count <- count[,-(2:6)]
rownames(count) <- count[,1]
count <- count[,-1]
tmp <- strsplit(colnames(count),"_")
colnames(count) <- unlist(lapply(tmp, function(x) strsplit(x,split="_")[1]))


######################################
# Deseq2 needs countData and colData
######################################
suppressMessages(library(DESeq2))
countData <- count[apply(count, 1, sum) > 0 ,]
(colData <- data.frame(row.names =colnames(countData), condition=factor(rep(c("day0","day1","day2","day3"),each=4))))
dds <- DESeqDataSetFromMatrix(
    countData = countData,
    colData = colData,
    design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition","day1", "day0"))
resOrdered <- res[order(res$padj),] 
head(resOrdered)
DEG <- as.data.frame(resOrdered)
View(DEG)

DEG <- na.omit(DEG)
Deseq_DEG <- DEG

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
    if(T){
        library(ggplot2)
        log2FoldChange_cutoff <- with(DEG,mean(abs(log2FoldChange)) + 2*sd(abs(log2FoldChange)) )
        
        DEG$change = as.factor(
            ifelse(DEG$pvalue < 0.05 & abs(DEG$log2FoldChange) > log2FoldChange_cutoff,
                   ifelse(DEG$log2FoldChange > log2FoldChange_cutoff ,'UP','DOWN'),'NOT')
        )
        this_tile <- paste0('Cutoff for log2FoldChange is ',round(log2FoldChange_cutoff,3),
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
        ggsave(g, filename = "vocalno.png", width = 10, height = 7 )
        
    }
}
