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
# edgeR
######################################
suppressPackageStartupMessages(library(edgeR))

countData <- count[apply(count, 1, sum) > 0 ,]
group_list <- rep(c("day0","day1","day2","day3"),each=4)

if(F){
    dge <- DGEList(counts=countData,group=factor(group_list))
    dge$samples$lib.size <- colSums(dge$counts)
    dge <- calcNormFactors(dge) 
    dge$samples
    
    
    design <- model.matrix(~0+factor(group_list))
    rownames(design)<-colnames(dge)
    colnames(design)<-levels(factor(group_list))
    
    dge <- estimateGLMCommonDisp(dge,design)
    dge <- estimateGLMTrendedDisp(dge, design)
    dge <- estimateGLMTagwiseDisp(dge, design)
    
    fit <- glmFit(dge, design)
    
    # adjust contrast by design
    lrt <- glmLRT(fit,  contrast=c(1,0,0,0)) 
    nrDEG=topTags(lrt, n=nrow(countData))
    nrDEG=as.data.frame(nrDEG)
    head(nrDEG)
    edgeR_DEG <- nrDEG
    
}

# Volcano
if(F){
    DEG <- edgeR_DEG
    library(ggplot2)
    logFC_cutoff <- with(DEG,mean(abs(logFC)) + 2*sd(abs(logFC)) )
    
    DEG$change = as.factor(
        ifelse(DEG$PValue < 0.05 & abs(DEG$logFC) > logFC_cutoff,
               ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
    )
    this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                        '\nThe number of up gene is ',nrow(DEG[DEG$change =='UP',]) ,
                        '\nThe number of down gene is ',nrow(DEG[DEG$change =='DOWN',])
    )
    g= ggplot(data=DEG, 
              aes(x=logFC, y=-log10(PValue), 
                  color=change)) +
        geom_point(alpha=0.4, size=1.75) +
        theme_set(theme_set(theme_bw(base_size=20)))+
        xlab("log2 fold change") + ylab("-log10 p-value") +
        ggtitle( this_tile ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
        scale_colour_manual(values = c('blue','black','red')) 
    ggsave(g, filename = "edgeR_vocalno.png", width = 10, height = 7 )
    
}