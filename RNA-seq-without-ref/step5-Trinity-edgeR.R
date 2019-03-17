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
# edgeR
##################
suppressMessages(library(edgeR))
expr <- countData
grp <- group_list

e <- DGEList(counts=expr,group=factor(grp))
keep <- rowSums(cpm(e)>1) >= 2
e <- e[keep, , keep.lib.sizes=FALSE]
e$samples$lib.size <- colSums(e$counts)
e <- calcNormFactors(e)
e$samples

DEG=e

design <- model.matrix(~0+factor(grp,levels = c("SS","LFC2","XXCAD","RS","DU")))
rownames(design)<-colnames(DEG)
colnames(design)<-levels(factor(grp,levels = c("SS","LFC2","XXCAD","RS","DU")))


DEG <- estimateGLMCommonDisp(DEG,design)
DEG <- estimateGLMTrendedDisp(DEG, design)
DEG <- estimateGLMTagwiseDisp(DEG, design)

fit <- glmFit(DEG, design)

# for salmon
if(F){
  # for DU vs SS
  DU_SS_lrt <- glmLRT(fit,  contrast=c(-1,1,0,0,0)) 
  DU_SS_DEG=topTags(DU_SS_lrt, n=nrow(DEG))
  salmon_edgeR_DU_SS=as.data.frame(DU_SS_DEG)

  
  # for LFC2 vs SS
  LFC2_SS_lrt <- glmLRT(fit,  contrast=c(-1,0,1,0,0)) 
  LFC2_SS_DEG=topTags(LFC2_SS_lrt, n=nrow(DEG))
  salmon_edgeR_LFC2_SS=as.data.frame(LFC2_SS_DEG)
  
  # for RS vs SS
  RS_SS_lrt <- glmLRT(fit,  contrast=c(-1,0,0,1,0)) 
  RS_SS_DEG=topTags(RS_SS_lrt, n=nrow(DEG))
  salmon_edgeR_RS_SS=as.data.frame(RS_SS_DEG)
  
  # for XXCAD vs SS
  XXCAD_SS_lrt <- glmLRT(fit,  contrast=c(-1,0,0,0,1)) 
  XXCAD_SS_DEG=topTags(XXCAD_SS_lrt, n=nrow(DEG))
  salmon_edgeR_XXCAD_SS=as.data.frame(XXCAD_SS_DEG)
  
  save("salmon_edgeR_DU_SS","salmon_edgeR_LFC2_SS",
       "salmon_edgeR_RS_SS","salmon_edgeR_XXCAD_SS", 
       file = "edgeR-salmon.Rdata")
  
}

# for RSEM
if(T){
  
  # for LFC2 vs SS
  LFC2_SS_lrt <- glmLRT(fit,  contrast=c(-1,1,0,0,0)) 
  LFC2_SS_DEG=topTags(LFC2_SS_lrt, n=nrow(DEG))
  RSEM_edgeR_LFC2_SS=as.data.frame(LFC2_SS_DEG)
  
  # for XXCAD vs SS
  XXCAD_SS_lrt <- glmLRT(fit,  contrast=c(-1,0,1,0,0)) 
  XXCAD_SS_DEG=topTags(XXCAD_SS_lrt, n=nrow(DEG))
  RSEM_edgeR_XXCAD_SS=as.data.frame(XXCAD_SS_DEG)

  # for RS vs SS
  RS_SS_lrt <- glmLRT(fit,  contrast=c(-1,0,0,1,0)) 
  RS_SS_DEG=topTags(RS_SS_lrt, n=nrow(DEG))
  RSEM_edgeR_RS_SS=as.data.frame(RS_SS_DEG)
  
  
  # for DU vs SS
  DU_SS_lrt <- glmLRT(fit,  contrast=c(-1,0,0,0,1)) 
  DU_SS_DEG=topTags(DU_SS_lrt, n=nrow(DEG))
  RSEM_edgeR_DU_SS=as.data.frame(DU_SS_DEG)
  
  
  save("RSEM_edgeR_DU_SS","RSEM_edgeR_LFC2_SS",
       "RSEM_edgeR_RS_SS","RSEM_edgeR_XXCAD_SS", 
       file = "edgeR-RSEM.Rdata")
  
}

###########################################
DEG = salmon_edgeR_DU_SS
DEG = salmon_edgeR_LFC2_SS
DEG = salmon_edgeR_RS_SS
DEG = salmon_edgeR_XXCAD_SS
if(T){
  library(ggplot2)
  logFC_cutoff <- with(DEG,mean(abs(logFC)) + 2*sd(abs(logFC)) )
  
  DEG$change = as.factor(
    ifelse(DEG$PValue < 0.05 & abs(DEG$logFC) > logFC_cutoff,
           ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
  )
  this_tile <- paste0('Salmon edgeR for XXCAD vs SS',
                      '\nCutoff for logFC is ',round(logFC_cutoff,3),
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
  ggsave(g, filename = "salmon_edgeR_XXCAD_SS.png", width = 10, height = 7 )
  
}

