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
load('DESeq2-salmon.Rdata')
load('edgeR-salmon.Rdata')
load('method_compare.Rdata')
# for RSEM
load('DESeq2-RSEM.Rdata')
load('edgeR-RSEM.Rdata')

####################################
## get all genes and FC to check correlation 
# 帮助判断不同方法是否正确（理论上两种方法得到的结果应该在一条直线上）
####################################
# for DU vs SS
DU_SS_deseq <- data.frame(gene=rownames(salmon_DESeq2_DU_SS),
                          DU_SS_Deseq=salmon_DESeq2_DU_SS$log2FoldChange)

DU_SS_edgeR <- data.frame(gene=rownames(salmon_edgeR_DU_SS),
                          DU_SS_edgeR=salmon_edgeR_DU_SS$logFC)

test1 <- merge(DU_SS_deseq,DU_SS_edgeR,by="gene")
plot(test1[2:3])

# for LFC2 vs SS
LFC2_SS_Deseq <- data.frame(gene=rownames(salmon_DESeq2_LFC2_SS),
                            LFC2_SS_Deseq=salmon_DESeq2_LFC2_SS$log2FoldChange)

LFC2_SS_edgeR <- data.frame(gene=rownames(salmon_edgeR_LFC2_SS),
                            LFC2_SS_edgeR=salmon_edgeR_LFC2_SS$logFC)

test2 <- merge(LFC2_SS_Deseq,LFC2_SS_edgeR,by="gene")
plot(test2[2:3])

# for RS vs SS
RS_SS_Deseq <- data.frame(gene=rownames(salmon_DESeq2_RS_SS),
                            RS_SS_Deseq=salmon_DESeq2_RS_SS$log2FoldChange)

RS_SS_edgeR <- data.frame(gene=rownames(salmon_edgeR_RS_SS),
                            RS_SS_edgeR=salmon_edgeR_RS_SS$logFC)

test3 <- merge(RS_SS_Deseq,RS_SS_edgeR,by="gene")
plot(test3[2:3])

# for XXCAD vs SS
XXCAD_SS_Deseq <- data.frame(gene=rownames(salmon_DESeq2_XXCAD_SS),
                            XXCAD_SS_Deseq=salmon_DESeq2_XXCAD_SS$log2FoldChange)

XXCAD_SS_edgeR <- data.frame(gene=rownames(salmon_edgeR_XXCAD_SS),
                          XXCAD_SS_edgeR=salmon_edgeR_XXCAD_SS$logFC)

test4 <- merge(XXCAD_SS_Deseq,XXCAD_SS_edgeR,by="gene")
plot(test4[2:3])

save(test1,test2,test3,test4,file="method_compare.Rdata")

# 组合准确性分析结果
par(mfrow=c(2,2))
plot(test1[2:3])
plot(test2[2:3])
plot(test3[2:3])
plot(test4[2:3])
####################################
## get UP and DOWN genes for each method to draw Venn plot
# 分别得到不同方法的上调、下调基因，然后分别做Venn图
####################################
  # for RSEM DESeq2
  DEG = RSEM_DESeq2_DU_SS
  DEG = RSEM_DESeq2_LFC2_SS
  DEG = RSEM_DESeq2_RS_SS
  DEG = RSEM_DESeq2_XXCAD_SS
  # for salmon DESeq2
  DEG = salmon_edgeR_DU_SS
  DEG = salmon_edgeR_LFC2_SS
  DEG = salmon_edgeR_RS_SS
  DEG = salmon_edgeR_XXCAD_SS
  
  logFC_cutoff <- with(DEG,mean(abs(logFC)) + 2*sd(abs(logFC)) )
  
  DEG$change = as.factor(ifelse(DEG$PValue < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                                ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
  )
  # for DU
  DU_up <- rownames(DEG[DEG$change =='UP',])
  DU_down <- rownames(DEG[DEG$change =='DOWN',])
  # for LFC2
  LFC2_up <- rownames(DEG[DEG$change =='UP',])
  LFC2_down <- rownames(DEG[DEG$change =='DOWN',])
  # for RS
  RS_up <- rownames(DEG[DEG$change =='UP',])
  RS_down <- rownames(DEG[DEG$change =='DOWN',])
  # for XXCAD
  XXCAD_up <- rownames(DEG[DEG$change =='UP',])
  XXCAD_down <- rownames(DEG[DEG$change =='DOWN',])
  
  save(DU_up,DU_down,LFC2_up,LFC2_down,XXCAD_up,
       XXCAD_down,RS_up,RS_down,file = "salmon_edgeR_DEGs.Rdata")

  ## Venn plot
  if(!require(VennDiagram))install.packages('VennDiagram')
  library (VennDiagram)
  if(T){
    # first look at UP comparison
    venn.diagram(x= list(DESeq2_DU_UP = DU_up,
                         DESeq2_LFC2_UP = LFC2_up,
                         DESeq2_RS_UP = RS_up,
                         DESeq2_XXCAD_UP = XXCAD_up),
                 filename = "DESeq2-RSEM-compare_UP.png",
                 height = 800, width = 1200,
                 resolution =500,
                 imagetype="png",
                 col="transparent",
                 fill=c("green","darkorchid1","yellow", "orange"),
                 alpha = 0.50,
                 cex=0.3,
                 cat.cex=0.25)
    # DOWN comparison
    venn.diagram(x= list(DESeq2_DU_DOWN = DU_down,
                         DESeq2_LFC2_DOWN = LFC2_down,
                         DESeq2_RS_DOWN = RS_down,
                         DESeq2_XXCAD_DOWN = XXCAD_down),
                 filename = "DESeq2-RSEM-compare_DOWN.png",
                 height = 900, width = 1400,
                 resolution =500,
                 imagetype="png",
                 col="transparent",
                 fill=c("green","darkorchid1","yellow", "orange"),
                 alpha = 0.50,
                 cex=0.3,
                 cat.cex=0.25)
  }
  
  ##############################################
  ### get edgeR UP and DOWN genes
  
  # for RSEM edgeR
  DEG = RSEM_edgeR_DU_SS
  DEG = RSEM_edgeR_LFC2_SS
  DEG = RSEM_edgeR_RS_SS
  DEG = RSEM_edgeR_XXCAD_SS
  
  logFC_cutoff <- with(DEG,mean(abs(logFC)) + 2*sd(abs(logFC)) )
  
  DEG$change = as.factor(ifelse(DEG$PValue < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                                ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
  )
  # for DU
  DU_up <- rownames(DEG[DEG$change =='UP',])
  DU_down <- rownames(DEG[DEG$change =='DOWN',])
  # for LFC2
  LFC2_up <- rownames(DEG[DEG$change =='UP',])
  LFC2_down <- rownames(DEG[DEG$change =='DOWN',])
  # for RS
  RS_up <- rownames(DEG[DEG$change =='UP',])
  RS_down <- rownames(DEG[DEG$change =='DOWN',])
  # for XXCAD
  XXCAD_up <- rownames(DEG[DEG$change =='UP',])
  XXCAD_down <- rownames(DEG[DEG$change =='DOWN',])


## Venn plot
if(!require(VennDiagram))install.packages('VennDiagram')
library (VennDiagram)
if(T){
  # first look at UP comparison
  venn.diagram(x= list(edegR_DU_UP = DU_up,
                       edegR_LFC2_UP = LFC2_up,
                       edegR_RS_UP = RS_up,
                       edegR_XXCAD_UP = XXCAD_up),
               filename = "edegR-RSEM-compare_UP.png",
               height = 800, width = 1200,
               resolution =500,
               imagetype="png",
               col="transparent",
               fill=c("green","darkorchid1","yellow", "orange"),
               alpha = 0.50,
               cex=0.3,
               cat.cex=0.25)
  # DOWN comparison
  venn.diagram(x= list(edegR_DU_DOWN = DU_down,
                       edegR_LFC2_DOWN = LFC2_down,
                       edegR_RS_DOWN = RS_down,
                       edegR_XXCAD_DOWN = XXCAD_down),
               filename = "edegR-RSEM-compare_DOWN.png",
               height = 900, width = 1400,
               resolution =500,
               imagetype="png",
               col="transparent",
               fill=c("green","darkorchid1","yellow", "orange"),
               alpha = 0.50,
               cex=0.3,
               cat.cex=0.25)
}
