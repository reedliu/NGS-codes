# This is a normal pipeline

########################################
# get up_gene & down_gene
########################################
#DEG = Deseq_DEG/edgeR_DEG/limma_voom_DEG
if(T){
  logFC_cutoff <- with(DEG,mean(abs(logFC)) + 2*sd(abs(logFC)) ) 
  DEG$result = as.factor(ifelse(DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                                ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT'))
}
gene <- rownames(DEG[DEG$result != 'NOT', ])
up_gene <- rownames(DEG[DEG$result =='UP',])
down_gene <- rownames(DEG[DEG$result =='DOWN',])

########################################
# trans ID
########################################
library(clusterProfiler) 
library(org.Hs.eg.db)
if(T){
  gene_tr <- bitr(gene, fromType = "SYMBOL",
                  toType = c("ENSEMBL", "ENTREZID"),
                  OrgDb = org.Hs.eg.db)
  up_gene_tr <- bitr(up_gene, fromType = "SYMBOL",
                     toType = c("ENSEMBL", "ENTREZID"),
                     OrgDb = org.Hs.eg.db)
  down_gene_tr <- bitr(down_gene, fromType = "SYMBOL",
                       toType = c("ENSEMBL", "ENTREZID"),
                       OrgDb = org.Hs.eg.db)
  gene <- gene_tr$ENTREZID
  up_gene <- up_gene_tr$ENTREZID
  down_gene <- down_gene_tr$ENTREZID
} 
head(gene_tr)

# CHECK:  e.g. Entrez ID: 23336 => https://www.ncbi.nlm.nih.gov/gene/23336

########################################
# KEGG
########################################
if(T){
  kk <- enrichKEGG(gene         = gene,
                   organism     = 'hsa',
                   pvalueCutoff = 0.05)
  kk_gene <- (kk)[,1:6]
}
dotplot(kk,showCategory = 14,color="pvalue",font.size=14)  

########################################
# gseKEGG GSEA
########################################
if(T){
  geneList <- DEG_mtx$logFC 
  names(geneList) <- rownames(DEG_mtx) 
  geneList_tr <- bitr(names(geneList), 
                      fromType = "SYMBOL",
                      toType = c("ENSEMBL","ENTREZID"),
                      OrgDb = org.Hs.eg.db) 
  
  new_list <- data.frame(SYMBOL=names(geneList), logFC = as.numeric(geneList)) 
  new_list <- merge(new_list, geneList_tr, by  = "SYMBOL")
  
  geneList <- new_list$logFC 
  names(geneList) <- geneList_tr$ENTREZID 
  
  geneList <- sort(geneList,decreasing = T) 
}

if(T){
  kk2 <- gseKEGG(geneList     = geneList,
                 organism     = 'hsa',
                 nPerm        = 1000,
                 minGSSize    = 120,
                 pvalueCutoff = 0.2,
                 verbose      = FALSE)
}
# CHECK
head(kk2)[,1:6]
gseaplot(kk2, geneSetID = "hsa04310")

########################################
# pathway
########################################
library(pathview)
hsa01212 <- pathview(gene.data  = geneList,
                     pathway.id = "hsa01212", 
                     species    = "hsa",
                     limit      = list(gene=max(abs(geneList)), cpd=1))

########################################
# GO
########################################
if(T){
  ego_CC <- enrichGO(gene          = gene,
                     universe      = names(geneList),
                     OrgDb         = org.Hs.eg.db,
                     ont           = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)
}
ego_gene <- head(ego_CC)[,1:6]

dotplot(ego_CC,title="EnrichmentGO_CC_dot")

barplot(ego_CC, showCategory=20,title="EnrichmentGO_CC")

library(Rgraphviz)
library(topGO)
plotGOgraph(ego_CC)





