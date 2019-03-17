### ---------------
###
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-3-3
### Email: jieandze1314@gmail.com
### Blog: https://reedliu.github.io/
### CAAS/AGIS/SDAU 
###
### ---------------

########################################
# 各个比较总差异基因
########################################
rm(list = ls())
options(stringsAsFactors = F)
load('salmon_edgeR_DEGs.Rdata')
DU <- c(DU_up,DU_down)
LFC2 <- c(LFC2_up,LFC2_down)
RS <- c(RS_up,RS_down)
XXCAD <- c(XXCAD_up,XXCAD_down)

########################################
# 共有转录本差异分析
########################################
up_gene <- Reduce(intersect, list(DU_up,LFC2_up,RS_up,XXCAD_up))
down_gene <- Reduce(intersect, list(DU_down,LFC2_down,RS_down,XXCAD_down))
all_gene <- c(up_gene,down_gene) #共有转录本

# 看下共有转录本表达差异
count <- read.csv("Salmon.gene.counts.matrix", sep = "\t", row.names = 1)
count <- apply(count, 1:2, round)
count <- as.data.frame(count)
countData <- count[apply(count, 1, sum) > 0 ,]
# 先看上调
up_matrix=countData[up_gene,]
up_matrix=t(scale(t(up_matrix)))
if(T){
  library(ComplexHeatmap)
  #install.packages("dendextend")
  library(dendextend)
  library(dendsort)
  e_mean = tail(sort(apply(countData,1,mean)),length(up_gene))
  #top30_cluster_cols = hclust(dist(t(top30_matrix)))
  sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...))) 
  top30_cluster_cols = sort_hclust(hclust(dist(t(up_matrix))))
  top30_cluster_rows <- sort_hclust(hclust(dist(up_matrix)))
  base_mean = e_mean
  ef_lable = HeatmapAnnotation(text = anno_text(colnames(countData), rot = 45,offset = unit(23, "mm"), gp = gpar(fontsize	= 8.5)))
  UP_DEGs <- Heatmap(up_matrix, name = "expresssion", 
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
  print(UP_DEGs)
}
# 再看下调
down_matrix=countData[down_gene,]
down_matrix=t(scale(t(down_matrix)))
if(T){
  library(ComplexHeatmap)
  #install.packages("dendextend")
  library(dendextend)
  library(dendsort)
  e_mean = tail(sort(apply(countData,1,mean)),length(down_gene))
  #top30_cluster_cols = hclust(dist(t(top30_matrix)))
  sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...))) 
  top30_cluster_cols = sort_hclust(hclust(dist(t(down_matrix))))
  top30_cluster_rows <- sort_hclust(hclust(dist(down_matrix)))
  base_mean = e_mean
  ef_lable = HeatmapAnnotation(text = anno_text(colnames(countData), rot = 45,offset = unit(23, "mm"), gp = gpar(fontsize	= 8.5)))
  DOWN_DEGs <- Heatmap(down_matrix, name = "expresssion", 
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
  print(DOWN_DEGs)
}
# 全部基因
all_matrix=countData[all_gene,]
all_matrix=t(scale(t(all_matrix)))
if(T){
  library(ComplexHeatmap)
  #install.packages("dendextend")
  library(dendextend)
  library(dendsort)
  e_mean = apply(countData[all_gene,],1,mean)
  #top30_cluster_cols = hclust(dist(t(top30_matrix)))
  sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...))) 
  top30_cluster_cols = sort_hclust(hclust(dist(t(all_matrix))))
  top30_cluster_rows <- sort_hclust(hclust(dist(all_matrix)))
  base_mean = e_mean
  ef_lable = HeatmapAnnotation(text = anno_text(colnames(countData), rot = 45,offset = unit(23, "mm"), gp = gpar(fontsize	= 8.5)))
  ALL_DEGs <- Heatmap(all_matrix, name = "expresssion", 
                       column_title = "Samples", 
                       column_title_gp = gpar(fontsize = 14, fontface = "bold"),
                       row_title = "Genes",
                       row_title_gp = gpar(fontsize = 14, fontface = "bold"),
                       row_names_gp = gpar(fontsize = 4),
                       row_names_side = "left",
                       cluster_rows = color_branches(top30_cluster_rows, k = 4),
                       cluster_columns = color_branches(top30_cluster_cols, k = 3),
                       show_column_names = F,
                       bottom_annotation = ef_lable, bottom_annotation_height = unit(0.5, "cm"))+
    Heatmap(base_mean, name = "base_mean", show_row_names = FALSE, width = unit(5, "mm"))
  print(ALL_DEGs)
}

########################################
# 转录本ID转换（Linux中利用longest_isoform.fasta比对结果）
########################################
write.table(up_gene,file="up_gene.txt",quote = F,row.names = F)
write.table(down_gene,file="down_gene.txt",quote = F,row.names = F)

write.table(DU,file = "DU_DEGs.txt",quote = F,row.names = F)
write.table(LFC2,file = "LFC2_DEGs.txt",quote = F,row.names = F)
write.table(RS,file = "RS_DEGs.txt",quote = F,row.names = F)
write.table(XXCAD,file = "XXCAD_DEGs.txt",quote = F,row.names = F)


# 根据转录本ID从longest_isoform.fasta中提取序列
# while read -r line; do awk -v pattern=$line -v RS=">" \
#  '$0 ~ pattern { printf(">%s", $0); }' longest_isoform.fasta;\
#   done < up_gene.txt > up_gene.fa

# 然后 blastx
# blastx -query $wkd/up_gene.fa -db nr \
#   -outfmt '7 qseqid sseqid sscinames scomnames pident length evalue bitscore' \
#   -out "up_fm7.blastn@nt" -evalue 1e-20 -num_threads 60 -num_alignments 3

#将结果开头的#行去掉
# sed '/# /d' up_fm7.blastx@nr >new_up_fm7.blastx@nr

# 提取genebank ID
# cat up_gene.txt | while read id;do if grep -q $id new_up_fm7.blastx@nr ;then (i=${id};grep ${i} new_up_fm7.blastx@nr |\
#    sort -k 1,1 -u);else echo "${id}=no";fi ;done | awk -F '\t' '{print $2}'| awk -F '|' '{print $3"|"$4"|"}'| grep ref

########################################
# 对于非模式生物
########################################
# 在线抓取
if(F){
  # 抓取 OrgDb
  require(AnnotationHub)
  hub <- AnnotationHub()
  query(hub, "Helicoverpa")
  Helicoverpa <- hub[['AH66951']]
  # 查看 OrgDb 里的Gene数量
  length(keys(Helicoverpa))
  # 查看 OrgDb 里的Gene ID类型
  columns(Helicoverpa)
  # 自己的数据需要转换ID，例如：
  require(clusterProfiler)
  bitr("ref|XP_021196409.1|", 'GID', c("GENENAME"), org.Harmigera.eg.db)
  
}
# 自己构建
if(F){
  library(tidyverse)
  library(stringr)
  library(KEGGREST)
  library(AnnotationForge)
  
  #' Title
  #'
  #' @param f_emapper_anno eggnog-mapper annotation result
  #' @param author Who is the creator of this package? like "xxx <xxx@xxx.xx>"
  #' @param tax_id The Taxonomy ID that represents your organism. (NCBI has a nice online browser for finding the one you need)
  #' @param genus Single string indicating the genus
  #' @param species Single string indicating the species
  #'
  #' @return OrgDb name
  #' @export
  #'
  #' @examples
  makeOrgPackageFromEmapper <- function(f_emapper_anno, 
                                        author, 
                                        tax_id = "0", 
                                        genus = "default", 
                                        species = "default") {
    
    # read emapper result
    emapper <- read_delim(f_emapper_anno,
                          "\t", escape_double = FALSE, trim_ws = TRUE)
    
    # extract gene name from emapper
    gene_info <- emapper %>%
      dplyr::select(GID = query_name, GENENAME = `eggNOG annot`) %>%
      na.omit()
    
    # extract go annotation from emapper
    gos <- emapper %>%
      dplyr::select(query_name, GO_terms) %>%
      na.omit()
    
    gene2go = data.frame(GID = character(),
                         GO = character(),
                         EVIDENCE = character())
    
    for (row in 1:nrow(gos)) {
      the_gid <- gos[row, "query_name"][[1]]
      the_gos <- str_split(gos[row,"GO_terms"], ",", simplify = FALSE)[[1]]
      
      df_temp <- data_frame(GID = rep(the_gid, length(the_gos)),
                            GO = the_gos,
                            EVIDENCE = rep("IEA", length(the_gos)))
      gene2go <- rbind(gene2go, df_temp)
    }
    
    # extract kegg pathway annotation from emapper
    gene2ko <- emapper %>%
      dplyr::select(GID = query_name, Ko = KEGG_KOs) %>%
      na.omit()
    
    load(file = "kegg_info.RData")
    gene2pathway <- gene2ko %>% left_join(ko2pathway, by = "Ko") %>% 
      dplyr::select(GID, Pathway) %>%
      na.omit()
    
    # make OrgDb
    makeOrgPackage(gene_info=gene_info,
                   go=gene2go,
                   ko=gene2ko,
                   pathway=gene2pathway,
                   # gene2pathway=gene2pathway,
                   version="0.0.2",
                   maintainer=author,
                   author=author,
                   outputDir = ".",
                   tax_id=tax_id,
                   genus=genus,
                   species=species,
                   goTable="go")
    
    my_orgdb <- str_c("org.", str_to_upper(str_sub(genus, 1, 1)) , species, ".eg.db", sep = "")
    return(my_orgdb)
  }
  
  my_orgdb <- makeOrgPackageFromEmapper("input/tmp1.annotations",
                                        "liuyunze <jieandze1314@gmail.com>",  
                                        tax_id = "29058", 
                                        genus = "Helicoverpa", 
                                        species = "armigera")
  
  
  
  
}
# 加载注释包
if(F){
  # install.packages('org.Harmigera.eg.db',repos = NULL, type="source")
  library(org.Harmigera.eg.db)
  # columns(org.Harmigera.eg.db)
  up_core_genes <- read.table('test_up_genes.txt')
  down_core_genes <- read.table('down_gene.txt')
  # for up_core_genes
  if(F){
    up.core.BP <- enrichGO(gene         = up_core_genes[,1],
                       OrgDb         = org.Harmigera.eg.db,
                       keyType       = 'GID',
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 1,
                       qvalueCutoff  = 1)
    dotplot(up.core.BP, showCategory=15)
    barplot(up.core.BP)
  }
  # for down_core_genes
  if(F){
    down_core_names <- bitr(down_core_genes[,1], 'GID', c("GENENAME"), org.Harmigera.eg.db)
    down.core.BP <- enrichGO(gene         = down_core_genes[,1],
                           OrgDb         = org.Harmigera.eg.db,
                           keyType       = 'GID',
                           ont           = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 1,
                           qvalueCutoff  = 1)
    dotplot(down.core.BP, showCategory=20)
    barplot(down.core.BP)
  }
  # for all DEGs
  if(F){
    DEGs <- read.table('XXCAD.txt')
    library(clusterProfiler)
    library(org.Harmigera.eg.db)
    DEG_names <- bitr(DEGs[,1], 'GID', c("GENENAME"), org.Harmigera.eg.db)
    # BP
    if(F){
      DEG.BP <- enrichGO(gene         = DEGs[,1],
                         OrgDb         = org.Harmigera.eg.db,
                         keyType       = 'GID',
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.8,
                         qvalueCutoff  = 0.8)
      dotplot(DEG.BP, showCategory=10)
      barplot(DEG.BP)
      write.csv(DEG.BP,"BP_output.csv",quote = F)
    }
    # MF
    if(F){
      DEG.MF <- enrichGO(gene         = DEGs[,1],
                         OrgDb         = org.Harmigera.eg.db,
                         keyType       = 'GID',
                         ont           = "MF",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.5,
                         qvalueCutoff  = 0.5)
      dotplot(DEG.MF, showCategory=20)
      barplot(DEG.MF)

      write.csv(DEG.MF,"MF_output.csv",quote = F)
    }
    # CC
    if(F){
      DEG.CC <- enrichGO(gene         = DEGs[,1],
                         OrgDb         = org.Harmigera.eg.db,
                         keyType       = 'GID',
                         ont           = "CC",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.5,
                         qvalueCutoff  = 0.5)
      dotplot(DEG.MF, showCategory=20)
      barplot(DEG.MF)
      write.csv(DEG.MF,"MF_output.csv",quote = F)
    }
    # 合并
    if(F){
      ego_CC <- DEG.CC
      ego_BP <- DEG.BP
      ego_MF <- DEG.MF
      library(tidyverse)
      df <- rbind(ego_CC@result,ego_BP@result,ego_MF@result)
      df$grouplist <- c(rep("CC",times=nrow(ego_CC@result)),
                        rep("BP",times=nrow(ego_BP@result)),
                        rep("MF",times=nrow(ego_MF@result)))
      df <- arrange(df,df$pvalue)
      df_s <- df[which(df$pvalue<0.05),]
      df_s$v <- -log10(df_s$pvalue)
      df_s$row <- 1:nrow(df_s)
      write.csv(df_s,file = '../RS_GOterms2Gene.csv',quote = F)
      out <- data_frame(ID=df_s$ID,Description=df_s$Description,
                        pvalue=df_s$pvalue,ontology=df_s$grouplist)
      write.csv(out,file="RS_GO.csv",quote = F)
      # 判断各个组分包含多少terms
      rownames(df_s)[df_s$grouplist == "BP"] %>% length()
      rownames(df_s)[df_s$grouplist == "CC"] %>% length()
      rownames(df_s)[df_s$grouplist == "MF"] %>% length()
      
      # 作图
      library("ggsci")
      library(yyplot)
      
      df_s <- df_s[1:20,]
      df_s <- df_s[order(df_s$grouplist),]
      rownames(df_s) <- 1:nrow(df_s)
      df_s$row <- nrow(df_s):1
      bar <- ggplot(data = df_s, mapping = aes(x = reorder(Description,row,order=T), 
                                               y = v,
                                               fill=grouplist))+ 
        geom_bar(stat = 'identity') + 
        #scale_fill_manual(values=c("#666666","#000000","#c0c0c0"))+
        coord_flip() +
        geom_text(aes(label=Count),hjust = -0.4, nudge_x = 0)+
        labs(y = "-log10(pvalue)", x = "GO Terms", 
             fill = "type") +
        theme_bw()+
        scale_fill_d3()
      g <- set_font(bar, size=3)
      plot(g)
    }
    
  }
  # 可选：合并相似GO term
  if(F){
    sim_down_core_BP <- simplify(down.core.BP, cutoff = 0.7, by = "p.adjust", select_fun = min)
  }

}

# cnetplot
# 先得到基因的ID与logFC值两列
ego <- up.core.BP
load('edgeR-salmon.Rdata')
all_core_gene <- salmon_edgeR_DU_SS
DEGs <- read.table('LFC2.txt')
DEG_names <- bitr(DEGs[,1], 'GID', c("GENENAME"), org.Harmigera.eg.db)
if(F){
  geneList <- all_core_gene$logFC
  names(geneList) <- DEG_names[,1]
  cnetplot(ego, 
           foldChange = geneList, 
           #foldChange = NULL, #不展示倍数
           circular = TRUE,
           #node_label = FALSE, #如果太多，就不要显示基因名了
           showCategory = 5, #显示富集的term数量，默认5
           colorEdge = TRUE)
  ggsave("DEG_circle.pdf", width = 8, height = 5)
}

# GO plot
ego <- read.csv("BP_output.csv", header = T)
ego[1,]
if(F){
  go <- data.frame(Category = "BP",
                   ID = ego$ID,
                   Term = ego$Description, 
                   Genes = gsub("/", ", ", ego$geneID), 
                   adj_pval = ego$p.adjust)
  all_gene <- read.table('DU_DEGs.txt')
  tmp <- salmon_edgeR_DU_SS[rownames(salmon_edgeR_DU_SS)%in%all_gene,]
  
  genelist <- data.frame(ID = down_core_names$GID, logFC = tmp$logFC)    
  circ <- circle_dat(go, genelist)
  head(circ)
  id.gsym <- bitr(gsub("REF","ref",circ$genes), 'GID', c("GENENAME"), org.Harmigera.eg.db)
  rownames(id.gsym) <- id.gsym$GID
  circ.gsym <- circ
  circ.gsym$genes <- id.gsym[gsub("REF","ref",circ$genes),]$GENENAME
  head(circ.gsym)
  
  #GOBar(subset(circ, category == 'BP'))
  #GOBubble(circ, labels = 3)
  #GOCircle(circ)
  
}

########################################
# KEGG
########################################
if(F){
  library(purrr)
  library(tidyverse)
  library(clusterProfiler)
  
  ########################################################################################
  # clusterProfiler 可能是目前最优秀的富集分析软件，参考网站：
  # https://www.bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html
  # 这里富集分析使用的是通用富集分析方法，参考：http://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html#universal-enrichment-analysis
  ########################################################################################
  
  ################################################
  # 导入自己构建的 OrgDb
  ################################################
  library(org.Harmigera.eg.db)
  columns(org.Harmigera.eg.db)
  
  ########################################################################################
  # 导入需要进行富集分析的基因列表，并转换为向量
  #########################################################################################
  
  gene_list <- DEGs[,1]
    
    ################################################
  # 从 OrgDB 提取 Pathway 和基因的对应关系
  ################################################
  
  pathway2gene <- AnnotationDbi::select(org.Harmigera.eg.db, 
                                        keys = keys(org.Harmigera.eg.db), 
                                        columns = c("Pathway","Ko")) %>%
    na.omit() %>%
    dplyr::select(Pathway, GID)
  
  ################################################
  # 导入 Pathway 与名称对应关系
  ################################################
  load("kegg_info.RData")
  
   #KEGG pathway 富集
  ekp <- enricher(gene_list, 
                  TERM2GENE = pathway2gene, 
                  TERM2NAME = pathway2name, 
                  pvalueCutoff = 1, 
                  qvalueCutoff = 1,
                  pAdjustMethod = "BH",
                  minGSSize = 1)
  
  ekp_results <- as.data.frame(ekp)
  
  barplot(ekp, showCategory=20,color="pvalue",
          font.size=10)
  dotplot(ekp)
  
  emapplot(ekp)
  
  
  
}





