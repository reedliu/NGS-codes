##############################
## install cran
##############################
if(length(getOption("CRAN"))==0) options(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")

if(!require("tidyr")) install.packages("tidyr",update = F,ask = F)

if(!require("dplyr")) install.packages("dplyr",update = F,ask = F)

if(!require("ggplot2")) install.packages("ggplot2",update = F,ask = F)

if(!require("devtools")) install.packages("devtools",update = F,ask = F)

if(!require("pheatmap")) install.packages("pheatmap",update = F,ask = F)

if(!require("ggfortify")) install.packages("ggfortify",update = F,ask = F)

if(!require("stringr")) install.packages("stringr",update = F,ask = F)

if(!require("ggpubr")) install.packages("ggpubr",update = F,ask = F)


##############################
## install bioconductor
##############################

# first prepare BioManager on CRAN
if(length(getOption("CRAN"))==0) options(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")

if(!require("BiocManager")) install.packages("BiocManager",update = F,ask = F)

if(length(getOption("BioC_mirror"))==0) options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")


# use BiocManager to install
if(!require("limma")) BiocManager::install("limma",update = F,ask = F)

if(!require("GO.db")) BiocManager::install("GO.db",update = F,ask = F)

if(!require("clusterProfiler")) BiocManager::install("clusterProfiler",update = F,ask = F)

if(!require("DESeq2")) BiocManager::install("DESeq2",update = F,ask = F)

if(!require("edgeR")) BiocManager::install("edgeR",update = F,ask = F)

# load multiple packages at one time
lapply(c('limma','DESeq2','edgeR'),function(x) suppressMessages(library(x,character.only = T)))
