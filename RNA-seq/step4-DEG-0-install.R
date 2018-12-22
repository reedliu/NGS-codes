##############################
## install cran
##############################


if(length(getOption("CRAN"))==0) options(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")

for (pkg in c("tidyr","dplyr","ggplot2","devtools",'pheatmap')){
  if (! require(pkg,character.only=T) ) {
    install.packages(pkg,ask = F,update = F)
    require(pkg,character.only=T) 
    }
}

##############################
## install bioconductor
##############################

# first prepare BioManager on CRAN
if(length(getOption("CRAN"))==0) options(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")

if(!require("BiocManager")) install.packages("BiocManager",update = F,ask = F)

if(length(getOption("BioC_mirror"))==0) options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")


# use BiocManager to install
for (pkg in c("limma","GO.db","clusterProfiler","DESeq2",'edgeR')){
  if (! require(pkg,character.only=T) ) {
    BiocManager::install(pkg,ask = F,update = F)
    require(pkg,character.only=T) 
    }
}
