##############################
## install cran
##############################


if(length(getOption("CRAN"))==0) options(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
cran_packages <- c('tidyverse',
                   'ggpubr',
                   'ggstatsplot')
for (pkg in cran_packages){
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
Biocductor_packages <- c('org.Hs.eg.db',
                         'hgu133a.db',
                         'CLL',
                         'hgu95av2.db',
                         'survminer',
                         'survival',
                         'hugene10sttranscriptcluster',
                         'limma')
for (pkg in Biocductor_packages){
  if (! require(pkg,character.only=T) ) {
    BiocManager::install(pkg,ask = F,update = F)
    require(pkg,character.only=T) 
    }
}
