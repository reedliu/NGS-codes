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
# limma
######################################
suppressPackageStartupMessages(library(limma))

countData <- count[apply(count, 1, sum) > 0 ,]
group_list <- rep(c("day0","day1","day2","day3"),each=4)

if(F){
    design <- model.matrix(~0+factor(group_list))
    colnames(design)=levels(factor(group_list))
    rownames(design)=colnames(countData)
    
    dge <- DGEList(counts=countData)
    dge <- calcNormFactors(dge)
    logCPM <- cpm(dge, log=TRUE, prior.count=3)
    
    v <- voom(dge,design,plot=TRUE, normalize="quantile")
    fit <- lmFit(v, design)
    
    group_list
    cont.matrix <- makeContrasts(contrasts=c("day0-day1", "day0-day2", "day0-day3"),levels = design) 
    fit2=contrasts.fit(fit,cont.matrix)
    fit2=eBayes(fit2)
    
    tempOutput = topTable(fit2, coef=c("day0-day1", "day0-day2", "day0-day3"), n=Inf)
    limma_voom_DEG = na.omit(tempOutput)
    
}

##################################################################
save(Deseq_DEG, edgeR_DEG, limma_voom_DEG, count,group_list, file = 'DEG_result.Rdata')

