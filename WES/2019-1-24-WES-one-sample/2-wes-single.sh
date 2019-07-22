#!/usr/bin/bash
###############################################
############# WES for one sample ##############
# Author: Reed Liu
# Mail: jieandze1314@gmail.com
# Date: 2019-01-24

## Usage: 
#cat $wkd/raw/conf | while read i;do
#	fqs=($i)
#	fq1=${fqs[0]}
#	fq2=${fqs[1]}
#	sh wes-single.sh $fq1 $fq2
#done
###############################################
## 针对一个样本PE fq文件

# 软件、工具路径
# 可以用conda，也可以自己指定

# Arguments
wkd=/vol2/agis/xiaoyutao_group/liuyunze/project/single-sap-wes
#mkdir $wkd/{raw,ref}
#rsync -av ~/reference/genome/hg38/ $wkd/ref # fast copy files(hg38.fa + hg38.fai)
ref_dir=$wkd/ref/hg38
ref_gnm=$ref_dir/hg38.fa
ref_idx=$ref_dir/hg38/hg38
GATK_bundle=/vol2/agis/xiaoyutao_group/liuyunze/biosoft/GATK/resources/bundle/hg38
conda_wes=/vol2/agis/xiaoyutao_group/liuyunze/biosoft/miniconda3/envs/wes

## shell执行参数
fq1=$1
fq2=$2
GID=Test_gid
library=Wes
sample=Test_sample
outdir=$wkd/outdir
core=16


## 按样本设置目录
outdir=${outdir}/${sample}

## 利用fq1得到fq的前缀，假设名字是*.1.fq.gz
fq_name=`basename $fq1`
fq_name=${fq_name%%.*}

## output dir
# 存储过滤好的数据
if [ ! -d $outdir/clean ]
	then mkdir -p $outdir/clean
fi
# 存储比对、索引直到BQSR前
if [ ! -d $outdir/bwa ]
	then mkdir -p $outdir/bwa
fi
# 存储变异检测结果
if [ ! -d $outdir/gatk ]
	then mkdir -p $outdir/gatk
fi
# 存储临时文件
if [ ! -d $outdir/tmp ]
	then mkdir -p $outdir/tmp
fi



###############################################
## Trimmomatic 质控过滤
###############################################
source activate wes
time trimmomatic PE \
	$fq1 $fq2 \
	$outdir/clean/${fq_name}.paired.1.fq.gz $outdir/tmp/${fq_name}.unpaired.1.fq.gz \
	$outdir/clean/${fq_name}.paired.2.fq.gz $outdir/tmp/${fq_name}.unpaired.2.fq.gz \
	ILLUMINACLIP:$conda_wes/share/trimmomatic-0.38-1/adapters/TruSeq3-PE-2.fa:2:30:10:8:True \
	SLIDINGWINDOW:5:15 LEADING:5 TRAILING:5 MINLEN:50 && echo "** fq QC done **"
source deactivate 


###############################################
## bwa mem 比对，对长度大于40bp小于2000bp的read非常有效
###############################################
## index
source activate wes
cd $ref_dir
bwa index -a bwtsw $ref_gnm -p hg38 && echo "** ref index done **"

## align

time bwa mem -t $core -M -R "@RG\tID:$GID\tSM:$sample\tLB:$library\tPL:Illumina" $ref_idx \
		$outdir/clean/${fq_name}.paired.1.fq.gz $outdir/clean/${fq_name}.paired.1.fq.gz |\
		samtools view -Sb - > $outdir/bwa/${sample}.bam && \
		echo "** BWA MEM done **"
time samtools sort -@ $core -m 100G -O bam -o $outdir/bwa/${sample}.sorted.bam $outdir/bwa/${sample}.bam && \
		echo "** sort raw bamfile done **"

# 可以对sorted后的bam进行index【可选】
# time samtools index $outdir/bwa/${sample}.sorted.bam && echo "** ${sample}.sorted.bam index done **"
source deactivate 

###############################################
## 标记重复序列
###############################################
source activate wes
gatk MarkDuplicates \
	-I $outdir/bwa/${sample}.sorted.bam \
	-M $outdir/bwa/${sample}.markup_metrics.txt \
	-O $outdir/bwa/${sample}.sorted.markdup.bam && \
	echo "** ${sample}.sorted.bam MarkDuplicates done **"

# 为${sample}.sorted.bam 构建index，后续需要
time samtools index $outdir/bwa/${sample}.sorted.markdup.bam && \
	echo "** ${sample}.sorted.markdup.bam index done **"
source deactivate 

###############################################
## BQSR
## 【注意】vcf文件需要有构建好的index
##  GATK4不支持对新的vcf进行index
###############################################
source activate wes
time gatk BaseRecalibrator \
	-R $ref_gnm \
	-I $outdir/bwa/${sample}.sorted.markdup.bam \
	--known-sites $GATK_bundle/1000G_phase1.indels.hg38.vcf \
	--known-sites $GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf \
	--known-sites $GATK_bundle/dbsnp_146.hg38.vcf \
	-O $outdir/bwa/${sample}.sorted.markdup.recal_data.table && \
	echo "** ${sample}.sorted.markdup.recal_data.table done **"

time gatk ApplyBQSR \
	--bqsr-recal-file $outdir/bwa/${sample}.sorted.markdup.recal_data.table \
	-R $ref_gnm \
	-I $outdir/bwa/${sample}.sorted.markdup.bam \
	-O $outdir/bwa/${sample}.sorted.markdup.BQSR.bam && \
	echo "** ApplyBQSR done **"

## 为${sample}.sorted.markdup.BQSR.bam构建索引，后续需要
time samtools index $outdir/bwa/${sample}.sorted.markdup.BQSR.bam && \
	echo "** ${sample}.sorted.markdup.BQSR.bam index done **"
source deactivate 


###############################################
## 开始变异检测
## 【注意】单个样本有四种检测方式（结果一样）
###############################################
## 第一种，直接调用HaplotypeCaller 输出样本vcf（较大文件比较慢）
## 这适合于单样本，或者那种固定样本数量的情况,也就是执行一次HaplotypeCaller之后就老死不相往来了
source activate wes
time gatk HaplotypeCaller \
	-R $ref_gnm \
	-I $outdir/bwa/${sample}.sorted.markdup.BQSR.bam \
	-O $outdir/gatk/${sample}.HC.vcf.gz && \
	echo "** ${sample}.HC.vcf.gz done **"
source deactivate 

## 第二种，输出每条染色体vcf，然后再合并结果。目的是提高速度，但仅仅是通过区分染色体提速
##  -L 参数，通过这个参数我们可以指定特定的染色体，写成chr1还是1需要看基因组fa格式
source activate wes
chrom=( chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 \
	chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM )
for i in ${chrom[@]};do
	time gatk HaplotypeCaller \
		-R $ref_gnm \
		-I $outdir/bwa/${sample}.sorted.markdup.BQSR.bam \
		-L $i \
		-O $outdir/gatk/${sample}.HC.${i}.vcf.gz && \
		echo "** ${sample}.HC.${i}.vcf.gz done **" &
done && wait
merge_vcfs=""
for i in ${chrom[@]};do
		merge_vcfs=${merge_vcfs}" -I $outdir/gatk/${sample}.HC.${i}.vcf.gz \\"\n
done && time gatk MegeVcfs ${merge_vcfs} -O $outdir/gatk/${sample}.HC.vcf.gz && \
	echo "** MergeVcfs done **"

source deactivate 

## 第三种，先输出样本的全gVCF，再进行GenotypeGCVFs。单样本非必需，多样本的标配，文件较大速度会慢
## gVCF全称是genome VCF，是每个样本用于变异检测的中间文件，格式类似于VCF，
## 它把joint-genotype过程中所需的所有信息都记录在这里面，文件无论是大小还是数据量都远远小于原来的BAM文件
## 一旦新增加样本也不需要再重新去读取所有人的BAM文件，只需为新样本生成一份gVCF，然后重新执行这个joint-genotype就行
source activate wes
time gatk HaplotypeCaller \
	--emit-ref-confidence GVCF \
	-R $ref_gnm \
	-I $outdir/bwa/${sample}.sorted.markdup.BQSR.bam \
	-O $outdir/gatk/${sample}.HC.g.vcf.gz && \
	echo "** GVCF ${sample}.HC.g.vcf.gz done **"

time gatk GenotypeGCVFs \
	-R $ref_gnm \
	-V $outdir/gatk/${sample}.HC.g.vcf.gz \
	-O $outdir/gatk/${sample}.HC.vcf.gz && \
	echo "** GVCF ${sample}.HC.vcf.gz done **"
source deactivate 

## 第四种，输出每个染色体的gvcf，然后对每个染色体单独进行GenotypeGVCFs。目的是提高速度，但仅仅是通过区分染色体提速
source activate wes
chrom=( chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 \
	chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM )
for i in ${chrom[@]};do
	time gatk HaplotypeCaller \
		--emit-ref-confidence GVCF \
		-R $ref_gnm \
		-I $outdir/bwa/${sample}.sorted.markdup.BQSR.bam \
		-L $i \
		-O $outdir/gatk/${sample}.HC.${i}.g.vcf.gz && \
	time gatk GenotypeGCVFs \
		-R $ref_gnm \
		-V $outdir/gatk/${sample}.HC.${i}.g.vcf.gz \
		-O $outdir/gatk/${sample}.HC.${i}.vcf.gz && \
		echo "** ${sample}.HC.${i}.vcf.gz done **" &
done && wait
merge_vcfs=""
for i in ${chrom[@]};do
	merge_vcfs=${merge_vcfs}" -I $outdir/gatk/${sample}.HC.${i}.vcf.gz \\"\n
done && time gatk MegeVcfs ${merge_vcfs} -O $outdir/gatk/${sample}.HC.vcf.gz && \
	echo "** MergeVcfs done **"

source deactivate 

###############################################
## 变异注释
## 使用VEP
###############################################
source activate wes
time VEP --fasta $ref_gnm \
 --vcf --merged --fork 10 --hgvs --force_overwrite --everything \
   --offline --dir_cache $outdir/tmp/ \
   -i $outdir/gatk/${sample}.HC.VQSR.vcf.gz \
   -o $outdir/gatk/${sample}.HC.VQSR.VEP.vcf.gz

source deactivate 













