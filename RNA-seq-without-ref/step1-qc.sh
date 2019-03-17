##################################################
###### RNA-seq for species without reference #####
##################################################
## Created by Yunze(Reed) Liu 
## Date: 2019-03-17
## Email: jieandze1314@gmail.com
## Wechat Public: 生信星球


##################################################
Step1: QC and filter
##################################################
wkd=/YOUR/WORK/PATH
mkdir $wkd/qc && cd $wkd/qc
fastqc -o ./ -t 10 $wkd/raw/*.fastq
multiqc ./

mkdir $wkd/clean && cd $wkd/clean
ls $wkd/raw/test*_1.fastq >1 
ls $wkd/raw/test*_2.fastq >2 
paste 1 2 > filter_conf
cat filter_conf | while read i;do
	fqs=($i)
	fq1=${fqs[0]}
	fq2=${fqs[1]}
	trim_galore -q 25 --phred33 --length 50 -e 0.1 --stringency 3 \
		--paired -o $wkd/clean $fq1 $fq2
done

# 再一次质控
cd $wkd/clean
source activate rna-seq
fastqc -o ./ -t 10 $wkd/clean/*.fq
multiqc ./
source deactivate




