#!/bin/bash

### ---------------
###
### Create: Yunze Liu (Reed Liu)
### Date: 2018-12-19
### CAAS/AGIS/SDAU
###
###
### ---------------

############
# qc
############
mkdir -p ~/RNA-seq/qc
cd ~/RNA-seq/qc
source activate tair-rna
fastqc ~/RNA-seq/raw/*.gz -o ./ -t 10
multiqc ./
source deactivate

############
# filter
############
# make index for raw reads
mkdir ~/RNA-seq/filter && cd ~/RNA-seq/filter
ls ~/RNA-seq/raw/*1.fastq.gz >1 
ls ~/RNA-seq/raw/*2.fastq.gz >2
paste 1 2 > conf

# first method: trim_galore
mkdir -p ~/RNA-seq/filter/trim-galore
cd -p ~/RNA-seq/filter/trim-galore

source activate tair-rna
cat >filter.sh
dir=~/RNA-seq/filter/trim-galore
conf=~/RNA-seq/filter/conf
cat $conf | while read i
do
	fqs=($i)
	fq1=${fqs[0]}
	fq2=${fqs[1]}
	nohup trim_galore -q 25 --phred33 --length 50 -e 0.1 --stringency 3 --paired -o $dir $fq1 $fq2 &
done

bash filter.sh
source deactivate

# second method: trimmomatic
mkdir -p ~/RNA-seq/filter/trimmomatic
cd ~/RNA-seq/filter/trimmomatic

source activate tair-rna
cat >filter-2.sh
conf=~/RNA-seq/filter/conf
cat $conf | while read i
do
	fqs=($i)
	fq1=${fqs[0]}
	fq2=${fqs[1]}
	tri1=`basename $fq1`
	tri2=`basename $fq1`
	nohup trimmomatic PE -phred33 \
	-trimlog trim.logfile\
	$fq1 $fq2 \
	clean.$tri1 unpaired.$tri1 \
	clean.$tri2 unpaired.$tri2 \
	SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:20 &
done

bash filter-2.sh
source deactivate






