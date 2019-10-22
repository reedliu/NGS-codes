#######################################
######### RNA-seq from SRA ############
# Reedliu created on 2019-01-17
# Email: jieandze1314@gmail.com
#######################################
#PBS -N rna-pipe
#PBS -o /vol2/agis/xiaoyutao_group/liuyunze/project/RNA-seq/public-ha-rna/src/rna-pipe.out
#PBS -e /vol2/agis/xiaoyutao_group/liuyunze/project/RNA-seq/public-ha-rna/src/rna-pipe.err
#PBS -l nodes=1:ppn=16

## Hard Code Params
wkd=/vol2/agis/xiaoyutao_group/liuyunze/project/RNA-seq/public-ha-rna

#########################################################################
##################### STEP1：Download data ##############################
#########################################################################
# raw data => https://www.ncbi.nlm.nih.gov/sra?term=SRP179439
# => Send to => Run selector => Download runinfo table
if [ ! -d $wkd/raw-sra ]
then mkdir -p $wkd/raw-sra
fi
less -SN $wkd/raw-sra/SraRunTable.txt | awk '{print $6}' | sed 1d > $wkd/raw-sra/sra.list
source activate rna-seq
cat $wkd/raw-sra/sra.list | while read id; do 
time ascp -v -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
	-k 1 -T -l160m anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRR844/$id/${id}.sra $wkd/raw-sra && \
	echo "** downloaded $id **"
done
source deactivate

#reference genome and gff
if [ ! -d $wkd/ref ]
then mkdir -p $wkd/ref
fi
cd $wkd/ref
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/Helicoverpa_armigera/CHR_Un/29058_ref_Harm_1.0_chrUn.fa.gz
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/Helicoverpa_armigera/GFF/ref_Harm_1.0_top_level.gff3.gz
gunzip *
mv *.fa ha.fa
mv *.gff ha.gff

#########################################################################
##################### STEP2:Pre-process ## ##############################
#########################################################################
## sra2fq
if [ ! -d $wkd/raw-fq ]
then mkdir -p $wkd/raw-fq
fi
cd $wkd/raw-fq
source active rna-seq
time fastq-dump --split-files --gzip \
	-O ./ $wkd/raw-sra/*.sra && echo "** sra2fq done **"
source deactivate

## qc-before
if [ ! -d $wkd/qc-before ]
then mkdir -p $wkd/qc-before
fi
cd $wkd/qc-before
source active rna-seq
time fastqc $wkd/raw-fq/*.gz -o ./ -t 10 && echo "** fq QC done **"
source deactivate

## filter
if [ ! -d $wkd/clean ]
then mkdir -p $wkd/clean
fi
cd $wkd/clean
source activate rna-seq
ls $wkd/raw-fq/*1.fastq.gz >1 
ls $wkd/raw-fq/*2.fastq.gz >2
paste 1 2 > conf
cat conf | while read i
do
	fqs=($i)
	fq1=${fqs[0]}
	fq2=${fqs[1]}
	time trim_galore -q 25 --phred33 --length 50 -e 0.1 --stringency 3 \
		--paired -o ./ $fq1 $fq2 && echo "** fq filter done **"
done
source deactivate

## qc-after
if [ ! -d $wkd/qc-after ]
then mkdir -p $wkd/qc-after
fi
cd $wkd/qc-after
source active rna-seq
time fastqc $wkd/clean/*.gz -o ./ -t 10 && echo "** filtered fq QC done **"
source deactivate

#########################################################################
##################### STEP3:align & quant  ##############################
#########################################################################
## SUBREAD
########################
if [ ! -d $wkd/quant/subread ]
then mkdir -p $wkd/quant/subread
fi
cd $wkd/quant/subread
source activate rna-seq
# index
time subread-buildindex -o subread_ha_index $wkd/ref/ha.fa && \
	echo "** subread ref index done **"
# align
cat $wkd/clean/conf | while read i
do
	fqs=($i)
	fq1=${fqs[0]}
	fq2=${fqs[1]}
	sample=`basename $fq1 | cut -d '_' -f1`
	time subjunc -i subread_ha_index -r $fq1 -R $fq2 -T 16 | \
		samtools sort -O bam -@ 16  -o ${sample}_subjunc.sorted.bam && \
		echo "** align & sort ${sample}_subjunc.bam done **"
	time samtools index -@ 16 -b ${sample}_subjunc.sorted.bam && \
		echo "** index ${sample}_subjunc.bam done **"
done
# quant
if [ ! -d $wkd/quant/subread/count ]
then mkdir -p $wkd/quant/subread/count
fi
subread_count=$wkd/quant/subread/count
time featureCounts  -T 16 -p -t exon -g gene_name -a $wkd/ref/ha.gtf \
	 -o $subread_count/subread_counts.txt *.bam && \
	echo "** subread-counts genename done **"
time featureCounts  -T 16 -p -t exon -g gene_id -a $wkd/ref/ha.gtf \
	 -o $subread_count/subread_counts_id.txt *.bam && \
	echo "** subread-counts geneid done **"

source deactivate


## HISAT2
########################
if [ ! -d $wkd/quant/hisat2 ]
then mkdir -p $wkd/quant/hisat2
fi
cd $wkd/quant/hisat2
source activate rna-seq
# index
time hisat2-build -p 16 $wkd/ref/ha.fa hisat2-ha-index && \
	echo "** hisat2 ref index done **"
# align
cat $wkd/clean/conf | while read i
do
	fqs=($i)
	fq1=${fqs[0]}
	fq2=${fqs[1]}
	sample=`basename $fq1 | cut -d '_' -f1`
	time hisat2 -p 16 -x hisat2-ha-index -1 $fq1 -2 $fq2  | \
		samtools sort -O bam -@ 16  -o ${sample}_hisat.sorted.bam && \
		echo "** align & sort ${sample}_hisat.bam done"
	time samtools index -@ 16 -b ${sample}_hisat.sorted.bam \
		echo "** index ${sample}_hisat.sorted.bam done **"
done
# quant
if [ ! -d $wkd/quant/hisat2/count ]
then mkdir -p $wkd/quant/hisat2/count
fi
hisat_count=$wkd/quant/hisat2/count
time featureCounts  -T 16 -p -t exon -g gene_name -a $wkd/ref/ha.gtf \
	 -o $hisat_count/hisat2_counts.txt *.bam && \
	echo "** hisat2-counts genename done **"
time featureCounts  -T 16 -p -t exon -g gene_id -a $wkd/ref/ha.gtf \
	 -o $hisat_count/hisat2_counts_id.txt *.bam && \
	echo "** hisat2-counts geneid done **"

source deactivate

## BOWTIE2
########################
if [ ! -d $wkd/quant/bowtie2 ]
then mkdir -p $wkd/quant/bowtie2
fi
cd $wkd/quant/bowtie2
source activate rna-seq
# index
time bowtie2-build --threads 16 $wkd/ref/ha.fa ha \
	echo "** bowtie2 ref index done **"
# align
cat $wkd/clean/conf | while read i
do
	fqs=($i)
	fq1=${fqs[0]}
	fq2=${fqs[1]}
	sample=`basename $fq1 | cut -d '_' -f1`
	time bowtie2 -p 16 -x $index -1 $fq1 -2 $fq2 | \
		samtools sort -O bam -@ 16 -o ${sample}_bowtie.sorted.bam && \
		echo "** align & sort ${sample}_bowtie.bam done"
	time samtools index -@ 16 -b ${sample}_bowtie.sorted.bam \
		echo "** index ${sample}_bowtie.sorted.bam done **"
done
# quant
if [ ! -d $wkd/quant/bowtie2/count ]
then mkdir -p $wkd/quant/bowtie2/count
fi
bowtie_count=$wkd/quant/bowtie2/count
time featureCounts  -T 16 -p -t exon -g gene_name -a $wkd/ref/ha.gtf \
	 -o $bowtie_count/bowtie2_counts.txt *.bam && \
	echo "** bowtie2-counts genename done **"
time featureCounts  -T 16 -p -t exon -g gene_id -a $wkd/ref/ha.gtf \
	 -o $bowtie_count/bowtie2_counts_id.txt *.bam && \
	echo "** bowtie2-counts geneid done **"
source deactivate
















