mkdir ~/RNA-seq/quant/hisat2
cd ~/RNA-seq/quant/hisat2

#####################################
# Index
#####################################
source activate tair-rna
cat >index.sh
ref=~/RNA-seq/ref/Arabidopsis_thaliana.TAIR10.28.dna.genome.fa
hisat2-build -p 10 $ref ath

bash index.sh


#####################################
# Alignment
#####################################
cat >align.sh
index=~/RNA-seq/quant/hisat2/ath
conf=~/RNA-seq/filter/conf
cat $conf | while read i
do
	fqs=($i)
	fq1=${fqs[0]}
	fq2=${fqs[1]}
	sample=`basename $fq1 | cut -d '_' -f1`
	hisat2 -p 10 -x $index -1 $fq1 -2 $fq2  | samtools sort -O bam \
		-@ 10  -o ${sample}_hisat.bam
	samtools index -@ 10 -b ${sample}_hisat.bam 
done

nohup bash align.sh &



#####################################
# Count
#####################################
cat >quant.sh
gtf=~/RNA-seq/ref/Arabidopsis_thaliana.TAIR10.28.gtf
count=~/RNA-seq/count
featureCounts  -T 10 -p -t exon -g gene_name \
	-a $gtf -o $count/hisat2_counts.txt *.bam 
featureCounts  -T 10 -p -t exon -g gene_id \
	-a $gtf -o $count/hisat2_counts_id.txt *.bam 

bash quant.sh

source deactivate