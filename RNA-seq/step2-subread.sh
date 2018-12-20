mkdir ~/RNA-seq/quant/subread
cd ~/RNA-seq/quant/subread

#####################################
# Index
#####################################
source activate tair-rna
cat >index.sh
ref=~/RNA-seq/ref
subread-buildindex -o subread_ath_index \
  $ref/Arabidopsis_thaliana.TAIR10.28.dna.genome.fa

bash index.sh

#####################################
# Alignment
#####################################
cat >align.sh
index=~/RNA-seq/quant/subread/subread_ath_index
conf=~/RNA-seq/filter/conf
cat $conf | while read i
do
	fqs=($i)
	fq1=${fqs[0]}
	fq2=${fqs[1]}
	sample=`basename $fq1 | cut -d '_' -f1`
	subjunc -i $index -r $fq1 -R $fq2 \
		-T 10 -o ./${sample}_subjunc.bam
	samtools index -@ 10 -b ${sample}_subjunc.bam
done

bash align.sh

#####################################
# Count
#####################################

mkdir ~/RNA-seq/count
cat >quant.sh
gtf=~/RNA-seq/ref/Arabidopsis_thaliana.TAIR10.28.gtf
count=~/RNA-seq/count
featureCounts  -T 10 -p -t exon -g gene_name \
	-a $gtf -o $count/subread_counts.txt *.bam 
featureCounts  -T 10 -p -t exon -g gene_id \
	-a $gtf -o $count/subread_counts_id.txt *.bam 

bash quant.sh
source deactivate
