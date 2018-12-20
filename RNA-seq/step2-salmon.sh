mkdir ~/RNA-seq/quant/salmon
cd ~/RNA-seq/quant/salmon

############################
# Index
############################
source activate tair-rna
cat >index.sh
ref=~/RNA-seq/ref
salmon index -t $ref/Arabidopsis_thaliana.TAIR10.28.cdna.all.fa -i salmon_ath_index

bash index.sh

############################
# Quant
# -l means library type and A means Antomatic
############################
cat >quant.sh
index=~/RNA-seq/quant/salmon/salmon_ath_index
conf=~/RNA-seq/filter/conf
cat $conf | while read i
do
	fqs=($i)
	fq1=${fqs[0]}
	fq2=${fqs[1]}
	sample=`basename $fq1 | cut -d '_' -f1`
	salmon quant -i $index -l A \
		-1 $fq1 -2 $fq2 \
		-p 10 -o ./${sample}_quant
done

bash quant.sh
source deactivate
