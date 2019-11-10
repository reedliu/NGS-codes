## step3-trim.sh

# conda create -n rna python=2 trim-galore -y
# conda activate rna
wkd=/home/liuyunze/mouse-c1

if [ ! -d $wkd/clean ]
then mkdir -p $wkd/clean
fi
cd $wkd/clean

ls $wkd/rawdata/*1.fastq.gz >1.raw 
ls $wkd/rawdata/*2.fastq.gz >2.raw
paste 1.raw 2.raw > conf.raw
cat conf.raw | while read i
do
	fqs=($i)
	fq1=${fqs[0]}
	fq2=${fqs[1]}
	sample=`basename $fq1 | cut -d '_' -f1`
	start=$(date +%s.%N)
	# echo trim_galore `date`
	trim_galore -q 25 --phred33 --length 50 -e 0.1 --stringency 3 \
		--paired -o $wkd/clean $fq1 $fq2 1>${sample}.trim.log 2>&1 && 
	# echo trim_galore `date`
	dur=$(echo "$(date +%s.%N) - $start" | bc)
	printf "Execution time for $sample trim_galore : %.6f seconds\n" $dur >>${sample}.trim.log
done 
# conda deactivate

cd $wkd/clean
fastqc *gz -t 10
multiqc ./
