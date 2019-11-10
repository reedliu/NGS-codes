## step5-salmon.sh

wkd=/home/liuyunze/mouse-c1

if [ ! -d $wkd/salmon ]
then mkdir -p $wkd/salmon
fi
cd $wkd/salmon

# salmon需要对转录本构建索引。因此只能使用参考转录组，而不能使用基因组
ref=$wkd/reference/gencode.vM23.transcripts.fa
salmon index -t $ref -i salmon.mm10

index=$wkd/salmon/salmon.mm10
conf=$wkd/clean/conf.raw
mkdir $wkd/salmon/quant && cd  $wkd/salmon/quant 
cat $conf | while read i
do
	fqs=($i)
	fq1=${fqs[0]}
	fq2=${fqs[1]}
	sample=`basename $fq1 | cut -d '_' -f1`
	start=$(date +%s.%N)
	salmon quant -i $index -l A --validateMappings \
		-1 $fq1 -2 $fq2 \
		-p 10 -o ${sample}_quant 1>${sample}.salmon.log 2>&1 && 
	dur=$(echo "$(date +%s.%N) - $start" | bc)
	printf "Execution time for $sample salmon : %.6f seconds\n" $dur >>${sample}.salmon.log
done

grep "Execution time" *log >salmon.summary


##################
# salmon数据整合
# in Rstudio
##################













