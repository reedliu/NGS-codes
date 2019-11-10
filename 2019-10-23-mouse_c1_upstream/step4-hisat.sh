## step4-hisat.sh

wkd=/home/liuyunze/mouse-c1

if [ ! -d $wkd/hisat ]
then mkdir -p $wkd/hisat
fi
cd $wkd/hisat

##################
# 准备index
##################
# 自己构建index
if [ ! -d $wkd/reference ]
then mkdir -p $wkd/reference
fi
cd $wkd/reference
# 从Gencode下载参考基因组
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.primary_assembly.genome.fa.gz
# 从Gencode下载参考转录组
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.transcripts.fa.gz

# 从UCSC下载参考基因组
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz
tar zxvf chromFa.tar.gz
for i in $(seq 1 19) X Y M;
do cat chr${i}.fa >> mm10.fasta
done
rm -fr chr*.fa

cd $wkd/hisat
ref=$wkd/reference/mm10.fasta
hisat2-build -p 40 $ref hisat.mm10
# 历时：00:14:39

# 下载hisat官网做好的index
# wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/mm10.tar.gz

##################
# 比对
##################
wkd=/home/liuyunze/mouse-c1
index=$wkd/hisat/hisat.mm10
conf=$wkd/clean/conf.raw
cat $conf | while read i
do
	fqs=($i)
	fq1=${fqs[0]}
	fq2=${fqs[1]}
	sample=`basename $fq1 | cut -d '_' -f1`
	start=$(date +%s.%N)
	hisat2 -p 10 -x $index -1 $fq1 -2 $fq2  --new-summary -S ${sample}_hisat.sam 1>${sample}.hisat.log 2>&1
	samtools sort -O bam -@ 10  -o ${sample}_hisat.bam ${sample}_hisat.sam 
	rm ${sample}_hisat.sam 
	samtools index -@ 10 -b ${sample}_hisat.bam  && 
	dur=$(echo "$(date +%s.%N) - $start" | bc)
	printf "Execution time for $sample hisat and sam2bam : %.6f seconds\n" $dur >>${sample}.hisat.log

done

grep "Overall alignment" * >hisat.summary
grep "Execution time" * >>hisat.summary


##################
# featureCount定量
##################
if [ ! -d $wkd/hisat/count ]
then mkdir -p $wkd/hisat/count
fi
cd $wkd/hisat/count

# 下载GTF文件
cd $wkd/reference
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.annotation.gtf.gz
gunzip gencode.vM23.annotation.gtf.gz

cd $wkd/hisat/count
gtf=$wkd/reference/gencode.vM23.annotation.gtf
featureCounts  -T 10 -p -t exon -g gene_name -a $gtf -o $wkd/hisat/count/hisat2_counts.txt $wkd/hisat/*.bam 1>featureCounts.log 2>&1


##################
# qualimap对转录组bam质控
##################
wkd=/home/liuyunze/mouse-c1
mkdir $wkd/hisat/qualimap && cd $wkd/hisat/qualimap
gtf=$wkd/reference/gencode.vM23.annotation.gtf
ls $wkd/hisat/*.bam> bam_path.txt
cat bam_path.txt | while read i
do
	sample=`basename $i | cut -d '_' -f1`
	qualimap rnaseq --java-mem-size=20G -gtf $gtf -bam $i -pe  -oc $sample 
done
# 结果会保存在bam所在目录下，每个bam一个质控结果

cd $wkd/hisat 
multiqc ./
















