##################################################
###### RNA-seq for species without reference #####
##################################################
## Created by Yunze(Reed) Liu 
## Date: 2019-03-17
## Email: jieandze1314@gmail.com
## Wechat Public: 生信星球


###############################################
## 转录本质控
###############################################
wkd=/YOUR/WORK/PATH

##0. basic stats
###############################################
mkdir -p $wkd/evaluate/0-basic && cd $wkd/evaluate/0-basic
# 拼接了多少转录本
grep '>' $wkd/assembly/trinity_out_dir/Trinity.fasta | wc -l

# N50
TrinityStats.pl $wkd/assembly/trinity_out_dir/Trinity.fasta

##1. Bowtie2比对
###############################################
mkdir -p $wkd/evaluate/1-bowtie2 && cd $wkd/evaluate/1-bowtie2

bowtie2-build  $wkd/assembly/trinity_out_dir/Trinity.fasta $wkd/evaluate/1-bowtie2/Trinity.fasta

cat $wkd/assembly/samples.txt | while read i;do
	fqs=($i)
	fq1=${fqs[2]}
	fq2=${fqs[3]}
	file=`basename $fq1`
	surname=${file%%_*}
	#echo $fq1 $fq2 $surname
	bowtie2 -p 10 -q --no-unal -k 20 -x $wkd/evaluate/1-bowtie2/Trinity.fasta \
	        -1 $fq1 -2 $fq2  \
     		2>$wkd/evaluate/1-bowtie2/${surname}_align_stats.txt| \
     		samtools view -@ 10 -Sb -o $wkd/evaluate/1-bowtie2/${surname}.bowtie2.bam 
done
 
cat *.txt | grep overall > aln_stats.txt ##比对结果统计
# 例如：
#47.61% overall alignment rate
#46.46% overall alignment rate
#47.44% overall alignment rate
#38.82% overall alignment rate
#40.19% overall alignment rate
#36.25% overall alignment rate

# 可以将比对结果放入IGV中查看，但首先需要
ls $wkd/evaluate/1-bowtie2/*.bam | while read i; do
	file=`basename $i`
	#echo $file
	surname=${file%%.bowtie2.bam}
	#echo $surname
	samtools sort $i -o  $wkd/evaluate/1-bowtie2/${surname}.sorted.bam
	samtools index $wkd/evaluate/1-bowtie2/${surname}.sorted.bam
done

samtools faidx $wkd/assembly/trinity_out_dir/Trinity.fasta


## 2.uniprot比对
###############################################
mkdir -p $wkd/evaluate/2-uniprot && cd $wkd/evaluate/2-uniprot
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

cp $wkd/assembly/trinity_out_dir/Trinity.fasta ./
makeblastdb -in uniprot_sprot.fasta -dbtype prot

blastx -query Trinity.fasta -db uniprot_sprot.fasta -out blastx.outfmt6 \
        -evalue 1e-20 -num_threads 6 -outfmt 6 -max_target_seqs 1

analyze_blastPlus_topHit_coverage.pl blastx.outfmt6 Trinity.fasta uniprot_sprot.fasta

cat blastx.outfmt6.hist ##结果如下

#hit_pct_cov_bin	count_in_bin	>bin_below
#100	71	71
#90	22	93
#80	11	104
#70	11	115
#60	19	134
#50	22	156
#40	38	194
#30	25	219
#20	62	281
#10	14	295


# 计算比对上的序列种类及数量
cat blastx.outfmt6 | cut -f 2 | cut -d '|' -f 3 | \
	cut -d '_' -f 2| sort | uniq -c| sort -nr 


## 3.nr比对
###############################################
mkdir -p $wkd/evaluate/3-nr && cd $wkd/evaluate/3-nr

blastx -query $wkd/assembly/trinity_out_dir/Trinity.fasta \
-db /vol2/agis/xiaoyutao_group/liuyunze/ncbi/blast/db_nr/nr -out blastx.outfmt6 \
-evalue 1e-20 -num_threads 6 -outfmt 6 

analyze_blastPlus_topHit_coverage.pl blastx.outfmt6 \
	$wkd/assembly/trinity_out_dir/Trinity.fasta \
	/vol2/agis/xiaoyutao_group/liuyunze/ncbi/blast/db_nr/nr

###blast太慢，可以换用diamond，并且利用最长转录本进行比对（获得最长转录本见BUSCO）
axel http://mirrors.vbi.vt.edu/mirrors/ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz ##解压后110G
wget http://mirrors.vbi.vt.edu/mirrors/ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz.md5
gunzip nr.gz
mv nr nr.faa
diamond makedb --in nr.faa -d nr -p 20
diamond blastx -d nr -q longest_isoform.fasta -o isoforms_nr_matches

# 长度分布统计
analyze_blastPlus_topHit_coverage.pl isoforms_nr_matches \
	$wkd/evaluate/5-busco/longest_isoform.fasta \
	nr.faa >diamond-nr.stats


# 得到需要转换的Genbank ID
cat isoforms_nr_matches | cut -f2 | sort | uniq -c | sort -nr | awk '{if ($1>50) print $2}' > genebank.name

# 得到对应的物种名，然后就可以统计比对各个物种的占比
cat genebank.name | while read i;do esearch -db protein -query $i| elink -target taxonomy |\
	efetch -format xml |xtract -pattern TaxaSet -element ScientificName|cut -f1 ;done


## 4.pfam比对
###############################################
mkdir -p $wkd/evaluate/4-pfam && cd $wkd/evaluate/4-pfam
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.fasta.gz

makeblastdb -in Pfam-A.fasta -dbtype prot
blastx -query $wkd/assembly/trinity_out_dir/Trinity.fasta \
-db Pfam-A.fasta -out blastx.outfmt6 \
-evalue 1e-20 -num_threads 6 -outfmt 6 

# 利用diamond的比对结果
analyze_blastPlus_topHit_coverage.pl isoforms_pfam_matches \
	$wkd/assembly/trinity_out_dir/Trinity.fasta \
	Pfam-A.fasta >pfam.stats



## 5.BUSCO比对
###############################################
mkdir -p $wkd/evaluate/5-busco && cd $wkd/evaluate/5-busco
# 下载数据库
wget http://busco.ezlab.org/v2/datasets/fungi_odb9.tar.gz
# 提取每个基因的最长转录本
get_longest_isoform_seq_per_trinity_gene.pl $wkd/assembly/trinity_out_dir/Trinity.fasta \
	>longest_isoform.fasta
# 运行 BUSCO 
# -l/--lineage: 数据库的路径
# -m:运行模式，geno|tran|prot
run_BUSCO.py -i longest_isoform.fasta \
        -l insecta_odb9 \
        -o busco \
        -m tran \
        --cpu 10

# 看结果的summary文件，主要看C(complete)值，即为num%有同源保守基因



















