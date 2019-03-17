##################################################
###### RNA-seq for species without reference #####
##################################################
## Created by Yunze(Reed) Liu 
## Date: 2019-03-17
## Email: jieandze1314@gmail.com
## Wechat Public: 生信星球

###############################################
## Trinity拼接 
###############################################

# 在正式分析前，先做一个样本表格samples.txt[【根据Trinity的要求：对于多个样本多个重复，这样方便识别】
# cond_A    cond_A_rep1    A_rep1_left.fq    A_rep1_right.fq
# cond_A    cond_A_rep2    A_rep2_left.fq    A_rep2_right.fq
# cond_B    cond_B_rep1    B_rep1_left.fq    B_rep1_right.fq
# cond_B    cond_B_rep2    B_rep2_left.fq    B_rep2_right.fq

wkd=/YOUR/WORK/PATH
mkdir $wkd/assembly && cd $wkd/assembly
echo -e "GSNO\nGSNO\nGSNO\nwt\nwt\nwt" >1
echo -e "GSNO_SRR1582646\nGSNO_SRR1582647\nGSNO_SRR1582648\nwt_SRR1582649\nwt_SRR1582650\nwt_SRR1582651" >2
ls $wkd/clean/test*1_val_1.fq >3
ls $wkd/clean/test*2_val_2.fq >4
paste 1 2 3 4 >samples.txt

Trinity --seqType fq  --samples_file $wkd/assembly/samples.txt \
      --CPU 10 --max_memory 10G --min_contig_length 150
# Running Trinity on this data set may take 10 to 15 minutes.
#  ~1 hour and ~1G RAM per ~1 million PE reads
# 运行过程：starting with Jellyfish to generate the k-mer catalog =》
# Inchworm to assemble 'draft' contigs =》
# Chrysalis to cluster the contigs and build de Bruijn graphs =》
# Butterfly for tracing paths through the graphs and reconstructing the final isoform sequences.

# conda安装samtools遇到动态库的问题=》最简单：降级，如最新版1.9，我们安装1.8
# conda install samtools=1.8 

# ~ 总体流程
#1. Trinity的第一个程序：seqtk-trinity =》 将fq转为fa
#2. 第二个程序：Jellyfish : jellyfish count + histo + dump => Normalization
#3. 第三个程序：Inchworm
#4. 第四个：Chrysalis ： inchworm_target -》 both.fa
#5. 第五个：Trinity Phase 2: Assembling Clusters of Reads =》ParaFly
#6. Trinity assemblies =》 Trinity.fasta

















