##################################################
###### RNA-seq for species without reference #####
##################################################
## Created by Yunze(Reed) Liu 
## Date: 2019-03-17
## Email: jieandze1314@gmail.com
## Wechat Public: 生信星球

###############################################
## 定量
###############################################

## 1.salmon
###############################################
wkd=/YOUR/WORK/PATH

mkdir -p $wkd/quant/salmon && cd $wkd/quant/salmon
align_and_estimate_abundance.pl --seqType fq  \
    --samples_file $wkd/assembly/samples.txt \
    --transcripts $wkd/assembly/trinity_out_dir/Trinity.fasta \
    --est_method salmon \
    --trinity_mode \
    --prep_reference \
    --output_dir $wkd/quant/salmon \
    --thread_count 10


## 2.rsem
###############################################
mkdir -p $wkd/quant/rsem && cd $wkd/quant/rsem
align_and_estimate_abundance.pl --seqType fq  \
    --samples_file $wkd/assembly/samples.txt \
    --transcripts $wkd/assembly/trinity_out_dir/Trinity.fasta \
    --est_method RSEM \
    --aln_method bowtie2 \
    --trinity_mode \
    --prep_reference \
    --output_dir $wkd/quant/rsem \
    --thread_count 10

## 3.kallisto
###############################################
mkdir -p $wkd/quant/kallisto && cd $wkd/quant/kallisto
align_and_estimate_abundance.pl --seqType fq  \
    --samples_file $wkd/assembly/samples.txt \
    --transcripts $wkd/assembly/trinity_out_dir/Trinity.fasta \
    --est_method kallisto \
    --trinity_mode \
    --prep_reference \
    --output_dir $wkd/quant/kallisto \
    --thread_count 10

###############################################
## 定量结果组合
###############################################
cd $wkd/quant/salmon
find $wkd/quant/salmon -name "*quant.sf.genes" >genes.quant_files.txt
abundance_estimates_to_matrix.pl 
	--est_method salmon \
    --out_prefix salmon-gene \
    --name_sample_by_basedir \
    --gene_trans_map none \
	--quant_files  $wkd/quant/salmon/genes.quant_files.txt

cd $wkd/quant/rsem
find $wkd/quant/rsem -name "RSEM.isoforms.results" >isoforms.quant_files.txt
abundance_estimates_to_matrix.pl 
	--est_method RSEM \
    --out_prefix RSEM-isoform \
    --name_sample_by_basedir \
    --gene_trans_map  $wkd/assembly/trinity_out_dir/Trinity.fasta.gene_trans_map \
	--quant_files  $wkd/quant/rsem/isoforms.quant_files.txt

###############################################
## ExN50
###############################################
mkdir -p $wkd/evaluate/0-basic && cd $wkd/evaluate/0-basic
contig_ExN50_statistic.pl \
	$wkd/quant/rsem/RSEM-isoform.isoform.TMM.EXPR.matrix \
	$wkd/assembly/trinity_out_dir/Trinity.fasta >ExN50.stats
plot_ExN50_statistic.Rscript  ExN50.stats
















