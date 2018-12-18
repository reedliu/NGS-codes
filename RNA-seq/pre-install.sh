#!/bin/bash

### ---------------
###
### Create: Yunze Liu
### Date: 2018-12-18
### CAAS/AGIS/SDAU
### 
###
### ---------------

conda  create -n tair-rna python=2
source activate tair-rna
conda install -y fastqc multiqc trim-galore trimmomatic salmon \
	subread hisat2 bowtie2


mkdir -p ~/RNA-seq/raw 
cd ~/RNA-seq/raw 

################################################
# Download sample info

# Data from "Temporal dynamics of gene expression and histone marks at the Arabidopsis shoot meristem during flowerin”
################################################
wget http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-5130/E-MTAB-5130.sdrf.txt

tail -n +2 E-MTAB-5130.sdrf.txt | awk '{print $(NF-6)}' >id

cat id | while read id; do wget -c $id; done

################################################
# Download reference data(genome/annotation)
################################################

mkdir -p ~/RNA-seq/ref
cd ~/RNA-seq/ref

axel ftp://ftp.ensemblgenomes.org/pub/plants/release-28/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.28.dna.genome.fa.gz 
axel ftp://ftp.ensemblgenomes.org/pub/plants/release-28/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.28.gff3.gz
axel ftp://ftp.ensemblgenomes.org/pub/plants/release-28/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.28.gtf.gz 











