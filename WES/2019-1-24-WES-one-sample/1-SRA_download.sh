#!/usr/bin/bash
###############################################
############# WES for one sample ##############
###############################################
# Author: Reed Liu
# Mail: jieandze1314@gmail.com
# Date: 2019-01-24
###############################################
## Download 数据下载
###############################################
# Arguments
wkd=/vol2/agis/xiaoyutao_group/liuyunze/project/single-sap-wes
if [ ! -d $wkd/raw ]
	then mkdir -p $wkd/raw
fi

source activate wes

cd $wkd/raw
cat $wkd/raw/SRR_Acc_List.txt |while read i
do
# SRP or
ascp -v -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -k 1 -T \
	-l200m anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByStudy/sra/SRP/SRP072/$i/${i}.sra ./ && \
	echo "** ${i}.sra done **"
# SRR
ascp -v -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -k 1 -T \
	-l200m anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRR327/$i/${i}.sra ./ && \
	echo "** ${i}.sra done **"
time fastq-dump --gzip --split-3 -A $i ${i}.sra && echo "** ${i}.sra to fastq done **"
done 

# 对于NCBI没有的，或者下载失败的，可以去EBI搜索
# https://www.ebi.ac.uk/ena/data/view/SRR7722939
ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
	era-fasp@fasp.sra.ebi.ac.uk:vol1/srr/SRR772/009/SRR7722939 ./


source deactivate 

# ls *.sra | while read i;do (time fastq-dump --gzip --split-3 -A ${i%%.*} $i);done

ls $wkd/raw/*_1* >1
ls $wkd/raw/*_2* >2
paste 1 2 >conf
