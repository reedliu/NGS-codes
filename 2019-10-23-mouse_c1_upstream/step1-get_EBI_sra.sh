wkd=/home/liuyunze/mouse-c1

if [ ! -d $wkd/sra ]
then mkdir -p $wkd/sra
fi
cd $wkd/sra

cat XX_SRR_Acc_List.txt  |while read i
do
	# echo fasp.sra.ebi.ac.uk:vol1/srr/SRR781/000/$i
# These SRAs in EBI not NCBI! ex.https://www.ebi.ac.uk/ena/browser/view/SRR7815790
	ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/srr/SRR781/000/$i ./ && \
	echo "** ${i}.sra done **";
done


# 需要注意的是：EBI的数据结构是这样的：
# 例如要下载SRR7815980.sra，那么它就是ftp://ftp.sra.ebi.ac.uk/vol1/srr/SRR781/000/SRR7815980；例如要下载SRR7815792.sra，它就是ftp://ftp.sra.ebi.ac.uk/vol1/srr/SRR781/002/SRR7815792。
# 000、002随着SRA ID的结尾数字发生变化

############################
# 【附】如果要下载NCBI SRA
############################
# cat $wkd/raw/SRR_Acc_List.txt |while read i
# do
# # 针对 SRR
# ascp -v -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -k 1 -T \
# 	-l200m anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRR327/$i/${i}.sra ./ && \
# 	echo "** ${i}.sra done **"
# time fastq-dump --gzip --split-3 -A $i ${i}.sra && echo "** ${i}.sra to fastq done **"
# done 
