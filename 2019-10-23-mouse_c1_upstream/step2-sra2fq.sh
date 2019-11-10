wkd=/home/liuyunze/mouse-c1

if [ ! -d $wkd/rawdata ]
then mkdir -p $wkd/rawdata
fi
cd $wkd/rawdata

# 如果是对于NCBI sra文件
# ls *.sra | while read i;do (fastq-dump --gzip --split-3 -A ${i%%.*} $i);done

# 这里的EBI文件本身就没有后缀
ls $wkd/sra/SRR* | while read i;do \
	(fastq-dump --gzip --split-3 -A `basename $i` -O $wkd/rawdata $i &);done