mkdir ~/RNA-seq/stat

# sambtools flagstat (take hisat2 bam result as an example)
cd ~/RNA-seq/quant/hisat2
ls *.bam | while read i; do (samtools flagstat -@ 10 $i > $(basename $i ".bam").flagstat);done
mv *flagstat ~/RNA-seq/stat

# copy stat to Excel
cd ~/RNA-seq/stat

# 一般flasgstat结果都有13行
cat * | awk '{print $1}'| paste - - - - - - - - - - - - -
# 把结果复制、粘贴到excel中=》转置=〉得到行为统计信息，列为全部样本
ls * | paste - - - - - - - - - - - - - - - - #得到全部样本名（作为列名）

# cat 其中任意一个.flagstat | cut -d'+' -f 2 | cut -d' ' -f 3-9 #得到统计信息（行名）

# 参考 https://upload-images.jianshu.io/upload_images/9376801-c3c4e38ec24404ad.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240