# Author: Reed Liu
# Mail: jieandze1314@gmail.com
# Date: 2019-01-24

#PBS -N wes-test
#PBS -o /public/home/liuyunze/result/wes-test.out
#PBS -e /public/home/liuyunze/result/wes-test.err
#PBS -l nodes=2:ppn=20

src=/vol2/agis/xiaoyutao_group/liuyunze/project/single-sap-wes/src
cd $src
cat conf | while read i;do
	fqs=($i)
	fq1=${fqs[0]}
	fq2=${fqs[1]}
	sh wes-single.sh $fq1 $fq2
done
