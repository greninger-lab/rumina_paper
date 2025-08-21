# for file in *.fastq.gz
# do 
# 	ext_file=${file%%.*}_ext.fastq.gz
# 	umi_tools extract -I $file --extract-method=string --bc-pattern=NNNXXXXNN --log2stderr | gzip > $ext_file
# done;
#
#
parallel -j 10 'umi_tools extract -I {} --extract-method=string --bc-pattern=NNNXXXXNN --log2stderr | pigz -c > {.}_ext.fastq.gz' ::: *.fastq.gz
