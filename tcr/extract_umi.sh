mkdir -p umis_in_name 
for file in *.fastq; do 
	umi_tools extract --bc-pattern NNNNNNNNNNNNCCCCCCCCCCCC -I $file -S umis_in_name/${file%%.*}_extracted.fastq
done
