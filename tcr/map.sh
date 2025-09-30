for file in ./umis_in_name/*.fastq; do 
	base_name=$(basename "$file")

	echo "writing to mapped/${base_name%%_*}.bam"

	bowtie2 -x hg38 -U "$file" --local -p 10 \
		| samtools view -bS --min-MQ 35 -o mapped/"${base_name%%_*}.bam"
done
