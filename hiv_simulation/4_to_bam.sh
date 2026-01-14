bwa index HIV_HXB2.fasta
for file in *_R1.fastq.gz; do
  prefix=${file%_R1.fastq.gz}
  r2=${prefix}_R2.fastq.gz

  echo $file $r2
  bwa mem -t 10 HIV_HXB2.fasta $file $r2 | samtools view -h -b > ${prefix}.bam
done
