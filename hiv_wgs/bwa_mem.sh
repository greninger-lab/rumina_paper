bwa index hiv.fasta
bwa mem -t 8 hiv.fasta SRR11207257_1_ex_trim_fil.fastq.gz SRR11207257_2_ex_trim_fil.fastq.gz | samtools view -bS -F 2052 > hiv_aln.bam
samtools sort -@ 10 hiv_aln.bam -o hiv_bwa.bam

