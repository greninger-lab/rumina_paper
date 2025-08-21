fastp -i SRR11207257_1_ex.fastq.gz -I SRR11207257_2_ex.fastq.gz \
  -o SRR11207257_1_ex_trim_fil.fastq.gz --out2 SRR11207257_2_ex_trim_fil.fastq.gz \
  --adapter_sequence AGATCGGAAGAG --adapter_sequence_r2 AGATCGGAAGAG \
  -e 20 -w 10 -q 0
