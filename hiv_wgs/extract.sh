pigz SRR11207257_1.fastq
pigz SRR11207257_2.fastq

rm SRR11207257_1.fastq
rm SRR11207257_2.fastq

umi_tools extract $file --bc-pattern NNNNNNNNNNNNCCCCCCCCCCCC --bc-pattern2 NNNNNNNNNNNNCCCCCCCCCCCC -I SRR11207257_1.fastq.gz --read2-in SRR11207257_2.fastq.gz --stdout SRR11207257_1_ex.fastq.gz --read2-out SRR11207257_2_ex.fastq.gz
                                                                                                                                            

