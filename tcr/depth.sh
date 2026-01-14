mkdir -p depth

printf 'sample\tmean_cov\tmean_depth\n' > depth_report.tsv

for file in *.bam; do
  prefix=${file%.bam}
  pandepth -i $file -o depth/${prefix} -t 10
  covdepth=$(gzcat depth/${prefix}.chr.stat.gz | tail -n 1 | cut -f 3,4)
  cov=$(echo $covdepth | cut -d ' ' -f 2)
  depth=$(echo $covdepth | cut -d ' ' -f 4)

  printf "${prefix}\t${cov}\t${depth}\n" >> depth_report.tsv
done

