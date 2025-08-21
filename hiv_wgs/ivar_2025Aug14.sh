get_date() {
  echo $(date '+D%m%d% %H%M%S' | sed 's/ /_/g')
}

exec_date=$(get_date)

outdir="ivar"/${exec_date}_"ivar_reports"
mkdir -p $outdir

run_ivar() {
  input=$1
  outdir=$2
  prefix=${input%.bam}

  outprefix=$(basename $prefix)
  outtop=$(dirname $prefix)

  outfile=${outdir}/${outtop}_${outprefix}
  echo $outfile

  # ivar variants
  samtools mpileup -aa -A -d 0 -B -Q 0 --reference hiv.fasta $input | ivar variants -p $outfile -t 0.001 -m 10 -r hiv.fasta

  # depth per position
  samtools depth -a $input > ${outfile}_depth.tsv

}

export -f run_ivar
{ ls */*/*.bam; } | parallel run_ivar {} $outdir

