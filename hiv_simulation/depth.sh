
OUTDIR=report_depth

mkdir -p $OUTDIR

get_depth() {
  file=$1
  outdir=$2

  prefix=${file%.bam}
  prefix=$(basename $prefix)
  outfile=${outdir}/${prefix}_depth.tsv

  set -x
  samtools depth -a -d 0 $file | awk -v prefix=$prefix -v OFS="\t" '{ print $0, prefix }' > ${outfile}
  set +x

}

export -f get_depth

ls *.bam | parallel get_depth {} $OUTDIR
