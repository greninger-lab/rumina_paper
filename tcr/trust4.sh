# for file in ./*/*.bam; do 
#

trust4 () {

  outdir="out_trust4"
  mkdir -p $outdir

  bam=$1
  prefix=${bam%.bam}
  fname=$(basename $prefix)

  desc=$(dirname $bam)
  out=${outdir}/${desc}

  mkdir -p $out
  echo $desc

  # if [ ! -f ${prefix}_report.tsv ]; then
  set -x 
  run-trust4 -b $bam -f hum_bcrtcr.fa --ref human_IMGT+C.fa -t 10 -o ${out}/${fname}
  set +x 
}

export -f trust4

# ls ./*[148]_threads*directional*None/*.bam | parallel trust4 {}
ls ./*UMI-tools*/*.bam | parallel trust4 {}
find out_trust4/*/ -type f ! -name "*_report.tsv" -delete
