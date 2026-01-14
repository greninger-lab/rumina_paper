OUTDIR=grouped
TEMPDIR=${OUTDIR}/temp

mkdir -p $OUTDIR
mkdir -p $TEMPDIR
mkdir -p cluster_reports

## first do grouping and tag clusters
function gen_cluster_report() {
  file=$1
  outdir=$2
  tempdir=$3
  prefix=${file%.bam}

  split_dir=$(basename $prefix)
  split_dir="${split_dir}_split"

  # tag groups, and sort by tag.
  rumina -t 2 $file -g directional --singletons --only-group --outdir $outdir -s _ &&
    samtools sort -@ 10 -t UG $outdir/${prefix}_RUMINA.bam -o $outdir/${prefix}_tagsorted.bam &&
    rm $outdir/${prefix}_RUMINA.bam

  # generate the cluster report
  tag_splitter -i $outdir/${prefix}_tagsorted.bam -d $split_dir -t UG -o ${prefix}_cluster_report.tsv > ${prefix}_tag_split_report.tsv
  rm -r $split_dir
}

export -f gen_cluster_report

ls *.bam | parallel gen_cluster_report {} $OUTDIR $TEMPDIR

mv *cluster_report.tsv cluster_reports/

## now make the report
REPORT=${OUTDIR}/all_files_cluster_report.tsv
[ -f $REPORT ] && rm $REPORT

header_written=0

for file in cluster_reports/*_cluster_report.tsv; do
    srr=$(echo "$file" | sed -E 's/.*_(reads.*).*_sorted.*$/\1/')
    echo $srr
    if [ $header_written -eq 0 ]; then
        # Write header from first file, skipping first column and prepending "SRR"
        awk -v srr="SRR" 'NR==1 {print srr "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9}' "$file" > "$REPORT"
        header_written=1
    fi
    # Append data rows, skipping first column and prepending SRR
    awk -v srr="$srr" 'NR>1 {print srr "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9}' "$file" >> "$REPORT"
done
