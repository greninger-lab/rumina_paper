### NOTE: RUN with Python 2.7
mkdir -p strict 
mkdir -p singletons
mkdir -p original

direc_dir="mapped/directional"
acyclic_dir="mapped/acyclic"
raw_dir="mapped/raw"

runs=($direc_dir $acyclic_dir $raw_dir)

##########################
### WITHOUT SINGLETONS ###
##########################

for dir in "${runs[@]}"; do 
	mkdir -p $dir
	rumina mapped/ --grouping_method $(basename $dir) --separator _ --threads 8 --split_window 0 --outdir $(basename $dir) 
	# source $HOME/stuff/pytwo/bin/activate
	python repro.py $dir _rumina.bam --csv-suffix _$(basename $dir) --csv-dir strict
done;

#######################
### WITH SINGLETONS ###
#######################

for dir in "${runs[@]}"; do 
	mkdir -p $dir
	rumina mapped/ --grouping_method $(basename $dir) --separator _ --singletons --threads 8 --split_window 0 --outdir $(basename $dir) 
	##source $HOME/stuff/pytwo/bin/activate
	python repro.py $dir _rumina.bam --csv-suffix _$(basename $dir) --csv-dir singletons
done;


runs=(
  "UMI_tools_0_None"
  "UMICollapse_0_None"
)


for dir in "${runs[@]}"; do 
	for file in $dir/*.bam; do 
    [ ! -f ${file}.bai ] && samtools index -@ 8 $file; 
  done
	# method=$(basename $dir)
	# mkdir -p $dir
	# source $HOME/stuff/pytwo/bin/activate
	python repro.py $dir _dedup.bam --csv-suffix _$(basename $dir) --csv-dir other_tools
  python3 repro.py $dir .bam
done;

## now look at original, unprocessed bams
python3 repro.py ./ .bam
