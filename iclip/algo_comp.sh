### NOTE: RUN with Python 2.7
mkdir -p strict 
mkdir -p singletons
mkdir -p original

#direc_dir="mapped/directional"
#acyclic_dir="mapped/acyclic"
#raw_dir="mapped/raw"

#runs=($direc_dir $acyclic_dir $raw_dir)

###########################
#### WITHOUT SINGLETONS ###
###########################

#for dir in "${runs[@]}"; do 
#	mkdir -p $dir
#	rumina mapped/ --grouping_method $(basename $dir) --separator _ --threads 8 --split_window 0 --outdir $(basename $dir) 
#	source $HOME/stuff/pytwo/bin/activate
#	python repro.py $dir _rumina.bam --csv-suffix _$(basename $dir) --csv-dir strict
#done;

## Just run already output files
########################
#### WITH SINGLETONS ###
########################

#for dir in "${runs[@]}"; do 
#	mkdir -p $dir
#	rumina mapped/ --grouping_method $(basename $dir) --separator _ --singletons --threads 8 --split_window 0 --outdir $(basename $dir) 
#	source $HOME/stuff/pytwo/bin/activate
#	python repro.py $dir _rumina.bam --csv-suffix _$(basename $dir) --csv-dir singletons
#done;


runs=(
  "RUMINA_8_threads.acyclic_0_--singletons"
  "RUMINA_8_threads.acyclic_0_None"
  "RUMINA_8_threads.directional_0_--singletons"
  "RUMINA_8_threads.directional_0_None"
)

# runs=($dir1 $dir2 $dir3)


for dir in "${runs[@]}"; do 
	for file in $dir/*.bam; do 
    [ ! -f ${file}.bai ] && samtools index -@ 8 $file; 
  done
	# method=$(basename $dir)
	# mkdir -p $dir
	# source $HOME/stuff/pytwo/bin/activate
	# python repro.py $dir _dedup.bam --csv-suffix _$(basename $dir) --csv-dir extra ## check suffix
  python3 repro.py $dir .bam
done;

## now look at original, unprocessed bams
python3 repro.py ./ .bam
