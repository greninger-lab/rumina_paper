"""
WARNING: THIS SCRIPT IS MEANT TO BE RUN WITH PYTHON v2.7

Retrieve the deduplicated files across the first 3 replicates for 
each sample, and run the the calculateReproducibility function from UMI-tools pipelines.
"""

import os 
import argparse
import pandas as pd
import subprocess
import sys

parser = argparse.ArgumentParser()
parser.add_argument("dir")
parser.add_argument("suffix")
parser.add_argument("--csv-suffix", required = False)
parser.add_argument("--csv-dir", required = False)
args = parser.parse_args()

replicates = pd.read_csv("filereport_read_run_PRJNA286202_tsv.txt", sep = "\t")

## get only replicates 1-3
replicates['rep_number'] = replicates['sample_title'].str.extract(r'_(Rep(\d+))')[1]
replicates['rep_number'] = pd.to_numeric(replicates["rep_number"])
replicates = replicates[replicates['rep_number'].astype(int) <= 3]

## get list of samples
samples = replicates["sample_title"].str.split("_").str[0].unique()

def run_repro(dir, suffix, csv_suffix, csv_dir):

    for sample in samples:
        print("sample:" + sample)
        ## get SRAs for all replicates of a sample
        replicate_files = replicates[
                replicates["sample_title"].str.contains(sample)]["run_accession"].tolist()

        ## gather replicate files for reproducibility analysis
        to_run = []
        for rep in replicate_files:
            for file in os.listdir(dir):
                if rep in file and file.endswith(suffix):
                    print("FILE:" + file)
                    to_run.append(os.path.join(dir, file))

        if csv_suffix:
            csv_name = sample + csv_suffix + ".csv" 
        else:
            csv_name = sample + ".csv"

        if csv_dir:
            csv_name = os.path.join(csv_dir, csv_name)
        else:
            csv_name = os.path.join(dir, csv_name)

        subprocess.call(["python3", "calc_repro.py", "--dtype=uint32", "-m", "2", "--stdout="+csv_name] + to_run)
        sys.exit(1)

if __name__ == "__main__":
    run_repro(args.dir, args.suffix, args.csv_suffix, args.csv_dir)
