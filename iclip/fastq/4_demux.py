import os
import re
import gzip
import glob
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import shutil
import pysam
import stdin

def demux_fastq(infiles):
    '''Demultiplex each fastq file into a separate file for each
    barcode/UMI combination.'''

    infile, meta, samples = infiles
    track = infile.split("_ext.fastq.gz")[0]

    # Create output directory if it doesn't exist
    output_dir = 'demux_fq'
    lint_dir = os.path.join(output_dir, "lint")
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(lint_dir, exist_ok=True)

    # Process the samples file
    with open(samples, 'r') as sample_file:
        for line in sample_file:
            line = line.strip().split("\t")
            if len(line) < 3:
                continue  # Skip lines that don't have enough columns

            bc, name, lanes = line[1:]
            name = name.strip()
            lanes_list = lanes.strip().split(",")

            if track in lanes_list:
                # Construct the output filename
                clean_filename = f"{output_dir}/{track}_{bc}.clean.gz"
                new_filename = f"{output_dir}/{name}_{track}.fastq.gz"

                # Check if the clean file exists before renaming
                if os.path.exists(clean_filename):
                    os.rename(clean_filename, new_filename)
                    print(f"Renamed {clean_filename} to {new_filename}")
                else:
                    print(f"Warning: {clean_filename} does not exist.")

    # Call reaper command (assuming it's still needed)
    reaper_command = [
        'reaper',
        '-geom', '5p-bc',
        '-meta', meta,
        '-i', f"<( gzcat {infile} | sed 's/ /_/g')",
        '--noqc',
        '-3p-head-to-tail', '2',
        '-3p-prefix', '6/2/1',
        '-basename', f"{output_dir}/{track}_",
        '-clean-length', '15'
    ]

    subprocess.run(' '.join(reaper_command),
                   shell = True,
                   executable='/bin/bash')


if __name__ == "__main__":

    fastqs = glob.glob("*_ext.fastq.gz")
    metasheets = glob.glob("*_reaper_metadata.tsv")

    fastqs.sort()
    metasheets.sort()

    with ThreadPoolExecutor(max_workers = os.cpu_count()) as executor:
        futures = {executor.submit(demux_fastq, [f, m, "sample_table.tsv"]): (f, m) for (f, m) in zip(fastqs, metasheets)}

    lint_files = glob.glob("*.lint.gz")

    for file in lint_files:
        shutil.move(file, lint_dir)

    clean_fastqs = glob.glob("*.gz")

                        


