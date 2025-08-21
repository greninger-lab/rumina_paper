"""
NOTE: parts of this code are inspired from the following repo:
https://github.com/JappyPing/Deduplication_ErrorCorrection


ALL TESTING PERFORMED ON M1 Mac Studio:
    CPU: M1 Pro Max (10-core)
    RAM: 64GB
"""

import subprocess
import psutil
import os
import time
import pandas as pd
import pysam
import datetime

class OveruseException(Exception):
    pass


"""
NOTE: max heap size was adjusted to 110GB:

this change was made to line 20 of the UMICollapse executable:
    default_jvm_mem_opts = ["-Xms2g", "-Xmx110g"]
                                        ^^^^^^^

UMICollapse was installed via conda as described here: https://github.com/Daniel-Liu-c0deb0t/UMICollapse

the UMICollapse executable modified was found here:
/opt/homebrew/Caskroom/miniforge/base/bin/umicollapse
"""

def run_umicollapse(dir, _group_method, _singletons, file, iteration, separator, _ref, _threads, _length, _paired_method):
    method = f"UMICollapse_{iteration}"

    outfile = os.path.join(dir, method, file.split(".bam")[0] + "_dedup.bam")

    cmd = [
        "umicollapse",
        "bam",
        "-i",
        os.path.join(dir, file),
        "--umi-sep",
        f"{separator}",
        "-o",
        # os.path.join(dir, method, outfile),
        outfile,
    ]

    return outfile, method, cmd



"""
UMI-tools was installed through pip as described here: https://github.com/CGATOxford/UMI-tools.
No changes were made.
"""

def run_umitools(dir, _group_method, _singletons, file, iteration, separator, _ref, _threads, length, _paired_method):
    method = f"UMI-tools_{iteration}"

    outfile = os.path.join(dir, method, file.split(".bam")[0] + "_dedup.bam")

    cmd = [
        "umi_tools",
        "dedup",
        "-I",
        os.path.join(dir, file),
        "--umi-separator",
        f"{separator}",
        "--method",
        "directional",
        "--stdout",
        outfile,
        "--random-seed=0",
    ]

    if length: cmd.extend(["--read-length"])

    return outfile, method, cmd


def run_rumina(dir, group_method, singletons, file, iteration, separator, _ref, threads, length, paired_method):
    method = f"RUMINA_{threads}_threads.{group_method}_{iteration}_{singletons}"

    outfile = os.path.join(dir, method, f"{file.split('.bam')[0]}_RUMINA.bam")

    cmd = [
        "rumina",
        "--grouping_method",
        group_method,
        os.path.join(dir, file),
        "--separator",
        f"{separator}",
        "--threads",
        str(threads),
        "--outdir",
        os.path.join(dir, method),
    ]
    if singletons:
        cmd.append("--singletons")

    if length:
        cmd.append("--length")

    return outfile, method, cmd


# tuple structure
### cmd to run, grouping algorithm to use, whether or not to keep singletons, threads, paired-end handling method
### 'None' means default value was used
METHODS = [
    ## other tools 

    # (run_umitools, None, None, None, None),
    # (run_umicollapse, None, None, None, None),

    ## directional method

    # (run_rumina, "directional", "--singletons", 1, None),
    # (run_rumina, "directional", "--singletons", 4, None),
    # (run_rumina, "directional", "--singletons", 8, None),
    
    # (run_rumina, "directional", None, 1, None),
    # (run_rumina, "directional", None, 4, None),
    # (run_rumina, "directional", None, 8, None),
    
    ## acyclic method

    # (run_rumina, "acyclic", "--singletons", 1, None),
    # (run_rumina, "acyclic", "--singletons", 4, None),
    # (run_rumina, "acyclic", "--singletons", 8, None),
    
    # (run_rumina, "acyclic", None, 1, None),
    # (run_rumina, "acyclic", None, 4, None),
    # (run_rumina, "acyclic", None, 8, None),

]

NUM_ITERATIONS = 3

# amount of time a tool is allowed to use before being terminated, in seconds
MAX_TIME = 14400

# dataset name, whether or not it's paired, and samtools-indexed reference fasta path (required by umierrorcorrect), extra
# whether or not to stratify by length
DATASETS = [
    ("iclip", False, "mm9.fa", False),
    ("tcr", False, "hg38.analysisSet.fa", False),
    ("hiv_wgs", True, "hiv.fasta", False),
]  # directories containing bamfiles to test


def update_report(report, columns, dat):
    minireport = dict()
    for col, d in zip(columns, dat):
        minireport[col] = [d]

    report = pd.concat([report, pd.DataFrame(minireport)], ignore_index=True)
    return report


def get_reads(file):
    total_reads = 0
    mapped_reads = 0

    try:
        total_reads += int(pysam.view("-@ 8", "-c", file))

        mapped_reads += int(pysam.view("-@ 8", "-F 4", "-c", file))
    except pysam.utils.SamtoolsError:
        total_reads = -1
        mapped_reads = -1

    return total_reads, mapped_reads


def get_mem_usage(pid):
    try:
        process = psutil.Process(pid)
        children = process.children(recursive=True)

        memory_gb = process.memory_info().rss / 1024**3

        ## get mem of child procs too (looking at you, JVM)
        for c in children:
            memory_gb += c.memory_info().rss / 1024**3

        return memory_gb

    except psutil.NoSuchProcess:
        return 0.0


def monitor_cmd(cmd, maxtime):
    process = subprocess.Popen(cmd)
    start_time = time.time()

    pid = process.pid
    peak_memory_usage = 0.0
    parent = psutil.Process(pid)

    try:
        while True:
            if process.poll() is not None:
                break

            memory_usage = get_mem_usage(pid)
            if memory_usage > peak_memory_usage:
                peak_memory_usage = memory_usage

            runtime = time.time() - start_time
            if runtime > maxtime:
                proc = psutil.Process(pid)
                for c in proc.children(recursive=True):
                    c.kill()

                proc.kill()
                raise OveruseException(
                    f"\033[31mMaximum runtime of {maxtime} seconds exceeded! Killing process...\033[0m"
                )

            print(
                f"\033[31mCurrent Peak Memory Usage: {peak_memory_usage:.2f} GB\033[0m",
                end="\r",
            )

            time.sleep(0.1)  # Sleep for a while before checking again

    except OveruseException as e:
        print(e)
        elapsed_time = time.time() - start_time

        return peak_memory_usage, elapsed_time, "Terminated"

    elapsed_time = time.time() - start_time

    print(
        f"\033[31mElapsed Time: {elapsed_time:.2f} seconds | Peak Memory Usage: {peak_memory_usage:.2f} GB\033[0m",
        flush=True,
    )

    return (peak_memory_usage, elapsed_time, "Completed")


def main():
    """
    Run each tool on each provided dataset, recording input and output reads, runtime and peak memory usage.
    """
    report = pd.DataFrame()

    columns = [
        "Tool",
        "File",
        "Peak memory used",
        "Runtime",
        "Grouping Method",
        "Singletons",
        "Dataset",
        "Completed?",
        "Input reads (total)",
        "Input reads (mapped)",
        "Ouput reads (total)",
        "Output reads (mapped)",
        "Iteration",
    ]

    separator = None
    group_method = None

    outdirs_to_check = set()

    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S").replace(" ", "_")

    for col in columns:
        report[col] = ""

    for dir, is_paired, ref, use_read_length in DATASETS:

        print(dir, is_paired)

        if dir == "iclip" or dir == "tcr" or dir == "hiv_wgs":
            separator = "_"
        else:
            separator = ":"

        for iteration in range(NUM_ITERATIONS):
            for method, group_method, singletons, threads, pair_method in METHODS:
                for file in os.listdir(dir):
                    if file.endswith(".bam"):
                        outfile, tool, cmd = method(
                            dir,
                            group_method,
                            singletons,
                            file,
                            iteration,
                            separator,
                            ref, 
                            threads,
                            use_read_length,
                            pair_method,
                        )

                        total_reads_in, mapped_reads_in = get_reads(
                            os.path.join(dir, file)
                        )

                        outdirs_to_check.add(tool)

                        if not os.path.exists(os.path.join(dir, tool)):
                            os.mkdir(os.path.join(dir, tool))

                        # modify params for dataset features
                        if is_paired:
                            if "UMI-tools" in tool:
                                cmd.extend(["--paired"])
                            if "UMICollapse" in tool:
                                cmd.extend(["--paired"])
                            elif "RUMINA" in tool:
                                if dir == "hiv_wgs" and not pair_method:
                                    cmd.extend(["--merge_pairs", ref])
                                else:
                                    cmd.append("--paired")
                                outfile = outfile.split(".bam")[0] + "_MERGED.bam"

                        if dir == "hiv_wgs":
                            if "RUMINA" in tool:
                                cmd.extend(["-x", "100"])
                                cmd.extend(["--length"])

                        print(cmd)
                    
                        mem, time, status = monitor_cmd(cmd, MAX_TIME)

                        # replace with tool default
                        if singletons is None and "RUMINA" in tool:
                            singleton_policy = "Strict"
                        else:
                            singleton_policy = "Singletons"

                        total_reads_out, mapped_reads_out = get_reads(outfile)

                        if status == "Terminated":
                            break  ## don't try with the other file(s) if time/mem limit is exceeded

                        report = update_report(
                            report,
                            columns,
                            [
                                tool,
                                file,
                                mem,
                                time,
                                group_method,
                                singleton_policy,
                                dir,
                                status,
                                total_reads_in,
                                mapped_reads_in,
                                total_reads_out,
                                mapped_reads_out,
                                iteration,
                            ],
                        )

                        ## writing report with every iteration just to check
                        report.to_csv(os.path.join("report", f"bench_report_{timestamp}.csv"), index=None)


if __name__ == "__main__":
    main()
