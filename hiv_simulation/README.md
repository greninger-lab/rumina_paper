##Scripts for Simulating and Analyzing a Diverse HIV Population##

The following scripts were used to simulate a diverse HIV population for downstream PCR error and sequencing error modeling:

1_simulate_pupulation_fasta_EN.py: This script generates a multi-FASTA file from a reference nucleotide sequence (downloaded from GenBank in FASTA format). It creates multiple copies of the reference and introduces non-overlapping mutations at predefined frequencies of 1%, 0.1%, and 0.01%.
Note: This template mutation process is not deterministic. The specific multi-FASTA file used for all downstream simulations (HIVpop_simulation_v2.fasta.gz) is included in this repository.

2_simulate_fragmentation_and_pcr_EN.py: This script takes the multi-FASTA generated in Step 1 and simulates the generation of full-length amplicons. Key features include:
Configurable number of PCR cycles
Customizable UMI length, which is added to the FASTA header for tracking
Application of UMI error rate (applied only to amplicons, not to the original templates)

3_art_illumina_chunks_EN.sh: This script uses the output from Step 2 and simulates Illumina sequencing errors using the ART Illumina software.

##Downstream Analysis##

The output FASTQ files from Step 3 are aligned to the original GenBank reference using BWA-MEM.

The resulting sorted BAM files are used as input for benchmarking.

Variant calling is performed using iVar. Output files can be found in the iVar_reports folder.

File names in the iVar_reports directory encode simulation parameters:
ERXX – UMI error rate
PCRXX – number of PCR cycles
UMIXX – UMI length
(Where XX indicates the value used in that specific run)

The main simulation used for all benchmarking and key figures in the paper is: umi10pcr7er05

Important: No variants were filtered using iVar’s Fisher statistical test. This test was designed for datasets without UMIs and assumes different error characteristics. Therefore, all detected variants were retained. However, only variants detected in at least two independent simulation replicates were considered valid.

Benchmarking Information
The file bench_hiv-scalability_2025-08-29_15_50_23.csv contains runtime and memory usage data from benchmark tests.

Benchmarks were performed using 12 PCR cycles with subsets of 3M, 6M, 9M, and 12M reads.

Data Visualization
All plots and analysis summaries are documented in the RMarkdown report:
📄 rumina_hiv_20250926.html