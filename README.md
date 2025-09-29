
This repository contains the code used for analysis and figures for the manuscript titled:   
**"RUMINA: high-throughput UMI deduplication for amplicon and whole-genome sequencing with enhanced error correction."**

NOTE: RUMINA 0.9.81 was used for all analysis.

To install it with Cargo: 
```bash 
cargo install rumina --version 0.9.81
```

#### Files:
- analysis: scripts used to create tables and figures for the manuscript.
- iclip, tcr: files used to prepare the iCLiP and TCR datasets.
- hiv_simulation: files used to prepare the simulated HIV dataset.
- memtest.py: the benchmarking script used to output runtime and memory usage. Comments on usage are included within the script.
- hum_bcrtcr.fa and human_IMGT+C.fa: auxiliary files needed to run TRUST4.
- tag_splitter: the UMI cluster inspection CLI tool used to gather clsuter composition metrics. Use the Makefile inside to build (requires Cargo).

#### Software needed:

- R and Rstudio (along with packages specified in .Rmd files)
- Python 2.7 (iCLIP calc_repro.py) and Python 3.11+ (all other analysis) along with packages specified in .py files.
- A recent version of Bash
- Cargo 1.86.0+
- [GNU Parallel](https://www.gnu.org/software/parallel/)
- [samtools 1.21 (htslib 1.22)](https://github.com/samtools/samtools)
- [iVar v1.3](https://github.com/andersen-lab/ivar)
- [UMI-tools v1.1.6](https://github.com/CGATOxford/UMI-tools)
- [UMICollapse v1.0.0](https://github.com/Daniel-Liu-c0deb0t/UMICollapse)
- [RUMINA v0.99](https://github.com/epiliper/rumina)
- [fastp v0.23.2](https://github.com/OpenGene/fastp)
- [Reaper](https://gensoft.pasteur.fr/docs/reaper/15-065/reaper.html)
- Bowtie2 v2.5.4
- BWA MEM 0.7.19-r1273

**Before starting: note that you will need to decompress `analysis/tcr.tar.gz` to get the TRUST4 files needed for analysis. Run `sh analysis/decompress.sh` to do this.**
