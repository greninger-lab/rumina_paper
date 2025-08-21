## RUMINA: 

This repository contains the code used for analysis and figures for the manuscript titled: 

#### Files: 
- analysis: scripts used to create tables and figures for the manuscript.
- hiv_wgs, iclip, tcr: files used to prepare the three datasets used
- memtest.py: the benchmarking script used to output runtime and memory usage. Comments on usage are included within the script.
- hum_bcrtcr.fa and human_IMGT+C.fa: auxiliary files needed to run trust4.

#### Software needed:

- R and Rstudio (along with packages specified in .Rmd files)
- Python 2.7 (iCLIP calc_repro.py) and Python 3.11+ (all other analysis) along with packages specified in .py files.
- A recent version of Bash
- [GNU Parallel](https://www.gnu.org/software/parallel/)
- [samtools 1.21 (htslib 1.22)](https://github.com/samtools/samtools)
- [iVar v1.3](https://github.com/andersen-lab/ivar)
- [UMI-tools v1.6](https://github.com/CGATOxford/UMI-tools)
- [UMICollapse v1.0.0](https://github.com/Daniel-Liu-c0deb0t/UMICollapse)
- [RUMINA v0.99](https://github.com/epiliper/rumina)
- [fastp v0.23.2](https://github.com/OpenGene/fastp)
- [Reaper](https://gensoft.pasteur.fr/docs/reaper/15-065/reaper.html)
- Bowtie2 v2.5.4
- BWA MEM 0.7.19-r1273

**Before starting: note that you will need to decompress `analysis/tcr.tar.gz` to get the TRUST4 files needed for analysis. Run `sh analysis/decompress.sh` to do this.**
