set -u
set -e

SEED1=31
SEED2=42
SEED3=53

#    UMI LENGTH   SEED            N PCR CYCLES      UMI_ERR_RATE        OUTPUT FILE PREFIX
conditions=(
  "--UMI_LEN 12 --SEED ${SEED1} --N_PCR_CYCLES 7  --UMI_ERR_RATE 0.005  reads_UMI12_s${SEED1}"
  "--UMI_LEN 12 --SEED ${SEED2} --N_PCR_CYCLES 7  --UMI_ERR_RATE 0.005  reads_UMI12_s${SEED2}"
  "--UMI_LEN 12 --SEED ${SEED3} --N_PCR_CYCLES 7  --UMI_ERR_RATE 0.005  reads_UMI12_s${SEED3}"

  "--UMI_LEN 10 --SEED ${SEED1} --N_PCR_CYCLES 7  --UMI_ERR_RATE 0.005  reads_UMI10_s${SEED1}"
  "--UMI_LEN 10 --SEED ${SEED2} --N_PCR_CYCLES 7  --UMI_ERR_RATE 0.005  reads_UMI10_s${SEED2}"
  "--UMI_LEN 10 --SEED ${SEED3} --N_PCR_CYCLES 7  --UMI_ERR_RATE 0.005  reads_UMI10_s${SEED3}"

  "--UMI_LEN 8  --SEED ${SEED1} --N_PCR_CYCLES 7  --UMI_ERR_RATE 0.005  reads_UMI8_s${SEED1}"
  "--UMI_LEN 8  --SEED ${SEED2} --N_PCR_CYCLES 7  --UMI_ERR_RATE 0.005  reads_UMI8_s${SEED2}"
  "--UMI_LEN 8  --SEED ${SEED3} --N_PCR_CYCLES 7  --UMI_ERR_RATE 0.005  reads_UMI8_s${SEED3}"

  "--UMI_LEN 10 --SEED ${SEED1} --N_PCR_CYCLES 10 --UMI_ERR_RATE 0.005  reads_PCR10_s${SEED1}"
  "--UMI_LEN 10 --SEED ${SEED2} --N_PCR_CYCLES 10 --UMI_ERR_RATE 0.005  reads_PCR10_s${SEED2}"
  "--UMI_LEN 10 --SEED ${SEED3} --N_PCR_CYCLES 10 --UMI_ERR_RATE 0.005  reads_PCR10_s${SEED3}"

  "--UMI_LEN 10 --SEED ${SEED1} --N_PCR_CYCLES 13 --UMI_ERR_RATE 0.005  reads_PCR13_s${SEED1}"
  "--UMI_LEN 10 --SEED ${SEED2} --N_PCR_CYCLES 13 --UMI_ERR_RATE 0.005  reads_PCR13_s${SEED2}"
  "--UMI_LEN 10 --SEED ${SEED3} --N_PCR_CYCLES 13 --UMI_ERR_RATE 0.005  reads_PCR13_s${SEED3}"

  "--UMI_LEN 8  --SEED ${SEED1} --N_PCR_CYCLES 7  --UMI_ERR_RATE 0.01   reads_ER01_s${SEED1}"
  "--UMI_LEN 8  --SEED ${SEED2} --N_PCR_CYCLES 7  --UMI_ERR_RATE 0.01   reads_ER01_s${SEED2}"
  "--UMI_LEN 8  --SEED ${SEED3} --N_PCR_CYCLES 7  --UMI_ERR_RATE 0.01   reads_ER01_s${SEED3}"

  "--UMI_LEN 8  --SEED ${SEED1} --N_PCR_CYCLES 7  --UMI_ERR_RATE 0.001  reads_ER1_s${SEED1}"
  "--UMI_LEN 8  --SEED ${SEED2} --N_PCR_CYCLES 7  --UMI_ERR_RATE 0.001  reads_ER1_s${SEED2}"
  "--UMI_LEN 8  --SEED ${SEED3} --N_PCR_CYCLES 7  --UMI_ERR_RATE 0.001  reads_ER1_s${SEED3}"
  )

# python3 1_simulate_population_fasta_EN.py

## run the simulation with all our parameters above
for condition in "${conditions[@]}"; do
  outfile=$(echo "$condition" | awk '{ print $NF }')
  args=$(echo "$condition" | awk 'NF--')

  if compgen -G ${outfile}_R1.fastq.gz > /dev/null; then
    echo "Already processed condition for $outfile. Skipping..."
    continue
  fi

  echo "${outfile} <- ${args}"

  python3 2_simulate_fragmentation_and_pcr_EN.py $args
  sh 3_art_illumina_chunks_EN.sh

  mv final_amplicon_reads_v2_R1.fastq.gz ${outfile}_R1.fastq.gz
  mv final_amplicon_reads_v2_R2.fastq.gz ${outfile}_R2.fastq.gz

done
