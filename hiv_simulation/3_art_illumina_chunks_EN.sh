#!/bin/bash

# === CONFIGURATION ===

INPUT_FASTA="HIVpop_simulation_v2_frg=.fasta.gz"
CHUNK_DIR="chunks"
OUTPUT_DIR="fastq_chunks"
FINAL_OUTPUT_PREFIX="final_amplicon_reads"
READ_LEN=150
COVERAGE=1
SEQ_SYSTEM="HS25"
JOBS=8  # change this depending on your available system (Mac M1: 8-10)
SEED=31 #used in art-illumina (changed to seed 42 and 53 for other two iterations)

# === MAKE DIRECTORIES ===

mkdir -p "$CHUNK_DIR"
mkdir -p "$OUTPUT_DIR"

# === SPLIT MULTIFASTA INTO CHUNKS OF 1M SEQS ===

echo "🔧 Splitting multiFASTA into chunks of 1M sequences..."

python3 <<EOF
from Bio import SeqIO
import os
import gzip

input_file = "$INPUT_FASTA"
chunk_dir = "$CHUNK_DIR"
chunk_size = 1000000

os.makedirs(chunk_dir, exist_ok=True)
if input_file.endswith(".gz"):
    handle = gzip.open(input_file, "rt")
else:
    handle = open(input_file, "r")

chunk_idx = 0
records = []
for i, record in enumerate(SeqIO.parse(handle, "fasta")):
    records.append(record)
    if (i + 1) % chunk_size == 0:
        out_path = os.path.join(chunk_dir, f"chunk_{chunk_idx:04d}.fasta")
        SeqIO.write(records, out_path, "fasta")
        print(f"✅ Chunk {chunk_idx:04d} escrito con {len(records)} secuencias")
        chunk_idx += 1
        records = []
if records:
    out_path = os.path.join(chunk_dir, f"chunk_{chunk_idx:04d}.fasta")
    SeqIO.write(records, out_path, "fasta")
    print(f"✅ Chunk {chunk_idx:04d} escrito con {len(records)} secuencias")
handle.close()
EOF

# === DEFINE FUNCTION TO PARALEL WORK === 

simulate_chunk() {
    chunk_path="$1"
    chunk_base=$(basename "$chunk_path" .fasta)

    echo "🚀 Processing $chunk_base..."

    art_illumina -ss "$SEQ_SYSTEM" -amp -p -na \
        -i "$chunk_path" \
        -l "$READ_LEN" -f "$COVERAGE" \
        -rs "$SEED" \
        -o "$OUTPUT_DIR/${chunk_base}_"

    # Compress FASTQs when exist
    if [[ -f "$OUTPUT_DIR/${chunk_base}_1.fq" && -f "$OUTPUT_DIR/${chunk_base}_2.fq" ]]; then
        gzip -f "$OUTPUT_DIR/${chunk_base}_1.fq"
        gzip -f "$OUTPUT_DIR/${chunk_base}_2.fq"
    else
        echo "⚠️ ART failed or didn't generated $chunk_base"
    fi
}

export -f simulate_chunk
export OUTPUT_DIR READ_LEN COVERAGE SEQ_SYSTEM SEED

# === RUN WITH GNU PARALLEL ===

echo "🎬 Starting simulation with ART-Illumina in paralel..."

find "$CHUNK_DIR" -name "*.fasta" | parallel --jobs "$JOBS" simulate_chunk {}

# === COMBINE ALL FASTQs ===

echo "🧬 Combining all FASTQ.gz..."

cat "$OUTPUT_DIR"/*_1.fq.gz > "${FINAL_OUTPUT_PREFIX}_R1.fastq.gz"
cat "$OUTPUT_DIR"/*_2.fq.gz > "${FINAL_OUTPUT_PREFIX}_R2.fastq.gz"

echo "✅ Ready:"
echo "  👉 ${FINAL_OUTPUT_PREFIX}_R1.fastq.gz"
echo "  👉 ${FINAL_OUTPUT_PREFIX}_R2.fastq.gz"
