import os
import random
import gzip
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ==== SEED FOR DETERMINISM ====
SEED = 31 #changed to seed 42 and 53 for other two iterations
random.seed(SEED)

# ==== PARAMETERS ====

INPUT_FASTA = "HIVpop_simulation_v2.fasta"
OUTPUT_FASTA_GZ = "HIVpop_simulation_v2_frg.fasta.gz"

UMI_LENGTH = 10 # changed in conditions UMI=8 and UMI=12
FRAGMENT_MEAN = 300
FRAGMENT_STD = 100
PCR_CYCLES = 7 # changed in conditions PCR=10 and PCR=13
PCR_EFFICIENCY = 0.9
CHUNK_SIZE = 10000
MAX_PCR_COPIES = 300 # changed in conditions PCR=10 (to 1200) and PCR=13 (to 8500)
ERROR_RATE_PCR = 1e-6 # doi: 10.1155/2014/287430
TEMP_DIR = "tmp_chunks"

# ==== FUNCTION ====

def random_umi(length=UMI_LENGTH):
    return ''.join(random.choices("ACGT", k=length))

def pcr_amplification(cycles, efficiency):
    expected = (1 + efficiency) ** cycles
    copies = max(int(random.gauss(expected, expected * 0.1)), 1)
    return min(copies, MAX_PCR_COPIES)

# remove fragments <150pb as limitation of art-illumina 2x150 simulation (script #3)
def fragment_sequence(seq, mean_len=FRAGMENT_MEAN, std=FRAGMENT_STD):
    fragments = []
    i = 0
    while i < len(seq) - 100:
        frag_len = int(random.gauss(mean_len, std))
        frag_len = max(150, min(frag_len, len(seq) - i))
        frag = seq[i:i + frag_len]
        if len(frag) >= 150:
            fragments.append(frag)
        i += frag_len
    return fragments
    
def introduce_pcr_errors(seq, error_rate=ERROR_RATE_PCR):
    bases = ['A', 'C', 'G', 'T']
    mutated = []
    for base in seq:
        if random.random() < error_rate:
            alt = random.choice([b for b in bases if b != base])
            mutated.append(alt)
        else:
            mutated.append(base)
    return ''.join(mutated)

def write_chunk(records, chunk_idx):
    os.makedirs(TEMP_DIR, exist_ok=True)
    tmp_path = os.path.join(TEMP_DIR, f"chunk_{chunk_idx:04d}.fasta")
    with open(tmp_path, "w") as f:
        SeqIO.write(records, f, "fasta")
    return tmp_path

def compress_final(output_uncompressed, output_gz):
    with open(output_uncompressed, 'rb') as f_in, gzip.open(output_gz, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(output_uncompressed)
    
def introduce_umi_errors(umi, error_rate=0.005): #changed in conditions UMI error rates 0.01 and 0.001
    bases = ['A', 'C', 'G', 'T']
    mutated = []
    for base in umi:
        if random.random() < error_rate:
            alt = random.choice([b for b in bases if b != base])
            mutated.append(alt)
        else:
            mutated.append(base)
    return ''.join(mutated)

# ==== PROCESSING ====

def process_templates(input_path, output_path_gz):
    total_templates = 0
    total_fragments = 0
    total_copies = 0
    chunk_idx = 0
    records = []
    chunk_paths = []
    templates_in_chunk = 0

    for record in SeqIO.parse(input_path, "fasta"):
        total_templates += 1
        templates_in_chunk += 1
        template_id = record.id
        fragments = fragment_sequence(str(record.seq))

        for frag_idx, frag in enumerate(fragments):
            umi = random_umi()
            n_copies = pcr_amplification(PCR_CYCLES, PCR_EFFICIENCY)

            for copy_idx in range(n_copies):
                mutated_frag = introduce_pcr_errors(frag)
                mutated_umi = introduce_umi_errors(umi)
                copy_id = f"{template_id}_frag{frag_idx+1}_copy{copy_idx+1}_UMI:{mutated_umi}"
                fragment_record = SeqRecord(Seq(mutated_frag), id=copy_id, description="")
                records.append(fragment_record)
                total_copies += 1

        total_fragments += len(fragments)
        
        if templates_in_chunk >= 10000:
            print(f"🗂 Processing chunk {chunk_idx} (processed templates: {total_templates})...")
            path = write_chunk(records, chunk_idx)
            chunk_paths.append(path)
            chunk_idx += 1
            records = []
            templates_in_chunk = 0

    if records:
        print(f"🗂 Processing chunk {chunk_idx} (processed templates: {total_templates})...")
        path = write_chunk(records, chunk_idx)
        chunk_paths.append(path)

    # Concatenate chunks without gz (temporary file)
    final_uncompressed = output_path_gz.replace(".gz", "_tmp.fasta")
    with open(final_uncompressed, "w") as outfile:
        for chunk_file in chunk_paths:
            with open(chunk_file, "r") as infile:
                shutil.copyfileobj(infile, outfile)

    # Compression final file
    print("💾 Compress final file...")
    compress_final(final_uncompressed, output_path_gz)

    # Clean temporary chunks
    for f in chunk_paths:
        os.remove(f)
    if os.path.exists(TEMP_DIR) and not os.listdir(TEMP_DIR):
        os.rmdir(TEMP_DIR)

    # Summary report
    print(f"✅ Processing complete")
    print(f"  - Process templates: {total_templates}")
    print(f"  - Fragments generated: {total_fragments}")
    print(f"  - Total copies PCR:   {total_copies}")
    print(f"  - File output:       {output_path_gz}")

# ==== MAIN ====

if __name__ == "__main__":
    process_templates(INPUT_FASTA, OUTPUT_FASTA_GZ)
