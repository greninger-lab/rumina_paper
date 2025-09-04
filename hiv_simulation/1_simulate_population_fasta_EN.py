import random
from Bio import SeqIO
from Bio.Seq import Seq

# parameters
REFERENCE_FASTA = "HIV_HXB2.fasta"
OUTPUT_FASTA = "HIVpop_simulation_v2.fasta"
TOTAL_TEMPLATES = 30000
CHUNK_SIZE = 10000 

# making SNV groups with defined frequencies
groups = {
    "high": {"freq": 0.01, "count": 3},   # SNV at 1%
    "medium": {"freq": 0.001, "count": 3},# SNV at 0.1%
    "low": {"freq": 0.0001, "count": 3}   # SNV at 0.01%
}

def generate_mutations(ref_seq, groups):
    length = len(ref_seq)
    used_positions = set()
    mutations = {}

    for group in groups:
        mutations[group] = []
        for _ in range(groups[group]["count"]):
            while True:
                pos = random.randint(100, length - 100)
                if pos not in used_positions:
                    ref_base = ref_seq[pos]
                    alt_base = random.choice([b for b in "ACGT" if b != ref_base])
                    used_positions.add(pos)
                    mutations[group].append((pos, ref_base, alt_base))
                    break
    return mutations

def write_chunk(out_handle, chunk_records):
    for rec_id, seq in chunk_records:
        out_handle.write(f">{rec_id}\n")
        out_handle.write(f"{seq}\n")

def main():
    print("Reading ref seq...")
    ref_record = SeqIO.read(REFERENCE_FASTA, "fasta")
    ref_seq = str(ref_record.seq)

    print("Making unique Single Nucleotide Variants...")
    mutations = generate_mutations(ref_seq, groups)

    # Calcular número de templates por mutación
    templates_per_mutation = {}
    for group in groups:
        freq = groups[group]["freq"]
        count = groups[group]["count"]
        # Total templates con mutación de este grupo = freq * TOTAL_TEMPLATES
        # Cada mutación tiene igual cantidad dentro del grupo
        num = int(freq * TOTAL_TEMPLATES)
        templates_per_mutation[group] = num

    # Calcular cantidad total mutada para restar del WT
    total_mutated = sum(templates_per_mutation[g] * groups[g]["count"] for g in groups)
    templates_wt = TOTAL_TEMPLATES - total_mutated
    print(f"Total mutated templates: {total_mutated}")
    print(f"Total WT templates: {templates_wt}")

    chunk = []
    count_written = 0

    with open(OUTPUT_FASTA, "w") as out_fasta:
        # making mutates
        for group in groups:
            mut_list = mutations[group]
            num_per_mut = templates_per_mutation[group]
            for pos, ref_base, alt_base in mut_list:
                for i in range(num_per_mut):
                    seq_list = list(ref_seq)
                    seq_list[pos] = alt_base
                    mutated_seq = "".join(seq_list)
                    rec_id = f"{group}_mut_pos{pos+1}_{i+1}"
                    chunk.append((rec_id, mutated_seq))

                    if len(chunk) >= CHUNK_SIZE:
                        write_chunk(out_fasta, chunk)
                        count_written += len(chunk)
                        print(f"Write {count_written} sequences...")
                        chunk = []

        # making WT
        for i in range(templates_wt):
            rec_id = f"WT_{i+1}"
            chunk.append((rec_id, ref_seq))

            if len(chunk) >= CHUNK_SIZE:
                write_chunk(out_fasta, chunk)
                count_written += len(chunk)
                print(f"Write {count_written} sequences...")
                chunk = []

        # Summary
        if chunk:
            write_chunk(out_fasta, chunk)
            count_written += len(chunk)
            print(f"Write {count_written} sequences...")

    print(f"✅ Finished. Generated file: {OUTPUT_FASTA}")

if __name__ == "__main__":
    main()
