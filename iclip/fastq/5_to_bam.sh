modify_qname() {
    awk '{
        if (NR % 4 == 1) {  
            qname = $1;  
            split(qname, parts, "_");  
            if (length(parts) > 1) {
                if (length(parts[2]) == 5) {
                    new_qname = parts[1];  
                    for (i = 3; i <= length(parts); i++) {
                        new_qname = new_qname "_" parts[i];  
                    }
                    new_qname = new_qname "_" parts[2];  
		    sub(/\/[12]/, "", new_qname);  
                    sub($1, new_qname);  
                }
            }
            print;  
        } else {
            print $0;  
        }
    }'
}

for file in demux_fq/*.gz; do 
    bn=$(basename "$file")
    bam_file="${bn%%_*}_unsorted.bam"
    echo "Processing: $bam_file from $file"

    gzcat "$file" | modify_qname | ../bowtie-1.1.2/bowtie ../mm9 - -v 2 -m 10 -a -S --threads 7 | samtools view -bS > "../$bam_file"


    if [[ ! -f "../mapped/$bam_file" ]]; then
        echo "Error: BAM file not created."
        continue
    fi

    sorted_file="${bam_file%%_*}.bam"
    samtools sort -@ 8 "../$bam_file" > "../$sorted_file"
    samtools index -@ 8 "../$sorted_file"
    rm $bam_file

    echo "Finished processing: $sorted_file"
done
