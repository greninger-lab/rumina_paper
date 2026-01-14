import os
import glob

def generateReaperMetaData(infile, sample_table_file):
    '''Take the sample_table and use it to generate a metadata table
    for REAPER in the correct format.'''

    adaptor_5prime = "AGATCGGAAGAGCGACGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
    adaptor_3prime = "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG"

    outlines = []
    lane = os.path.basename(infile).split('_ext.fastq.gz')[0]

    with open(sample_table_file, 'r') as sample_table:
        for line in sample_table:
            fields = line.strip().split("\t")
            barcode = fields[1]
            lanes = fields[-1].strip().split(",")
            if lane in lanes:
                outlines.append([barcode, adaptor_3prime, adaptor_5prime, "-"])

    output_filename = f"{lane}_reaper_metadata.tsv"

    with open(output_filename, 'w') as outf:
        outf.write("\t".join(["barcode", "3p-ad", "tabu", "5p-si"]) + "\n")
        for outline in outlines:
            outf.write("\t".join(outline) + "\n")

if __name__ == "__main__":
    fastq_files = glob.glob("*_ext.fastq.gz")
    sample_table_path = "sample_table.tsv"

    for fastq_file in fastq_files:
        generateReaperMetaData(fastq_file, sample_table_path)
