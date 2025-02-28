#!/bin/bash

# Make a directory for our referene transcriptome and genome
if [ ! -d reference_seqs ]; then
    mkdir reference_seqs
fi 

# make sure salmon runs in the correct directory
cd reference_seqs

# reference transcriptome and genome from GENCODE (M36, MGRCm39)

wget -O GRCm39_transcripts.fa.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/gencode.vM36.transcripts.fa.gz

# whole genome

wget -O GRCm39.fa.gz  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/GRCm39.genome.fa.gz

GENOME="GRCm39.fa.gz"
TRANSCRIPTOME="GRCm39_transcripts.fa.gz"

# Get all headers from the genome fasta file 
# cut command used to remove anyhting after the first space in the header
gunzip -c "$GENOME" | grep "^>" | cut -d " " -f1 > decoys.txt

# remove the leading ">" presenet in fasta headers, create backup of original file
sed -i.bak 's/>//g' decoys.txt

# concatenate the genome to the transcriptome to create a combined reference
cat "$TRANSCRIPTOME" "$GENOME" > gentrome.fa.gz

# salmon index

salmon index -t gentrome.fa.gz -d decoys.txt -p 18 -i salmon_index --gencode
