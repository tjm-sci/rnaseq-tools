#!/bin/bash

# Create directory if needed
mkdir -p reference_seqs
cd reference_seqs

# Download transcriptome and genome from GENCODE M36 (GRCm39)
wget -O GRCm39_transcripts.fa.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/gencode.vM36.transcripts.fa.gz
wget -O GRCm39.fa.gz  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/GRCm39.genome.fa.gz

GENOME="GRCm39.fa.gz"
TRANSCRIPTOME="GRCm39_transcripts.fa.gz"

# Extract decoy sequences
gunzip -c "$GENOME" | grep "^>" | cut -d " " -f1 | sed 's/>//g' > decoys.txt

# Concatenate correctly (fix!)
gunzip -c "$TRANSCRIPTOME" "$GENOME" > gentrome.fa

# Run Salmon indexing
salmon index -t gentrome.fa -d decoys.txt -p 18 -i salmon_index --gencode

# Optional: compress the gentrome.fa file after indexing (save space)
gzip gentrome.fa
