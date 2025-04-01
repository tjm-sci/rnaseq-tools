#!/bin/bash

# Script to generate STAR index.

# Specify paths to required files, as well as output directory
BASE_DIR="/home/tmurphy/phd_work/FANS_PILOT_DEV/250331_reference_seqs"

GENOME="${BASE_DIR}/GRCm39.fa"
GTF="${BASE_DIR}/GRCm39.gtf"

STAR_DIR="/home/tmurphy/phd_work/FANS_PILOT_DEV/star"
mkdir -p "${STAR_DIR}/star_index"


# Run STAR
STAR --runThreadN 18 \
  --runMode genomeGenerate \
  --genomeDir ${STAR_DIR}/star_index \
  --genomeFastaFiles ${GENOME} \
  --sjdbGTFfile ${GTF} \
  --sjdbOverhang 156