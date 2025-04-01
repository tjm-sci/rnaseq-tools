#!/bin/bash

help() {
 echo $0
 echo "A wrapper to run Salmon quantification using a config file."
 echo ""
 echo "Usage:"
 echo " ./run_salmon.sh --config /path/to/config.sh"
}

# this stops the script if any command fails 
set -euo pipefail

# --- Parse arguments ---
if [ "$#" -eq 0 ]; then
  help
  exit 1
fi

while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --config)
      CONFIG_FILE="$2"
      shift 2
      ;;
    *)
      echo "Unknown parameter: $1"
      exit 1
      ;;
  esac
done

if [ ! -f "$CONFIG_FILE" ]; then
  echo "Config file not found: $CONFIG_FILE"
  exit 1
fi

# --- Source configuration file ---
# The config file should export variables such as:
#   FASTQ_DIR, OUTPUT_DIR, INDEX, and optionally THREADS.
source "$CONFIG_FILE"

# --- Check that required variables are set ---
: "${FASTQ_DIR:?Please set FASTQ_DIR in the config file}"
: "${OUTPUT_DIR:?Please set OUTPUT_DIR in the config file}"
: "${INDEX:?Please set INDEX (path to your Salmon index) in the config file}"

# check that the fastq directory exists and is not empty
if [ ! -d "$FASTQ_DIR" ] || [ -z "$(ls -A "$FASTQ_DIR")" ]; then
  echo "Error: FASTQ_DIR is empty or does not exist: $FASTQ_DIR"
  exit 1
fi

# check index directory exists and is not empty
if [ ! -d "$INDEX" ] || [ -z "$(ls -A "$INDEX")" ]; then
  echo "Error: INDEX is empty or does not exist: $INDEX"
  exit 1
fi

# --- Create the output directory if needed ---
mkdir -p "$OUTPUT_DIR"

# --- Loop over read1 files and determine matching read2 ---
# The script finds files matching the pattern *_R1*.fastq* and then constructs the
# corresponding read2 filename by replacing _R1 with _R2.
for r1 in "${FASTQ_DIR}"/*_R1*.fastq*; do
  # Construct the corresponding R2 filename
  r2="${r1/_R1/_R2}"
  
  if [ ! -f "$r2" ]; then
    echo "Warning: Corresponding R2 file not found for $r1 (expected: $r2). Skipping."
    continue
  fi
  
  # Extract a sample name by stripping everything from _R1 onwards.
  sample=$(basename "$r1")
  sample="${sample%%_R1*}"
  sample_output="${OUTPUT_DIR}/${sample}"
  mkdir -p "$sample_output"
  
  echo "Processing sample: $sample"
  
  # Run Salmon with selective alignment, sequence and GC bias correction, using the specified threads.
  salmon quant \
    -i "$INDEX" \
    -l A \
    -1 "$r1" \
    -2 "$r2" \
    --validateMappings \
    --seqBias \
    --gcBias \
    --threads "$THREADS" \
    --softclip \
    -o "$sample_output"


  # Rename the quant.sf file to include the sample ID
  if [ -f "${sample_output}/quant.sf" ]; then
    mv "${sample_output}/quant.sf" "${sample_output}/${sample}_quant.sf"
    echo "Renamed quant.sf to ${sample}_quant.sf in $sample_output"
  else
    echo "Warning: quant.sf not found in ${sample_output}"
  fi
    done