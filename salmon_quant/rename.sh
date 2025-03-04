#!/bin/bash
# rename_quant.sh
# This script renames the quant.sf file in each sample subdirectory within a given output directory.
#
# Usage:
#   ./rename_quant.sh /path/to/output_directory

set -euo pipefail

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 /path/to/output_directory"
  exit 1
fi

OUTPUT_DIR="$1"

# Loop over each subdirectory in the OUTPUT_DIR
for sample_dir in "$OUTPUT_DIR"/*/ ; do
  # Remove the trailing slash and get the sample name
  sample_dir="${sample_dir%/}"
  sample=$(basename "$sample_dir")
  old_file="${sample_dir}/quant.sf"
  new_file="${sample_dir}/${sample}_quant.sf"
  
  if [ -f "$old_file" ]; then
    mv "$old_file" "$new_file"
    echo "Renamed $old_file to $new_file"
  else
    echo "Warning: quant.sf not found in ${sample_dir}"
  fi
done
