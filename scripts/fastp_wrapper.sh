#!/bin/bash
# fastp_config_wrapper.sh
# This script runs fastp in a Docker container on paired-end FASTQ files.
# It requires a configuration file (supplied with --config) that exports:
#   FASTQ_DIR, FASTP_OUTPUT_DIR, FASTP_REPORT_DIR, THREADS, and OUTPUT_SUFFIX.
#
# Usage:
#   ./fastp_config_wrapper.sh --config <config_file>
#
# The script processes paired FASTQ files, outputs trimmed files with the OUTPUT_SUFFIX
# appended to the base filename, and generates HTML and JSON reports with a date prefix.

help() {
    echo ""
    echo "Usage: $0 --config <config_file>"
    echo ""
    echo "The config file must export the following variables:"
    echo "  FASTQ_DIR           Path to raw FASTQ input directory"
    echo "  FASTP_OUTPUT_DIR    Directory for trimmed FASTQ output files"
    echo "  FASTP_REPORT_DIR    Directory for fastp HTML and JSON reports"
    echo "  THREADS             Number of threads to use"
    echo "  OUTPUT_SUFFIX       Suffix to append to output FASTQ filenames (e.g. _trimmed)"
    echo ""
}

# Validate command-line arguments.
if [ "$#" -ne 2 ]; then
    help
    exit 1
fi

if [ "$1" != "--config" ]; then
    help
    exit 1
fi

CONFIG_FILE="$2"
if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Config file '$CONFIG_FILE' not found."
    exit 1
fi

echo "Sourcing configuration from $CONFIG_FILE..."
# shellcheck source=/home/tmurphy/phd_work/FANS_PILOT_DEV/fastp_config_tom.sh
source "$CONFIG_FILE"

# Validate that required environment variables are set.
if [ -z "$FASTQ_DIR" ] || [ -z "$FASTP_OUTPUT_DIR" ] || [ -z "$FASTP_REPORT_DIR" ] || [ -z "$THREADS" ] || [ -z "$OUTPUT_SUFFIX" ]; then
    echo "Error: Config file must export FASTQ_DIR, FASTP_OUTPUT_DIR, FASTP_REPORT_DIR, THREADS, and OUTPUT_SUFFIX."
    exit 1
fi

# Check that the input FASTQ directory exists.
if [ ! -d "$FASTQ_DIR" ]; then
    echo "Error: FASTQ_DIR '$FASTQ_DIR' does not exist."
    exit 1
fi

# Create the output directories 
# fastp processed fastq
if [ ! -d "$FASTP_OUTPUT_DIR" ]; then
    mkdir -p "$FASTP_OUTPUT_DIR"
fi
# fastp reports
if [ ! -d "$FASTP_REPORT_DIR" ]; then
    mkdir -p "$FASTP_REPORT_DIR"
fi
# fastp json reports
if [ ! -d "${FASTP_REPORT_DIR}/json" ]; then
    mkdir -p "${FASTP_REPORT_DIR}/json"
fi

# Pull the fastp Docker image.
FASTP_IMAGE="staphb/fastp:latest"
docker pull $FASTP_IMAGE

# Define a date prefix (yy-mm-dd) for our filenames.
date_prefix=$(date +%y-%m-%d)

# Process each paired FASTQ file in the input directory.
for r1 in "${FASTQ_DIR}"/*_R1*.fastq*; do
    # Construct the corresponding R2 filename.
    r2="${r1/_R1/_R2}"
    if [ ! -f "$r2" ]; then
        echo "Warning: Paired file for $r1 not found. Skipping."
        continue
    fi

    echo "Processing paired files:"
    echo "  R1: $r1"
    echo "  R2: $r2"
    # extract the base filename from path
    base1=$(basename "$r1")
    base2=$(basename "$r2")

    # Build output filenames by removing the .fastq* extension and appending OUTPUT_SUFFIX.
    out1="${FASTP_OUTPUT_DIR}/${base1%%.fastq*}${OUTPUT_SUFFIX}.fastq.gz"
    out2="${FASTP_OUTPUT_DIR}/${base2%%.fastq*}${OUTPUT_SUFFIX}.fastq.gz"

    # Build report filenames with the date prefix.
    html_report="${FASTP_REPORT_DIR}/${date_prefix}_${base1%%.fastq*}_fastp.html"
    json_report="${FASTP_REPORT_DIR}/json/${date_prefix}_${base1%%.fastq*}_fastp.json"

    docker run --rm \
      -u "$(id -u):$(id -g)" \
      -v "$(realpath "$FASTQ_DIR")":/data/input \
      -v "$(realpath "$FASTP_OUTPUT_DIR")":/data/fastq_output \
      -v "$(realpath "$FASTP_REPORT_DIR")":/data/reports \
      $FASTP_IMAGE fastp \
      -i /data/input/"$base1" -o /data/fastq_output/"$(basename "$out1")" \
      -I /data/input/"$base2" -O /data/fastq_output/"$(basename "$out2")" \
      -q 30 \
      -x \
      -p \
      -w "$THREADS" \
      --detect_adapter_for_pe \
      --html /data/reports/"$(basename "$html_report")" \
      --json /data/reports/json/"$(basename "$json_report")"
done

echo "fastp processing completed."