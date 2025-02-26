#!/bin/bash

# Usage function
usage() {
    echo "This is a wrapper script for fastp processing of paired end Illumina reads."
    echo "The script takes an input directory containing FASTQ files and runs fastp on each pair."
    echo "Results from fastp are stored in a user-specified output directory."
    echo ""
    echo "fastp is run in a Docker container."
    echo "Warning: Docker must be installed and set up correctly to run this program."
    echo ""
    echo "Usage: $0 --fastq <input_fastq_directory> --output <output_directory> --threads <number_of_threads>"
    echo
    
}

# Check for no arguments
if [ "$#" -eq 0 ]; then
    usage
    exit 1
fi

# Initialize variables
fastq_dir=""
output_dir=""
threads=""

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --fastq|-f)
            fastq_dir="$2"
            shift 2
            ;;
        --output|-o)
            output_dir="$2"
            shift 2
            ;;
        --threads |-t)
            threads="$2"
            shift 2
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# Validate required arguments
if [ -z "$fastq_dir" ] || [ -z "$output_dir" ] || [ -z "$threads" ]; then
    echo "Error: --fastq, --output, and --threads must all be specified."
    usage
    exit 1
fi

# Check that input directory exists
if [ ! -d "$fastq_dir" ]; then
    echo "Error: Input directory '$fastq_dir' does not exist."
    exit 1
fi

# Create output directory if it doesn't exist
if [ ! -d "$output_dir" ]; then
    echo "Output directory '$output_dir' does not exist."
    echo "Creating output directory..."
    mkdir -p "$output_dir"
fi

# Loop through FASTQ files in the input directory matching a flexible _R1 pattern
for r1 in "${fastq_dir}"/*_R1*.fastq*; do
    # Construct the corresponding R2 filename by replacing the first occurrence of _R1 with _R2
    r2="${r1/_R1/_R2}"
    
    # Check if the paired R2 file exists; if not, skip this pair
    if [ ! -f "$r2" ]; then
        echo "Warning: Paired file for $r1 not found. Skipping."
        continue
    fi

    echo "Processing paired files:"
    echo "  R1: $r1"
    echo "  R2: $r2"

    # Extract just the base filenames (so that inside the container we refer to /data/input/<filename>)
    base1=$(basename "$r1")
    base2=$(basename "$r2")

    # Create output filenames by removing the .fastq* extension and appending _trimmed.fastq.gz.
    out1="${output_dir}/${base1%%.fastq*}_trimmed.fastq.gz"
    out2="${output_dir}/${base2%%.fastq*}_trimmed.fastq.gz"

    # Set report names (HTML and JSON) based on the R1 file name.
    html_report="${output_dir}/${base1%%.fastq*}_fastp.html"


    # Run fastp in Docker.
    FASTP_IMAGE="staphb/fastp:latest"
    # - Mount the input directory to /data/input and the output directory to /data/output.
    # - Run as the current user to avoid permission issues.
    # - adapters should be automatically detected.
    # - Pass the fastp options:
    #    - -i and -I point to the input files (inside /data/input)
    #    - -o and -O are the output files (inside /data/output)
    #    - -q 30 sets the qualified quality threshold to Q30
    #    - -x enables polyX trimming
    #    - -w uses the user-specified number of threads
    #    - --html generates report in the output directory
    docker run --rm \
      -u "$(id -u):$(id -g)" \
      -v "$(realpath "$fastq_dir")":/data/input \
      -v "$(realpath "$output_dir")":/data/output \
      $FASTP_IMAGE fastp \
      -i /data/input/"$base1" -o /data/output/"$(basename "$out1")" \
      -I /data/input/"$base2" -O /data/output/"$(basename "$out2")" \
      -q 30 \
      -x \
      -w "$threads" \
      --html /data/output/"$(basename "$html_report")"
      

done

echo "fastp processing completed."
