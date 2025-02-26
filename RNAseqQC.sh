#!/bin/bash
# This script takes an input directory containing FASTQ files and runs FastQC on each file.
# Results from FastQC are then stored in a user-specified output directory.
# A MutliQC report is generated from the FastQC results and stored in the output directory.

help(){
    echo ""
    echo "Wrapper script for running FastQC on fastq.gz files genrated by Illumina paired-end seqeuncing"
    echo "Generates a MultiQC report from the FastQC results"
    echo ""
    echo "Warning: Docker must be installed and set up correctly to run this script"
    echo ""
    echo "Usage: $0 --fastq <fastq_directory> --output <output_directory>"
    echo ""
    echo "Options:"
    echo "  --fastq         Path to directory containing FASTQ files (paired-end; _R1.fastq and _R2.fastq)"
    echo "  --output, -o    Path to directory where QC reports will be stored"
    echo "  --help, -h      Display this help message and exit"
}

# if no parameters are provided to the command line, print the help dialgoue
if [ $# -eq 0 ]; then
  help
  exit 1
fi

# initialise variables
fastq_dir=""
output_dir=""

# Parse command-line arguments:
# iterate over argument whilst there's more than 0 arguments remaining
# when the argument found by in the current iteration of the while loop
# is one of the options, assign the value of the next argument to the
# appropriate variable. If the argument is not recognised, print the help.
# shift 2 moves the indexes two places left effectively move to the next argument.

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --fastq)
            fastq_dir="$2"
            echo "fastq_dir: $fastq_dir"
            shift 2
            ;;
        --output|-o)
            output_dir="$2"
            echo "output_dir: $output_dir"
            shift 2
            ;;
        --help|-h)
            help
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            help
            exit 1
            ;;
    esac
done

# Check that both required directories have been provided.
# The user must provide both, but it doesn't matter if the output exists or not.
if [ -z "$fastq_dir" ] || [ -z "$output_dir" ]; then
    echo "Error: Both --fastq and --output options must be specified."
    usage
    exit 1
fi

# Check the input directory exists on the filesystem.
if [ ! -d "$fastq_dir" ]; then
    echo "Error: FASTQ directory '$fastq_dir' does not exist."
    exit 1
fi

# Check the output directory exists on the filesystem.
# if it doesn't, create it
if [ ! -d "$output_dir" ]; then
    echo "Output directory '$output_dir' does not exist."
    echo "Creating output directory..."
    mkdir -p "$output_dir"
fi

# Running fastqc using a docker image.
# Define the Docker image name

FASTQC_DOCKER_IMAGE="staphb/fastqc:latest"

# pull image from dockerhub
docker pull $FASTQC_DOCKER_IMAGE

# Loop through FASTQ files in the input directory matching the flexible _R1 pattern.

for r1 in "${fastq_dir}"/*_R1*.fastq*; do

    # Construct the corresponding R2 filename by replacing the first occurrence of _R1 with _R2.
    r2="${r1/_R1/_R2}"

    # Check if the paired R2 file exists.
    if [ -f "$r2" ]; then
        echo "Processing paired files:"
        echo "  Read 1: $r1"
        echo "  Read 2: $r2"

        # Extract only the base filenames for use inside the container.
        r1_base=$(basename "$r1")
        r2_base=$(basename "$r2")

        # Run FastQC using Docker:
        docker run --rm \
          -u "$(id -u):$(id -g)" \  # ensures output files are owned by the current user
          -v "$(realpath "$fastq_dir")":/data/input \
          -v "$(realpath "$output_dir")":/data/output \
          $FASTQC_DOCKER_IMAGE \
          fastqc /data/input/"$r1_base" /data/input/"$r2_base" --outdir /data/output
    else
        echo "Warning: Paired file for $r1 not found. Skipping."
    fi
done

# Generate a MultiQC report from the FastQC results.
# Define the MultiQC Docker image name.

MULTIQC_DOCKER_IMAGE="multiqc/multiqc"
docker pull $MULTIQC_DOCKER_IMAGE

# Run MultiQC using Docker:
# Generate a MultiQC report from the FastQC outputs
report_name="$(date +%Y%m%d)_mutliqc_report.html"

docker run --rm -u "$(id -u):$(id -g)" \
  -v "$(realpath "$output_dir")":/data \
  $MULTIQC_DOCKER_IMAGE \
  multiqc -n "$report_name" /data