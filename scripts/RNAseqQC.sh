#!/bin/bash
# fastqc_multiqc_wrapper.sh
# This script runs FastQC in a Docker container on paired-end FASTQ files
# and then runs MultiQC to aggregate the reports.
# If RUN_MULTI_QC_ONLY is set to "true" in the config file, it skips the FastQC run
# and only generates a new MultiQC report based solely on the directories specified
# in MULTIQC_SEARCH_DIRS.
#
# The config file (supplied via --config <config_file>) must export:
#   FASTQ_DIR, FASTQC_REPORT_DIR, MULTIQC_SEARCH_DIR, MULTIQC_OUTPUT_DIR,
#   RUN_MULTI_QC_ONLY, and MULTIQC_REPORT_TITLE.
#
# Usage: ./fastqc_multiqc_wrapper.sh --config <config_file>

help() {
    echo ""
    echo "Usage: $0 --config <config_file>"
    echo ""
    echo "The config file must export the following variables:"
    echo "  FASTQ_DIR           Path to raw FASTQ input directory (ignored if RUN_MULTI_QC_ONLY is \"true\")"
    echo "  FASTQC_REPORT_DIR   Directory for FastQC reports (ignored if RUN_MULTI_QC_ONLY is \"true\")"
    echo "  MULTIQC_SEARCH_DIR  Space-separated list of directories for MultiQC to search"
    echo "                      (if empty, defaults to FASTQC_REPORT_DIR)"
    echo "  MULTIQC_OUTPUT_DIR  Directory for the final MultiQC report"
    echo "  RUN_MULTI_QC_ONLY   Set to \"true\" to skip FastQC run and only run MultiQC."
    echo "                      If empty or \"false\", the normal mode is used."
    echo "  MULTIQC_REPORT_TITLE  Title for the MultiQC report (a date prefix will be added)"
    echo ""
}

# this stops the script if any command fails 
set -euo pipefail

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
source "$CONFIG_FILE"

# Ensure RUN_MULTI_QC_ONLY is set; if empty, default to "false"
if [ -z "$RUN_MULTI_QC_ONLY" ]; then
    RUN_MULTI_QC_ONLY="false"
fi
# Shellcheck source=/home/tmurphy/phd_work/FANS_PILOT_DEV/RNAseqQC_config_tom.sh
# If MULTIQC_SEARCH_DIR is empty, default to FASTQC_REPORT_DIR.
if [ -z "$MULTIQC_SEARCH_DIRS" ]; then
    MULTIQC_SEARCH_DIRS="$FASTQC_REPORT_DIR"
    echo "MULTIQC_SEARCH_DIRS not set"
    echo "If a fastqc run is specified, $MULTIQC_SEARCH_DIRS will be used"
fi

# Create the MultiQC output directory if it doesn't exist.
mkdir -p "$MULTIQC_OUTPUT_DIR"

if [ "$RUN_MULTI_QC_ONLY" = "true" ]; then
    echo ""
    echo "RUN_MULTI_QC_ONLY is true; skipping FastQC processing."
    echo ""
else
    # Normal mode: FastQC processing is required.
    if [ -z "$FASTQ_DIR" ] || [ -z "$FASTQC_REPORT_DIR" ]; then
        echo ""
        echo "Error: When RUN_MULTI_QC_ONLY is not true, FASTQ_DIR and FASTQC_REPORT_DIR must be specified."
        echo ""
        exit 1
    fi
    if [ ! -d "$FASTQ_DIR" ]; then
        echo ""
        echo "Error: FASTQ_DIR '$FASTQ_DIR' does not exist."
        echo ""
        exit 1
    fi
    mkdir -p "$FASTQC_REPORT_DIR"

    # Pull the FastQC Docker image.
    FASTQC_DOCKER_IMAGE="staphb/fastqc:latest"
    docker pull $FASTQC_DOCKER_IMAGE

    echo "Running FastQC on paired FASTQ files..."
    for r1 in "${FASTQ_DIR}"/*_R1*.fastq*; do
        r2="${r1/_R1/_R2}"
        if [ ! -f "$r2" ]; then
            echo ""
            echo "Warning: Paired file for $r1 not found. Skipping."
            echo ""
            continue
        fi
        echo ""
        echo "Processing paired files:"
        echo "  R1: $r1"
        echo "  R2: $r2"
        echo ""

        base1=$(basename "$r1")
        base2=$(basename "$r2")

        docker run --rm \
          -u "$(id -u):$(id -g)" \
          -v "$(realpath "$FASTQ_DIR")":/data/input \
          -v "$(realpath "$FASTQC_REPORT_DIR")":/data/output \
          $FASTQC_DOCKER_IMAGE \
          fastqc /data/input/"$base1" /data/input/"$base2" -t $THREADS --outdir /data/output 
    done
    echo ""
    echo "FastQC processing completed."
    echo ""
fi

# Build volume mounts for MultiQC.
# initialise two empty arrays for the volumes and local directories 
# and the correpsonding mount points in the container. 
volumes=()
mount_points=()
i=1 # used to number the directories in the container
for dir in $MULTIQC_SEARCH_DIRS; do
    if [ ! -d "$dir" ]; then
        echo ""
        echo "Warning: MultiQC search directory '$dir' does not exist. Skipping."
        echo ""
        continue
    fi
    real_dir=$(realpath "$dir")
    mount_point="/data/dir$i"
    volumes+=("-v" "$real_dir:$mount_point") # this builds multiple bind mounts, one for each place multiQC should check.
    mount_points+=("$mount_point")
    i=$((i+1))
done

if [ ${#mount_points[@]} -eq 0 ]; then
    echo ""
    echo "Error: No valid directories provided for MultiQC search."
    echo ""
    exit 1
fi

# Pull the MultiQC Docker image.
MULTIQC_DOCKER_IMAGE="multiqc/multiqc"
docker pull $MULTIQC_DOCKER_IMAGE

# Generate a report name (and title) with a date prefix.
report_title="$(date +%y-%m-%d)_${MULTIQC_REPORT_TITLE}.html"

# Build the MultiQC command. The loop appends each mount point (directory) to the cmd array.
# This passes all user-specified directories to MultiQC.
cmd=(multiqc -n "$report_title")
for mp in "${mount_points[@]}"; do
    cmd+=("$mp")
done

echo "Running MultiQC..."
# -u "$(id -u):$(id -g)" ensures that the output files are owned by the current user
# "${volumes[@]}" expands the volumes array into separate -v arguments
# "${cmd[@]}" expands to creat a command with each mounted directory as an argument

docker run --rm \
  -u "$(id -u):$(id -g)" \
  "${volumes[@]}" \
  -v "$(realpath "$MULTIQC_OUTPUT_DIR")":/data_output \
  $MULTIQC_DOCKER_IMAGE \
  "${cmd[@]}" -o /data_output

echo ""
echo "MultiQC report generated in $MULTIQC_OUTPUT_DIR"
