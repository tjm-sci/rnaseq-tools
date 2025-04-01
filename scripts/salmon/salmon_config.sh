#!/bin/bash
# Description: Configuration file for running Salmon quantification wrapper script

# Directory containing your FASTQ files (the script will auto-detect paired-end files)
export FASTQ_DIR="/path/to/fastq_files"

# Output directory for Salmon quantification results
export OUTPUT_DIR="/path/to/output_directory"

# Path to your Salmon index (a directory created when building the index)
export INDEX="/path/to/salmon_index"

# Number of threads to use (optional; defaults to 18 if not set)
export THREADS=18
