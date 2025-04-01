#!/bin/bash
# This script defines the input and output options for the fastp_wrapper.sh script
# The user must set the following variables.
# Only the FASTQ_DIR MUST exist, the other directories will be created if they do not exist.

# Path to directory containing FASTQ files 
export FASTQ_DIR="/path/to/raw_fastq"

# Path to output directory for trimmed FASTQ files
export FASTP_OUTPUT_DIR="/path/to/trimmed_fastq"    

# Path to directory for fastp reports
export FASTP_REPORT_DIR="/path/to/reports"     

# Number of threads to use (max 16)
export THREADS=4

# Suffix to append to output FASTQ files. underscore not automatically added
export OUTPUT_SUFFIX="_trimmed"