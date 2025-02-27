#!/bin/bash
# This script defines the input and output options for the RNAseqQC.sh script
# The user must set the following variables.
# Only the FASTQ_DIR MUST exist, the other directories will be created if they do not exist.


# direcotry containing the raw fastq files
export FASTQ_DIR="/path/to/raw_fastq"

# directory for the fastqc html reports
export FASTQC_REPORT_DIR="/path/to/fastqc_reports"

# directories multiqc will search for reports to aggregate. 
# If left empty, FASTQC_REPORT_DIR is used
export MULTIQC_SEARCH_DIRS="/path/to/fastqc_reports /another/path/to/reports"

# directory to save mutliqc output
export MULTIQC_OUTPUT_DIR="/path/to/multiqc_output"

# whether to skip fastqc and run multiqc only
# set to false by default, set to true to skip fastqc.
# this will overide FASTQ_DIR and FASTQC_REPORT_DIR
export RUN_MULTI_QC_ONLY="false"

# title for the multiqc report
export MULTIQC_REPORT_TITLE="My Project MultiQC Report"