#!/bin/bash

help(){
    echo "This script provdes a wrapper for Salmon for the quantification Illumnia paired-end RNA-seq data."
    echo ""
    echo "Usage: $0 --config <config_file>"
    echo ""
    echo "The config file must export the following variables:"
    echo "  FASTQ_DIR           Path to raw FASTQ input directory"
    echo "PATH_TO_SALMON_INDEX  Path to the Salmon index"  
    
    exit 1
}



validate_args() {

if [[ "$#" -ne 2 ]]; then
    echo ""
    echo "Usage: $0 --config <config_file>"
    exit 1
fi

if [[ "$1" != "--config" ]]; then
    echo ""
    echo "Usage: $0 --config <config_file>"
    exit 1
fi
}

validate_config() {
    
    CONFIG_FILE="$2"
    if [[ ! -f "$CONFIG_FILE" ]]; then
        echo "Error: Config file '$CONFIG_FILE' not found."
        exit 1
    else
        echo "Sourcing configuration from $CONFIG_FILE..."
        source "$CONFIG_FILE"
    fi
}

