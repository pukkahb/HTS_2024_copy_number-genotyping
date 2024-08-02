#!/bin/bash

# Directory and parameter definitions
FASTQC_DIR="$HOME/project/results/fastqc"
MULTIQC_OUT="$HOME/project/results/multiqc"
#PREFIX="$3"

# Check if the input and output directories exist
if [[ ! -d "$FASTQC_DIR" ]]; then
    echo "Input directory $FASTQC_DIR does not exist!"
    exit 1
fi

if [[ ! -d "$MULTIQC_OUT" ]]; then
    mkdir -p "$MULTIQC_OUT"
fi

# Run MultiQC to aggregate FastQC results
echo "Running MultiQC to aggregate fastqc output"
multiqc $FASTQC_DIR/ -o $MULTIQC_OUT

echo "MultiQC process completed"

