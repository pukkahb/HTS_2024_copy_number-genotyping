#!/bin/bash

# Directory and parameter definitions
REFERENCE_DIR="$HOME/project/data/reference" # Directory containing the reference genome file
REFERENCE_FILE=$(find "$REFERENCE_DIR" -name "*.fasta" -print -quit) # Reference genome file
OUTPUT_DIR="$HOME/project/results/indexed" # Output directory for alignment files
THREADS=2 # Number of threads to use

# Check if the reference directory exist
if [[ ! -d "$REFERENCE_DIR" ]]; then
    echo "Input directory $REFERENCE_DIR does not exist!"
    exit 1
fi

# Check if the reference genome file exists
if [[ ! -f "$REFERENCE_FILE" ]]; then
    echo "Reference genome file $REFERENCE_FILE does not exist!"
    exit 1
fi

# Create the output directory if it does not exist
if [[ ! -d "$OUTPUT_DIR" ]]; then
    mkdir -p "$OUTPUT_DIR"
fi

# Run BWA index on the reference genome
echo "Running bwa index on $REFERENCE_FILE"
bwa index "$REFERENCE_FILE"
if [[ $? -ne 0 ]]; then
    echo "BWA indexing failed"
    exit 1
fi
echo "BWA indexing completed"
