#!/bin/bash

# Directory and parameter definitions
WORKDIR="$HOME/project"
RAW_DIR="$WORKDIR/data/raw"
OUT_DIR="$WORKDIR/results/fastqc"
THREADS=2

# Check if the input and output directories exist
if [[ ! -d "$INPUT_DIR" ]]; then
    echo "Input directory $INPUT_DIR does not exist!"
    exit 1
fi

if [[ ! -d "$OUTPUT_DIR" ]]; then
    mkdir -p "$OUTPUT_DIR"
fi

# Loop through each fastq || fastq.gz file and run fastqc
for file in $RAW_DIR/*.{fastq,fastq.gz}
do
    if [[ -e "$file" ]]; then
        echo "Running fastqc on $file"
        fastqc -o $OUT_DIR $file -t $THREADS
    else
        echo "No matching files found in $RAW_DIR"
    fi
done

