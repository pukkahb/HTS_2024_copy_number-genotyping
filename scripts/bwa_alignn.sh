#!/bin/bash

# Directory and parameter definitions
FASTQ_DIR="$HOME/project/results/trimmed"
REFERENCE_DIR="$HOME/project/results/indexed"
OUT_DIR="$HOME/project/results/aligned"

# Check if the input and output directories exist
if [[ ! -d "$FASTQ_DIR" ]]; then
    echo "Input directory $FASTQ_DIR does not exist!"
    exit 1
fi

if [[ ! -d "$REFERENCE_DIR" ]]; then
    echo "Reference directory $REFERENCE_DIR does not exist!"
    exit 1
fi

if [[ ! -d "$OUT_DIR" ]]; then
    mkdir -p "$OUT_DIR"
fi

# Loop through each file and run BWA mem
for file in $FASTQ_DIR/*.{fastq,fastq.gz,fq.gz}
do
    if [[ -e "$file" ]]; then
        BASE_NAME=$(basename $file .fq.gz)
        echo "Running bwa mem on $file"
        bwa mem $REFERENCE_DIR/$(basename $REFERENCE_DIR) $file > $OUT_DIR/$BASE_NAME.sam
    else
        echo "No matching files found in $FASTQ_DIR"
    fi
done

