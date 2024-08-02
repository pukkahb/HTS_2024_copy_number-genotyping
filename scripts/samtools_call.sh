#!/bin/bash

# Directory and parameter definitions
SAM_DIR="$HOME/project/results/aligned"
OUT_DIR="$HOME/project/results/variants"

# Check if the input and output directories exist
if [[ ! -d "$SAM_DIR" ]]; then
    echo "Input directory $SAM_DIR does not exist!"
    exit 1
fi

if [[ ! -d "$OUT_DIR" ]]; then
    mkdir -p "$OUT_DIR"
fi

# Loop through each SAM file and run Samtools
for file in $SAM_DIR/*.sam
do
    if [[ -e "$file" ]]; then
        BASE_NAME=$(basename $file .sam)
        echo "Converting $file to BAM"
        
        # Convert SAM to BAM
        samtools view -S -b $file -o $OUT_DIR/$BASE_NAME.bam
        if [[ $? -ne 0 ]]; then
            echo "Failed to convert $file to BAM"
            continue
        fi
        
        echo "Sorting $OUT_DIR/$BASE_NAME.bam"
        
        # Sort BAM file
        samtools sort -o $OUT_DIR/$BASE_NAME.sorted.bam $OUT_DIR/$BASE_NAME.bam
        if [[ $? -ne 0 ]]; then
            echo "Failed to sort $OUT_DIR/$BASE_NAME.bam"
            continue
        fi
        
        echo "Indexing $OUT_DIR/$BASE_NAME.sorted.bam"
        
        # Index BAM file
        samtools index $OUT_DIR/$BASE_NAME.sorted.bam
        if [[ $? -ne 0 ]]; then
            echo "Failed to index $OUT_DIR/$BASE_NAME.sorted.bam"
            continue
        fi
        
        # Clean up intermediate BAM file
        rm $OUT_DIR/$BASE_NAME.bam
    else
        echo "No matching files found in $SAM_DIR"
    fi
done

echo "SAM to BAM conversion, sorting, and indexing process completed"
