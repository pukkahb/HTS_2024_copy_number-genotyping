#!/bin/bash

# Directory and parameter definitions
BAM_DIR="$HOME/project/results/variants"
REFERENCE="$HOME/project/data/reference"
OUT_DIR="$HOME/project/results/cnv"
HIS_SIZE=10

# Check if the input and output directories exist
if [[ ! -d "$BAM_DIR" ]]; then
    echo "Input directory $BAM_DIR does not exist!"
    exit 1
fi

if [[ ! -d "$OUT_DIR" ]]; then
    mkdir -p "$OUT_DIR"
fi

# Loop through each BAM file and run CNVnator
for file in $BAM_DIR/*.sorted.bam
do
    if [[ -e "$file" ]]; then
        BASE_NAME=$(basename $file .sorted.bam)
        echo "Running CNV analysis on $file"
        
        # Clean up existing output files
        rm -f $OUT_DIR/$BASE_NAME.root
        rm -f $OUT_DIR/$BASE_NAME.cnv
        
        # Check if BAM file is indexed
        if [[ ! -e "$file.bai" ]]; then
            echo "Indexing BAM file: $file"
            samtools index $file
        fi
        
        # Run CNVnator commands
        cnvnator -root $OUT_DIR/$BASE_NAME.root -tree $file
        if [[ $? -ne 0 ]]; then
            echo "Error in cnvnator -tree command for $file"
            continue
        fi
        
        cnvnator -root $OUT_DIR/$BASE_NAME.root -his $HIS_SIZE -d $REFERENCE
        if [[ $? -ne 0 ]]; then
            echo "Error in cnvnator -his command for $file"
            continue
        fi
        
        cnvnator -root $OUT_DIR/$BASE_NAME.root -stat $HIS_SIZE
        if [[ $? -ne 0 ]]; then
            echo "Error in cnvnator -stat command for $file"
            continue
        fi
        
        cnvnator -root $OUT_DIR/$BASE_NAME.root -partition $HIS_SIZE
        if [[ $? -ne 0 ]]; then
            echo "Error in cnvnator -partition command for $file"
            continue
        fi
        
        cnvnator -root $OUT_DIR/$BASE_NAME.root -call $HIS_SIZE > $OUT_DIR/$BASE_NAME.cnv
        if [[ $? -ne 0 ]]; then
            echo "Error in cnvnator -call command for $file"
            continue
        fi
        
        echo "CNV analysis completed for $file"
    else
        echo "No matching files found in $BAM_DIR"
    fi
done
