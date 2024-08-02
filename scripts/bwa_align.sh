#!/bin/bash

# Directory and parameter definitions
REFERENCE_DIR="$HOME/project/data/reference" # Directory containing the reference genome file
TRIMMED_DIR="$HOME/project/results/trimmed" # Directory containing trimmed FASTQ files
OUTPUT_DIR="$HOME/project/results/aligned" # Output directory for alignment files
THREADS=2 # Number of threads to use

# Find the reference genome file in the REFERENCE_DIR
REFERENCE_FILE=$(find "$REFERENCE_DIR" -name "*.fasta" -print -quit)

# Check if the reference genome file exists
if [[ -z "$REFERENCE_FILE" ]]; then
    echo "Reference genome file not found in $REFERENCE_DIR"
    exit 1
fi

# Check if the trimmed reads directory exists
if [[ ! -d "$TRIMMED_DIR" ]]; then
    echo "Trimmed reads directory $TRIMMED_DIR does not exist!"
    exit 1
fi

# Create the output directory if it does not exist
if [[ ! -d "$OUTPUT_DIR" ]]; then
    mkdir -p "$OUTPUT_DIR"
fi

# Check if the reference genome has already been indexed
index_files=( "$REFERENCE_FILE.bwt" "$REFERENCE_FILE.pac" "$REFERENCE_FILE.ann" "$REFERENCE_FILE.sa" "$REFERENCE_FILE.fai" )
indexed=true
for index in "${index_files[@]}"; do
    if [[ ! -e "$index" ]]; then
        indexed=false
        break
    fi
done

if [[ "$indexed" == false ]]; then
    echo "Reference genome not indexed. Running bwa index on $REFERENCE_FILE"
    bwa index "$REFERENCE_FILE"
    if [[ $? -ne 0 ]]; then
        echo "BWA indexing failed"
        exit 1
    fi
    echo "BWA indexing completed"
else
    echo "Reference genome already indexed. Skipping indexing."
fi

# Loop through each FASTQ file and run bwa mem
for file1 in "$TRIMMED_DIR"/*_1_val_1.fq.gz
do
    # Check if it is paired-end read (ends with _1_val_1.fq.gz)
    if [[ "$file1" == *_1_val_1.fq.gz ]]; then
        file2="${file1/_1_val_1/_2_val_2}"
        if [[ -e "$file2" ]]; then
            base_name=$(basename "$file1" "_1_val_1.fq.gz")
            echo "Running bwa mem on paired-end files $file1 and $file2"
            
            bwa mem -t "$THREADS" "$REFERENCE_FILE" "$file1" "$file2" > "$OUTPUT_DIR/${base_name}.sam"
            if [[ $? -ne 0 ]]; then
                echo "BWA MEM failed for $file1 and $file2"
                exit 1
            fi
        else
            echo "Paired file for $file1 not found, skipping"
        fi
    fi
done

# Loop through each single-end FASTQ file and run bwa mem
for file in "$TRIMMED_DIR"/*_val_2.fq.gz
do
    # Check if it is a single-end read (ends with _val_2.fq.gz)
    if [[ "$file" == *_val_2.fq.gz ]]; then
        base_name=$(basename "$file" "_val_2.fq.gz")
        echo "Running bwa mem on single-end file $file"
        
        bwa mem -t "$THREADS" "$REFERENCE_FILE" "$file" > "$OUTPUT_DIR/${base_name}.sam"
        if [[ $? -ne 0 ]]; then
            echo "BWA MEM failed for $file"
            exit 1
        fi
    fi
done

echo "BWA MEM process completed"
