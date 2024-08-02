#!/bin/bash

# Directory and parameter definitions
INPUT_DIR="$HOME/project/data/raw" # Input directory containing raw FASTQ files
OUTPUT_DIR="$HOME/project/results/trimmed" # Output directory for trimmed FASTQ files
THREADS=2 # Number of threads to use

# Check if the input and output directories exist
if [[ ! -d "$INPUT_DIR" ]]; then
    echo "Input directory $INPUT_DIR does not exist!"
    exit 1
fi

if [[ ! -d "$OUTPUT_DIR" ]]; then
    mkdir -p "$OUTPUT_DIR"
fi

# Find all FASTQ files
fastq_files=($INPUT_DIR/*.fastq.gz)

# Check if there are files in the directory
if [ ${#fastq_files[@]} -eq 0 ]; then
    echo "No FASTQ files found in $INPUT_DIR"
    exit 1
fi

# Process files for paired-end or single-end
if [ ${#fastq_files[@]} -gt 1 ]; then
    # Sort files to handle paired-end files
    IFS=$'\n' sorted_files=($(sort <<<"${fastq_files[*]}"))
    unset IFS

    # Check for paired-end files
    if [ $((${#sorted_files[@]} % 2)) -ne 0 ]; then
        echo "Please provide an even number of input files for paired-end FastQ trimming!"
        exit 1
    fi

    echo "Running trim_galore on paired-end files"
    for ((i=0; i<${#sorted_files[@]}; i+=2))
    do
        file1=${sorted_files[$i]}
        file2=${sorted_files[$i+1]}
        
        base_name1=$(basename "$file1" "_1.fastq.gz")
        base_name2=$(basename "$file2" "_2.fastq.gz")
        
        if [[ "$base_name1" != "$base_name2" ]]; then
            echo "Paired files do not match: $file1 and $file2"
            exit 1
        fi
        
        echo "Trimming paired-end files $file1 and $file2"
        
        trim_galore --paired -o "$OUTPUT_DIR" --cores "$THREADS" "$file1" "$file2"
    done
else
    echo "Running trim_galore on single-end files"
    for file in "${fastq_files[@]}"
    do
        base_name=$(basename "$file" ".fastq.gz")
        echo "Trimming single-end file $file"
        
        trim_galore -o "$OUTPUT_DIR" --cores "$THREADS" "$file"
    done
fi

echo "Trimming process completed"
