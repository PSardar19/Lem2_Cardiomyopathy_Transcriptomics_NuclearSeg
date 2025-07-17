#!/bin/bash
#SBATCH --job-name=trim_fastq  # Job name
#SBATCH --output=../logs/trim_fastq_%j.log  # Output log file
#SBATCH --error=../logs/trim_fastq_%j.err  # Error log file
#SBATCH --mem=24G
#SBATCH --gres=gpu

# Check if user provided both input and output directories
if [ $# -ne 2 ]; then
    echo "Usage: $0 <input_fastq_directory> <output_directory>"
    exit 1
fi

# Load cutadapt module
module load py-cutadapt/4.4-gcc-13.2.0-python-3.11.6

# Set input and output directories
INPUT_DIR=$1
OUTPUT_DIR=$2

# Create output directiry if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through each fastq.gz file in the input directory
for read_file in "$INPUT_DIR"/*.fastq.gz; do

    # Extract filename without extension
    filename=$(basename "${read_file%_merged.fastq.gz}")
    tmp_output_dir="$OUTPUT_DIR/${filename}"

    mkdir -p "$tmp_output_dir"

    # Run cutadapt to trim the reads
    cutadapt --poly-a -o "$tmp_output_dir/${filename}_trimmed.fastq.gz" "$read_file"
    echo "Trimming complete for sample: ${filename}"


done

echo "All trimming complete!"
