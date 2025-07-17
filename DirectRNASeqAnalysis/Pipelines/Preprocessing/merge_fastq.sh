#!/bin/bash
#SBATCH --job-name=merge_fastq  # Job name
#SBATCH --output=../logs/merge_fastq_%j.log  # Output log file
#SBATCH --error=../logs/merge_fastq_%j.err  # Error log file
#SBATCH --mem=8G  # Memory allocation
#SBATCH --cpus-per-task=1  # Number of CPU cores

# Check if user provided both input and output directories
if [ $# -ne 2 ]; then
    echo "Usage: $0 <parent_directory> <output_directory>"
    exit 1
fi

PARENT_DIR=$1
OUTPUT_DIR=$2

# Check if parent directory exists
if [ ! -d "$PARENT_DIR" ]; then
    echo "Error: Parent directory '$PARENT_DIR' does not exist."
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all subfolders in the parent directory
for folder_path in "$PARENT_DIR"/*/; do
    # Remove trailing slash and extract folder name
    folder_name=$(basename "$folder_path")

    # Check if there are any fastq.gz files in the subfolder
    if compgen -G "${folder_path}"*.fastq.gz > /dev/null; then
        echo "Merging files in $folder_name..."
        ls "${folder_path}"*.fastq.gz | sort -V | xargs cat > "$OUTPUT_DIR/${folder_name}_merged.fastq.gz"
    else
        echo "No FASTQ files found in $folder_name, skipping..."
    fi
done
