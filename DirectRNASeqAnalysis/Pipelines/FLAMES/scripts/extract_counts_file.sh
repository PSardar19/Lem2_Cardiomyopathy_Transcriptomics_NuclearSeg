#!/bin/bash
#SBATCH --output=../logs/slurm_%j.out
#SBATCH --ntasks=2
#SBATCH --nodes=1
#SBATCH --mem=8G

# Check if correct number of arguments is provided
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <source_directory> <destination_directory>"
    exit 1
fi

# Assign input arguments to variables
SOURCE_DIR="$1"
DEST_DIR="$2"

# Create the destination directory if it doesn't exist
mkdir -p "$DEST_DIR"

# Loop through each sample folder in the source directory
for sample_folder in "$SOURCE_DIR"/*; do
    if [[ -d "$sample_folder" ]]; then  # Ensure it's a directory
        sample_name=$(basename "$sample_folder")  # Extract sample name
        input_file="$sample_folder/transcript_count.csv.gz"
        output_file="$DEST_DIR/${sample_name}_counts.csv"

        # Check if the compressed file exists
        if [[ -f "$input_file" ]]; then
            # Unzip only if the unzipped version doesn't already exist
            if [[ ! -f "$output_file" ]]; then
                echo "Unzipping and moving $input_file..."
                gunzip -c "$input_file" > "$output_file"
            else
                echo "File already processed: $output_file"
            fi
        else
            echo "File not found: $input_file"
        fi
    fi
done

echo "All files processed and moved to $DEST_DIR."
