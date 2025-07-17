#!/bin/bash
#SBATCH --job-name=FLAMES  # Job name
#SBATCH --output=../logs/FLAMES_%j.log  # Output log file
#SBATCH --error=../logs/FLAMES_%j.err  # Error log file
#SBATCH --mem=24G
#SBATCH --gres=gpu

#Note: This script has to be run under appropriate environment

# Checking if correct number of arguments are passed
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_directory> <output_directory>"
    exit 1
fi

# Assign arguments to variables
INPUT_DIR="$1"
OUTPUT_DIR="$2"
REFERENCE_PATH="/scratch/prj/cmm_stroud/Nanopore_DirectRNA/ensembl_ref"
FLAMES_PATH="/scratch/prj/cmm_stroud/FLAMES"

# Creating output directory if not created already
mkdir -p "$OUTPUT_DIR"

# Running the Python script with dynamic paths
"$FLAMES_PATH/python/bulk_long_pipeline.py" \
    --gff3 "$REFERENCE_PATH/Mus_musculus.GRCm39.113.gtf" \
    --genomefa "$REFERENCE_PATH/Mus_musculus.GRCm39.dna.primary_assembly.fa" \
    --outdir "$OUTPUT_DIR" \
    --config_file "$FLAMES_PATH/examples/SIRV/data/SIRV_config.json" \
    --fq_dir "$INPUT_DIR"
