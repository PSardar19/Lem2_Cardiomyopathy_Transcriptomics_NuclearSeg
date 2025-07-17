#!/bin/bash
#SBATCH --job-name=nanoseq_run
#SBATCH --output=logs/nanoseq_%j.out
#SBATCH --error=logs/nanoseq_%j.err
#SBATCH --mem=64G
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --partition=interruptible_gpu

module load nextflow/24.10.2-gcc-13.2.0
module load python/3.11.6-gcc-13.2.0

# Set Singularity cache and work directories to local paths
export NXF_SINGULARITY_CACHEDIR=/scratch/prj/cmm_stroud/NanoSeq/singularity_cache
export NXF_WORK=/scratch/prj/cmm_stroud/NanoSeq/work

# Make sure the folders exist
mkdir -p $NXF_SINGULARITY_CACHEDIR
mkdir -p $NXF_WORK

# Run Nextflow
nextflow run nf-core/nanoseq \
  --input metadata/master_sample_sheet.csv \
  --protocol directRNA \
  --outdir outputs/master_outputs \
  --skip_demultiplexing \
  --aligner minimap2 \
  --quantification_method bambu \
  --skip_modification_analysis \
  --skip_xpore \
  --skip_m6anet \
  --skip_fusion_analysis \
  --max_cpus 16 \
  --gpu_device auto \
  --gpu_cluster_options '--partition=interruptible_gpu --gres=gpu:1' \
  -w $NXF_WORK \
  -profile singularity \
  -resume


