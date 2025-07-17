# Direct RNA Sequencing Analysis

This module contains the complete pipeline and analysis scripts for processing and analyzing **Nanopore direct RNA-seq (dRNA-seq)** data in the context of Lem2-associated dilated cardiomyopathy (DCM) in a murine model.

It includes:
- Quantification workflows using **FLAMES** and **NanoSeq**
- Differential analysis (DGE, DTE, DTU)
- Gene set enrichment analysis (GSEA) and overrepresentation analysis (ORA)

---

## 🧬 Overview

This pipeline performs:

- **Base-level quantification** with two complementary tools:
  - [`FLAMES`](https://github.com/LuyiTian/FLAMES): Full-length transcriptome splicing and mutation analysis, R-based framework
  - [`NanoSeq`](https://github.com/nf-core/nanoseq): Nextflow pipeline for efficient transcript-level quantification for dRNA-seq

- **Downstream analyses** including:
  - DGE: Differential Gene Expression (DESeq2)
  - DTE: Differential Transcript Expression (DESeq2)
  - DTU: Differential Transcript Usage (IsoformSwitchAnalyzeR)
  - GSEA and ORA: Functional enrichment (fgsea and clusterProfiler)
 
---

## ⚙️ Reproducibility & Setup

### 1. Install Environments

Use the provided conda environment YAML files:

```bash
conda env create -f Pipelines/FLAMES/flames_environment.yml
conda env create -f Pipelines/NanoSeq/nanoseq_environment.yml
```

### 2. Run Preprocessing (if needed)

These optional preprocessing steps help prepare the raw sequencing data before quantification:

```bash
sbatch Pipelines/Preprocessing/merge_fastq.sh
sbatch Pipelines/Preprocessing/trim_reads.sh
```
merge_fastq.sh: Merges multiple .fastq files corresponding to the same biological sample, which may have been sequenced in separate batches or lanes.
trim_reads.sh: Performs poly-A tail trimming on raw reads before downstream quantification, improving mapping and transcript boundary resolution.

ℹ️ These scripts are tailored for use on a SLURM-based HPC system. If running locally, adapt them accordingly (e.g., ./scriptname.sh).

### 3. Run Quantification Pipelines

Before running the FLAMES pipeline, clone the official FLAMES GitHub repository:
```bash
git clone https://github.com/LuyiTian/FLAMES
```
Refer to the bash scripts in:

1. Pipelines/FLAMES/scripts/
2. Pipelines/NanoSeq/scripts/

Please update paths to reflect your local directory and the location of the cloned FLAMES repo.

NOTE:
1. NanoSeq stores the gene-level counts and transcript-level counts in output_fol/bambu
2. FLAMES stores the count files in the specific sample folder. To get the collated counts matrix, run extract_counts_file.sh -> merge_counts.py in the FLAMES/scripts

### 4. Downstream Analysis

1. Open the .qmd or .R scripts in the Analysis/ folder in RStudio.
2. Modify file paths or parameters as needed.
3. Knit or run interactively.

Note: The GSEA and ORA analysis scripts can only be run after DGE (DGE_DESeq2_master.qmd)

## Additional Notes:
1. Data not included: Raw fastq files, GFF/GTF annotations, and counts tables are excluded.
2. Alignment and transcript quantification were performed against the GRCm39 primary assembly (ENSEMBL) using the Mus_musculus.GRCm39.113.gtf annotation (ENSEMBL release 113).
3. FLAMES configuration is stored in metadata/FLAMES_config.json.
4. Both Quantification pipelines were run on a Linux HPC environment (HPC CREATE)
