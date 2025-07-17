################################################################################
# Title: Over-Representation Analysis (ORA) Pipeline using clusterProfiler + msigdbr
# Author: Payel Sardar
#
# Description:
# This R script implements a reproducible pipeline for performing Over-Representation 
# Analysis (ORA) using gene sets from MSigDB (via `msigdbr`) and `clusterProfiler`. 
# The pipeline reads in differential gene expression results (CSV files), maps gene 
# identifiers to symbols if needed, selects significantly differentially expressed (DE) 
# genes based on user-defined thresholds, and performs ORA across multiple gene set 
# collections (Hallmark, KEGG, Reactome, GO).
#
# Supported species: Mus musculus
# Input: Differential expression CSV files with log2FoldChange and adjusted p-values
# Output: ORA results (CSV) saved by analysis type (stagewise vs temporal)
#
# Input Expectations:
# -------------------
# - Each CSV must contain:
#     - Gene ID column (default: "X", e.g., ENSEMBL or SYMBOL)
#     - log2FoldChange column
#     - padj (adjusted p-value) column
#
# Output Directory Structure:
# ---------------------------
# Outputs/ORA/
# ├── Stagewise/
# │   ├── *_ORA_combined.csv  <- ORA results for each contrast
# └── Temporal/
#     ├── *_ORA_combined.csv
#
# Required R Packages:
# --------------------
# dplyr, msigdbr, clusterProfiler, AnnotationDbi, org.Mm.eg.db
#
# Notes:
# ------
# - Gene sets are filtered to include only those with size between `min_size` and `max_size`.
# - The gene ID type (ENSEMBL or SYMBOL) is auto-detected based on the first 10 IDs.
# - Only significant pathways (FDR < 0.05) are saved in final results.
#
# Example Usage:
# --------------
# - Place your DE CSV files in the `DGE/` folder with appropriate naming based on the analysis type:
#     e.g., "StageWiseDGE_stage_P1.csv", "TemporalDGE_condition_WT.csv"
# - Run the script to perform ORA and save output to `Outputs/ORA/`
#
################################################################################

# Load required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(msigdbr)
  library(clusterProfiler)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
})

################################################################################

# 1. run_ora_analysis():
#    - Reads a differential expression file.
#    - Maps ENSEMBL to SYMBOL if required.
#    - Filters DE genes (padj < 0.05 and |log2FC| ≥ log2(1.2)).
#    - Loads gene sets from MSigDB using `msigdbr`, filtered by size.
#    - Performs ORA using `clusterProfiler::enricher()`.
#    - Combines and filters results (FDR < 0.05).
#    - Saves results to output directory.

run_ora_analysis <- function(
    csv_path,
    gene_id_col = "X",
    output_dir = NULL,
    kegg_version = "CP:KEGG_LEGACY",
    min_size = 5,
    max_size = 1000
) {
  # Prepare output
  if (!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }
  prefix <- tools::file_path_sans_ext(basename(csv_path))
  message("\nProcessing: ", prefix)
  
  # Read data and filter NAs
  df <- read.csv(csv_path, header = TRUE) %>%
    dplyr::filter(!is.na(log2FoldChange), !is.na(padj)) %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(gene_id_col)), .keep_all = TRUE)
  
  # Detect ID type
  sample_ids <- head(df[[gene_id_col]], 10)
  id_type <- if (all(grepl("^ENS", sample_ids))) "ENSEMBL" else "SYMBOL"
  
  # Map to SYMBOL if needed
  if (id_type == "ENSEMBL") {
    # Explicitly use AnnotationDbi::select
    mapped <- AnnotationDbi::select(
      org.Mm.eg.db,
      keys = df[[gene_id_col]],
      keytype = "ENSEMBL",
      columns = "SYMBOL"
    ) %>%
      dplyr::as_tibble() %>%
      dplyr::filter(!is.na(SYMBOL)) %>%
      dplyr::distinct(ENSEMBL, .keep_all = TRUE)
    
    df <- df %>%
      dplyr::inner_join(mapped, by = c(setNames("ENSEMBL", gene_id_col)))
    message("Mapped ", nrow(df), " genes to SYMBOL")
  } else {
    df$SYMBOL <- df[[gene_id_col]]
    message("Using provided SYMBOL identifiers")
  }
  
  # Define DE genes (criteria: padj < 0.05 & |log2FC| >= log2(1.2))
  de_genes <- df %>%
    dplyr::filter(padj < 0.05, abs(log2FoldChange) >= log2(1.2)) %>%
    dplyr::pull(SYMBOL) %>%
    unique()
  
  # Define background genes
  background_genes <- unique(df$SYMBOL)
  message("DE genes: ", length(de_genes), 
          " | Background genes: ", length(background_genes))
  
  # Function to get filtered gene sets
  get_genesets <- function(collection, subcollection = NULL) {
    gs_df <- msigdbr(
      species = "Mus musculus",
      category = collection,
      subcollection = subcollection
    ) %>%
      dplyr::as_tibble()
    
    # Filter by size using dplyr
    gs_counts <- gs_df %>%
      dplyr::group_by(gs_name) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
      dplyr::filter(dplyr::between(n, min_size, max_size))
    
    gs_df %>% 
      dplyr::inner_join(gs_counts, by = "gs_name") %>%
      dplyr::select(gs_name, gene_symbol)  # Explicit dplyr::select
  }
  
  # Define collections
  collections <- list(
    HALLMARK = list(collection = "H"),
    KEGG = list(collection = "C2", subcollection = kegg_version),
    REACTOME = list(collection = "C2", subcollection = "CP:REACTOME"),
    GO_BP = list(collection = "C5", subcollection = "GO:BP"),
    GO_MF = list(collection = "C5", subcollection = "GO:MF"),
    GO_CC = list(collection = "C5", subcollection = "GO:CC")
  )
  
  # Run ORA for each collection
  combined_results <- list()
  
  for (name in names(collections)) {
    message("  Analyzing: ", name)
    gs_df <- tryCatch({
      do.call(get_genesets, collections[[name]])
    }, error = function(e) {
      message("    Error retrieving gene sets: ", e$message)
      return(NULL)
    })
    
    if (is.null(gs_df)) next
    
    # Skip if no DE genes or small gene sets
    if (length(de_genes) == 0) {
      message("    No DE genes - skipping")
      next
    }
    if (nrow(gs_df) == 0) {
      message("    No gene sets after filtering - skipping")
      next
    }
    
    # Run ORA
    ora_result <- tryCatch({
      clusterProfiler::enricher(
        gene = de_genes,
        universe = background_genes,
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        minGSSize = min_size,
        maxGSSize = max_size,
        TERM2GENE = data.frame(
          term = gs_df$gs_name,
          gene = gs_df$gene_symbol
        )
      )
    }, error = function(e) {
      message("    ORA failed: ", e$message)
      return(NULL)
    })
    
    # Process results
    if (!is.null(ora_result) && nrow(ora_result) > 0) {
      res_df <- ora_result@result %>%
        dplyr::mutate(
          Collection = name,
          Comparison = prefix,
          DE_Genes = length(de_genes),
          Background_Genes = length(background_genes),
          .before = 1
        ) %>%
        dplyr::filter(p.adjust < 0.05)
      
      if (nrow(res_df) > 0) {
        combined_results[[name]] <- res_df
        message("    Found ", nrow(res_df), " significant pathways")
      } else {
        message("    No significant pathways")
      }
    }
  }
  
  # Combine results
  final_results <- dplyr::bind_rows(combined_results)
  
  # Save results if any
  if (nrow(final_results) > 0) {
    output_file <- file.path(output_dir, paste0(prefix, "_ORA_combined.csv"))
    write.csv(final_results, output_file, row.names = FALSE)
    message("Saved results to: ", output_file)
  } else {
    message("No significant results found for any collection")
  }
  
  return(final_results)
}

################################################################################

# 2. Batch Processing:
#    - Lists all CSV files in the `DGE/` directory.
#    - Separates files into `stagewise` and `temporal` analyses based on filename.
#    - Runs `run_ora_analysis()` on each file and saves results to:
#         Outputs/ORA/Stagewise/
#         Outputs/ORA/Temporal/

# List all CSV files in DGE/ folder
all_files <- list.files("DGE", pattern = "\\.csv$", full.names = TRUE)

# Separate files for stagewise and temporal DGE analysis
stagewise <- grep("StageWiseDGE_stage", all_files, value = TRUE)
temporal <- grep("TemporalDGE_condition", all_files, value = TRUE)

# Define output directory
output_dir_gsea <- "Outputs/ORA/"

# Create directory if it doesn't exist
if (!dir.exists(output_dir_gsea)) {
  dir.create(output_dir_gsea, recursive = TRUE)
}

# Stagewise analysis
stagewise_results <- lapply(stagewise, function(file) {
  run_ora_analysis(
    csv_path = file,
    output_dir = "Outputs/ORA/Stagewise"
  )
})

# Temporal analysis
temporal_results <- lapply(temporal, function(file) {
  run_ora_analysis(
    csv_path = file,
    output_dir = "Outputs/ORA/Temporal"
  )
})
