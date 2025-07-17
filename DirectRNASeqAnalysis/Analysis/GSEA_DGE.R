################################################################################

# Title: GSEA Analysis Pipeline using fgsea and msigdbr
# Author: Payel Sardar
#
# Description:
# This script performs an enhanced Gene Set Enrichment Analysis (GSEA) pipeline 
# using the `fgsea` package and gene sets from `msigdbr`. It reads differential 
# expression results (DESeq2-like CSV files), ranks genes, maps ENSEMBL IDs to 
# gene symbols, loads multiple gene set collections (Hallmark, KEGG, Reactome, GO),
# runs GSEA, saves results, and visualizes enrichment for the Oxidative Phosphorylation 
# (OxPhos) pathway if significant.
#
# Supported species: Mus musculus
# Input: Differential expression CSV files
# Output: Enrichment results (CSV, RDS), diagnostic messages, optional plots
#
# Output Directory Structure:
# ---------------------------
# Outputs/GSEA/
# ├── Stagewise/
# │   ├── *_FGSEA_*.csv       <- Per-collection GSEA results
# │   ├── *_FGSEA_results.rds <- Full result object
# │   └── Plots/              <- OxPhos enrichment plots
# └── Temporal/
#     ├── ...
#
# Required Packages:
# ------------------
# dplyr, msigdbr, fgsea, ggplot2, AnnotationDbi, org.Mm.eg.db, tibble
#
################################################################################

# Load all necessary packages
suppressPackageStartupMessages({
  library(dplyr)
  library(msigdbr)
  library(fgsea)
  library(ggplot2)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(tibble)
})

################################################################################

# 1. run_gsea_analysis_fgsea_enhanced():
#    - Core pipeline to perform GSEA using ranked gene lists.
#    - Supports signed ranking metric (-log10(padj) * sign(log2FC)).
#    - Loads and filters gene sets from multiple MSigDB collections.
#    - Performs ENSEMBL to SYMBOL mapping.
#    - Saves FGSEA results (CSV and RDS) to output directory.

run_gsea_analysis_fgsea_enhanced <- function(
    csv_path,
    species = "Mus musculus",
    gene_id_col = "X",
    output_dir = NULL,
    kegg_version = "CP:KEGG_LEGACY",
    min_size = 5,             # Reduced min size
    max_size = 1000,          # Increased max size
    nperm = 20000,            # Increased permutations
    rank_metric = "signed_p"  # Better ranking metric
) {
  
  # Prepare output
  if (!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }
  prefix <- tools::file_path_sans_ext(basename(csv_path))
  message("Processing: ", prefix)
  
  # Read and prepare data
  df <- read.csv(csv_path, header = TRUE) %>%
    filter(!is.na(log2FoldChange)) %>%  # Fixed syntax error here
    distinct(across(all_of(gene_id_col)), .keep_all = TRUE)
  
  # Create ranking metric
  if (rank_metric == "signed_p") {
    df <- df %>%
      filter(!is.na(padj)) %>%
      mutate(rank_val = sign(log2FoldChange) * -log10(padj)) %>%
      # Cap extreme values
      mutate(rank_val = ifelse(rank_val > 20, 20, 
                               ifelse(rank_val < -20, -20, rank_val)))
  } else {
    df$rank_val <- df$log2FoldChange
  }
  
  # ENSEMBL to SYMBOL conversion with diagnostics
  ens2symbol <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys = df[[gene_id_col]],
    keytype = "ENSEMBL",
    columns = "SYMBOL"
  ) %>%
    filter(!is.na(SYMBOL)) %>%
    distinct(ENSEMBL, .keep_all = TRUE)
  
  mapped_genes <- nrow(ens2symbol)
  total_genes <- nrow(df)
  message(
    "Gene mapping: ", mapped_genes, "/", total_genes, " (", 
    round(mapped_genes/total_genes*100, 1), "%) genes mapped"
  )
  
  gene_df <- df %>%
    inner_join(ens2symbol, by = c(setNames("ENSEMBL", gene_id_col))) %>%
    group_by(SYMBOL) %>%
    slice_max(order_by = abs(rank_val), n = 1, with_ties = FALSE) %>%
    ungroup()
  
  gene_list_named <- gene_df$rank_val
  names(gene_list_named) <- gene_df$SYMBOL
  gene_list_named <- sort(gene_list_named, decreasing = TRUE)
  
  # Load gene sets with ENSEMBL fallback option
  get_genesets <- function(collection, subcollection = NULL) {
    # Try SYMBOL first
    gs <- tryCatch({
      gs_df <- msigdbr(
        species = species, 
        category = collection,
        subcategory = subcollection
      )
      split(gs_df$gene_symbol, gs_df$gs_name)
    }, error = function(e) {
      # Fallback to ENSEMBL if SYMBOL fails
      message("Using ENSEMBL IDs for ", collection)
      gs_df <- msigdbr(
        species = species, 
        category = collection,
        subcategory = subcollection,
        gene_id = "ensembl_gene"
      )
      split(gs_df$ensembl_gene, gs_df$gs_name)
    })
    # Filter gene sets by size
    gs_lengths <- lengths(gs)
    gs[gs_lengths >= min_size & gs_lengths <= max_size]
  }
  
  message("Loading gene sets...")
  msigdb_sets <- list(
    HALLMARK = get_genesets("H"),
    KEGG = get_genesets("C2", kegg_version),
    REACTOME = get_genesets("C2", "CP:REACTOME"),
    GO_BP = get_genesets("C5", "GO:BP"),
    GO_MF = get_genesets("C5", "GO:MF"),
    GO_CC = get_genesets("C5", "GO:CC")
  )
  
  # Print gene set counts
  message(sprintf(
    "Filtered gene sets: Hallmark (%d), KEGG (%d), Reactome (%d), GO:BP (%d), GO:MF (%d), GO:CC (%d)",
    length(msigdb_sets$HALLMARK),
    length(msigdb_sets$KEGG),
    length(msigdb_sets$REACTOME),
    length(msigdb_sets$GO_BP),
    length(msigdb_sets$GO_MF),
    length(msigdb_sets$GO_CC)
  ))
  
  # Run fgsea
  fgsea_results <- lapply(names(msigdb_sets), function(name) {
    message("Running FGSEA for: ", name, " (", length(msigdb_sets[[name]]), " gene sets)")
    res <- fgsea(
      pathways = msigdb_sets[[name]],
      stats = gene_list_named,
      nperm = nperm,
      minSize = min_size,
      maxSize = max_size
    )
    res$collection <- name
    res <- res[order(res$padj), ]
    message("  - Significant pathways (padj < 0.05): ", sum(res$padj < 0.05, na.rm = TRUE))
    res
  }) %>% setNames(names(msigdb_sets))
  
  # Save results
  if (!is.null(output_dir)) {
    for (name in names(fgsea_results)) {
      result_for_csv <- fgsea_results[[name]]
      if ("leadingEdge" %in% colnames(result_for_csv)) {
        result_for_csv$leadingEdge <- sapply(result_for_csv$leadingEdge, paste, collapse = ";")
      }
      write.csv(
        result_for_csv,
        file.path(output_dir, paste0(prefix, "_FGSEA_", name, ".csv")),
        row.names = FALSE
      )
    }
    full_results <- list(
      fgsea_results = fgsea_results,
      ranked_gene_list = gene_list_named
    )
    saveRDS(full_results, file.path(output_dir, paste0(prefix, "_FGSEA_results.rds")))
    #saveRDS(fgsea_results, file.path(output_dir, paste0(prefix, "_FGSEA_results.rds")))
    message("Results saved to: ", output_dir)
  }
  
  return(list(
    fgsea_results = fgsea_results,
    ranked_gene_list = gene_list_named
  ))
}

################################################################################

# 2. examine_gsea_results():
#    - Diagnostic helper to print summary stats for each gene set collection.
#    - Lists top up/downregulated pathways.

examine_gsea_results <- function(fgsea_results, top_n = 10) {
  cat("=== GSEA Results Summary ===\n")
  
  for (collection in names(fgsea_results)) {
    cat("\n", collection, ":\n")
    res <- fgsea_results[[collection]]
    
    cat("  Total pathways tested:", nrow(res), "\n")
    cat("  Significant (padj < 0.05):", sum(res$padj < 0.05, na.rm = TRUE), "\n")
    cat("  Significant (padj < 0.25):", sum(res$padj < 0.25, na.rm = TRUE), "\n")
    
    # Show top upregulated pathways
    up_paths <- res[res$NES > 0, ][1:min(top_n, sum(res$NES > 0)), ]
    if (nrow(up_paths) > 0) {
      cat("  Top upregulated pathways:\n")
      for (i in 1:nrow(up_paths)) {
        cat(sprintf("    %s (NES=%.2f, padj=%.3f)\n", 
                    up_paths$pathway[i], up_paths$NES[i], up_paths$padj[i]))
      }
    }
    
    # Show top downregulated pathways  
    down_paths <- res[res$NES < 0, ][1:min(top_n, sum(res$NES < 0)), ]
    if (nrow(down_paths) > 0) {
      cat("  Top downregulated pathways:\n")
      for (i in 1:nrow(down_paths)) {
        cat(sprintf("    %s (NES=%.2f, padj=%.3f)\n",
                    down_paths$pathway[i], down_paths$NES[i], down_paths$padj[i]))
      }
    }
  }
}

################################################################################

# 3. check_gene_list_stats():
#    - Helper function to print basic stats about the input gene list.
#    - Useful for quality control and exploring DE patterns.

check_gene_list_stats <- function(csv_path, gene_id_col = "X") {
  df <- read.csv(csv_path, header = TRUE)
  
  cat("=== Gene List Statistics ===\n")
  cat("Total genes in file:", nrow(df), "\n")
  cat("Genes with log2FoldChange:", sum(!is.na(df$log2FoldChange)), "\n")
  cat("Genes with padj < 0.05:", sum(df$padj < 0.05, na.rm = TRUE), "\n")
  cat("Genes with padj < 0.1:", sum(df$padj < 0.1, na.rm = TRUE), "\n")
  
  # Log2FC distribution
  lfc <- df$log2FoldChange[!is.na(df$log2FoldChange)]
  cat("Log2FC range:", round(min(lfc), 2), "to", round(max(lfc), 2), "\n")
  cat("Log2FC > 1:", sum(lfc > 1), "\n")
  cat("Log2FC < -1:", sum(lfc < -1), "\n")
  
  # Show most up/down regulated genes
  df_sorted <- df[!is.na(df$log2FoldChange), ]
  df_sorted <- df_sorted[order(-df_sorted$log2FoldChange), ]
  
  cat("\nTop 5 upregulated genes:\n")
  print(df_sorted[1:5, c(gene_id_col, "log2FoldChange", "padj")])
  
  cat("\nTop 5 downregulated genes:\n")
  print(tail(df_sorted[, c(gene_id_col, "log2FoldChange", "padj")], 5))
}

################################################################################

# 4. check_msigdb_collections():
#    - Helper to list available MSigDB collections and subcollections.

check_msigdb_collections <- function(species = "Mus musculus") {
  cat("Available collections for", species, ":\n")
  collections <- msigdbr_collections()
  print(collections)
  
  # Show specific subcollections for C2 and C5
  cat("\nC2 subcollections:\n")
  c2_data <- msigdbr(species = species, collection = "C2")
  print(unique(c2_data$gs_subcollection))
  
  cat("\nC5 subcollections:\n")
  c5_data <- msigdbr(species = species, collection = "C5")
  print(unique(c5_data$gs_subcollection))
}

################################################################################

# List all CSV files in DGE/ folder
all_files <- list.files("DGE", pattern = "\\.csv$", full.names = TRUE)

# Separate files for stagewise and temporal DGE analysis
stagewise <- grep("StageWiseDGE_stage", all_files, value = TRUE)
temporal <- grep("TemporalDGE_condition", all_files, value = TRUE)

# Define output directory
output_dir_gsea <- "Outputs/GSEA/"

# Create directory if it doesn't exist
if (!dir.exists(output_dir_gsea)) {
  dir.create(output_dir_gsea, recursive = TRUE)
}

# Preload KEGG gene sets for OxPhos plotting
kegg_genesets <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG_LEGACY") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  base::split(x = .$gene_symbol, f = .$gs_name)

################################################################################

# 5. plot_oxphos_if_enriched():
#    - Automatically detects if OxPhos is significantly enriched in any collection.
#    - If so, generates and saves an enrichment plot.

# Function to plot OxPhos if enriched in any relevant collection
plot_oxphos_if_enriched <- function(
    fgsea_result, 
    ranked_genes, 
    output_dir, 
    file_prefix,
    kegg_version = "CP:KEGG_MEDICUS"  # Add parameter here
) {
  # Define all possible OxPhos pathway identifiers
  oxphos_pathways <- c(
    KEGG = "KEGG_OXIDATIVE_PHOSPHORYLATION",
    HALLMARK = "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    GO_BP = "GOBP_OXIDATIVE_PHOSPHORYLATION"
  )
  
  # Preload gene sets with consistent KEGG version
  gene_sets <- list(
    KEGG = msigdbr(species = "Mus musculus", category = "C2", subcategory = kegg_version) %>% 
      filter(gs_name == oxphos_pathways["KEGG"]) %>%
      pull(gene_symbol),
    HALLMARK = msigdbr(species = "Mus musculus", category = "H") %>%
      filter(gs_name == oxphos_pathways["HALLMARK"]) %>%
      pull(gene_symbol),
    GO_BP = msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP") %>%
      filter(gs_name == oxphos_pathways["GO_BP"]) %>%
      pull(gene_symbol)
  )
  
  # Check for significant OxPhos in any collection
  oxphos_res <- list()
  for (collection in names(oxphos_pathways)) {
    if (collection %in% names(fgsea_result)) {
      res <- fgsea_result[[collection]] %>%
        filter(pathway == oxphos_pathways[collection] & padj < 0.05)
      if (nrow(res) > 0) {
        res$collection <- collection
        oxphos_res[[collection]] <- res
      }
    }
  }
  
  oxphos_res <- bind_rows(oxphos_res)
  
  if (nrow(oxphos_res) > 0) {
    # Use the most significant result
    best_res <- oxphos_res %>%
      slice_min(padj, n = 1) %>%
      as.list()
    
    # Create plot using the appropriate gene set
    p <- fgsea::plotEnrichment(
      gene_sets[[best_res$collection]],
      ranked_genes
    ) + 
      labs(
        title = paste("Oxidative Phosphorylation\n", file_prefix),
        subtitle = paste(
          best_res$collection, "pathway\n",
          "NES =", round(best_res$NES, 2), 
          "padj =", format.pval(best_res$padj, digits = 2)
        )
      )
    
    # Save plot
    plot_dir <- file.path(output_dir, "OxPhos_Plots")
    dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
    plot_path <- file.path(
      plot_dir, 
      paste0(file_prefix, "_", best_res$collection, "_OxPhos_enrichment.png")
    )
    ggsave(plot_path, p, width = 8, height = 6, dpi = 300)
    message("Saved OxPhos plot to: ", plot_path)
    return(plot_path)
  }
  return(NULL)
}

################################################################################

# 6. Batch processing
#    - Iterates through CSV files in `DGE/` folder.
#    - Separates `Stagewise` and `Temporal` analyses based on filename.
#    - Runs full GSEA + OxPhos plotting pipeline for each.

# Run analysis for stagewise files
for (file in stagewise) {
  file_prefix <- tools::file_path_sans_ext(basename(file))
  message("\nProcessing stagewise file: ", file_prefix)
  
  # Run enhanced GSEA
  gsea_result <- run_gsea_analysis_fgsea_enhanced(
    csv_path = file,
    output_dir = file.path(output_dir_gsea, "Stagewise"),
    kegg_version = "CP:KEGG_MEDICUS",
    rank_metric = "signed_p",
    gene_id_col = "X"
  )
  
  # Plot OxPhos if enriched
  plot_oxphos_if_enriched(
    fgsea_result = gsea_result$fgsea_results,
    ranked_genes = gsea_result$ranked_gene_list,
    output_dir = file.path(output_dir_gsea, "Stagewise/Plots"),
    file_prefix = file_prefix,
    kegg_version = "CP:KEGG_MEDICUS"  # Add this parameter
  )
}

# Run analysis for temporal files
for (file in temporal) {
  file_prefix <- tools::file_path_sans_ext(basename(file))
  message("\nProcessing temporal file: ", file_prefix)
  
  # Run enhanced GSEA
  gsea_result <- run_gsea_analysis_fgsea_enhanced(
    csv_path = file,
    output_dir = file.path(output_dir_gsea, "Temporal"),
    kegg_version = "CP:KEGG_MEDICUS",
    rank_metric = "signed_p",
    gene_id_col = "X" 
  )
  
  # Plot OxPhos if enriched
  plot_oxphos_if_enriched(
    fgsea_result = gsea_result$fgsea_results,
    ranked_genes = gsea_result$ranked_gene_list,
    output_dir = file.path(output_dir_gsea, "Temporal/Plots"),
    file_prefix = file_prefix,
    kegg_version = "CP:KEGG_MEDICUS"
  )
}