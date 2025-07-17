library(dplyr)
library(stringr)
library(here)
library(purrr)
library(msigdbr)
library(tibble)

# Define directories
gsea_dirs <- c(
  FLAMES = here("../FLAMES/Outputs/GSEA"),
  NanoSeq = here("../NanoSeq/Outputs/GSEA")
)

# Enhanced metadata parser
parse_gsea_filename <- function(filename) {
  if (str_detect(filename, "StageWiseDGE")) {
    stage <- str_extract(filename, "E18\\.5|P1|P3")
    return(list(
      analysis_type = "Stagewise",
      stage = stage,
      condition = "Mut_vs_WT",
      transition = NA_character_
    ))
  } else if (str_detect(filename, "TemporalDGE")) {
    condition <- str_extract(filename, "(?<=condition_)(WT|Mut)")
    stages <- str_match(filename, "(E18\\.5|P1|P3)_vs_(E18\\.5|P1|P3)")
    transition <- paste(stages[1,3], stages[1,2], sep = " -> ")
    return(list(
      analysis_type = "Temporal",
      stage = NA_character_,
      condition = condition,
      transition = transition
    ))
  } else {
    return(list(
      analysis_type = "Unknown",
      stage = NA_character_,
      condition = NA_character_,
      transition = NA_character_
    ))
  }
}

# Get pathway descriptions
pathway_meta <- msigdbr(species = "Mus musculus") %>%
  distinct(gs_name, gs_description, gs_collection, gs_subcollection) %>%
  mutate(
    collection = case_when(
      gs_collection == "H" ~ "HALLMARK",
      gs_collection == "C2" & grepl("KEGG", gs_subcollection) ~ "KEGG",
      gs_collection == "C2" & grepl("REACTOME", gs_subcollection) ~ "REACTOME",
      gs_collection == "C5" & grepl("GO:BP", gs_subcollection) ~ "GO_BP",
      gs_collection == "C5" & grepl("GO:MF", gs_subcollection) ~ "GO_MF",
      gs_collection == "C5" & grepl("GO:CC", gs_subcollection) ~ "GO_CC",
      TRUE ~ "Other"
    )
  ) %>%
  select(-gs_collection, -gs_subcollection) %>%
  filter(collection %in% c("HALLMARK", "KEGG", "REACTOME", "GO_BP", "GO_MF", "GO_CC"))

# Helper function to process a single GSEA file
process_gsea_file <- function(rds_file, pipeline_name) {
  tryCatch({
    filename <- basename(rds_file)
    prefix <- str_remove(filename, "_FGSEA_results\\.rds")
    meta <- parse_gsea_filename(prefix)
    
    message("Processing: ", rds_file)
    message("Parsed metadata: ", 
            "analysis=", meta$analysis_type,
            " stage=", meta$stage, 
            " condition=", meta$condition)
    
    # Read full results
    full_results <- readRDS(rds_file)
    
    # Create unique key for this analysis
    file_key <- paste(pipeline_name, prefix, sep = "::")
    
    # Process each collection
    df <- map_dfr(names(full_results$fgsea_results), function(collection) {
      result_df <- full_results$fgsea_results[[collection]]
      
      if (nrow(result_df) == 0) return(tibble())
      
      result_df %>%
        mutate(
          pipeline = pipeline_name,
          collection = collection,
          analysis_type = meta$analysis_type,
          stage = meta$stage,
          condition = meta$condition,
          transition = meta$transition,
          file_prefix = prefix,
          file_key = file_key
        ) %>%
        left_join(
          pathway_meta, 
          by = c("pathway" = "gs_name", "collection" = "collection")
        ) %>%
        mutate(
          gs_description = ifelse(is.na(gs_description), pathway, gs_description),
          Significance = case_when(
            padj < 0.01 ~ "Highly Significant",
            padj < 0.05 ~ "Significant",
            TRUE ~ "Not Significant"
          ),
          Direction = ifelse(NES > 0, "Activated", "Suppressed")
        ) %>%
        # Store gene counts and preview
        mutate(
          n_genes = map_int(leadingEdge, length),
          gene_preview = map_chr(leadingEdge, ~paste(head(., 5), collapse = ";"))
        ) %>%
        # Keep only essential columns
        select(pathway, gs_description, collection, analysis_type, pipeline,
               stage, condition, transition, file_prefix, file_key,
               NES, pval, padj, size, Significance, Direction, n_genes, gene_preview)
    })
    
    return(list(
      df = df,
      ranked_genes = full_results$ranked_gene_list,
      key = file_key,
      # Store leading edge genes only for significant pathways
      leading_edges = full_results$fgsea_results %>% 
        map(~.x %>% filter(padj < 0.05) %>% select(pathway, collection, leadingEdge))
    ))
  }, error = function(e) {
    warning("Error processing ", rds_file, ": ", e$message)
    return(list(df = tibble(), ranked_genes = NULL, key = NULL, leading_edges = NULL))
  })
}

# Main processing loop
gsea_data <- tibble()
ranked_genes_list <- list()
leading_edge_list <- list()

for (pipeline in names(gsea_dirs)) {
  dir_path <- gsea_dirs[pipeline]
  rds_files <- list.files(
    path = dir_path,
    pattern = "_FGSEA_results\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  message("Found ", length(rds_files), " files in ", pipeline)
  
  for (rds_file in rds_files) {
    result <- process_gsea_file(rds_file, pipeline)
    
    # Only add if we have valid results
    if (nrow(result$df) > 0 && !is.null(result$key)) {
      gsea_data <- bind_rows(gsea_data, result$df)
      ranked_genes_list[[result$key]] <- result$ranked_genes
      
      # Store leading edges with unique keys
      for (collection in names(result$leading_edges)) {
        if (nrow(result$leading_edges[[collection]]) > 0) {
          result$leading_edges[[collection]] %>%
            rowwise() %>%
            mutate(edge_key = paste(
              result$key, collection, pathway, sep = "::"
            )) %>%
            ungroup() -> edges_df
          
          leading_edge_list <- c(
            leading_edge_list,
            setNames(edges_df$leadingEdge, edges_df$edge_key)
          )
        }
      }
    }
  }
}

# Filter to get significant results
gsea_data_sig <- gsea_data %>% 
  filter(padj < 0.05)

# Save for Shiny
fgsea_data <- list(
  gsea_results = gsea_data_sig,
  ranked_genes_list = ranked_genes_list,
  leading_edge_genes = leading_edge_list
)

saveRDS(fgsea_data, "data/gsea_full_data.rds")

# Print summary
message("\nGSEA data processing complete!")
message("Total significant pathways: ", nrow(gsea_data_sig))
message("Unique analyses: ", length(ranked_genes_list))
message("Pathways with stored genes: ", length(leading_edge_list))
message("Memory usage: ", format(object.size(fgsea_data), units = "auto"))