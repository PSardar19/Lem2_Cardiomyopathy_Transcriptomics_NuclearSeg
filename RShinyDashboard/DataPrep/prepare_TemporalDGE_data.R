library(biomaRt)
library(dplyr)
library(stringr)

# List files from FLAMES and NanoSeq folders
flames_paths <- list.files("../FLAMES/DGE", full.names = TRUE, pattern = "TemporalDGE_condition")
nanoseq_paths <- list.files("../NanoSeq/DGE", full.names = TRUE, pattern = "TemporalDGE_condition")

# Extract metadata from filename (condition, target stage, reference stage)
extract_metadata <- function(filename) {
  condition <- str_extract(filename, "(WT|Mut)")
  stages_match <- str_match(filename, "_(P?[0-9\\.]+|E[0-9\\.]+)_vs_(P?[0-9\\.]+|E[0-9\\.]+)")
  target_stage <- stages_match[, 2]
  reference_stage <- stages_match[, 3]
  
  list(condition = condition,
       transition = paste(target_stage, "vs", reference_stage, sep = "_"))
}

# Create nested list structure: condition -> transition -> data.frame
load_and_structure_dge <- function(paths) {
  # Read all files into a named list by combined "condition_transition"
  all_data <- lapply(paths, read.csv)
  names(all_data) <- sapply(paths, function(p) {
    meta <- extract_metadata(basename(p))
    paste(meta$condition, meta$transition, sep = "_")
  })
  
  # Split into nested list: condition -> transition -> data.frame
  nested_list <- list()
  for (nm in names(all_data)) {
    parts <- str_split(nm, "_", n = 2)[[1]]  # splits into c(condition, transition)
    condition <- parts[1]
    transition <- parts[2]
    
    if (!condition %in% names(nested_list)) nested_list[[condition]] <- list()
    nested_list[[condition]][[transition]] <- all_data[[nm]]
  }
  
  nested_list
}

temporal_dge_data <- list(
  FLAMES = load_and_structure_dge(flames_paths),
  NanoSeq = load_and_structure_dge(nanoseq_paths)
)

# Combine all gene IDs across all datasets & conditions & transitions
flames_ids <- unique(unlist(lapply(temporal_dge_data$FLAMES, function(cond) {
  unlist(lapply(cond, function(df) df$X))
})))

nano_ids <- unique(unlist(lapply(temporal_dge_data$NanoSeq, function(cond) {
  unlist(lapply(cond, function(df) df$X))
})))

all_gene_ids <- unique(c(flames_ids, nano_ids))

# Annotate gene_ids using biomaRt
ensembl_mouse <- useEnsembl(biomart = 'genes', 
                            dataset = 'mmusculus_gene_ensembl',
                            mirror = 'www')

annotations <- getBM(
  attributes = c('ensembl_gene_id', 'external_gene_name', 'description'),
  filters = 'ensembl_gene_id',
  values = all_gene_ids,
  mart = ensembl_mouse
)
colnames(annotations)[1:2] <- c("gene_id", "gene_name")

annotate_dge <- function(df) {
  if ("X" %in% colnames(df)) colnames(df)[colnames(df) == "X"] <- "gene_id"
  merge(df, annotations, by = "gene_id", all.x = TRUE)
}

# Annotate and clean stat column for nested lists
annotate_and_clean <- function(nested_list) {
  lapply(nested_list, function(cond_list) {
    lapply(cond_list, function(df) {
      df2 <- annotate_dge(df)
      if ("stat" %in% colnames(df2)) df2 <- df2[, !(colnames(df2) %in% "stat")]
      df2
    })
  })
}

temporal_dge_data$FLAMES <- annotate_and_clean(temporal_dge_data$FLAMES)
temporal_dge_data$NanoSeq <- annotate_and_clean(temporal_dge_data$NanoSeq)

# Save annotated data (optional, adjust file paths if needed)
save_annotated_files <- function(nested_list, dir_path) {
  for (condition in names(nested_list)) {
    for (transition in names(nested_list[[condition]])) {
      file_name <- paste0("TemporalDGE_", condition, "_", transition, "_annotated.csv")
      write.csv(nested_list[[condition]][[transition]],
                file = file.path(dir_path, file_name), row.names = FALSE)
    }
  }
}

save_annotated_files(temporal_dge_data$FLAMES, "../FLAMES/DGE")
save_annotated_files(temporal_dge_data$NanoSeq, "../NanoSeq/DGE")

# Save entire RData object
save(temporal_dge_data, file = "data/temporal_dge_data.RData")

# Add this section to your existing script

# Load expression data
flames_norm_expr <- read.csv("../FLAMES/DGE/Grouped_norm_gene_counts_FLAMES.csv") %>%
  mutate(pipeline = "FLAMES")
nanoseq_norm_expr <- read.csv("../NanoSeq/DGE/Grouped_norm_gene_counts_NanoSeq.csv") %>%
  mutate(pipeline = "NanoSeq")

# Combine and prepare expression data
expression_data <- bind_rows(flames_norm_expr, nanoseq_norm_expr) %>%
  mutate(
    stage = factor(stage, levels = c("E18.5", "P1", "P3")),
    condition = factor(condition, levels = c("WT", "Mut"))
  )

# Save both objects together
save(temporal_dge_data, expression_data, file = "data/temporal_data.RData")