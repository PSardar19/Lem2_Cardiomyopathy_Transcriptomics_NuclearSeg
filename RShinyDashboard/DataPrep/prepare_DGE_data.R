## Importing all necessary packages
library(biomaRt)
library(dplyr)

## List files from FLAMES and NanoSeq
flames_paths <- list.files("../FLAMES/DGE", full.names = TRUE, pattern = "StageWiseDGE_stage_")
nanoseq_paths <- list.files("../NanoSeq/DGE", full.names = TRUE, pattern = "StageWiseDGE_stage_")

# Extract stage from filenames using regex
extract_stage <- function(filename) {
  if (grepl("E18.5", filename)) return("E18.5")
  if (grepl("P1", filename)) return("P1")
  if (grepl("P3", filename)) return("P3")
  return(NA)
}

# Create named lists by stage
flames_files <- setNames(flames_paths, sapply(flames_paths, extract_stage))
nanoseq_files <- setNames(nanoseq_paths, sapply(nanoseq_paths, extract_stage))

# Load the data
dge_data <- list(
  FLAMES = lapply(flames_files, read.csv),
  NanoSeq = lapply(nanoseq_files, read.csv)
)

# Get unique gene IDs from both datasets
flames_ids <- unique(unlist(lapply(dge_data$FLAMES, function(df) df$X)))
nano_ids   <- unique(unlist(lapply(dge_data$NanoSeq, function(df) df$X)))

# Combine and deduplicate
all_gene_ids <- unique(c(flames_ids, nano_ids))

# Annotating the gene_ids
ensembl <- useEnsembl(biomart = "genes")
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

# Saving the annotated data

# Rename "X" to "gene_id" and merge annotations in both NanoSeq and FLAMES
dge_data$NanoSeq <- lapply(dge_data$NanoSeq, function(df) {
  if ("X" %in% colnames(df)) colnames(df)[colnames(df) == "X"] <- "gene_id"
  merged_df <- merge(df, annotations, by = "gene_id", all.x = TRUE)
  return(merged_df)
})

dge_data$FLAMES <- lapply(dge_data$FLAMES, function(df) {
  if ("X" %in% colnames(df)) colnames(df)[colnames(df) == "X"] <- "gene_id"
  merged_df <- merge(df, annotations, by = "gene_id", all.x = TRUE
                     )
  return(merged_df)
})

# Clean 'stat' column if present
clean_stat_column <- function(df) {
  if ("stat" %in% colnames(df)) {
    df <- df[, !(colnames(df) %in% "stat")]
  }
  return(df)
}

dge_data$NanoSeq <- lapply(dge_data$NanoSeq, clean_stat_column)
dge_data$FLAMES <- lapply(dge_data$FLAMES, clean_stat_column)

# Assign names again after the lapply operations
stages <- c("E18.5", "P1", "P3")
names(dge_data$NanoSeq) <- stages
names(dge_data$FLAMES) <- stages

# Save annotated data files
lapply(names(dge_data$NanoSeq), function(stage) {
  file_path <- paste0("../NanoSeq/DGE/StageWiseDGE_", stage, "_Mut_vs_WT_annotated.csv")
  write.csv(dge_data$NanoSeq[[stage]], file = file_path, row.names = FALSE)
})

lapply(names(dge_data$FLAMES), function(stage) {
  file_path <- paste0("../FLAMES/DGE/StageWiseDGE_", stage, "_Mut_vs_WT_annotated.csv")
  write.csv(dge_data$FLAMES[[stage]], file = file_path, row.names = FALSE)
})

# Save the full dge_data object
save(dge_data, file = "data/dge_data2.RData")
