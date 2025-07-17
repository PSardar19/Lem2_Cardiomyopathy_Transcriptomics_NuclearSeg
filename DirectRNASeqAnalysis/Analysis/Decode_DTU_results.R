# DTU Results Extraction and Analysis

library(biomaRt)

# Loading the results for stage P1 (For results + discussion)
# Change path for other stages (Future scope: Convert into a function that can be reused )
switchList_P1 <- readRDS("DTU/switchList_P1.rds")

# =====================================
# 1. BASIC OVERVIEW AND SUMMARY STATS
# =====================================

# Checking the structure of key components
print("=== BASIC INFO ===")
print(paste("Conditions analyzed:", paste(switchList_P1$conditions, collapse = ", ")))
print(paste("Number of samples:", nrow(switchList_P1$designMatrix)))
print(paste("Number of isoforms:", nrow(switchList_P1$isoformFeatures)))

# Checking the design matrix
print("=== DESIGN MATRIX ===")
print(switchList_P1$designMatrix)
print(table(switchList_P1$designMatrix$condition))

# =====================================
# 2. ISOFORM SWITCHING RESULTS
# =====================================

print("=== ISOFORM SWITCHING ANALYSIS ===")

# Checking if switching analysis was successful
if("isoformSwitchAnalysis" %in% names(switchList_P1) && 
   !is.null(switchList_P1$isoformSwitchAnalysis)) {
  
  # Extracting switching results
  switching_results <- switchList_P1$isoformSwitchAnalysis
  print("Switching analysis results:")
  print(head(switching_results))
  
  # Summary of significant switches
  if("padj" %in% names(switching_results)) {
    sig_switches <- switching_results[switching_results$padj < 0.05, ]
    print(paste("Number of significant switching isoforms (p.adj < 0.05):", nrow(sig_switches)))
    
    if(nrow(sig_switches) > 0) {
      print("Top 10 significant switches:")
      print(sig_switches[order(sig_switches$padj)[1:min(10, nrow(sig_switches))], 
                         c( "isoform_id", "condition_1", "condition_2", 
                           "dIF", "padj")])
    }
  }
} else {
  print("No switching analysis results found - may need to run isoformSwitchTestDEXSeq()")
}

#Annotating isoform features:
# Creating an ensembl object
ensembl <- useEnsembl(biomart = "genes")

# For selecting the dataset
datasets <- listDatasets(ensembl)
head(datasets)

ensembl_mouse <- useEnsembl(biomart = 'genes', 
                            dataset = 'mmusculus_gene_ensembl',
                            mirror = 'www')

isoform_list <- unique(sig_switches$isoform_id)

annot_tx <- getBM(
  attributes = c('ensembl_transcript_id', 'ensembl_gene_id',
                 'external_gene_name', 'transcript_biotype', 'description'),
  filters = 'ensembl_transcript_id',
  values = isoform_list,
  mart = ensembl_mouse
)
head(annot_tx)

annot_sig_switch <- merge(sig_switches, annot_tx,
                   by.x = "isoform_id", by.y = "ensembl_transcript_id",
                   all.x = TRUE)

# =====================================
# 3. ISOFORM FEATURES AND EXPRESSION
# =====================================

print("=== ISOFORM FEATURES ===")
isoform_features <- switchList_P1$isoformFeatures
print("Column names in isoformFeatures:")
print(names(isoform_features))

# Basic stats about isoforms
print(paste("Total genes:", length(unique(isoform_features$gene_id))))
print(paste("Total isoforms:", nrow(isoform_features)))

# Expression summary
if("isoformRepExpression" %in% names(switchList_P1)) {
  expr_data <- switchList_P1$isoformRepExpression
  print("=== EXPRESSION SUMMARY ===")
  print("Expression data dimensions:")
  print(dim(expr_data))
  print("Expression summary (first few rows):")
  print(head(expr_data))
}

# Isoform fractions
if("isoformRepIF" %in% names(switchList_P1)) {
  if_data <- switchList_P1$isoformRepIF
  print("=== ISOFORM FRACTION SUMMARY ===")
  print("IF data dimensions:")
  print(dim(if_data))
  print("IF summary (first few rows):")
  print(head(if_data))
}

# =====================================
# 4. FUNCTIONAL CONSEQUENCES (WITH GENE ANNOTATION)
# =====================================

print("=== SWITCH CONSEQUENCES ===")
if("switchConsequence" %in% names(switchList_WT) && 
   !is.null(switchList_WT$switchConsequence)) {
  
  consequences <- switchList_WT$switchConsequence
  print("Available consequence columns:")
  print(names(consequences))
  
  # Load biomaRt for gene annotation
  if (!require("biomaRt")) install.packages("biomaRt")
  library(biomaRt)
  
  # Creating mouse gene annotation function
  annotate_genes <- function(gene_ids) {
    # Use Ensembl mouse database
    mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    
    # Get gene symbols
    annotation <- getBM(
      attributes = c("ensembl_gene_id", "external_gene_name"),
      filters = "ensembl_gene_id",
      values = unique(gene_ids),
      mart = mart
    )
    
    # Return as named vector
    setNames(annotation$external_gene_name, annotation$ensembl_gene_id)
  }
  
  # Annotate gene symbols
  gene_symbols <- annotate_genes(consequences$gene_id)
  
  # Add gene symbols to consequences table
  consequences$gene_symbol <- gene_symbols[consequences$gene_id]
  
  # REVISION 1: Handle detailed consequence types
  print("=== FUNCTIONAL CONSEQUENCE SUMMARY ===")
  
  # Create frequency table of consequence types
  consequence_counts <- as.data.frame(
    table(Consequence = consequences$switchConsequence)
  )
  
  # Filter for actual consequences (Freq > 0)
  significant_consequences <- consequence_counts[consequence_counts$Freq > 0, ]
  
  # Print formatted summary
  if(nrow(significant_consequences) > 0) {
    print("Functional consequences detected:")
    print(significant_consequences[order(-significant_consequences$Freq), ])
  } else {
    print("No functional consequences found")
  }
  
  # Count by consequence type
  if("switchConsequence" %in% names(consequences)) {
    ir_switches <- sum(
      grepl("intron retention", consequences$switchConsequence, ignore.case = TRUE),
      na.rm = TRUE
    )
    print(paste("Genes with intron retention consequences:", ir_switches))
    
    nmd_switches <- sum(
      grepl("NMD", consequences$switchConsequence, ignore.case = TRUE),
      na.rm = TRUE
    )
    print(paste("Genes with NMD-related consequences:", nmd_switches))
    
    orf_switches <- sum(
      grepl("ORF", consequences$switchConsequence, ignore.case = TRUE),
      na.rm = TRUE
    )
    print(paste("Genes with ORF changes:", orf_switches))
  }
  
  # Adding detailed consequence view with gene symbols
  if(nrow(consequences) > 0) {
    print("Top consequences by gene:")
    
    # Create display names (use gene symbol if available, otherwise gene ID)
    consequences$display_name <- ifelse(
      !is.na(consequences$gene_symbol) & consequences$gene_symbol != "",
      consequences$gene_symbol,
      consequences$gene_id
    )
    
    # Get top consequences
    top_conseq <- consequences[
      order(consequences$display_name), 
      c("display_name", "switchConsequence")
    ]
    
    # Remove NAs and empty consequences
    top_conseq <- top_conseq[!is.na(top_conseq$switchConsequence) & 
                               top_conseq$switchConsequence != "", ]
    
    # Print first 10 unique genes
    print(head(unique(top_conseq), 10))
  }
  
  # Saving annotated consequences
  switchList_WT$switchConsequence <- consequences
  print("Added gene symbols to switchConsequence table")
  
} else {
  print("No switch consequences found - may need to run analyzeSwitchConsequences()")
}

# =====================================
# 5. ORF ANALYSIS RESULTS
# =====================================

print("=== ORF ANALYSIS ===")
if("orfAnalysis" %in% names(switchList_WT) && 
   !is.null(switchList_WT$orfAnalysis)) {
  
  orf_data <- switchList_WT$orfAnalysis
  print("ORF analysis columns:")
  print(names(orf_data))
  
  orf_length_col <- "orfTransciptLength"
  if(orf_length_col %in% names(orf_data)) {
    print(paste("Using ORF length column:", orf_length_col))
    orf_lengths <- orf_data[[orf_length_col]]
    
    # Calculating statistics
    total_isoforms <- length(orf_lengths)
    orf_detected <- sum(!is.na(orf_lengths))
    orf_percentage <- round(orf_detected/total_isoforms * 100, 1)
    
    print(paste("ORF detected in", orf_detected, "of", total_isoforms, 
                "isoforms (", orf_percentage, "%)"))
    
    # Summary statistics
    print("ORF length summary (bp):")
    print(summary(na.omit(orf_lengths)))
    
    # Quantiles for better distribution understanding
    quantiles <- quantile(na.omit(orf_lengths), probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
    print("ORF length quantiles (5%, 25%, 50%, 75%, 95%):")
    print(quantiles)
    
    # Additional useful metrics
    mean_length <- mean(na.omit(orf_lengths))
    median_length <- median(na.omit(orf_lengths))
    print(paste("Mean ORF length:", round(mean_length, 1), "bp"))
    print(paste("Median ORF length:", median_length, "bp"))
    
  } else {
    print("No ORF length column found")
  }
  
} else {
  print("No ORF analysis found - may need to run analyzeORF()")
}

# =====================================
# 6. INTRON RETENTION ANALYSIS
# =====================================

print("=== INTRON RETENTION ANALYSIS ===")

# Check if switch analysis exists
if("isoformSwitchAnalysis" %in% names(switchList) && 
   !is.null(switchList$isoformSwitchAnalysis)) {
  
  switch_analysis <- switchList$isoformSwitchAnalysis
  print("Switch analysis columns:")
  print(names(switch_analysis))
  
  # Subset for intron retention events only
  ir_rows <- switch_analysis[switch_analysis$featureCompared == "intron_retention", ]
  
  if(nrow(ir_rows) > 0) {
    print(paste("Found", nrow(ir_rows), "differential intron retention events"))
    
    # Extract dIF values
    dIF_vals <- ir_rows$dIF
    
    # Filter significant IR events (|dIF| > 0.1)
    sig_ir <- sum(abs(dIF_vals) > 0.1, na.rm = TRUE)
    
    print("Intron retention dIF summary:")
    print(summary(dIF_vals))
    print(paste("Significant IR events (|dIF| > 0.1):", sig_ir))
    
    # Merge with genomic coordinates from intronRetentionAnalysis(Optional, only if genomic coordinates are available)
    if("intronRetentionAnalysis" %in% names(switchList)) {
      ir_data <- switchList$intronRetentionAnalysis
      merged_ir <- merge(ir_rows, ir_data, by = "isoform_id")
      print("Merged IR data with genomic coordinates:")
      print(head(merged_ir))
    }
    
  } else {
    print("No differential intron retention events found")
  }
  
} else {
  print("No switch analysis results found")
}

# =====================================
# 7. CREATE SUMMARY TABLE FOR TOP RESULTS
# =====================================

print("=== CREATING SUMMARY TABLE ===")

# Try to create a comprehensive summary table
if(exists("switching_results") && nrow(switching_results) > 0) {
  
  # Merge with isoform features for gene names
  summary_table <- merge(switching_results, 
                         isoform_features[, c("isoform_id", "gene_name", "gene_id")], 
                         by = "isoform_id", all.x = TRUE)
  summary_table <- summary_table %>% distinct()
  # Add consequences if available
  if("switchConsequence" %in% names(switchList_WT) && 
     !is.null(switchList_WT$switchConsequence)) {
    summary_table <- merge(summary_table, 
                           switchList_WT$switchConsequence, 
                           by = "isoform_id", all.x = TRUE)
  }
  
  # Select key columns and sort by significance
  key_cols <- intersect(names(summary_table), 
                        c("gene_name", "isoform_id", "condition_1", "condition_2",
                          "dIF", "isoform_switch_q_value", "intron_retention", 
                          "NMD_status", "ORF_seq_similarity"))
  
  final_summary <- summary_table[, key_cols]
  final_summary <- final_summary[order(final_summary$isoform_switch_q_value), ]
  
  print("=== TOP 20 SWITCHING EVENTS ===")
  print(head(final_summary, 20))
  
  # Save results
  write.csv(final_summary, "dtu_switching_results_summary.csv", row.names = FALSE)
  print("Results saved to: dtu_switching_results_summary.csv")
  
} else {
  print("No switching results available for summary table creation")
}
