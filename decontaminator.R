# ============================================
# Decontamination Analysis Script
# ============================================

# Load required libraries
library(tidyverse)
library(microDecon)
library(readxl)

# ============================================
# FILE READING FUNCTIONS
# ============================================

read_otu_file <- function(filepath, filename) {
  ext <- tools::file_ext(filename) |> tolower()
  
  if(ext %in% c("csv")) {
    read.csv(filepath, header = TRUE, check.names = FALSE)
  } else if(ext %in% c("txt", "tsv")) {
    read.delim(filepath, header = TRUE, check.names = FALSE)
  } else if(ext %in% c("xls", "xlsx")) {
    readxl::read_excel(filepath)
  } else {
    stop("Unsupported file type: ", ext)
  }
}

# ============================================
# DATA PROCESSING FUNCTIONS
# ============================================

extract_otu_and_metadata <- function(df, otu_id_col, otu_start_col, otu_end_col) {
  col_names <- colnames(df)
  
  # Find column indices
  otu_id_idx <- which(col_names == otu_id_col)
  otu_start_idx <- which(col_names == otu_start_col)
  otu_end_idx <- which(col_names == otu_end_col)
  
  if(length(otu_id_idx) == 0) stop("OTU ID column not found: ", otu_id_col)
  if(length(otu_start_idx) == 0) stop("Start column not found: ", otu_start_col)
  if(length(otu_end_idx) == 0) stop("End column not found: ", otu_end_col)
  if(otu_start_idx > otu_end_idx) stop("Start column must come before end column")
  
  # OTU count columns
  otu_col_indices <- otu_start_idx:otu_end_idx
  
  # Metadata columns (everything except OTU ID and OTU count columns)
  all_indices <- 1:ncol(df)
  metadata_indices <- setdiff(all_indices, c(otu_id_idx, otu_col_indices))
  
  # Build OTU table (OTU ID + counts)
  otu_table <- df[, c(otu_id_idx, otu_col_indices), drop = FALSE]
  colnames(otu_table)[1] <- "OTU_ID"
  
  # Build metadata (if any)
  metadata <- if(length(metadata_indices) > 0) {
    df[, c(otu_id_idx, metadata_indices), drop = FALSE]
  } else {
    NULL
  }
  if(!is.null(metadata)) colnames(metadata)[1] <- "OTU_ID"
  
  list(
    otu_table = otu_table,
    metadata = metadata,
    otu_id_col_name = otu_id_col,
    otu_col_names = col_names[otu_col_indices],
    metadata_col_names = if(length(metadata_indices) > 0) col_names[metadata_indices] else character(0),
    num_otu_cols = length(otu_col_indices),
    num_metadata_cols = length(metadata_indices)
  )
}

prepare_for_decontamination <- function(otutab, blank_prefixes = c("CB"), ntc_prefixes = c("NTC")) {
  # Assumes first column is OTU ID, all other columns are counts
  
  otu_col <- 1
  
  # Find extraction blank columns
  blank_pattern <- paste0("(", paste(blank_prefixes, collapse="|"), ")")
  extraction_blanks <- which(grepl(blank_pattern, colnames(otutab), ignore.case = TRUE))
  
  # Find PCR blank columns (NTCs) - these will ALSO be treated as blanks
  ntc_pattern <- paste0("(", paste(ntc_prefixes, collapse="|"), ")")
  pcr_blanks <- which(grepl(ntc_pattern, colnames(otutab), ignore.case = TRUE))
  
  # COMBINE all blanks for decontamination
  all_blanks <- c(extraction_blanks, pcr_blanks)
  
  # All other columns are samples
  other_cols <- setdiff(2:ncol(otutab), all_blanks)
  
  # Order: OTU ID, then ALL blanks, then samples
  otutab_ordered <- otutab[, c(otu_col, all_blanks, other_cols)]
  
  OTU_ID <- as.character(otutab_ordered[[1]])
  counts <- otutab_ordered[, -1]
  
  counts <- counts %>% mutate(across(everything(), ~as.numeric(as.character(.))))
  counts[is.na(counts)] <- 0
  
  # Keep all blanks; remove zero-sum samples
  blank_keep <- if(length(all_blanks) > 0) colnames(counts)[1:length(all_blanks)] else character(0)
  sample_keep <- if(length(other_cols) > 0) colnames(counts)[-(1:length(all_blanks))] else character(0)
  
  if(length(sample_keep) > 0) {
    sample_keep <- sample_keep[colSums(counts[, sample_keep, drop=FALSE]) != 0]
  }
  
  counts <- counts[, c(blank_keep, sample_keep), drop=FALSE]
  
  # Remove zero-sum OTUs (rows with all zeros)
  nonzero_idx <- rowSums(counts) > 0
  OTU_ID <- OTU_ID[nonzero_idx]
  counts <- counts[nonzero_idx, , drop=FALSE]
  
  otu_table_only <- cbind(OTU_ID, counts)
  colnames(otu_table_only) <- make.unique(colnames(otu_table_only))
  
  list(
    otu_table_only = otu_table_only,
    num_blanks = length(blank_keep),
    num_samples = length(sample_keep),
    blank_cols = blank_keep,
    sample_cols = sample_keep,
    extraction_blank_matches = colnames(otutab)[extraction_blanks],
    pcr_blank_matches = colnames(otutab)[pcr_blanks],
    all_blank_matches = colnames(otutab)[all_blanks]
  )
}

# Simple subtraction method (sum only) with zero-row removal
simple_blank_subtract <- function(df, num_blanks) {
  
  OTU_ID <- df[,1]
  counts <- df[,-1, drop=FALSE]
  
  blank_cols <- counts[, 1:num_blanks, drop=FALSE]
  sample_cols <- counts[, (num_blanks+1):ncol(counts), drop=FALSE]
  
  # Compute blank signal per OTU using SUM (total reads in blanks)
  blank_signal <- rowSums(blank_cols)
  
  # Subtract from each sample
  corrected_samples <- sweep(sample_cols, 1, blank_signal, FUN = "-")
  
  # Prevent negative values
  corrected_samples[corrected_samples < 0] <- 0
  
  # Track removed reads
  removed_reads <- sample_cols - corrected_samples
  
  # Ensure integers
  corrected_samples <- round(corrected_samples)
  removed_reads <- round(removed_reads)
  
  # Remove rows where all corrected sample counts are zero
  nonzero_rows <- rowSums(corrected_samples) > 0
  OTU_ID <- OTU_ID[nonzero_rows]
  corrected_samples <- corrected_samples[nonzero_rows, , drop=FALSE]
  removed_reads <- removed_reads[nonzero_rows, , drop=FALSE]
  
  # Build outputs
  corrected_table <- cbind(OTU_ID, corrected_samples)
  removed_table <- cbind(OTU_ID, removed_reads)
  
  list(
    decon.table = corrected_table,
    reads.removed = removed_table
  )
}

# ============================================
# MAIN DECONTAMINATION FUNCTION
# ============================================

run_decontamination <- function(otu_file,
                                otu_id_col,
                                otu_start_col,
                                otu_end_col,
                                blank_prefixes = c("CB", "EXT", "EB"),
                                ntc_prefixes = c("NTC", "PCR", "PB"),
                                decon_method = "microdecon",  # "microdecon" or "simple"
                                # microDecon parameters
                                runs = 2,
                                thresh = 0.7,
                                prop_thresh = 0.00005,
                                regression_type = 0,  # 0=auto, 1=aggressive, 2=conservative
                                # Output options
                                output_dir = ".",
                                output_prefix = "decontaminated") {
  
  # Create output directory if it doesn't exist
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Step 1: Read OTU table
  cat("Reading OTU table...\n")
  filename <- basename(otu_file)
  uploaded_data <- read_otu_file(otu_file, filename)
  
  # Step 2: Extract OTU table and metadata
  cat("Extracting OTU table and metadata based on column selections...\n")
  extracted <- extract_otu_and_metadata(uploaded_data, otu_id_col, otu_start_col, otu_end_col)
  
  cat(sprintf("  OTU count columns: %d\n", extracted$num_otu_cols))
  cat(sprintf("  Metadata columns: %d\n", extracted$num_metadata_cols))
  
  # Step 3: Prepare for decontamination
  cat("\nPreparing for decontamination...\n")
  cat(sprintf("  Blank prefixes: %s\n", paste(blank_prefixes, collapse=", ")))
  cat(sprintf("  NTC prefixes: %s\n", paste(ntc_prefixes, collapse=", ")))
  
  prep <- prepare_for_decontamination(extracted$otu_table, blank_prefixes, ntc_prefixes)
  
  cat(sprintf("\nControl Detection:\n"))
  cat(sprintf("  Extraction blanks detected: %d\n", length(prep$extraction_blank_matches)))
  if(length(prep$extraction_blank_matches) > 0) {
    for(col in prep$extraction_blank_matches) {
      cat(sprintf("    - %s\n", col))
    }
  }
  
  cat(sprintf("  PCR blanks detected: %d\n", length(prep$pcr_blank_matches)))
  if(length(prep$pcr_blank_matches) > 0) {
    for(col in prep$pcr_blank_matches) {
      cat(sprintf("    - %s\n", col))
    }
  }
  
  cat(sprintf("\n  TOTAL BLANKS for decontamination: %d\n", prep$num_blanks))
  cat(sprintf("  Total samples: %d\n", prep$num_samples))
  
  if(prep$num_blanks == 0) {
    stop("No blank samples detected! Check your naming settings.")
  }
  
  data_for_decon <- prep$otu_table_only
  numb_ind <- ncol(data_for_decon) - 1 - prep$num_blanks
  
  # Step 4: Run selected decontamination method
  cat(sprintf("\nRunning %s decontamination...\n", 
              if(decon_method == "simple") "simple blank subtraction" else "microDecon"))
  
  if (decon_method == "simple") {
    decon_res <- simple_blank_subtract(
      df = data_for_decon,
      num_blanks = prep$num_blanks
    )
  } else {
    decon_res <- microDecon::decon(
      data = data_for_decon,
      numb.blanks = prep$num_blanks,
      numb.ind = numb_ind,
      taxa = FALSE,
      runs = as.integer(runs),
      thresh = as.numeric(thresh),
      prop.thresh = as.numeric(prop_thresh),
      regression = as.numeric(regression_type)
    )
  }
  
  cleaned <- decon_res$decon.table
  
  # Remove blank columns from final output (only for microDecon)
  if (decon_method == "microdecon") {
    cleaned <- cleaned[, !colnames(cleaned) %in% c(prep$blank_cols, "Mean.blank"), drop=FALSE]
  }
  
  # Store removed reads and merge with metadata
  removed_reads_df <- decon_res$reads.removed
  if(!is.null(removed_reads_df) && !is.null(extracted$metadata)) {
    removed_reads_df <- left_join(removed_reads_df, extracted$metadata, by = "OTU_ID")
  }
  
  colnames(cleaned)[1] <- "OTU_ID"
  
  # Reattach metadata
  final_table <- if(!is.null(extracted$metadata)) {
    left_join(cleaned, extracted$metadata, by="OTU_ID")
  } else {
    cleaned
  }
  
  # Step 5: Save outputs
  cat("\nSaving outputs...\n")
  
  # Generate output filenames
  method_tag <- if(decon_method == "simple") "simple_subtraction" else "microdecon"
  decon_file <- file.path(output_dir, paste0(output_prefix, "_", method_tag, "_decontaminated.csv"))
  removed_file <- file.path(output_dir, paste0(output_prefix, "_", method_tag, "_removed_reads.csv"))
  
  # Write files
  write.csv(final_table, decon_file, row.names = FALSE)
  write.csv(removed_reads_df, removed_file, row.names = FALSE)
  
  # Step 6: Print summary
  cat("\n========== DECONTAMINATION SUMMARY ==========\n")
  cat(sprintf("Method: %s\n", if(decon_method == "simple") "Simple blank subtraction" else "microDecon"))
  
  if(decon_method == "microdecon") {
    cat(sprintf("  Runs: %d\n", runs))
    cat(sprintf("  Threshold (proportion zeros): %.2f\n", thresh))
    cat(sprintf("  Proportion threshold: %.6f\n", prop_thresh))
    regression_names <- c("Auto", "Aggressive", "Conservative")
    cat(sprintf("  Regression type: %s\n", regression_names[regression_type + 1]))
  }
  
  cat(sprintf("\nInput statistics:\n"))
  cat(sprintf("  Original OTUs: %d\n", nrow(uploaded_data)))
  cat(sprintf("  Original samples: %d\n", ncol(uploaded_data) - 1))
  cat(sprintf("  Blanks used: %d\n", prep$num_blanks))
  cat(sprintf("  Samples retained: %d\n", prep$num_samples))
  
  cat(sprintf("\nOutput statistics:\n"))
  cat(sprintf("  OTUs retained: %d\n", nrow(cleaned)))
  cat(sprintf("  OTUs removed (all zeros after decontamination): %d\n", 
              nrow(data_for_decon) - nrow(cleaned)))
  
  cat(sprintf("\nOutput files saved:\n"))
  cat(sprintf("  Decontaminated table: %s\n", decon_file))
  cat(sprintf("  Removed reads table: %s\n", removed_file))
  cat("============================================\n")
  
  # Return results as a list
  return(invisible(list(
    decontaminated_table = final_table,
    removed_reads = removed_reads_df,
    decon_file = decon_file,
    removed_file = removed_file,
    prep_info = prep,
    extracted_info = extracted
  )))
}

# ============================================
# USER SETTINGS - MODIFY THESE BEFORE RUNNING
# ============================================

# File path (use file.choose() for interactive selection)
cat("Please select your OTU table file...\n")
otu_file <- file.choose()

# Column names in your OTU table (modify these to match your file)
otu_id_col <- "OTU_ID"        # Column containing OTU/ASV identifiers
otu_start_col <- "Sample01"    # First column with count data (replace with actual column name)
otu_end_col <- "Sample64"      # Last column with count data (replace with actual column name)

# Control sample detection (case-insensitive, matches anywhere in column name)
blank_prefixes <- c("CB", "EXT", "EB")     # Extraction blank prefixes
ntc_prefixes <- c("NTC", "PCR", "PB")      # PCR blank/negative prefixes

# Decontamination method
decon_method <- "microdecon"   # Options: "microdecon" or "simple" for simple subtraction of blank reads from sample reads (skips microdecon)

# microDecon parameters (only used if decon_method = "microdecon")
runs <- 2                      # Number of decontamination iterations
thresh <- 0.7                  # Zero-prevalence threshold (0-1)
prop_thresh <- 0.00005         # Proportion threshold
regression_type <- 0           # 0=Auto, 1=Aggressive, 2=Conservative

# Simple subtraction parameters (only used if decon_method = "simple")
# No additional parameters needed - uses sum of all blanks

# Output settings
output_dir <- "."              # Directory to save results (current directory)
output_prefix <- "decontaminated"     # Prefix for output files

# ============================================
# RUN DECONTAMINATION
# ============================================

results <- run_decontamination(
  otu_file = otu_file,
  otu_id_col = otu_id_col,
  otu_start_col = otu_start_col,
  otu_end_col = otu_end_col,
  blank_prefixes = blank_prefixes,
  ntc_prefixes = ntc_prefixes,
  decon_method = decon_method,
  runs = runs,
  thresh = thresh,
  prop_thresh = prop_thresh,
  regression_type = regression_type,
  output_dir = output_dir,
  output_prefix = output_prefix
)

# If you want to work with the results in R:
# decontaminated_data <- results$decontaminated_table
# removed_reads_data <- results$removed_reads