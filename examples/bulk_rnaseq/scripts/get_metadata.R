#!/usr/bin/env Rscript

# =========================================================
# get_metadata.R
#
# Purpose:
#   Generates a DESeq2-compatible sample metadata table by
#   matching SRA metadata (SraRunTable) with RNA-seq count
#   matrix sample names.
#
# Inputs:
#   - counts/gene_counts.txt (featureCounts output)
#   - counts/SraRunTable.csv (NCBI SRA metadata export)
#
# Outputs:
#   - counts/sample_info.txt (clean metadata for DESeq2)
#
# Key Steps:
#   1. Extract SRR sample IDs from count matrix
#   2. Load SRA metadata table
#   3. Standardize condition labels (control, treatment, etc.)
#   4. Join metadata to expression samples
#   5. Validate completeness and ordering
#
# Notes:
#   - Assumes SRR IDs are consistent between count matrix
#     and SRA metadata
#   - Automatically normalizes treatment labels
#   - Required for downstream DESeq2 analysis
#
# =========================================================

# --- PERSONAL LIBRARY ---
user_lib <- file.path(Sys.getenv("HOME"), "R", "library")
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)

.libPaths(c(user_lib, .libPaths()))

# --- INSTALL & LOAD ---
install_load <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE, lib.loc = user_lib)) {
      install.packages(pkg, lib = user_lib, repos = "https://cloud.r-project.org")
    }
    library(pkg, character.only = TRUE, lib.loc = user_lib)
  }
}

install_load(c("dplyr", "readr", "stringr"))

# --- READ COUNT MATRIX ---
count_data <- read.delim("counts/gene_counts.txt",
                         comment.char = "#",
                         check.names = FALSE)

# Ensure Geneid exists
if (!"Geneid" %in% colnames(count_data)) {
  stop("Geneid column not found in count file")
}

rownames(count_data) <- count_data$Geneid

# Extract sample columns (featureCounts format)
sample_cols <- colnames(count_data)[7:ncol(count_data)]

# Clean sample names (remove .bam if present)
sample_names <- stringr::str_extract(sample_cols, "SRR[0-9]+")

cat("Detected samples:\n")
print(sample_names)

# --- READ SRA METADATA ---
metadata_file <- "counts/SraRunTable.csv"
sra_meta <- read.csv(metadata_file, stringsAsFactors = FALSE)

# Check required columns
required_cols <- c("Run", "treatment")
missing_cols <- setdiff(required_cols, colnames(sra_meta))
if (length(missing_cols) > 0) {
  stop(paste("Missing columns in metadata:", paste(missing_cols, collapse=", ")))
}

# --- CLEAN CONDITION LABELS ---
# Automatically simplify treatment labels
sra_meta$condition <- tolower(sra_meta$treatment)

# Normalize common patterns
sra_meta$condition <- gsub(" shRNA", "", sra_meta$condition, ignore.case = TRUE)
sra_meta$condition <- gsub(" knockdown", "", sra_meta$condition, ignore.case = TRUE)

# Standardize control naming
sra_meta$condition[grepl("control", sra_meta$condition)] <- "control"

# --- MATCH TO COUNT MATRIX ---
sample_info <- data.frame(sample = sample_names, stringsAsFactors = FALSE)

sample_info <- dplyr::left_join(
  sample_info,
  sra_meta[, c("Run", "condition", "treatment", "batch", "cell_line")],
  by = c("sample" = "Run")
)

# --- SAFETY CHECKS ---
if (any(is.na(sample_info$condition))) {
  warning("Some samples missing condition info:")
  print(sample_info[is.na(sample_info$condition), ])
}

# Ensure order matches count matrix
sample_info <- sample_info[match(sample_names, sample_info$sample), ]

# Set rownames (required for DESeq2)
rownames(sample_info) <- sample_info$sample

# Convert to factor (important for DESeq2)
sample_info$condition <- factor(sample_info$condition)

# --- OUTPUT ---
write.table(sample_info,
            "counts/sample_info.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

cat("\nFinal sample info:\n")
print(sample_info)