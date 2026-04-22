#!/usr/bin/env Rscript

# =========================================================
# run_deseq2.R
#
# Purpose:
#   Performs differential gene expression analysis using
#   DESeq2 and pathway enrichment using FGSEA.
#
# Inputs:
#   - counts/gene_counts.txt (featureCounts output)
#   - counts/sample_info.txt (metadata table)
#
# Outputs:
#   - results/tables/ (DE results, normalized counts, DEGs)
#   - results/plots/ (PCA, MA plot, volcano plot)
#   - results/gsea/ (FGSEA results + enrichment plots)
#
# Methods:
#   - DESeq2 with lfcShrink (apeglm)
#   - Variance stabilizing transform (VST)
#   - FGSEA using MSigDB Hallmark gene sets
#
# Assumptions:
#   - Unstranded RNA-seq
#   - Gene IDs compatible with MSigDB Ensembl mapping
#
# =========================================================

# --- LIBRARY SETUP ---
user_lib <- file.path(Sys.getenv("HOME"), "R", "library")
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))

install_load <- function(pkgs, bioc = FALSE) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (bioc) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager", lib = user_lib,
                           repos = "https://cloud.r-project.org")
        }
        BiocManager::install(pkg, lib = user_lib,
                             ask = FALSE, update = FALSE)
      } else {
        install.packages(pkg, lib = user_lib,
                         repos = "https://cloud.r-project.org")
      }
    }
    library(pkg, character.only = TRUE, lib.loc = user_lib)
  }
}

# --- PACKAGES ---
install_load(c("DESeq2"), bioc = TRUE)
install_load(c("readr", "dplyr", "ggplot2",
               "pheatmap", "RColorBrewer", "stringr"))
install_load(c("fgsea"), bioc = TRUE)
install_load(c("msigdbr"), bioc = TRUE)

# --- INPUTS ---
counts_file <- "counts/gene_counts.txt"
sample_info_file <- "counts/sample_info.txt"

# --- OUTPUT ---
output_dir <- "results"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "tables"), showWarnings = FALSE)
dir.create(file.path(output_dir, "plots"), showWarnings = FALSE)
dir.create(file.path(output_dir, "gsea"), showWarnings = FALSE)

# --- LOAD COUNTS ---
raw_counts <- read.delim(counts_file, comment.char = "#", check.names = FALSE)

count_data <- raw_counts[, -(1:6)]
rownames(count_data) <- raw_counts$Geneid

# --- CLEAN SAMPLE NAMES ---

sample_names <- basename(colnames(count_data))
sample_names <- sub("\\.bam$", "", sample_names)
sample_names <- sub("\\.sorted.*", "", sample_names)
sample_names <- trimws(sample_names)

colnames(count_data) <- sample_names

if (any(!grepl("^SRR[0-9]+$", sample_names))) {
  stop("Sample names malformed")
}

# --- SAMPLE INFO ---
sample_info <- read.delim(sample_info_file, stringsAsFactors = FALSE)
rownames(sample_info) <- sample_info$sample

if (!all(sample_names %in% rownames(sample_info))) {
  stop("Sample mismatch between counts and metadata")
}

sample_info <- sample_info[sample_names, , drop = FALSE]
sample_info$condition <- factor(sample_info$condition)
sample_info$condition <- relevel(sample_info$condition, ref = "control")

# --- DESEQ2 ---
dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = sample_info,
  design = ~ condition
)

dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

coef_name <- resultsNames(dds)[2]
res <- lfcShrink(dds, coef = coef_name, type = "apeglm")

# --- RESULTS ---
res_df <- as.data.frame(res)
res_df <- res_df[!is.na(res_df$padj), ]
res_df <- res_df[order(res_df$padj), ]

write.csv(res_df, file.path(output_dir, "tables/deseq2_results.csv"))

norm_counts <- counts(dds, normalized = TRUE)
write.csv(as.data.frame(norm_counts),
          file.path(output_dir, "tables/normalized_counts.csv"))

# --- PCA ---
vsd <- vst(dds, blind = FALSE)
pca <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca, "percentVar"))

p <- ggplot(pca, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal()

ggsave(file.path(output_dir, "plots/PCA_plot.png"), p, width = 6, height = 5)

# --- MA PLOT ---
png(file.path(output_dir, "plots/MA_plot.png"))
plotMA(res, ylim = c(-5, 5))
dev.off()

# --- DEG FILTER ---
sig_all <- res_df[
  res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, ]

write.csv(sig_all,
          file.path(output_dir, "tables/DEGs_significant.csv"))

# --- VOLCANO ---
res_df$significance <- "NS"
res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Up"
res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Down"

padj_safe <- ifelse(res_df$padj == 0, 1e-300, res_df$padj)

volcano <- ggplot(res_df, aes(log2FoldChange, -log10(padj_safe), color = significance)) +
  geom_point(alpha = 0.6) +
  theme_minimal()

ggsave(file.path(output_dir, "plots/volcano_plot.png"), volcano)

# --- SUMMARY ---
cat("\n========================\nDESEQ2 SUMMARY\n========================\n")
cat("Genes tested:", nrow(res_df), "\n")
cat("Significant:", nrow(sig_all), "\n")
cat("Up:", sum(res_df$significance == "Up"), "\n")
cat("Down:", sum(res_df$significance == "Down"), "\n")

# --- GSEA ---
cat("\nRunning FGSEA...\n")

gene_ids <- gsub("\\..*", "", rownames(res_df))
gene_list <- res_df$log2FoldChange

valid <- !is.na(gene_list) & !is.na(gene_ids)

gene_list <- gene_list[valid]
gene_ids <- gene_ids[valid]
names(gene_list) <- gene_ids

gene_list <- gene_list[!duplicated(names(gene_list))]
gene_list <- gene_list * (-log10(res_df$padj[valid] + 1e-10))
gene_list <- gene_list[!is.na(gene_list)]
gene_list <- sort(gene_list, decreasing = TRUE)

msig <- msigdbr(species = "Homo sapiens", collection = "H")
pathways <- split(msig$ensembl_gene, msig$gs_name)

fgsea_res <- fgsea(
  pathways = pathways,
  stats = gene_list,
  minSize = 10,
  maxSize = 500
)

fgsea_res <- fgsea_res[order(fgsea_res$padj), ]

fgsea_out <- as.data.frame(fgsea_res)

if ("leadingEdge" %in% colnames(fgsea_out)) {
  fgsea_out$leadingEdge <- sapply(fgsea_out$leadingEdge, paste, collapse = ";")
}

write.csv(
  fgsea_out,
  file.path(output_dir, "gsea/fgsea_results.csv"),
  row.names = FALSE
)

# --- GSEA PLOTS ---
top10 <- head(fgsea_res, 10)

png(file.path(output_dir, "gsea/fgsea_top10_barplot.png"), width = 900, height = 600)

par(mar = c(5, 14, 4, 2))  # <-- FIX: more space for labels

barplot(
  rev(top10$NES),
  names.arg = rev(top10$pathway),
  las = 2,
  horiz = TRUE,
  col = ifelse(rev(top10$NES) > 0, "steelblue", "tomato"),
  main = "Top FGSEA Pathways (NES)"
)

dev.off()

cat("\nAnalysis complete. Results in:", output_dir, "\n")