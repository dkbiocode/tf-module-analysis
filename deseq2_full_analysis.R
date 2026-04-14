#!/usr/bin/env Rscript

# Complete DESeq2 analysis for Dineen data
# Performs all 2-way comparisons: elt2d, elt7d, and double mutant vs WT

library(DESeq2)
library(tibble)
library(dplyr)

cat(paste0(strrep("=", 70), "\n"))
cat("DESeq2 Differential Expression Analysis - Dineen Data\n")
cat(paste0(strrep("=", 70), "\n\n"))

# Read raw counts
cat("Reading Dineen_raw_counts.csv...\n")
dineen_raw <- read.csv("Dineen_raw_counts.csv", row.names=NULL)

# Extract count matrix
count_cols <- c("wt_sorted_1", "wt_sorted_2", "wt_sorted_3", "wt_sorted_4",
                "elt7d_sorted_1", "elt7d_sorted_2", "elt7d_sorted_3",
                "elt2d_sorted_1", "elt2d_sorted_2", "elt2d_sorted_3", "elt2d_sorted_4",
                "elt2delt7d_sorted_1", "elt2delt7d_sorted_2", "elt2delt7d_sorted_3")

count_matrix <- dineen_raw[, count_cols]
rownames(count_matrix) <- dineen_raw$wbid

cat(sprintf("Count matrix: %d genes x %d samples\n\n",
            nrow(count_matrix), ncol(count_matrix)))

# Create sample metadata
sample_info <- data.frame(
  sample = colnames(count_matrix),
  genotype = c(rep("wt", 4),
               rep("elt7d", 3),
               rep("elt2d", 4),
               rep("elt2d_elt7d", 3))
)
rownames(sample_info) <- sample_info$sample
sample_info$genotype <- factor(sample_info$genotype,
                               levels = c("wt", "elt2d", "elt7d", "elt2d_elt7d"))

# Create DESeq2 dataset
cat("Creating DESeq2 dataset...\n")
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info,
  design = ~ genotype
)

# Filter low-count genes
cat("Filtering low-count genes (>= 10 reads total)...\n")
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat(sprintf("Retained %d genes\n\n", sum(keep)))

# Run DESeq2
cat("Running DESeq2 analysis...\n")
cat("(This may take a few minutes...)\n\n")
dds <- DESeq(dds)

# Function to extract and annotate results
extract_results <- function(dds, contrast_name, numerator, denominator,
                           alpha = 0.05, lfc_threshold = 1) {
  cat(sprintf("Extracting results for: %s vs %s\n", numerator, denominator))

  res <- results(dds, contrast = c("genotype", numerator, denominator), alpha = alpha)

  # Convert to data frame
  res_df <- as.data.frame(res) %>%
    rownames_to_column(var = "wbid") %>%
    mutate(
      contrast = contrast_name,
      # Classify genes
      status = case_when(
        is.na(padj) ~ "NA",
        padj >= alpha ~ "not_sig",
        padj < alpha & abs(log2FoldChange) < lfc_threshold ~ "sig_small_FC",
        padj < alpha & log2FoldChange >= lfc_threshold ~ "up",
        padj < alpha & log2FoldChange <= -lfc_threshold ~ "down",
        TRUE ~ "other"
      )
    ) %>%
    arrange(padj)

  # Print summary
  cat(sprintf("  Total genes: %d\n", nrow(res_df)))
  cat(sprintf("  Significant (padj < %.2f): %d\n", alpha, sum(res_df$padj < alpha, na.rm=TRUE)))
  cat(sprintf("  Up (LFC >= %d): %d\n", lfc_threshold, sum(res_df$status == "up", na.rm=TRUE)))
  cat(sprintf("  Down (LFC <= -%d): %d\n", lfc_threshold, sum(res_df$status == "down", na.rm=TRUE)))
  cat(sprintf("  Sig but |LFC| < %d: %d\n", lfc_threshold, sum(res_df$status == "sig_small_FC", na.rm=TRUE)))
  cat("\n")

  return(res_df)
}

# Perform all three comparisons
cat(paste0(strrep("=", 70), "\n"))
cat("CONTRAST 1: ELT-2Δ vs WT\n")
cat(paste0(strrep("=", 70), "\n"))
res_elt2d <- extract_results(dds, "elt2d_vs_wt", "elt2d", "wt")

cat(paste0(strrep("=", 70), "\n"))
cat("CONTRAST 2: ELT-7Δ vs WT\n")
cat(paste0(strrep("=", 70), "\n"))
res_elt7d <- extract_results(dds, "elt7d_vs_wt", "elt7d", "wt")

cat(paste0(strrep("=", 70), "\n"))
cat("CONTRAST 3: ELT-2Δ/ELT-7Δ (double) vs WT\n")
cat(paste0(strrep("=", 70), "\n"))
res_double <- extract_results(dds, "double_vs_wt", "elt2d_elt7d", "wt")

# Save individual results
write.csv(res_elt2d, "Dineen_DESeq2_elt2d_vs_wt.csv", row.names = FALSE)
write.csv(res_elt7d, "Dineen_DESeq2_elt7d_vs_wt.csv", row.names = FALSE)
write.csv(res_double, "Dineen_DESeq2_double_vs_wt.csv", row.names = FALSE)

cat("Saved individual contrast files:\n")
cat("  - Dineen_DESeq2_elt2d_vs_wt.csv\n")
cat("  - Dineen_DESeq2_elt7d_vs_wt.csv\n")
cat("  - Dineen_DESeq2_double_vs_wt.csv\n\n")

# Combine all results into one file
all_results <- bind_rows(res_elt2d, res_elt7d, res_double)
write.csv(all_results, "Dineen_DESeq2_all_contrasts.csv", row.names = FALSE)
cat("Saved combined file: Dineen_DESeq2_all_contrasts.csv\n\n")

# Compare with existing results
cat(paste0(strrep("=", 70), "\n"))
cat("COMPARISON WITH EXISTING RESULTS\n")
cat(paste0(strrep("=", 70), "\n"))

if (file.exists("Dineen_DESeq2_results.csv")) {
  cat("Reading existing Dineen_DESeq2_results.csv...\n")
  existing <- read.csv("Dineen_DESeq2_results.csv")

  cat(sprintf("Existing results: %d genes\n", nrow(existing)))
  cat("Descriptions in existing file:\n")
  print(table(existing$description))
  cat("\n")

  # Merge with new elt2d results for comparison
  comparison <- existing %>%
    select(wbid,
           old_log2FC = log2foldchange,
           old_padj = padj,
           old_status = status,
           old_description = description) %>%
    left_join(
      res_elt2d %>% select(wbid, new_log2FC = log2FoldChange, new_padj = padj, new_status = status),
      by = "wbid"
    )

  # Calculate correlation
  cor_lfc <- cor(comparison$old_log2FC, comparison$new_log2FC, use = "complete.obs")
  cor_padj <- cor(-log10(comparison$old_padj), -log10(comparison$new_padj), use = "complete.obs")

  cat(sprintf("\nCorrelation between old and new ELT2Δ results:\n"))
  cat(sprintf("  Log2FC correlation: %.4f\n", cor_lfc))
  cat(sprintf("  -log10(padj) correlation: %.4f\n", cor_padj))

  # Check agreement on significant genes
  comparison <- comparison %>%
    mutate(
      old_sig = old_padj < 0.05,
      new_sig = new_padj < 0.05,
      agreement = old_sig == new_sig
    )

  cat(sprintf("\nAgreement on significance (padj < 0.05):\n"))
  cat(sprintf("  Both significant: %d\n", sum(comparison$old_sig & comparison$new_sig, na.rm=TRUE)))
  cat(sprintf("  Both not significant: %d\n", sum(!comparison$old_sig & !comparison$new_sig, na.rm=TRUE)))
  cat(sprintf("  Disagreement: %d\n", sum(!comparison$agreement, na.rm=TRUE)))
  cat(sprintf("  Agreement rate: %.1f%%\n",
              100 * sum(comparison$agreement, na.rm=TRUE) / sum(!is.na(comparison$agreement))))

  write.csv(comparison, "Dineen_DESeq2_comparison_old_vs_new.csv", row.names = FALSE)
  cat("\nSaved comparison: Dineen_DESeq2_comparison_old_vs_new.csv\n")
} else {
  cat("No existing results file found for comparison.\n")
}

cat("\n")
cat(paste0(strrep("=", 70), "\n"))
cat("ANALYSIS COMPLETE\n")
cat(paste0(strrep("=", 70), "\n"))
