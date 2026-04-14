#!/usr/bin/env Rscript

# DESeq2 Interaction Analysis - Double mutant vs ELT-2Δ
# Tests whether ELT-7 loss enhances the ELT-2Δ phenotype (buffering effect)

library(DESeq2)
library(tibble)
library(dplyr)

cat(paste0(strrep("=", 70), "\n"))
cat("DESeq2 Interaction Analysis - ELT-7 Buffering Effect\n")
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
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat(sprintf("Retained %d genes\n\n", sum(keep)))

# Run DESeq2
cat("Running DESeq2 analysis...\n\n")
dds <- DESeq(dds)

# Extract the key contrast: Double vs ELT-2Δ
cat(paste0(strrep("=", 70), "\n"))
cat("CONTRAST: Double Mutant vs ELT-2Δ (ELT-7 buffering effect)\n")
cat(paste0(strrep("=", 70), "\n"))

res_interaction <- results(dds, contrast = c("genotype", "elt2d_elt7d", "elt2d"), alpha = 0.05)

# Convert to data frame
res_df <- as.data.frame(res_interaction) %>%
  rownames_to_column(var = "wbid") %>%
  mutate(
    status = case_when(
      is.na(padj) ~ "NA",
      padj >= 0.05 ~ "not_sig",
      padj < 0.05 & abs(log2FoldChange) < 0.5 ~ "sig_small_FC",
      padj < 0.05 & log2FoldChange >= 0.5 ~ "enhanced_in_double",
      padj < 0.05 & log2FoldChange <= -0.5 ~ "suppressed_in_double",
      TRUE ~ "other"
    )
  ) %>%
  arrange(padj)

# Print summary
cat(sprintf("\nTotal genes: %d\n", nrow(res_df)))
cat(sprintf("Significant (padj < 0.05): %d\n", sum(res_df$padj < 0.05, na.rm=TRUE)))
cat(sprintf("Enhanced in double (LFC >= 0.5): %d\n", sum(res_df$status == "enhanced_in_double", na.rm=TRUE)))
cat(sprintf("Suppressed in double (LFC <= -0.5): %d\n", sum(res_df$status == "suppressed_in_double", na.rm=TRUE)))
cat(sprintf("Sig but |LFC| < 0.5: %d\n\n", sum(res_df$status == "sig_small_FC", na.rm=TRUE)))

# Save results
write.csv(res_df, "Dineen_DESeq2_double_vs_elt2d.csv", row.names = FALSE)
cat("Saved: Dineen_DESeq2_double_vs_elt2d.csv\n\n")

# Now load the previous contrasts to create a comprehensive view
cat(paste0(strrep("=", 70), "\n"))
cat("COMPREHENSIVE ANALYSIS: ELT-7 Buffering Patterns\n")
cat(paste0(strrep("=", 70), "\n\n"))

# Load previous results
elt2d_results <- read.csv("Dineen_DESeq2_elt2d_vs_wt.csv")
double_results <- read.csv("Dineen_DESeq2_double_vs_wt.csv")

# Merge all three for comparison
buffering_analysis <- elt2d_results %>%
  select(wbid, elt2d_log2FC = log2FoldChange, elt2d_padj = padj) %>%
  left_join(
    double_results %>% select(wbid, double_log2FC = log2FoldChange, double_padj = padj),
    by = "wbid"
  ) %>%
  left_join(
    res_df %>% select(wbid, interaction_log2FC = log2FoldChange, interaction_padj = padj),
    by = "wbid"
  ) %>%
  mutate(
    # Classify buffering patterns
    elt2d_sig = elt2d_padj < 0.05,
    double_sig = double_padj < 0.05,
    interaction_sig = interaction_padj < 0.05,

    # Calculate fold-change magnitude differences
    fc_magnitude_elt2d = abs(elt2d_log2FC),
    fc_magnitude_double = abs(double_log2FC),
    magnitude_increase = fc_magnitude_double - fc_magnitude_elt2d,

    # Classify genes by buffering pattern
    buffering_pattern = case_when(
      # Strong ELT-7 buffering: enhanced effect in double
      interaction_sig & interaction_log2FC > 0.5 & elt2d_sig ~ "ELT7_buffers_upregulation",
      interaction_sig & interaction_log2FC < -0.5 & elt2d_sig ~ "ELT7_buffers_downregulation",

      # Opposite effects (compensation)
      interaction_sig & abs(interaction_log2FC) > 0.5 &
        sign(elt2d_log2FC) != sign(interaction_log2FC) ~ "ELT7_compensatory",

      # No additional effect
      !interaction_sig & elt2d_sig ~ "ELT2_dominant_no_ELT7_effect",

      # Only significant in double
      !elt2d_sig & double_sig & interaction_sig ~ "ELT7_reveals_new_targets",

      TRUE ~ "other"
    )
  ) %>%
  arrange(desc(abs(interaction_log2FC)))

# Summary statistics
cat("Buffering Pattern Summary:\n")
pattern_table <- table(buffering_analysis$buffering_pattern)
print(pattern_table)
cat("\n")

# Identify top buffered genes (greatest enhancement in double)
top_buffered_up <- buffering_analysis %>%
  filter(buffering_pattern == "ELT7_buffers_upregulation") %>%
  arrange(desc(interaction_log2FC)) %>%
  head(20)

top_buffered_down <- buffering_analysis %>%
  filter(buffering_pattern == "ELT7_buffers_downregulation") %>%
  arrange(interaction_log2FC) %>%
  head(20)

cat("Top 10 genes where ELT-7 buffers UP-regulation (greater up in double):\n")
if (nrow(top_buffered_up) > 0) {
  print(top_buffered_up %>%
    select(wbid, elt2d_log2FC, double_log2FC, interaction_log2FC, interaction_padj) %>%
    head(10), row.names = FALSE)
} else {
  cat("  None found\n")
}
cat("\n")

cat("Top 10 genes where ELT-7 buffers DOWN-regulation (greater down in double):\n")
if (nrow(top_buffered_down) > 0) {
  print(top_buffered_down %>%
    select(wbid, elt2d_log2FC, double_log2FC, interaction_log2FC, interaction_padj) %>%
    head(10), row.names = FALSE)
} else {
  cat("  None found\n")
}
cat("\n")

# Save comprehensive results
write.csv(buffering_analysis, "Dineen_DESeq2_buffering_analysis.csv", row.names = FALSE)
cat("Saved comprehensive analysis: Dineen_DESeq2_buffering_analysis.csv\n\n")

# Final summary
cat(paste0(strrep("=", 70), "\n"))
cat("SUMMARY\n")
cat(paste0(strrep("=", 70), "\n"))
cat(sprintf("Total genes with ELT-7 buffering effect (sig in double vs elt2d): %d\n",
            sum(buffering_analysis$interaction_sig, na.rm=TRUE)))
cat(sprintf("  - ELT-7 buffers upregulation: %d\n",
            sum(buffering_analysis$buffering_pattern == "ELT7_buffers_upregulation", na.rm=TRUE)))
cat(sprintf("  - ELT-7 buffers downregulation: %d\n",
            sum(buffering_analysis$buffering_pattern == "ELT7_buffers_downregulation", na.rm=TRUE)))
cat(sprintf("  - ELT-7 reveals new targets: %d\n",
            sum(buffering_analysis$buffering_pattern == "ELT7_reveals_new_targets", na.rm=TRUE)))
cat(sprintf("\nGenes with ELT-2 dominant effect (no additional change in double): %d\n",
            sum(buffering_analysis$buffering_pattern == "ELT2_dominant_no_ELT7_effect", na.rm=TRUE)))
cat("\n")
cat(paste0(strrep("=", 70), "\n"))
cat("ANALYSIS COMPLETE\n")
cat(paste0(strrep("=", 70), "\n"))
