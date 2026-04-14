#!/usr/bin/env Rscript

# DESeq2 rlog transformation of Dineen raw counts
# Prepares data for WGCNA and GENIE3 analyses

library(DESeq2)
library(tibble)

# Read raw counts
cat("Reading Dineen_raw_counts.csv...\n")
dineen_raw <- read.csv("Dineen_raw_counts.csv", row.names=NULL)

# Remove metadata columns (wbint, update_time) and set wbid as rownames
cat("Preparing count matrix...\n")
count_cols <- c("wt_sorted_1", "wt_sorted_2", "wt_sorted_3", "wt_sorted_4",
                "elt7d_sorted_1", "elt7d_sorted_2", "elt7d_sorted_3",
                "elt2d_sorted_1", "elt2d_sorted_2", "elt2d_sorted_3", "elt2d_sorted_4",
                "elt2delt7d_sorted_1", "elt2delt7d_sorted_2", "elt2delt7d_sorted_3")

# Extract count matrix and set gene IDs as rownames
count_matrix <- dineen_raw[, count_cols]
rownames(count_matrix) <- dineen_raw$wbid

cat(sprintf("Count matrix dimensions: %d genes x %d samples\n",
            nrow(count_matrix), ncol(count_matrix)))

# Create sample metadata for DESeq2
sample_info <- data.frame(
  sample = colnames(count_matrix),
  genotype = c(rep("wt", 4),
               rep("elt7d", 3),
               rep("elt2d", 4),
               rep("elt2d_elt7d", 3))
)
rownames(sample_info) <- sample_info$sample

# Create DESeq2 dataset
cat("Creating DESeq2 dataset...\n")
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info,
  design = ~ genotype
)

# Filter low-count genes (recommended for rlog)
# Keep genes with at least 10 reads total across all samples
cat("Filtering low-count genes...\n")
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat(sprintf("Retained %d genes after filtering\n", sum(keep)))

# Perform rlog transformation
cat("Performing rlog transformation...\n")
cat("(This may take a few minutes...)\n")
rld <- rlog(dds, blind=FALSE)

# Extract rlog-transformed counts
rlog_counts <- assay(rld)

# Save to CSV file
output_file <- "Dineen_rlog_counts.csv"
cat(sprintf("Saving to %s...\n", output_file))

# Convert to data frame with gene IDs
rlog_df <- as.data.frame(rlog_counts)
rlog_df <- rownames_to_column(rlog_df, var = "wbid")

write.csv(rlog_df, output_file, row.names = FALSE)

cat(sprintf("\nSuccess! Saved %d genes x %d samples to %s\n",
            nrow(rlog_counts), ncol(rlog_counts), output_file))
cat("\nThis file is ready for WGCNA and GENIE3 analyses.\n")
