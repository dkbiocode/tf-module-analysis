#!/usr/bin/env Rscript

# Visualizations of ELT-7 buffering effects with TF highlighting

library(ggplot2)
library(dplyr)
library(gridExtra)

cat(paste0(strrep("=", 70), "\n"))
cat("ELT-7 Buffering Visualizations with TF Highlighting\n")
cat(paste0(strrep("=", 70), "\n\n"))

# Load data
cat("Loading data...\n")
buffering <- read.csv("Dineen_DESeq2_buffering_analysis.csv")
tf_list <- read.csv("TF_list.csv")

# Get baseMean from one of the DESeq2 results
elt2d_full <- read.csv("Dineen_DESeq2_elt2d_vs_wt.csv")

cat(sprintf("  Buffering analysis: %d genes\n", nrow(buffering)))
cat(sprintf("  TF list: %d transcription factors\n", nrow(tf_list)))

# Merge baseMean and mark TFs in buffering data
buffering <- buffering %>%
  left_join(elt2d_full %>% select(wbid, baseMean), by = "wbid") %>%
  mutate(is_TF = wbid %in% tf_list$wbid)

cat(sprintf("  TFs in buffering analysis: %d\n\n", sum(buffering$is_TF)))

# Filter for plotting (remove NA values and set significance threshold)
plot_data <- buffering %>%
  filter(!is.na(elt2d_padj) & !is.na(double_padj) & !is.na(interaction_padj)) %>%
  mutate(
    # Simplified buffering categories for plotting
    buffering_category = case_when(
      buffering_pattern == "ELT7_buffers_upregulation" ~ "Buffers UP",
      buffering_pattern == "ELT7_buffers_downregulation" ~ "Buffers DOWN",
      buffering_pattern == "ELT7_compensatory" ~ "Compensatory",
      buffering_pattern == "ELT7_reveals_new_targets" ~ "New targets",
      buffering_pattern == "ELT2_dominant_no_ELT7_effect" ~ "No ELT-7 effect",
      TRUE ~ "Other"
    ),
    # Point alpha based on significance
    point_alpha = ifelse(interaction_sig, 0.6, 0.3),
    # Calculate average expression for MA plot
    baseMean_log = log10(baseMean + 1)
  )

# Define color schemes
tf_colors <- c("TF" = "#E31A1C", "non-TF" = "#1F78B4")
buffering_colors <- c(
  "Buffers UP" = "#E31A1C",
  "Buffers DOWN" = "#1F78B4",
  "Compensatory" = "#33A02C",
  "New targets" = "#FF7F00",
  "No ELT-7 effect" = "#999999",
  "Other" = "#CCCCCC"
)

cat("Creating visualizations...\n\n")

# ============================================================================
# PLOT 1: MA Plot (Interaction)
# baseMean vs interaction_log2FC
# ============================================================================
cat("1. MA plot (interaction: Double vs ELT-2Δ)...\n")

# Prepare data for MA plot - highlight significant interactions
ma_data <- plot_data %>%
  arrange(is_TF) # TFs plotted last so they're on top

p1 <- ggplot(ma_data, aes(x = baseMean_log, y = interaction_log2FC)) +
  # All points
  geom_point(aes(color = is_TF), size = 0.8, alpha = 0.3) +
  # Highlight significant interactions
  geom_point(data = filter(ma_data, interaction_sig),
             aes(color = is_TF), size = 1.2, alpha = 0.6) +
  # Reference lines
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray40", linewidth = 0.3) +
  # Styling
  scale_color_manual(values = tf_colors, name = "Gene Type",
                     labels = c("FALSE" = "Non-TF", "TRUE" = "TF")) +
  labs(
    title = "MA Plot: Interaction Effect (Double vs ELT-2Δ)",
    subtitle = sprintf("Significant interactions (padj<0.05): %d genes (%d TFs)",
                      sum(ma_data$interaction_sig),
                      sum(ma_data$interaction_sig & ma_data$is_TF)),
    x = "log10(Mean Expression)",
    y = "log2FC (Double vs ELT-2Δ)"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave("buffering_MA_plot.pdf", p1, width = 8, height = 7)
ggsave("buffering_MA_plot.png", p1, width = 8, height = 7, dpi = 300)

# ============================================================================
# PLOT 2: Scatter Plot - ELT-2Δ vs Double (both vs WT)
# Shows buffering as deviation from diagonal
# ============================================================================
cat("2. Scatter plot (ELT-2Δ vs Double, both vs WT)...\n")

# Prepare data - focus on genes significant in at least one condition
scatter_data <- plot_data %>%
  filter(elt2d_sig | double_sig) %>%
  arrange(is_TF)

p2 <- ggplot(scatter_data, aes(x = elt2d_log2FC, y = double_log2FC)) +
  # Quadrant shading
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
           fill = "lightblue", alpha = 0.1) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
           fill = "pink", alpha = 0.1) +
  # Points colored by buffering pattern
  geom_point(aes(color = buffering_category, shape = is_TF),
             size = 1.5, alpha = 0.5) +
  # Diagonal line (y = x, no interaction)
  geom_abline(slope = 1, intercept = 0, linetype = "solid",
              color = "black", linewidth = 0.8) +
  # Reference lines
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.3) +
  # Styling
  scale_color_manual(values = buffering_colors, name = "Buffering Pattern") +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17),
                     name = "Gene Type",
                     labels = c("FALSE" = "Non-TF", "TRUE" = "TF")) +
  labs(
    title = "ELT-7 Buffering: Double vs ELT-2Δ Effects",
    subtitle = "Above diagonal = stronger in double (buffering); Below = weaker in double",
    x = "log2FC ELT-2Δ vs WT",
    y = "log2FC Double vs WT"
  ) +
  coord_fixed() +
  theme_bw() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

ggsave("buffering_scatter_plot.pdf", p2, width = 10, height = 8)
ggsave("buffering_scatter_plot.png", p2, width = 10, height = 8, dpi = 300)

# ============================================================================
# PLOT 3: Effect Decomposition - ELT-2Δ effect vs Interaction
# Shows how buffering relates to initial effect
# ============================================================================
cat("3. Effect decomposition plot (ELT-2Δ vs Interaction)...\n")

decomp_data <- plot_data %>%
  filter(elt2d_sig | interaction_sig) %>%
  arrange(is_TF)

p3 <- ggplot(decomp_data, aes(x = elt2d_log2FC, y = interaction_log2FC)) +
  # Quadrant shading
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
           fill = "lightblue", alpha = 0.1) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0,
           fill = "pink", alpha = 0.1) +
  # Points
  geom_point(aes(color = buffering_category, shape = is_TF),
             size = 1.5, alpha = 0.5) +
  # Reference lines
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray40", linewidth = 0.3) +
  # Styling
  scale_color_manual(values = buffering_colors, name = "Buffering Pattern") +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17),
                     name = "Gene Type",
                     labels = c("FALSE" = "Non-TF", "TRUE" = "TF")) +
  labs(
    title = "Buffering Effect Decomposition",
    subtitle = "Same sign = enhancing buffering; Opposite sign = compensatory",
    x = "log2FC ELT-2Δ vs WT (Initial Effect)",
    y = "log2FC Interaction (Double - ELT-2Δ)"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

ggsave("buffering_decomposition_plot.pdf", p3, width = 10, height = 8)
ggsave("buffering_decomposition_plot.png", p3, width = 10, height = 8, dpi = 300)

# ============================================================================
# PLOT 4: Combined panel figure
# ============================================================================
cat("4. Combined panel figure...\n")

# Simplified versions for panel
p1_panel <- p1 + theme(legend.position = "none", plot.title = element_text(size = 10))
p2_panel <- p2 + theme(legend.position = "none", plot.title = element_text(size = 10))
p3_panel <- p3 + theme(legend.position = "none", plot.title = element_text(size = 10))

combined <- grid.arrange(p1_panel, p2_panel, p3_panel, ncol = 2)
ggsave("buffering_combined_panels.pdf", combined, width = 14, height = 12)
ggsave("buffering_combined_panels.png", combined, width = 14, height = 12, dpi = 300)

cat("\nPlots saved:\n")
cat("  - buffering_MA_plot.pdf/png\n")
cat("  - buffering_scatter_plot.pdf/png\n")
cat("  - buffering_decomposition_plot.pdf/png\n")
cat("  - buffering_combined_panels.pdf/png\n\n")

# ============================================================================
# TF ENRICHMENT ANALYSIS
# ============================================================================
cat(paste0(strrep("=", 70), "\n"))
cat("TF ENRICHMENT ANALYSIS\n")
cat(paste0(strrep("=", 70), "\n\n"))

# Overall TF frequency
total_genes <- nrow(buffering)
total_tfs <- sum(buffering$is_TF)
tf_freq_overall <- total_tfs / total_genes

cat(sprintf("Overall TF frequency: %d / %d = %.2f%%\n\n",
            total_tfs, total_genes, 100 * tf_freq_overall))

# TF enrichment by buffering pattern
enrichment <- buffering %>%
  group_by(buffering_pattern) %>%
  summarise(
    total_genes = n(),
    num_TFs = sum(is_TF),
    pct_TFs = 100 * num_TFs / total_genes,
    .groups = "drop"
  ) %>%
  arrange(desc(pct_TFs))

cat("TF enrichment by buffering pattern:\n")
print(enrichment, n = Inf)
cat("\n")

# Fisher's exact test for each pattern
cat("Fisher's exact test for TF enrichment:\n")
patterns_to_test <- c("ELT7_buffers_upregulation", "ELT7_buffers_downregulation",
                      "ELT7_compensatory", "ELT7_reveals_new_targets",
                      "ELT2_dominant_no_ELT7_effect")

for (pattern in patterns_to_test) {
  in_pattern <- buffering$buffering_pattern == pattern

  # 2x2 contingency table
  tf_in_pattern <- sum(buffering$is_TF & in_pattern)
  nontf_in_pattern <- sum(!buffering$is_TF & in_pattern)
  tf_not_pattern <- sum(buffering$is_TF & !in_pattern)
  nontf_not_pattern <- sum(!buffering$is_TF & !in_pattern)

  contingency <- matrix(c(tf_in_pattern, nontf_in_pattern,
                         tf_not_pattern, nontf_not_pattern),
                       nrow = 2, byrow = TRUE)

  test <- fisher.test(contingency)

  cat(sprintf("  %s:\n", pattern))
  cat(sprintf("    TFs: %d / %d (%.1f%%), OR = %.2f, p = %.2e\n",
              tf_in_pattern, sum(in_pattern),
              100 * tf_in_pattern / sum(in_pattern),
              test$estimate, test$p.value))
}

cat("\n")

# Plot TF enrichment
enrichment_plot <- ggplot(enrichment, aes(x = reorder(buffering_pattern, pct_TFs),
                                         y = pct_TFs)) +
  geom_col(fill = "#E31A1C", alpha = 0.7) +
  geom_hline(yintercept = 100 * tf_freq_overall, linetype = "dashed",
             color = "black", linewidth = 0.8) +
  coord_flip() +
  labs(
    title = "TF Enrichment by Buffering Pattern",
    subtitle = sprintf("Dashed line = overall TF frequency (%.1f%%)", 100 * tf_freq_overall),
    x = "Buffering Pattern",
    y = "% Transcription Factors"
  ) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave("TF_enrichment_by_pattern.pdf", enrichment_plot, width = 10, height = 6)
ggsave("TF_enrichment_by_pattern.png", enrichment_plot, width = 10, height = 6, dpi = 300)

cat("Saved TF enrichment plot: TF_enrichment_by_pattern.pdf/png\n\n")

# Export TFs in each category for inspection
tfs_buffering <- buffering %>%
  filter(is_TF) %>%
  select(wbid, elt2d_log2FC, elt2d_padj, double_log2FC, double_padj,
         interaction_log2FC, interaction_padj, buffering_pattern) %>%
  arrange(buffering_pattern, desc(abs(interaction_log2FC)))

write.csv(tfs_buffering, "TFs_buffering_analysis.csv", row.names = FALSE)
cat(sprintf("Exported %d TFs to TFs_buffering_analysis.csv\n\n", nrow(tfs_buffering)))

cat(paste0(strrep("=", 70), "\n"))
cat("ANALYSIS COMPLETE\n")
cat(paste0(strrep("=", 70), "\n"))
