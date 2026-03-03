#!/usr/bin/env Rscript

# [Environment]: Write R code that runs in the [bch709_vibe_coding] conda environment.
# view the file "copilot-instructions.md" for more information about installed packages and the environment. Please save the code in the "scripts" folder
#
# **Input specification:**
# - File: results/chr_feature_counts.tsv, tab delimiter 
# - Structure: [column descriptions, special parsing rules]
# - Additional inputs: data/gasch2000.txt
#
# **Analysis conditions:**
# Analysis procedure:
# 1. From gasch2000.txt, extract only the top 200 genes (by gene_id list from results/chr_feature_counts.tsv,)
# 2. Skip metadata columns (NAME, description, GWEIGHT), keep only numeric condition columns
# 3. Subset to the first 30 condition columns only
# 4. Z-score normalize: for each row (gene), compute (value - mean) / sd
# 5. Hierarchical clustering: dist(euclidean) → hclust(ward.D2)
# 6. cutree(k=4) to assign 4 clusters
#
# For the heatmap:
# - Use pheatmap
# - Show cluster assignment as annotation_row color bar
# - Use only the first 30 condition columns (original order, cluster_cols = FALSE)
#
# **Output 1 — Table:**
# - Filename: please save in the results folder as "cluster_assignment.tsv"
# - Columns: "gene_id "(contains the gene identifier example: YAL001C); "cluster" (contains the cluster number 1-4) 
# - Sorting: Sort by cluster ascending
#
# **Output 2 — Plot:**
# - Filename: please save in the results folder as "csv_top200_cluster_heatmap.pdf"
# - Specification (value): Data (log2 expression ratios from results/chr_feature_counts.tsv,); Z-score ([value-row_mean]/row_sd); Rows (gene_id [heiarchal clustering, ward.D2, euclidean); Columns (first 30 stress condition columns only [original order, clister_cols=FALSE]); Annotation (k=4 cutree as color bar)
# - Size: 8x12 inches
# - Colors: Please use blue for cluster 1, purple for cluster 2, green for cluster 3, and pink for cluster 4.
# - Axes/labels/legend: Title:" Yeast Cluster Heatmap"

# Load required packages
library(tidyverse)
library(pheatmap)
library(stats)

# Set working directory
setwd("/home/gfornesi/bch709_vibe_coding")

# Read the gasch2000.txt file
cat("Loading expression data from gasch2000.txt...\n")
expr_data <- read.delim("data/gasch2000.txt", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

# Extract the top 200 genes (rows 1-200 after excluding header)
cat("Extracting top 200 genes...\n")
top_200_genes <- expr_data[1:200, ]

# Get gene IDs
gene_ids <- rownames(top_200_genes)

# Extract only numeric columns (skip NAME, GWEIGHT, and other metadata)
# Columns 2 onwards are numeric (column 1 is NAME, column 2 is GWEIGHT, columns 3+ are conditions)
numeric_cols <- 3:ncol(top_200_genes)  # Skip UID, NAME, GWEIGHT

# Subset to first 30 condition columns
expr_subset <- top_200_genes[, numeric_cols[1:30]]

cat("Expression data dimensions: ", nrow(expr_subset), "genes x", ncol(expr_subset), "conditions\n")

# Z-score normalize each row (gene)
cat("Performing Z-score normalization...\n")
expr_normalized <- t(apply(expr_subset, 1, function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}))

# Handle NaN values (from genes with 0 variance)
expr_normalized[is.nan(expr_normalized)] <- 0

# Hierarchical clustering using euclidean distance and ward.D2 linkage
cat("Performing hierarchical clustering...\n")
dist_matrix <- dist(expr_normalized, method = "euclidean")
hc <- hclust(dist_matrix, method = "ward.D2")

# Cut tree to assign 4 clusters
cluster_assignment <- cutree(hc, k = 4)

# Create output table 1: cluster assignment
cat("Creating cluster assignment table...\n")
cluster_df <- data.frame(
  gene_id = names(cluster_assignment),
  cluster = as.numeric(cluster_assignment),
  stringsAsFactors = FALSE
)

# Sort by cluster ascending
cluster_df <- cluster_df %>% arrange(cluster)

# Save cluster assignment table
write.table(cluster_df, "results/cluster_assignment.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved cluster assignment to: results/cluster_assignment.tsv\n")

# Print summary statistics
cat("\n--- Cluster Summary ---\n")
print(table(cluster_assignment))

# Create annotation data frame for heatmap (cluster assignments)
annotation_row <- data.frame(
  Cluster = as.factor(cluster_assignment),
  row.names = rownames(expr_normalized)
)

# Define cluster colors
cluster_colors <- list(
  Cluster = c("1" = "blue", "2" = "purple", "3" = "green", "4" = "pink")
)

# Create the heatmap
cat("Creating heatmap...\n")
pdf("results/csv_top200_cluster_heatmap.pdf", width = 8, height = 12)

pheatmap(
  expr_normalized,
  annotation_row = annotation_row,
  annotation_colors = cluster_colors,
  clustering_distance_rows = "euclidean",
  clustering_method = "ward.D2",
  cluster_cols = FALSE,
  main = "Yeast Cluster Heatmap",
  fontsize = 8,
  fontsize_row = 6,
  fontsize_col = 8,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-3, 3, length.out = 101)
)

dev.off()
cat("Saved heatmap to: results/csv_top200_cluster_heatmap.pdf\n")

cat("\n--- Analysis Complete ---\n")
