# R script to analyze gasch2000 expression data and generate a heatmap for top 10 variable genes

library(tidyverse)

# read data
file <- "data/gasch2000.txt"
expr <- read.delim(file, header = TRUE, sep = "\t", comment.char = "#", stringsAsFactors = FALSE)

# determine numeric expression columns
is_num <- sapply(expr, is.numeric)
num_cols <- which(is_num)
# keep only first 30 stress condition columns
if(length(num_cols) > 30) num_cols <- num_cols[1:30]
cat("Numeric columns indices used (first 30):", paste(num_cols, collapse = ","), "\n")

# compute CV per gene
id_col <- names(expr)[1]
data_mat <- as.matrix(expr[num_cols])
row_means <- rowMeans(data_mat, na.rm = TRUE)
row_sds <- apply(data_mat, 1, sd, na.rm = TRUE)
cv <- ifelse(row_means == 0, NA, row_sds / row_means)
expr$CV <- cv

# pick top 10
library(dplyr)
top10 <- expr %>% arrange(desc(CV)) %>% slice(1:10)

# reshape for plotting
long_df <- top10 %>%
  select(all_of(id_col), all_of(names(expr)[num_cols])) %>%
  pivot_longer(-all_of(id_col), names_to = "condition", values_to = "expression")

# plot
library(ggplot2)
heat <- ggplot(long_df, aes(x = .data[[id_col]], y = condition, fill = expression)) +
  geom_tile(color = "black") +
  scale_fill_distiller(palette = "PuGn") +
  theme_minimal(base_family = "Times New Roman", base_size = 11) +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  labs(title = "gene top 10", x = "gene", y = "conditions")
heat <- heat + theme(
  legend.position = c(-0.15, -0.15),
  legend.justification = c(0, 0),
  plot.margin = margin(40, 40, 40, 60, "pt")
)

# save image
outfile <- "results/version2_heat_map.png"
ggsave(outfile, plot = heat, width = 5, height = 5, dpi = 300, bg = "white")

cat("Saved heatmap to", outfile, "\n")
