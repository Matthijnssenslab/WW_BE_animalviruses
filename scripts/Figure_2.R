# Figure 2. Heatmap Overall

# Load required libraries
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)
library(readr)
library(tibble)

# Load data and ensure sample_ID is character
data <- read_csv("data.csv")
data$sample_ID <- as.character(data$sample_ID)

# (Optional) Keep all samples – note that M6 and M8 are retained
data_fam <- data

# Summarize total reads aligned by genus and sample_ID
heatmap_data <- data_fam %>%
  group_by(genus, sample_ID) %>%
  summarise(total_reads = sum(reads_aligned, na.rm = TRUE)) %>%
  ungroup()

# Pivot data to wide format and set genus as row names
heatmap_matrix <- heatmap_data %>%
  pivot_wider(names_from = sample_ID, values_from = total_reads, values_fill = 0) %>%
  column_to_rownames(var = "genus") %>%
  as.matrix()

# Log-transform the matrix (adding 1 to avoid log10(0))
heatmap_matrix_log <- log10(heatmap_matrix + 1)

# Build a binary annotation table for host associations per genus
genus_annotations <- data_fam %>%
  select(genus, grouped_association) %>%
  distinct() %>%
  group_by(genus) %>%
  summarise(grouped_associations = list(unique(grouped_association))) %>%
  ungroup()

# Extract all unique associations
all_associations <- unique(unlist(genus_annotations$grouped_associations))

# Create a binary annotation matrix (rows: genus; columns: associations)
annotation_matrix <- sapply(all_associations, function(x) {
  sapply(genus_annotations$grouped_associations, function(y) x %in% y)
})
annotation_df <- as.data.frame(annotation_matrix)
rownames(annotation_df) <- genus_annotations$genus

# Define colors for each association and prepare annotation color mapping
association_colors <- c(
  "Strictly animal" = "darkred",
  "Human and other animals" = "green4",
  "insects" = "gray20",
  "other" = "blue4"
)
annotation_colors_list <- lapply(association_colors, function(col) c("TRUE" = col, "FALSE" = "white"))

# Create the row annotation object for the heatmap
row_anno <- rowAnnotation(
  df = annotation_df,
  col = annotation_colors_list,
  annotation_name_side = "top",
  show_annotation_name = TRUE
)

# Define a custom color scale using colorRamp2 for the heatmap
col_fun <- colorRamp2(
  c(min(heatmap_matrix_log, na.rm = TRUE),
    mean(heatmap_matrix_log, na.rm = TRUE),
    max(heatmap_matrix_log, na.rm = TRUE)),
  c("white", "yellow4", "darkorange4")
)

# Plot the heatmap
Heatmap(
  heatmap_matrix_log,
  name = "log10(reads_aligned + 1)",
  col = col_fun,
  row_names_side = "left",
  column_names_side = "top",
  row_names_gp = gpar(fontsize = 12),
  column_names_gp = gpar(fontsize = 12),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  left_annotation = row_anno,
  border = TRUE,
  heatmap_legend_param = list(title = "Log10 Reads")
)
