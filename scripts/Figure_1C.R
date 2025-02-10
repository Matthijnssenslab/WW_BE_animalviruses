# Load required libraries
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(readr)

# Load and filter data (excluding samples M6 and M8)
data <- read_csv("data.csv")
data_fam <- data %>% filter(!sample_ID %in% c("M6", "M8"))

# Summarize total reads per host category and compute percentages/labels
reads_summary <- data_fam %>%
  group_by(grouped_association) %>%
  summarise(total_reads = sum(reads_aligned, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    percentage_reads = total_reads / sum(total_reads) * 100,
    label_reads = paste0(total_reads, " (", round(percentage_reads, 1), "%)")
  )

# Summarize unique genera count per host category and compute percentages/labels
genera_summary <- data_fam %>%
  group_by(grouped_association) %>%
  summarise(genera_count = n_distinct(species)) %>%
  ungroup() %>%
  mutate(
    percentage_genera = genera_count / sum(genera_count) * 100,
    label_genera = paste0(genera_count, " (", round(percentage_genera, 1), "%)")
  )

# Define a custom, consistent color palette for host categories
color_palette <- c(
  "Human and other animals" = "#EA1A0A",
  "insects" = "#FAAA00",
  "other" = "#AAB850",
  "Strictly animal" = "#377EB8"
)

# Create the donut chart for total reads
donut_chart_reads <- ggplot(reads_summary, aes(x = 2, y = total_reads, fill = grouped_association)) +
  geom_col(width = 0.75, color = "white") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +  # Controls the donut hole size
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 10),
    plot.title   = element_text(size = 16, face = "bold", hjust = 0.5, color = "#4b4b4b")
  ) +
  scale_fill_manual(values = color_palette) +
  labs(
    title = "Reads Assigned to Different Host Infecting Genera",
    fill = "Host Category"
  ) +
  geom_text(aes(label = label_reads), position = position_stack(vjust = 0.5), color = "black", size = 4)

# Create the donut chart for genera counts
donut_chart_genera <- ggplot(genera_summary, aes(x = 2, y = genera_count, fill = grouped_association)) +
  geom_col(width = 0.75, color = "white") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +  # Controls the donut hole size
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 10),
    plot.title   = element_text(size = 16, face = "bold", hjust = 0.5, color = "#4b4b4b")
  ) +
  scale_fill_manual(values = color_palette) +
  labs(
    title = "Number of Genera Infecting Different Hosts",
    fill = "Host Category"
  ) +
  geom_text(aes(label = label_genera), position = position_stack(vjust = 0.5), color = "black", size = 4)

# Display the donut charts
print(donut_chart_reads)
print(donut_chart_genera)
