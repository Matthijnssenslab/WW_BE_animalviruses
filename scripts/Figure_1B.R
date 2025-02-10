# Figure 1.B

# Load required libraries
library(ggplot2)
library(dplyr)
library(readr)

# Set working directory and load data
setwd("~")
data <- read_csv("data.csv") %>%
  mutate(
    completeness = (covered_bases / reference_length) * 100,
    completeness_group = case_when(
      completeness < 25 ~ "0-25%",
      completeness < 75 ~ "25-75%",
      TRUE              ~ ">75%"
    )
  )

# Filter out rows with family "Picobirnaviridae" and ensure family info is not missing
data_75 <- data %>%
  filter(family != "Picobirnaviridae") %>%
  mutate(family = if_else(is.na(family), "Unknown", family))

# Count distinct species per family for each sample_ID
species_per_family <- data_75 %>%
  group_by(sample_ID, family) %>%
  summarise(distinct_species_count = n_distinct(species), .groups = "drop")

# Define a custom color palette for the families
num_families <- n_distinct(species_per_family$family)
custom_family_colors <- if (num_families > 12) {
  colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(num_families)
} else {
  RColorBrewer::brewer.pal(num_families, "Paired")
}

# Create a stacked bar plot
ggplot(species_per_family, aes(x = sample_ID, y = distinct_species_count, fill = family)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = custom_family_colors) +
  labs(
    title = "Distinct Species per Family for Samples with >75% Completeness",
    x = "Sample ID",
    y = "Distinct Species Count",
    fill = "Family"
  ) +
  theme_minimal() +
  theme(
    plot.title    = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title.x  = element_text(size = 14, face = "bold"),
    axis.title.y  = element_text(size = 14, face = "bold"),
    axis.text.x   = element_text(size = 12, angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border  = element_rect(color = "grey80", fill = NA, linewidth = 1),
    plot.background = element_rect(fill = "white")
  )
