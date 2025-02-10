### Load required libraries
library(tidyverse)  # loads tidyr, readr, dplyr, ggplot2

### Read and pre‐process data
data <- read_csv("data.csv") %>%
  filter(!sample_ID %in% c("M6", "M8"))

### Prepare dataset for non-human animal viruses (Figure 3)
virus_data <- data %>%
  filter(grouped_association %in% c("Strictly animal", "insects")) %>%
  mutate(
    completeness      = covered_bases / reference_length,
    complete_status   = ifelse(completeness > 0.75, "Higher completeness (>75%)", "Lower completeness (<75%)"),
    # For insects use 'family', for others use 'species'
    identification_level = ifelse(grouped_association == "insects", family, species)
  )


### Figure 3 – Plot : Rotated axes (swap x and y)
plot <- ggplot(virus_data, aes(x = completeness, y = identification_level, fill = complete_status)) +
  geom_point(size = 4, alpha = 0.8, shape = 21, stroke = 0.3, color = "black") +
  facet_grid(`associated_with_(not_categorized)` ~ ., scales = "free_y", space = "free") +
  scale_fill_manual(values = c("Higher completeness (>75%)" = "#CA3D66",
                               "Lower completeness (<75%)" = "gray30")) +
  labs(
    title = "Completeness by Genus (Insects) and Species (Others) Grouped by Association",
    x = "Genome Completeness", y = "Family or Species", fill = "Genome completeness"
  ) +
  theme(
    axis.text.y     = element_text(hjust = 1, size = 12, color = "black"),
    axis.text.x     = element_text(size = 12, color = "black"),
    strip.background = element_rect(fill = "gray90", color = "black"),
    strip.text      = element_text(face = "bold", size = 11, color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray50", size = 0.1),
    panel.border    = element_rect(color = "black", fill = NA),
    plot.title      = element_text(face = "bold", size = 18, hjust = 0.5),
    legend.position = "bottom",
    legend.title    = element_text(face = "bold", size = 12),
    legend.text     = element_text(size = 10),
    plot.background = element_rect(fill = "white", color = "white"),
    plot.margin     = margin(10, 10, 10, 10)
  )

### Supplementary Figure – Insect-specific viruses at species level
# Filter for insects only and use species as the identification level
insect_data <- data %>%
  filter(grouped_association == "insects") %>%
  mutate(
    completeness      = covered_bases / reference_length,
    complete_status   = ifelse(completeness > 0.75, "Higher completeness (>75%)", "Lower completeness (<75%)"),
    identification_level = species
  ) %>%
  select(sample_ID, species, completeness, complete_status, grouped_association, identification_level)

supplementary_plot <- ggplot(insect_data, aes(x = identification_level, y = completeness, fill = complete_status)) +
  geom_point(size = 4, alpha = 0.8, shape = 21, stroke = 0.3, color = "black") +
  facet_grid(~ grouped_association, scales = "free_x", space = "free") +
  scale_fill_manual(values = c("Higher completeness (>75%)" = "darkred", "Lower completeness (<75%)" = "gray30")) +
  labs(
    title = "Completeness by Species (Insects)",
    x = "Species", y = "Completeness", fill = "Completeness Status"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y     = element_text(size = 12, color = "black"),
    strip.background = element_rect(fill = "gray90", color = "black"),
    strip.text      = element_text(face = "bold", size = 11, color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(color = "gray50", size = 0.1),
    panel.border    = element_rect(color = "black", fill = NA),
    plot.title      = element_text(face = "bold", size = 18, hjust = 0.5),
    legend.position = "bottom",
    legend.title    = element_text(face = "bold", size = 12),
    legend.text     = element_text(size = 10),
    plot.background = element_rect(fill = "white", color = "white"),
    plot.margin     = margin(10, 10, 10, 10)
  )

### Display or save the plots as needed
print(plot)           # Standard orientation
print(supplementary_plot)  # Supplementary figure (insect-specific)
