### Load required libraries
library(dplyr)
library(ggplot2)

### Data Loading & Filtering
data <- read.csv("data.csv") %>%
# Calculate genome completeness
  mutate(completeness = covered_bases / reference_length)

### Section 1: Completeness Bins & Bar Plot
# Create completeness bins and count sequences per bin
bins <- cut(data$completeness, breaks = c(0, 0.25, 0.5, 0.75, 1),
            labels = c("0-25%", "25-50%", "50-75%", "75-100%"))
bin_counts <- table(bins)
total_sequences <- nrow(data)
bin_percentages <- 100 * bin_counts / total_sequences

# Build a data frame for plotting
plot_data <- data.frame(
  bin = names(bin_counts),
  count = as.numeric(bin_counts),
  percentage = as.numeric(bin_percentages)
)

# Bar plot of completeness distribution
p1 <- ggplot(plot_data, aes(x = "", y = percentage, fill = bin)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(percentage, 1), "% (", count, ")")),
            color = "white", position = position_stack(vjust = 0.5)) +
  labs(x = "", y = "Percentage of Sequences (%)") +
  scale_fill_manual(values = c("0-25%" = "#1a9e77", "25-50%" = "#d76227",
                               "50-75%" = "#7672b2", "75-100%" = "#e52c8b")) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA),
        panel.grid.minor = element_blank())

print(p1)

### Section 2: Species Count per Sample by Completeness Category
# Count distinct species for sequences with 50-75% completeness
count_50_75 <- data %>%
  filter(completeness >= 0.50 & completeness <= 0.75) %>%
  group_by(sample_ID) %>%
  summarise(num_species = n_distinct(species))

# Count distinct species for sequences with >75% completeness
count_above_75 <- data %>%
  filter(completeness > 0.75) %>%
  group_by(sample_ID) %>%
  summarise(num_species = n_distinct(species))

# Merge the two results and rename columns
completeness_counts <- merge(count_50_75, count_above_75, by = "sample_ID", all = TRUE)
colnames(completeness_counts) <- c("Sample_ID", "50-75% Completeness", ">75% Completeness")
completeness_counts[is.na(completeness_counts)] <- 0

# Display the summarized species counts
print(completeness_counts)
