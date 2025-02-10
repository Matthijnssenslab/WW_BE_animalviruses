library(rtracklayer)
library(ggplot2)
library(GenomicRanges)
library(dplyr)
library(Biostrings)




plotCoverage <- function(bed_file_path,
                         seqname_filter,
                         bin_size,
                         cap_value,
                         xlim_values,  # These are now based on the FASTA file length
                         width_filter = 100,
                         output_path = NULL) {

  # Import BED file
  bed_data <- rtracklayer::import(bed_file_path)
  bed_df <- as.data.frame(bed_data)

  # Filter for the sequence of interest and width greater than width_filter
  filtered_bed <- bed_df %>%
    filter(seqnames == seqname_filter, width > width_filter)

  # Convert the filtered BED data to GRanges object
  gr <- GRanges(
    seqnames = Rle(filtered_bed$seqnames),
    ranges = IRanges(start = filtered_bed$start, end = filtered_bed$end),
    strand = Rle(filtered_bed$strand)
  )

  # Calculate coverage
  cov <- coverage(gr)

  # Convert coverage Rle object to data frame
  cov_df <- as.data.frame(cov)
  cov_df$position <- as.numeric(rownames(cov_df))
  rownames(cov_df) <- NULL

  # Bin the positions into bin_size nucleotide intervals
  cov_df$bin <- with(cov_df, cut(position, breaks = seq(from = min(xlim_values), to = max(xlim_values), by = bin_size), include.lowest = TRUE, labels = FALSE))

  # Aggregate depth by bin
  agg_cov <- aggregate(value ~ bin, data = cov_df, FUN = mean)

  # Create midpoints for each bin to use as the x-axis in the plot
  agg_cov$midpoint <- (agg_cov$bin - 0.5) * bin_size

  # Cap the highest value to cap_value
  agg_cov$value <- pmin(agg_cov$value, cap_value)

  # Create the plot with a fixed x-axis (based on FASTA length)
  p <- ggplot(agg_cov, aes(x = midpoint, y = value)) +
    geom_area(stat = "identity", fill = "#192841", alpha = 0.90) +
    labs(x = paste("Genomic position (binned in", bin_size, "bp intervals)"), y = "Depth of coverage") +
    scale_x_continuous(limits = xlim_values) +  # Use xlim_values based on FASTA length
    coord_cartesian(xlim = xlim_values) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 65, hjust = 1),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          panel.grid.major = element_line(color = "gray", size = 0.3),
          panel.grid.minor = element_blank())

  # Save plot to file if output_path is provided
  if (!is.null(output_path)) {
    ggsave(output_path, plot = p, width = 10, height = 6)
  }

  return(p)
}


#plot & parameters
plot_path <- "bam_mapped.bed"
# Load the BED file (not necessary)
bed_data <- rtracklayer::import(plot_path)
seqname_filter <- "NC_001427.1"
bin_size <- 100
cap_value <- 1000
xlim_values <- c(0, 2300)

# Plot and save
coverage_plot <- plotCoverage(plot_path, 
                                         seqname_filter = seqname_filter, 
                                         bin_size = bin_size, 
                                         cap_value = cap_value, 
                                         xlim_values = xlim_values, 
                                         output_path = "coverage_plot.png")
print(coverage_plot)


#####script to plot all of them.#####
library(rtracklayer)
library(ape)  # For reading FASTA files

# Path to your BAM and BED files directory
bam_dir <- "bam_fasta/"

# Import the BED file to get sequence names
plot_path <- paste0(bam_dir, "bam_mapped.bed")
bed_data <- rtracklayer::import(plot_path)
bed_df <- as.data.frame(bed_data)

# Get all unique sequence names
unique_seqnames <- unique(bed_df$seqnames)

# Define a loop to plot for each sequence
for (seqname in unique_seqnames) {

  # Define parameters for each sequence
  seqname_filter <- as.character(seqname)

  # Find the corresponding FASTA file based on the sequence name
  fasta_file <- paste0(bam_dir, seqname_filter, ".fasta")

  # Read the FASTA file and calculate the sequence length
  fasta_data <- read.FASTA(fasta_file)
  fasta_length <- length(fasta_data[[1]])

  # Use the sequence length to set xlim values
  xlim_values <- c(0, fasta_length)

  bin_size <- 100
  cap_value <- 1000

  # Set output file name dynamically based on the sequence name
  output_file <- paste0("coverage_plot_", seqname_filter, ".png")

  # Call your existing plotCoverage function
  coverage_plot <- plotCoverage(plot_path,
                                seqname_filter = seqname_filter,
                                bin_size = bin_size,
                                cap_value = cap_value,
                                xlim_values = xlim_values,
                                output_path = output_file)

  # Print a message indicating the plot has been generated
  print(paste("Plot saved for sequence:", seqname_filter))

  # Display the plot in RStudio viewer
  print(coverage_plot)
}









library(ggplot2)



###### Function to plot GFF file ######
library(rtracklayer)
library(gggenes)
library(ggplot2)
library(cowplot)  # For combining plots
library(Biostrings)

# Function to plot GFF with CDS features as arrows, colored and arranged in tracks without numbers or gridlines
plotGFF_CDS <- function(gff_file, xlim_values = NULL, output_path = NULL, arrow_height = 0.1) {
  
  # Load the GFF file
  gff_data <- rtracklayer::import(gff_file)
  gff_df <- as.data.frame(gff_data)
  
  # Filter the GFF file to only include CDS features
  cds_df <- gff_df[gff_df$type == "CDS", ]
  
  # Create a new column "track" to separate overlapping features
  cds_df <- cds_df %>%
    arrange(start) %>%
    mutate(track = cumsum(c(0, diff(start) > 0) * (start < lag(end, default = 0))))
  
  # Create the plot using gggenes, with dark pastel blue color and adjustable arrow height
  p <- ggplot(cds_df, aes(xmin = start, xmax = end, y = as.factor(track), forward = (strand == "+"))) +
    geom_gene_arrow(fill = "#779ECB", color = "#779ECB", size = 0.5, height = arrow_height) +  # Adjust the arrow height here
    labs(x = "Genomic Position", y = "CDS") +  # Add x-axis label
    theme_genes() +
    theme(
      legend.position = "none",
      axis.text.y = element_blank(),    # Remove y-axis labels (numbers for tracks)
      axis.ticks.y = element_blank(),   # Remove y-axis ticks
      axis.line.y = element_blank(),    # Remove y-axis line
      panel.grid.major = element_blank(), # Ensure major gridlines are removed
      panel.grid.minor = element_blank(), # Ensure minor gridlines are removed
      axis.title.y = element_blank(),   # Remove y-axis title
      panel.background = element_rect(fill = "white"),  # White background
      plot.background = element_rect(fill = "white"))   # White background for the entire plot
  
  # Set x-axis limits if provided
  if (!is.null(xlim_values)) {
    p <- p + scale_x_continuous(limits = xlim_values)
  }
  
  # Save the plot if an output path is provided
  if (!is.null(output_path)) {
    ggsave(output_path, plot = p, width = 10, height = 6)
  }
  
  return(p)
}

# Function to plot the coverage plot without x-axis labels
plotCoverage <- function(bed_file_path, 
                         seqname_filter, 
                         bin_size, 
                         cap_value, 
                         xlim_values,  # These are now based on the FASTA file length
                         width_filter = 100, 
                         output_path = NULL) {
  
  # Import BED file
  bed_data <- rtracklayer::import(bed_file_path)
  bed_df <- as.data.frame(bed_data)
  
  # Filter for the sequence of interest and width greater than width_filter
  filtered_bed <- bed_df %>%
    filter(seqnames == seqname_filter, width > width_filter)
  
  # Convert the filtered BED data to GRanges object
  gr <- GRanges(
    seqnames = Rle(filtered_bed$seqnames),
    ranges = IRanges(start = filtered_bed$start, end = filtered_bed$end),
    strand = Rle(filtered_bed$strand)
  )
  
  # Calculate coverage
  cov <- coverage(gr)
  
  # Convert coverage Rle object to data frame
  cov_df <- as.data.frame(cov)
  cov_df$position <- as.numeric(rownames(cov_df))
  rownames(cov_df) <- NULL
  
  # Bin the positions into bin_size nucleotide intervals
  cov_df$bin <- with(cov_df, cut(position, breaks = seq(from = min(xlim_values), to = max(xlim_values), by = bin_size), include.lowest = TRUE, labels = FALSE))
  
  # Aggregate depth by bin
  agg_cov <- aggregate(value ~ bin, data = cov_df, FUN = mean)
  
  # Create midpoints for each bin to use as the x-axis in the plot
  agg_cov$midpoint <- (agg_cov$bin - 0.5) * bin_size
  
  # Cap the highest value to cap_value
  agg_cov$value <- pmin(agg_cov$value, cap_value)
  
  # Create the plot without x-axis labels and no grid lines
  p <- ggplot(agg_cov, aes(x = midpoint, y = value)) +
    geom_area(stat = "identity", fill = "#192841", alpha = 0.90) +
    labs(x = NULL, y = "Depth of coverage") +  # Remove x-axis labels
    scale_x_continuous(limits = xlim_values) +  # Use xlim_values based on FASTA length
    coord_cartesian(xlim = xlim_values) + 
    theme_minimal() +
    theme(axis.text.x = element_blank(),  # Remove x-axis labels and ticks
          axis.ticks.x = element_blank(),  # Remove x-axis ticks
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          panel.grid.major = element_line(color = "gray", size = 0.3),
          panel.grid.minor = element_blank())
  
  # Save plot to file if output_path is provided
  if (!is.null(output_path)) {
    ggsave(output_path, plot = p, width = 10, height = 6)
  }
  
  return(p)
}

# Path to your BAM and BED files directory
bam_dir <- "bam_fasta"

# Import the BED file to get sequence names
plot_path <- paste0(bam_dir, "bam_mapped.bed")
bed_data <- rtracklayer::import(plot_path)
bed_df <- as.data.frame(bed_data)

# Get all unique sequence names
unique_seqnames <- unique(bed_df$seqnames)

# Define a loop to plot for each sequence (both coverage and GFF)
for (seqname in unique_seqnames) {
  
  # Define parameters for each sequence
  seqname_filter <- as.character(seqname)
  
  # Find the corresponding FASTA file based on the sequence name
  fasta_file <- paste0(bam_dir, seqname_filter, ".fasta")
  
  # Read the FASTA file and calculate the sequence length
  fasta_data <- read.FASTA(fasta_file)
  fasta_length <- length(fasta_data[[1]])
  
  # Use the sequence length to set xlim values
  xlim_values <- c(0, fasta_length)
  
  bin_size <- 100
  cap_value <- 1000
  
  # Set output file name dynamically based on the sequence name
  coverage_output_file <- paste0("coverage_plot_", seqname_filter, ".png")
  gff_output_file <- paste0("gff_plot_", seqname_filter, ".png")
  
  # Call your existing plotCoverage function to generate the coverage plot
  coverage_plot <- plotCoverage(plot_path, 
                                seqname_filter = seqname_filter, 
                                bin_size = bin_size, 
                                cap_value = cap_value, 
                                xlim_values = xlim_values, 
                                output_path = coverage_output_file)
  
  # Find the corresponding GFF file based on the sequence name
  gff_file <- paste0(bam_dir, seqname_filter, ".gff3")
  
  # Call the GFF plotting function with a custom arrow height
  gff_plot <- plotGFF_CDS(gff_file, xlim_values = xlim_values, output_path = gff_output_file, arrow_height = 0.1)  # Adjust arrow_height
  
  # Combine the coverage plot and GFF plot
  combined_plot <- plot_grid(coverage_plot, gff_plot, align = "v", ncol = 1, rel_heights = c(2, 1))
  
  # Add a small reference label at the bottom (replacing the title)
  combined_plot <- combined_plot + 
    labs(caption = paste("Reference sequence:", seqname_filter)) + 
    theme(plot.caption = element_text(size = 4, hjust = 0.5))  # Make the reference small and centered
  
  # Save the combined plot as PDF (custom size)
  combined_output_file_pdf <- paste0("combined_coverage_gff_plot_", seqname_filter, ".pdf")
  ggsave(combined_output_file_pdf, plot = combined_plot, width = 9, height = 3.5, device = "pdf")
  
  # Print a message indicating the plot has been generated
  print(paste("Combined plot saved for sequence:", seqname_filter, "as PDF"))
  
  # Display the combined plot in RStudio viewer
  print(combined_plot)
}
