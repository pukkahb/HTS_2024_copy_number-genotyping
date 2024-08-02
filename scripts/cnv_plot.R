#!/usr/bin/env Rscript

# Install & Load necessary libraries
install.packages(c("cn.mops", "ggplot2"))
library(cn.mops)
library(ggplot2)

# Set parameters for CNV and output directories
CNV_DIR <- "/home/project/results/cnv"  # Replace with your CNV directory path
OUTPUT_DIR <- "/home/project/results/plots"  # Replace with your output directory path

# Get CNV file names from the CNV directory
cnv_files <- list.files(CNV_DIR, pattern = "*.cnv", full.names = TRUE)

# Loop through each CNV file to read and plot
for (cnv_file in cnv_files) {
  # Create output file name based on input file
  output_file <- file.path(OUTPUT_DIR, paste0(basename(sub(".cnv", "", cnv_file)), ".png"))

  # Read CNV data
  cnv_data <- read.table(cnv_file, header = FALSE)
  colnames(cnv_data) <- c("chr", "start", "end", "CN", "quality")

  # Plot CNV data
  p <- ggplot(cnv_data, aes(x = start, y = CN)) +
    geom_point(alpha = 0.5) +
    theme_minimal() +
    labs(title = paste("CNV Plot for", basename(cnv_file)), x = "Position", y = "Copy Number")

  # Save plot
  ggsave(output_file, plot = p, width = 10, height = 5)
}

