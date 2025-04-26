#!/usr/local/lib/R-3.5.3/lib/R/bin/Rscript

# ------------------------------------------------------------------------------
# Loop Intensity Fold Change by Gene Expression Response to Nutlin
#
# This script evaluates how changes in gene expression following nutlin treatment
# relate to Hi-C loop intensity changes at gene-associated loops. Gene expression
# changes are binned, and loop fold changes are visualized across these bins.
#
# Steps:
# 1. Load gene expression and loop signal data for multiple cell lines.
# 2. Calculate fold change in expression and loop signal (nutlin vs control).
# 3. Bin genes by expression ratio and group associated loops accordingly.
# 4. Visualize loop fold changes across gene bins using boxplots.
# 5. Perform Pearson correlation between gene expression and loop IF.
#
# Inputs:
# - paths.tsv: file with cell line and reference paths
# - tssLoops_values.tsv: normalized Hi-C interaction frequencies
# - gene_expression_CNVcorrected.tsv: CNV-corrected gene expression matrix
#
# Outputs:
# - PNG boxplots for each cell line with annotated correlation results
#
# Author: Gony Shanel
# Date: March 2025
# ------------------------------------------------------------------------------

## Set the scripts directory as working directory.
## This script is designed with relative paths for the specific data structure
## with which it was uploaded to https://github.com/ElkonLab/micro-C.

# Arguments ----

# Define file paths for loop interaction and expression data
tssLoops_path <- "../../data/loops/tssLoops_values.tsv"
exp_path <- "../../data/gene_expression_CNVcorrected.tsv"

# Setup ----

library(ggplot2)
library(dplyr)
library(data.table)
library(ggsci)

# Load Data ----

loops_raw <- fread(tssLoops_path)              # Hi-C loop intensity data
exp_data <- fread(exp_path)                    # Gene expression data
cells <- readLines("../../data/cell_lines.txt")  # List of cell lines

# Loop Over Cell Lines ----

cor <- numeric()  # Initialize vector to store correlation p-values

cell <- cells[2]  # Example for dev/testing
for (cell in cells) {
  
  ## Define paired samples ----
  sample1 <- paste(cell, "control", sep = "-")
  sample2 <- paste(cell, "nutlin", sep = "-")
  
  # Filter genes expressed in the control condition
  exp_raw <- exp_data[exp_data[[sample1]] >= 1]
  
  ## Prepare Expression Data ----
  exp <- exp_raw[gene %in% loops_raw$gene, c("gene", ..sample1, ..sample2)]
  exp$R <- (exp[[3]] + 1) / (exp[[2]] + 1)        # Add pseudocount, compute ratio
  exp$logR <- log2(exp$R)                         # Log2 fold change
  
  # Bin genes into 10 groups based on fold change
  exp <- exp %>% arrange(R)
  group_size <- floor(nrow(exp) / 10)
  exp$bin <- cut(exp$R,
                 breaks = quantile(exp$R, probs = seq(0, 1, 0.1)),
                 include.lowest = TRUE, labels = 1:10)
  exp$bin <- factor(exp$bin, levels = 1:10)
  
  table(exp$bin)  # Optional: print bin sizes
  
  ## Prepare Loop Data ----
  cols <- c(colnames(loops_raw)[1:7], "gene", sample1, sample2)
  loops <- loops_raw[gene %in% exp$gene, ..cols]
  loops$R <- loops[[sample2]] / loops[[sample1]]     # Loop fold change
  loops$logR <- log2(loops$R)
  
  ## Merge Expression & Loop Data ----
  exp_loops <- merge(
    exp[, c("gene", "R", "logR", "bin")],
    loops[, c(1:8, 11, 12)],
    by = "gene",
    suffixes = c(".cpm", ".if")
  )
  
  exp_loops <- exp_loops %>% arrange(R.cpm)
  
  ## Correlation Test ----
  cor.IF <- cor.test(exp_loops$R.cpm, exp_loops$R.if, method = "pearson")
  r <- format(cor.IF$estimate, digits = 2)
  pval <- as.numeric(format(cor.IF$p.value, digits = 2))
  stats.IF <- paste0("r=", r, "\np=", pval)
  
  cor[[cell]] <- pval  # Store p-value
  
  ## Prepare for Plotting ----
  exp_loops$round <- round(exp_loops$logR.cpm * 2) / 2
  data <- exp_loops
  data$round <- factor(data$round)
  
  ## Bin Summary Counts ----
  bin_counts <- exp_loops %>%
    group_by(bin) %>%
    summarise(
      num_genes = length(gene[!duplicated(gene)]),
      num_loops = sum(!is.na(logR.if))
    )
  
  x.ticks <- paste0("(N=", bin_counts$num_loops, ")")
  
  ## Plot Boxplot ----
  ggplot(data, aes(x = bin, y = logR.if, fill = bin, group = bin)) +
    geom_boxplot(outlier.shape = NA) +
    coord_cartesian(ylim = c(-2, 2)) +
    theme_bw() +
    scale_fill_brewer(palette = "RdBu") +
    labs(
      x = 'loops associated with genes\nbinned by expression fold change',
      y = 'loop intensity log2 fold change'
    ) +
    annotate("text", x = 0, y = 1.5, label = stats.IF,
             hjust = -0.2, vjust = 0.15, size = 2) +
    scale_x_discrete(labels = x.ticks) +
    theme(
      text = element_text(size = 8),
      axis.text.x = element_text(angle = 45, vjust = 0.8),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5)
    )
  
  ## Save Plot ----
  ggsave(paste0(cell, "_pv-", pval, ".png"), width = 3, height = 2.5)
}
