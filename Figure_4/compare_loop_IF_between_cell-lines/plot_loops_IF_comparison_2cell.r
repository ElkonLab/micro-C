#!/usr/local/lib/R-3.5.3/lib/R/bin/Rscript

# ------------------------------------------------------------------------------
# Plot Pairwise Comparison of Loop Interaction Frequencies (IF) between Samples
# ------------------------------------------------------------------------------
# This script compares Hi-C loop interaction frequencies (IF) between two samples,
# such as different cell lines or treatment conditions, using Mustache loop calls.
#
# It performs the following steps:
# 1. Loads Mustache loop sets for each sample and a reference loop index.
# 2. Identifies overlapping and unique loops between the two samples.
# 3. Retrieves normalized interaction frequencies (IFs) for each loop.
# 4. Creates a scatterplot of log2(IF) values between the two samples,
#    color-coded by loop category (common, unique, or none).
# 5. Saves the plot as a PNG file.
#
# Inputs:
# - Loop sets: merged_loops.tsv for each sample
# - Normalized IF values table: path given in paths.tsv
# - List of cell lines: for automated batch plotting
#
# Output:
# - Scatterplots showing IF comparison for each sample pair
#   (e.g., HepG2-control_GM12878-control.png)
#
# Usage:
# Run via Rscript with two sample names as arguments:
# Rscript this_script.R <sample1> <sample2>
#
# Author: Gony Shanel
# Date: March 2025
# ------------------------------------------------------------------------------

## Set the scripts directory as working directory.
## This script is designed with relative paths for the specific data structure
## with which it was uploaded to https://github.com/ElkonLab/micro-C.

# Functions ---------------------------------------------------------------

plotIF <- function(sample1, sample2, loop_index=fread("../../data/loops/pooled_loops_values_normalized.tsv")){
  loopSet1_path <- paste0("../../data/loops/", sample1, "_merged_loops.tsv")
  loopSet2_path <- paste0("../../data/loops/", sample2, "_merged_loops.tsv")
  
  loopSet1_raw <- fread(loopSet1_path, 
                        drop = 7:8,
                        col.names = colnames(loop_index)[1:6])
  loopSet2_raw <- fread(loopSet2_path, 
                        drop = 7:8,
                        col.names = colnames(loop_index)[1:6])
  
  loopSet1 <- left_join(loopSet1_raw, 
                        loop_index[, 1:7])[["peak_name"]]
  loopSet2 <- left_join(loopSet2_raw, 
                        loop_index[, 1:7])[["peak_name"]]
  
  common_set <- intersect(loopSet1, loopSet2)
  uniq1_set  <- setdiff(loopSet1, loopSet2)
  uniq2_set  <- setdiff(loopSet2, loopSet1)
  non_set    <- setdiff(loop_index$peak_name, union(loopSet1, loopSet2))
  
  grouped_loops <- rbind(data.table("peak_name" = non_set,    "Identification" = "None"),
                         data.table("peak_name" = uniq2_set,  "Identification" = sample2),
                         data.table("peak_name" = uniq1_set,  "Identification" = sample1),
                         data.table("peak_name" = common_set, "Identification" = "Common"))
  
  data_to_plot <- left_join(grouped_loops,
                            loop_index[, c("peak_name", sample1, sample2), with = F])
  data_to_plot$Identification <- factor(data_to_plot$Identification, 
                                        levels = c("Common", sample1, sample2, "None"),
                                        ordered = T)
  data_to_plot$log1 <- log2(data_to_plot[[sample1]])
  data_to_plot$log2 <- log2(data_to_plot[[sample2]])
  
  if (sub(".*-", "", sample1) == sub(".*-", "", sample2)){
    axes_labs <- c(sub("-.*", "", sample1), sub("-.*", "", sample2))
  } else {
    axes_labs <- c(sub(".*-", "", sample1), sub(".*-", "", sample2))
  }
  
  p <- ggplot(data_to_plot, aes(log1, log2, color = Identification)) +
    geom_point(size = 0.05, alpha = 0.5) + theme_bw() +
    scale_color_manual(values=c("#FF8000", "#00E15A", "#0080FF", "#A0A0A0")) +
    labs(x = paste0(axes_labs[1], " [log2 interaction frequency]"),
         y = paste0(axes_labs[2], " [log2 interaction frequency]"),
         caption = paste0("N(total)=", 
                          nrow(loop_index), ", N(",
                          axes_labs[1], ")=", 
                          length(loopSet1), ",\nN(",
                          axes_labs[2], ")=", 
                          length(loopSet2), ", N(common)=",
                          length(common_set))) +
    theme(text = element_text(size = 8),
          plot.caption = element_text(size = 5, hjust = 0),
          legend.position = "none")
  
  fn <- paste0(sample1, "_", sample2, ".png")
  ggsave(fn, p, device = "png", width = 2.5, height = 2.5)
}

# Arguments ---------------------------------------------------------------
args <- c(sample1 <- "HepG2-control",
          sample2 <- "GM12878-control")

# args <- commandArgs(trailingOnly = T)
sample1 <- args[1]
sample2 <- args[2]

# Command -----------------------------------------------------------------
library(tidyverse)
library(parallel)
library(data.table)

loop_index <- fread("../../data/loops/pooled_loops_values_normalized.tsv")

plotIF(sample1, sample2)


# ## Plot all control-nutlin plots
# setwd("/specific/elkon/gonyshanel/2020-03_hsieh/hic/dci/control_treatment/compare_IF/plot_IF/")
# mclapply(readLines(paths["cell_lines", "path"]),
#          function(x){
#            plotIF(paste0(x,"-control"), paste0(x, "-nutlin"))},
#          mc.cores = 10)
# 
# ## Plot all cell line comparisons
# setwd("/specific/elkon/gonyshanel/2020-03_hsieh/hic/dci/inter_cell_line/manual/if_pairwise_comparison/")
# dt <- as.data.table(t(combn(readLines(paths["cell_lines", "path"]), 2)))
# mclapply(1:nrow(dt),
#          function(i){
#            plotIF(paste0(dt[i,1],"-control"), paste0(dt[i,2],"-control"))},
#          mc.cores = 15)


