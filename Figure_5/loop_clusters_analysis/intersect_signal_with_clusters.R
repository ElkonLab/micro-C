
# ------------------------------------------------------------------------------
# Heatmap Visualization of Epigenomic Signals Over CLICK Cluster Anchors
#
# This script generates heatmaps of normalized epigenomic signal intensities
# (e.g., CTCF, H3K27ac, H3K4me1) over genomic anchors associated with loop
# clusters defined by Expander's CLICK. Anchors are annotated as promoter
# or distal, and plotted separately for each cell line.
#
# Steps:
# 1. Load loop and cluster assignments, as well as normalized epigenomic signals.
# 2. Identify anchor points associated with loops per cell line.
# 3. Annotate anchors (promoter vs. distal).
# 4. Generate z-score-scaled heatmaps for both promoter and distal anchors.
# 5. Save results and plots for each cell line.
#
# Inputs:
# - CLICK cluster assignments: tssLoop_clusters.tsv
# - Cluster-gene dictionary: cluster_dictionary.tsv
# - Normalized signal over anchors: mean_signal_over_pooled_anchors_normalized.tsv
# - Loop list: mustache_5k_TSS
# - Cell line list: cell_lines
#
# Outputs:
# - Annotated anchor tables per cell line (promoter/distal)
# - Heatmap PNGs for each marker and anchor type
#
# Author: Gony Shanel
# Date: March 2025
# ------------------------------------------------------------------------------

## Set the scripts directory as working directory.
## This script is designed with relative paths for the specific data structure
## with which it was uploaded to https://github.com/ElkonLab/micro-C.

## This script requires bedtools to be installed (https://bedtools.readthedocs.io/en/latest/index.html)
path.bedtools <- "/path/to/bedtools"

# Functions ---------------------------------------------------------------

intersect_with_tss <- function(intervals,
                               tss.file="/specific/elkon/gonyshanel/data/public/ensembl_hg38_protein_coding_first_TSS_5000win.bed"){
  require(bedtoolsr)
  require(data.table)
  options(bedtools.path = path.bedtools)
  
  df <- data.frame(intervals[,1:3])
  tss <- fread(tss.file)
  
  inter <- bt.intersect(a = df, b = tss, c = T)
  out <- inter[,4] != 0
}

annotate_loops <- function(bedpe,
                           tss.file="../../data/ensembl_hg38_protein_coding_first_TSS_5000win.bed"){
  annot1 <- intersect_with_tss(bedpe[,1:3], tss.file)
  annot1 <- ifelse(annot1 == TRUE, "P", "D")
  annot2 <- intersect_with_tss(bedpe[,4:6], tss.file)
  annot2 <- ifelse(annot2 == TRUE, "P", "D")
  paste0(annot1, annot2)
}

plot.heatmap <- function(dataToPlot){
  ggplot(dataToPlot, aes(variable, anchor_name, fill = value)) +
    geom_tile() + theme_minimal() +
    scale_fill_distiller(palette = "RdBu",
                         direction = -1) +
    theme(text = element_text(size = 8),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(angle = 90,
                                     vjust = .5,
                                     hjust = 1),
          legend.position = "none",
          panel.grid.major = element_blank())
}

# Arguments ---------------------------------------------------------------

path.tssLoopClusters <- "../../data/loops/tssLoop_clusters.tsv"
path.clusterDictioary <- "../../data/loops/cluster_dictionary.tsv"
path.normalizedSignalOverAnchors <- "../../data/loops/H3K27ac_mean_signal_over_pooled_anchors_normalized.tsv"

marker <- "H3K27ac"

# Commands ----------------------------------------------------------------

## Set environment ----
library(data.table)
library(parallel)
library(tidyverse)
source("../../_functions_hic.r")
source("../../_functions.r")

dir.create("results")
dir.create("plots")

## Load data ----
cells <- readLines("../../data/cell_lines.txt")
raw.loops <- fread("../../data/loops/tssLoops_values.tsv")
raw.tssLoopClusters <- fread(path.tssLoopClusters)
clusterDictioary <- fread(path.clusterDictioary)
normalizedSignalOverAnchors <- fread(path.normalizedSignalOverAnchors)
cells <- cells[cells %in% colnames(normalizedSignalOverAnchors)]

## Find cell associated anchor sets ----
loopsSets <- right_join(right_join(raw.loops[, 1:7], 
                                   raw.tssLoopClusters, 
                                   by = "peak_name"),
                        clusterDictioary)
## Remove PP-loops
loopsSets$annotation <- annotate_loops(loopsSets)
loopsSets <- loopsSets[annotation != "PP"]
## Split into anchors
anchorSets <- as.data.table(loops2Anchors(loopsSets,
                                          remove.duplicates = F,
                                          nameing = "suffix"))

setsSignalOverAnchors <- left_join(anchorSets[, .(chr, start, end, anchor_name, cell_affiliation )],
                                   normalizedSignalOverAnchors)

cell <- cells[1]
allDataToPlots <- data.table()
dataToPlot <- as.data.table(do.call(rbind, lapply(cells, function(cell){
  
  setSignalOverAnchors <- setsSignalOverAnchors[cell_affiliation == cell]
  setSignalOverAnchors$promoter <- annotate_anchors(setSignalOverAnchors)
  
  a.promoter <- setSignalOverAnchors[promoter == T]
  a.distal <- setSignalOverAnchors[promoter == F]
  
  ## Save data
  write.table(a.promoter, paste0("results/", cell, "_promoter_anchors.tsv"),
              append = F, quote = F, sep = '\t', row.names = F, col.names = T)
  write.table(a.distal, paste0("results/", cell, "_distal_anchors.tsv"),
              append = F, quote = F, sep = '\t', row.names = F, col.names = T)
  
  ## Plot ----
  
  ### Promoter anchors
  tab <- a.promoter[, c("anchor_name", cells), with = F]
  tab[,2:ncol(tab)] <- as.data.table(t(scale(t(tab[,2:ncol(tab)]))))
  dataToPlot <- melt(tab)
  
  plot.heatmap(dataToPlot)
  
  ggsave(paste0("plots/p.heatmap.", cell, "_promoter_anchors.png"),
         height = 2, width = 1.5)
  
  ### Distal anchors
  tab <- a.distal[, c("anchor_name", cells), with = F]
  tab[,2:ncol(tab)] <- as.data.table(t(scale(t(tab[,2:ncol(tab)]))))
  dataToPlot <- melt(tab)
  
  plot.heatmap(dataToPlot)
  
  ggsave(paste0("plots/p.heatmap.", cell, "_distal_anchors.png"),
         height = 2, width = 1.5)
  
  cbind(dataToPlot, "cell" = cell)
  
})))

