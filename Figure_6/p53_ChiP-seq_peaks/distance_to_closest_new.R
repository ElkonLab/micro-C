# ------------------------------------------------------------------------------
# Distance from Induced p53 ChIP-seq Peaks to Nearest TSS (Filtered Cell Lines)
#
# This script calculates the genomic distance between induced p53 ChIP-seq peaks 
# (MACS2-called) and the closest annotated protein-coding TSS (hg38) for a selected 
# subset of cell lines (excluding HEK293 and HeLa). It produces distance tables 
# and comparative visualizations.
#
# Steps:
# 1. Load ChIP-seq peak files and TSS annotations.
# 2. For each cell line:
#    - Identify the closest TSS for each peak using bedtools.
#    - Compile distances across all selected cell lines.
# 3. Save the results in a table.
# 4. Generate plots:
#    - Density plots of log10-transformed distances across cell lines.
#    - Boxplots of distances per cell line, including a plot restricted to peaks 
#      located within 50kb from TSS.
#
# Inputs:
# - Induced p53 ChIP-seq peaks: ../../data/induced_p53-ChIP-seq_peaks/
# - TSS annotation: ../../data/ensembl_mart_hg38_1st_TSS_protein_coding.txt
# - Cell lines list: ../../data/cell_lines.txt
#
# Outputs:
# - peaks_distance_to_closest_gene.tsv: Table with closest gene and distance per peak.
# - plots/p.boxplot.distance_to_nearest_gene.png: Boxplot of log10 distances (all peaks).
# - plots/p.boxplot.distance_to_nearest_gene_50kbp-lim.png: Boxplot (peaks <50kbp).
# - p.density.png: Density plot of log10 peak-to-TSS distances by cell line.
#
# Author: Gony Shanel
# Date: April 2025
# ------------------------------------------------------------------------------

## Set the scripts directory as working directory.
## This script is designed with relative paths for the specific data structure
## with which it was uploaded to https://github.com/ElkonLab/micro-C.

## This script requires bedtools to be installed (https://bedtools.readthedocs.io/en/latest/index.html)
path.bedtools <- "/path/to/bedtools"

# Arguments ---------------------------------------------------------------

dir.chipseqPeaks <- "../../data/induced_p53-ChIP-seq_peaks/"
path.tss <- "../../data/ensembl_mart_hg38_1st_TSS_protein_coding.txt"

cells <- readLines("../../data/cell_lines.txt")
chipCells <- grep("HEK293|HeLa", cells, value = T, invert = T)

wd <- "."

# commands ----------------------------------------------------------------

## Set environment

library(data.table)
library(parallel)
library(ggplot2)
library(bedtoolsr)
options(bedtools.path = path.bedtools)

## Import data ----
paths.chipseq <- list.files(dir.chipseqPeaks, "narrow", recursive = T, full.names = T)
names(paths.chipseq) <- cells
tss <- fread(path.tss)

## Find closest ----
c <- chipCells[8]
distanceToClosest <- as.data.table(do.call(rbind, lapply(chipCells, function(c){
    
  chipseq <- fread(paths.chipseq[c])
  genes <- tss
  
  distanceToClosest <- as.data.table(bt.closest(bt.sort(chipseq),
                                                bt.sort(genes[, 1:4]),
                                                d = T))
  distanceToClosest <- distanceToClosest[,c(14,15)]
  names(distanceToClosest) <- c("gene", "distance")
  cbind(distanceToClosest[distance > 0], "cell" = c)
})))

write.table(distanceToClosest, "peaks_distance_to_closest_gene.tsv",
            append = F, quote = F, sep = '\t', row.names = F, col.names = T)

## Plot ----
distanceToClosest$distance <- as.numeric(distanceToClosest$distance)
ggplot(distanceToClosest, aes(log10(distance))) +
  geom_density() + theme_bw() +
  labs(x = "log10 distance [bp]") +
  theme(text = element_text(size = 8)) +
  facet_wrap(~cell, nrow = 2)
ggsave("p.density.png",
       height = 2, width = 4)


## Plot boxplot ----
dir.create("plot")

ggplot(distanceToClosest, aes(cell, log10(distance))) +
  geom_boxplot(outlier.shape = NA) + theme_bw() +
  coord_flip() +
  labs(y = "log10 distance to nearest gene") +
  theme(text = element_text(size = 8),
        axis.title.y = element_blank())
ggsave("plots/p.boxplot.distance_to_nearest_gene.png",
       height = 3, width = 3)  

# Limit distance to 50 kbp
dataToPlot <- distanceToClosest[distance < 50000]

ggplot(dataToPlot, aes(cell, log10(distance))) +
  geom_boxplot(outlier.shape = NA) + theme_bw() +
  coord_flip() +
  scale_x_discrete(limits = rev) +
  labs(y = "log10 distance to nearest gene") +
  theme(text = element_text(),
        axis.title.y = element_blank())
ggsave("plots/p.boxplot.distance_to_nearest_gene_50kbp-lim.png",
       height = 3, width = 3.2)  

       