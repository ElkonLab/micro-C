#!/usr/local/lib/R-3.5.3/lib/R/bin/Rscript

# ------------------------------------------------------------------------------
# Expression vs Loop Signal Dynamic Range: Cell Identity vs Nutlin-3a Response
#
# This script compares the dynamic range of expression and loop signal 
# fold changes in two contexts: (1) between cell lines (cell identity), 
# and (2) within a cell line following nutlin treatment (treatment response).
#
# Steps:
# 1. Load CNV-corrected expression and Hi-C loop interaction data.
# 2. Compute fold changes between control conditions of two cell lines,
#    and between control vs nutlin within the same cell.
# 3. Visualize log2 fold change distributions using density, violin,
#    boxplot, and ECDF plots.
# 4. Compare dynamic ranges across expression and loop signal levels.
#
# Inputs:
# - gene_expression_CNVcorrected.tsv
# - mustache_5k_pos_TSS (loop matrix)
# - cell_pairs list (e.g., HCT116-U2OS)
# - paths.tsv
#
# Outputs:
# - Density, violin, boxplot, ECDF plots for EXP and IF
#
# Author: Gony Shanel
# Date: March 2025
# ------------------------------------------------------------------------------

## Set the scripts directory as working directory.
## This script is designed with relative paths for the specific data structure
## with which it was uploaded to https://github.com/ElkonLab/micro-C.

# Arguments ----

path.exp <- "../../data/gene_expression_CNVcorrected.tsv"

# Setup ----

## Set environment ----
library(data.table)
library(ggplot2)

# Load Data ----

## Import data ----
raw.exp <- fread(path.exp)
exp <- raw.exp[, -c(1:2)] + 1  # Drop metadata columns, add pseudocount

loops <- fread("../../data/loops/positive_correlating_0.5_tssLoops.tsv")
pairs <- readLines("../../data/cell_lines_pairs.txt")

pair <- pairs[[1]]  # Use first cell pair (e.g., HCT116-U2OS)

# Expression Fold Change Comparisons ----

## mRNA expression ----

splitPair <- strsplit(pair, "-")[[1]]
cell1 <- splitPair[[1]]
cell2 <- splitPair[[2]]

# Fold change: between cell lines (control vs control)
FC.interCells <- exp[[paste0(cell2, "-control")]] / exp[[paste0(cell1, "-control")]]

# Fold change: within cell1 (nutlin vs control)
FC.controlTreatment <- exp[[paste0(cell1, "-nutlin")]] / exp[[paste0(cell1, "-control")]]

# Combine into one long-format table
cpm.dynamicRanges <- rbind(
  data.table(value = FC.interCells,      type = "cell idnetity"),
  data.table(value = FC.controlTreatment, type = "treatment response")
)
cpm.dynamicRanges$type <- factor(cpm.dynamicRanges$type,
                                 levels = c("cell idnetity", "treatment response"),
                                 labels = c(paste0(cell2, "\nvs. ", cell1), paste0("Response to\nnutlin (", cell1, ")")),
                                 ordered = T)

# Plot: density plot of expression log2 fold change
ggplot(cpm.dynamicRanges, aes(log2(value), color = type)) +
  geom_density() + theme_bw() +
  scale_color_brewer(palette = "Set1") +
  coord_cartesian(xlim = c(-5, 5)) +
  labs(x = "expression fold change (log2)") +
  theme(text = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"),
        legend.box.spacing = unit(0.2, "cm"),
        legend.position = "top",
        legend.title = element_blank())
ggsave("p.density.exp_dynmaic_range.png", height = 2.5, width = 2.5)

# Plot: violin + boxplot
ggplot(cpm.dynamicRanges, aes(type, log2(value))) +
  geom_violin() + theme_bw() +
  geom_boxplot(outlier.shape = NA, width = 0.1) +
  scale_color_brewer(palette = "Set1") +
  labs(y = "expression fold change (log2)") +
  theme(text = element_text(size = 8),
        axis.title.x = element_blank())
ggsave("p.violin.exp_dynamic_range.png", height = 2.5, width = 2.5)

# Plot: boxplot only
ggplot(cpm.dynamicRanges, aes(type, log2(value))) +
  theme_bw() +
  geom_boxplot(outlier.shape = NA) +
  scale_color_brewer(palette = "Set1") +
  coord_cartesian(ylim = c(-5, 5)) +
  labs(y = "expression fold change (log2)") +
  theme(text = element_text(size = 8),
        axis.title.x = element_blank())
ggsave("p.boxplot.exp_dynamic_range.png", height = 2.5, width = 2.5)

# Plot: ECDF
ggplot(cpm.dynamicRanges, aes(log2(value), color = type)) +
  stat_ecdf(geom = "step") + theme_bw() +
  labs(x = "| log2 expression fold change |", y = "cumulative density") +
  scale_color_brewer(palette = "Set1") +
  theme(text = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"),
        legend.box.spacing = unit(0.2, "cm"),
        legend.position = "top",
        legend.title = element_blank())
ggsave("p.cumulative.exp_dynmaic_range.png", height = 2.2, width = 2)

# Plot: absolute log2 fold change ECDF
ggplot(cpm.dynamicRanges, aes(abs(log2(value)), color = type)) +
  stat_ecdf(geom = "step") + theme_bw() +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
  labs(x = "| log2 expression fold change |", y = "cumulative density") +
  scale_color_brewer(palette = "Set1") +
  theme(text = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"),
        legend.box.spacing = unit(0.2, "cm"),
        legend.position = "top",
        legend.title = element_blank())
ggsave("p.cumulative.abs_exp_dynmaic_range.png", height = 2.2, width = 2)

# Loop Fold Change Comparisons ----

splitPair <- strsplit(pair, "-")[[1]]
cell1 <- splitPair[[1]]
cell2 <- splitPair[[2]]

# Fold change: loop intensity between cells
FC.interCells <- loops[[paste0(cell2, "-control")]] / loops[[paste0(cell1, "-control")]]

# Fold change: loop intensity under nutlin
FC.controlTreatment <- loops[[paste0(cell1, "-nutlin")]] / loops[[paste0(cell1, "-control")]]

IF.dynamicRanges <- rbind(
  data.table(value = FC.interCells, type = "cell idnetity"),
  data.table(value = FC.controlTreatment, type = "treatment response")
)

# Plot: loop IF dynamic range density
ggplot(IF.dynamicRanges, aes(log2(value), color = type)) +
  geom_density() + theme_bw() +
  coord_cartesian(xlim = c(-5, 5)) +
  labs(x = "log2 intensity fold change") +
  theme(text = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"),
        legend.box.spacing = unit(0.2, "cm"),
        legend.position = "top",
        legend.title = element_blank())
ggsave("p.density.IF_dynmaic_range.png", height = 2.2, width = 2)

# Compare EXP vs IF for Cell Identity ----

dynamicRanges <- rbind(
  cbind(cpm.dynamicRanges, datatype = "EXP"),
  cbind(IF.dynamicRanges, datatype = "IF")
)

# Plot: density comparison for EXP vs IF (cell identity only)
ggplot(dynamicRanges[type == "cell idnetity"], aes(log2(value), color = datatype)) +
  geom_density() + theme_bw() +
  coord_cartesian(xlim = c(-5, 5)) +
  labs(x = "log2 fold change") +
  theme(text = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"),
        legend.box.spacing = unit(0.2, "cm"),
        legend.position = "top",
        legend.title = element_blank())
ggsave("p.density.cell_identity_dynmaic_range.png", height = 2.2, width = 2)

