#!/usr/local/lib/R-3.5.3/lib/R/bin/Rscript

# ------------------------------------------------------------------------------
# Association Between Loop Clusters and Gene Expression
# ------------------------------------------------------------------------------
# This script evaluates how clusters of Hi-C loops, as identified by Expander's 
# CLICK algorithm, correspond to gene expression patterns across cell lines.
#
# The analysis includes:
# 1. Loading CLICK clustering results.
# 2. Mapping clustered loops to associated genes.
# 3. Standardizing interaction frequency (IF) and gene expression (GE) across cell lines.
# 4. Calculating the mean IF and GE per cluster for each cell line.
# 5. Testing the correlation between IF and GE (Pearson's r) per cluster.
# 6. Plotting cluster-level IF and GE trends and summarizing correlation statistics.
#
# Input:
# - CLICK clustering output: click_clusters.txt
# - Normalized Hi-C interaction frequencies: tssLoops_values.tsv
# - CNV-corrected gene expression: gene_expression_CNVcorrected.tsv
#
# Output:
# - tssLoop_clusters.tsv: loops and their assigned clusters
# - mean_pattern.png: faceted trends of IF and GE across clusters
# - plots/: individual cluster plots
# - clusters_IF_GE_correlation.png: summary of correlations across clusters
#
# Author: Gony Shanel
# Date: March 2025
# ------------------------------------------------------------------------------

## Set the scripts directory as working directory.
## This script is designed with relative paths for the specific data structure
## with which it was uploaded to https://github.com/ElkonLab/micro-C.

# Setup and Input Paths ----

# args <- comif (length(args) != 3){
#   stop("Usage: analyze_CLICK_loop_clusters_associated gene_expressoin.r {\"clustering\" dir} {CLICK clusters} {tssLoops file)")
# }

# Example hardcoded paths (used instead of command line args)
click_path <- paste0("../../data/loops/click_clusters.txt")
tssLoops_raw_path <- "../../data/loops//tssLoops_values.tsv"
exp_path <- "../../data/gene_expression_CNVcorrected.tsv"

# Alternative command-arg version (disabled)
# wd <- R.utils::getAbsolutePath(args[1])
# click_path <- R.utils::getAbsolutePath(args[2])
# tssLoops_raw_path <- R.utils::getAbsolutePath(args[3])

# Load Libraries and Utilities ----
library(data.table, quietly = T)
library(tidyverse, quietly = T)
library(patchwork, quietly = T)
library(ggpubr, quietly = T)
source("../../_functions.r")

# Load and Standardize Hi-C and Expression Data ----

## Standardize IF data ----
tssLoops_raw <- fread(tssLoops_raw_path)
tssLoops <- tssLoops_raw[, !grepl("nutlin", colnames(tssLoops_raw)), with = F]
colnames(tssLoops) <- sub("-control", "", colnames(tssLoops))
IF_std <- data.frame(
  t(scale(t(tssLoops[, 8:17]))),
  tssLoops[, c("peak_name", "gene", "sym")],
  check.names = F
)

## Standardize gene expression data ----
exp_raw <- fread(exp_path)
exp <- exp_raw[, !grepl("nutlin", colnames(exp_raw)), with = F]
colnames(exp) <- sub("-control", "", colnames(exp))
exp_std <- data.frame(
  t(scale(t(exp[, 3:12]))),
  exp[, c("gene", "sym")],
  check.names = F
)

# Parse CLICK Clusters ----

# clst_dic <- read.delim(clst_dic_path, header = F, row.names = 2) 
click <- fread(click_path, header = F)
click_cluster_names <- unique(click$V2)[grep("Cluster", unique(click$V2))]

# Initialize storage for loop-wise and cluster-wise results
loopSets <- data.frame()
geneSets <- data.frame()
clusters_mean_IF <- data.frame()
clusters_mean_exp <- data.frame()
corr <- data.table()

# Per-Cluster Analysis ----

cluster_name = click_cluster_names[2]  # Example index for dev/debug
for (cluster_name in click_cluster_names) {
  loopSet <- click[click$V2 == cluster_name,]
  tssLoopSet <- data.frame(IF_std[IF_std$peak_name %in% loopSet$V1,], "cluster" = cluster_name, check.names = F)
  geneSet <- data.frame(exp_std[exp_std$gene %in% tssLoopSet$gene,], "cluster" = cluster_name, check.names = F)
  
  ## Correlate IF with GE ----
  mean_IF <- colMeans(tssLoopSet[, 1:10], na.rm = T)
  mean_exp <- colMeans(geneSet[, 1:10])
  c <- cor.test(mean_IF, mean_exp, method = "pearson")
  corr <- rbind(corr, data.frame("cluster" = cluster_name, "r" = c$estimate, "p" = c$p.value))
  
  ## Aggregate loop and expression data ----
  loopSets <- rbind(loopSets, tssLoopSet)
  geneSets <- rbind(geneSets, geneSet)
  clusters_mean_IF <- rbind(clusters_mean_IF, data.frame("mean" = mean_IF, "cluster" = cluster_name, "cell_line" = names(mean_IF), check.names = F))
  clusters_mean_exp <- rbind(clusters_mean_exp, data.frame("mean" = mean_exp, "cluster" = cluster_name, "cell_line" = names(mean_exp), check.names = F))
}

# Save Results ----
write.tsv(loopSets[, 11:14], "tssLoop_clusters.tsv", header = T)

# Prepare Data for Plotting ----

clusters_mean_IF$data_type <- "IF"
clusters_mean_exp$data_type <- "GE"
data_to_plot_means <- rbind(clusters_mean_exp, clusters_mean_IF)

loopSets$cluster <- factor(loopSets$cluster, levels = unique(loopSets$cluster), ordered = T)
geneSets$cluster <- factor(geneSets$cluster, levels = unique(loopSets$cluster), ordered = T)

titles <- paste0(
  unique(loopSets$cluster),
  "\n(Nl=", table(loopSets[!duplicated(loopSets$peak_name), "cluster"]),
  " Ng=", table(geneSets$cluster), ")"
)

clst_title <- as.data.frame(cbind("cluster" = levels(loopSets$cluster), "cluster_title" = titles))
clst_title$cluster <- factor(clst_title$cluster, unique(clst_title$cluster), ordered = T)
clst_title$cluster_title <- factor(clst_title$cluster_title, unique(clst_title$cluster_title), ordered = T)

data_to_plot_means <- as.data.table(left_join(data_to_plot_means, clst_title, by = "cluster"))
data_to_plot_means$cluster <- factor(data_to_plot_means$cluster, unique(data_to_plot_means$cluster), ordered = T)

corr$label <- paste0("r=", formatC(corr$r, 2), ", p=", formatC(corr$p, format = "e", 2))
corr$cluster_title <- unique(data_to_plot_means$cluster_title)
corr$data_type <- "IF"

# Plot Faceted Cluster Trends ----

p <- ggplot(data = data_to_plot_means[cluster %in% levels(data_to_plot_means$cluster)[1:20]],
            mapping = aes(cell_line, mean, color = data_type, group = data_type)) +
  geom_line() + geom_point() + theme_bw() +
  labs(color = "Data type") +
  scale_color_brewer(palette = "Set1") +
  theme(text = element_text(size = 8),
        axis.text.x = element_text(size = 7, angle = 90, hjust = 0.95, vjust = 0.2),
        axis.title = element_blank()) +
  coord_cartesian(ylim = c(-2, 2.5)) +
  facet_wrap(~cluster_title, ncol = 4)

p + geom_text(data = corr[1:20,], aes(x = 1, y = -1.5, label = label), hjust = 0, size = 2.5, color = "black")

ggsave("mean_pattern.png", width = 5, height = 6)

# Plot Individual Cluster Panels ----

clst <- click_cluster_names[2]  # Example for dev
dir.create("plots", showWarnings = F)

for (clst in click_cluster_names) {
  data_to_plot <- data_to_plot_means[cluster == clst,]
  ggplot(data_to_plot, aes(cell_line, mean, color = data_type, group = data_type)) +
    geom_line() + geom_point() + theme_bw() +
    coord_cartesian(ylim = c(-1, 2)) +
    scale_color_brewer(palette = "Set1") +
    annotate("text", x = 3, y = -0.7, label = corr[cluster == clst, label], size = 2.5) +
    labs(color = data_to_plot[1, cluster_title]) +
    theme(text = element_text(size = 8),
          axis.text.x = element_text(size = 7, angle = 90, hjust = 0.95, vjust = 0.2),
          axis.title = element_blank())
  
  ggsave(paste0("plots/mean_pattern_", clst, ".png"), width = 3.5, height = 2)
}

# Correlation Summary Plots ----

corr$logP <- -log10(corr$p)

p_r <- ggplot(corr, aes(y = r)) +
  geom_boxplot() + theme_bw() +
  ylab("r") +
  ggtitle("Clusters IF~GE", paste0("(n=", length(unique(click$V2)) - 1, ")")) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 8))

p_p <- ggplot(corr, aes(y = logP)) +
  geom_boxplot() + theme_bw() +
  ylab(expression(-log[10]*P)) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 8))

p_corr <- p_r + p_p

ggsave("clusters_IF_GE_correlation.png", p_corr, height = 3, width = 3)
