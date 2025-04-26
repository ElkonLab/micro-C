#!/usr/local/lib/R-3.5.3/lib/R/bin/Rscript

# ------------------------------------------------------------------------------
# Expression Fold Change of Genes Near p53 ChIP-seq Peaks
#
# This script compares the expression response (log2 fold change) of genes 
# located within 50 kb of p53-induced p53 ChIP-seq peaks to all other genes 
# following nutlin treatment. It evaluates statistical enrichment of gene 
# induction near peaks and categorizes genes by expression fold change groups.
#
# Steps:
# 1. Load expression data, ChIP-seq peak distances, and significant DE genes.
# 2. For each cell line:
#    - Identify genes near peaks (< 50 kb).
#    - Compare fold change of these genes to background using boxplots.
#    - Test statistical significance and annotate plots.
#    - Output plots and per-cell results.
# 3. Analyze enrichment of significantly induced genes in fold change bins.
#
# Inputs:
# - paths.tsv
# - peaks_distance_to_closest_gene.tsv
# - cell_lines.txt
# - mean_gene_expression (from paths)
# - all_significant_differentail_genes.tsv
#
# Outputs:
# - Boxplots comparing nearby vs. background gene fold changes
# - Summary statistics per cell line
# - Induction enrichment per fold change bin
#
# Author: Gony Shanel
# Date: March 2025
# ------------------------------------------------------------------------------

## Set the scripts directory as working directory.
## This script is designed with relative paths for the specific data structure
## with which it was uploaded to https://github.com/ElkonLab/micro-C.

# Arguments ----

path.distanceToClosest <- "../../data/p53_ChIP-seq_peaks_distance_to_closest_gene.tsv"
path.cells <- "../../data/p53_cell_lines.txt"
path.deseq_sig_diff_genes <- "../../data/deseq2/all_significant_differentail_genes.tsv"

maxDistanceToGene <- 50000

# Setup ----

## Set environment ----
library(data.table)
library(parallel)
library(ggplot2)
library(ggpubr)
library(ggsignif)

dir.create("plots/closest_gene", recursive = T, showWarnings = F)
dir.create("results/gene_induction_by_FC_groups/", recursive = T, showWarnings = F)

# Load Data ----

raw.expression <- fread("../../data/gene_expression_CNVcorrected.tsv")
raw.distanceToClosest <- fread(path.distanceToClosest)
cells <- readLines(path.cells)

# Compare Expression Response of Nearby Genes ----

c <- cells[[2]]  # For debugging

stats <- as.data.table(do.call(rbind, lapply(cells, function(c){
  
  print(c)
  
  # Identify genes near p53 peaks (< 50kb)
  genes <- raw.distanceToClosest[distance < 50000 & cell == c, gene]
  
  # Extract expression columns for control and nutlin conditions
  columns <- grep(c, names(raw.expression), value = T)
  
  # Compute fold change for all genes, mark if near peaks
  fc <- data.table(
    gene = raw.expression[["gene"]],
    log2foldChange = log2((raw.expression[[columns[2]]] + 1) /
                            (raw.expression[[columns[1]]] + 1)),
    close = raw.expression[["gene"]] %in% genes
  )
  
  fc <- fc[complete.cases(fc)]
  
  # Wilcoxon rank test (greater) for logFC near peaks vs. others
  C <- compare_means(log2foldChange ~ close, fc, alternative = "greater")
  
  ## Plot Fold Change Comparison ----
  ggplot(fc, aes(close, log2foldChange)) +
    geom_boxplot(outlier.shape = NA) +
    theme_bw() +
    coord_cartesian(ylim = c(
      quantile(fc$log2foldChange, 0.01),
      quantile(fc$log2foldChange, 0.99)
    )) +
    scale_fill_brewer(palette = "Set1") +
    geom_signif(
      comparisons = list(c("FALSE", "TRUE")),
      annotations = C$p.adj,
      y_position = quantile(fc$log2foldChange, 0.9),
      tip_length = 0,
      extend_line = -0.02,
      textsize = 2.5
    ) +
    scale_x_discrete(labels = paste0(
      c("background\ngenes", "nearest\ngenes"),
      "\n(N=", table(fc$close), ")"
    )) +
    labs(
      y = "log2 expression fold change",
      title = c
    ) +
    theme(
      text = element_text(size = 8),
      legend.position = "none",
      axis.title = element_blank(),
      axis.title.x = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  
  ggsave(paste0("plots/closest_gene/", c, ".png"),
         height = 2, width = 1.5)
  
  # Return test results
  c(cell = c, n_close = sum(fc$close), C)
  
})))

# Save Statistics ----

fwrite(stats,
       "results/genes_close_to_peak_compared_to_all_others.tsv",
       append = F, quote = F, sep = '\t', row.names = F, col.names = T)

# Summary Plot of P-values ----

ggplot(stats, aes(unlist(cell), -log10(unlist(p.adj)), group = 1)) +
  geom_line() + geom_point() +
  theme_bw() +
  labs(y = "-log10 P-value") +
  geom_hline(yintercept = -log10(0.05)) +
  scale_x_discrete(labels = paste0(
    stats$cell,
    "\n(N=", stats$n_close, ")"
  )) +
  theme(
    text = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 0.8)
  )

ggsave("plots/p.sum_pval.png", height = 2, width = 3)

# Analyze Enrichment of Significantly Induced Genes per Fold Change Group ----

sig_diff_genes <- fread(path.deseq_sig_diff_genes)

c <- cells[[6]]  # For dev/testing

stats <- as.data.table(do.call(rbind, lapply(cells, function(c){
  
  print(c)
  
  # Identify genes close to p53 peaks
  genes <- raw.distanceToClosest[distance < 50000 & cell == c, gene]
  columns <- grep(c, names(raw.expression), value = T)
  
  # Compute expression fold change
  fc <- data.table(
    gene = raw.expression[["gene"]],
    log2foldChange = log2((raw.expression[[columns[2]]] + 1) /
                            (raw.expression[[columns[1]]] + 1)),
    significantly_induced = raw.expression[["gene"]] %in% sig_diff_genes[cell == c, gene]
  )
  
  fc <- fc[complete.cases(fc)]
  
  # Classify genes by fold change
  fc$fc_group <- fc[, ifelse(
    log2foldChange < log2(1.2), "non",
    ifelse(log2foldChange > log2(1.7), "high", "mild")
  )]
  
  # Compute % of significant genes per FC group
  true_induction <- fc[, .(
    N_genes = .N,
    `significantly_induced [%]` = mean(significantly_induced) * 100
  ), by = fc_group]
  
  true_induction[, fc_group := factor(fc_group, levels = c("non", "mild", "high"))]
  setorder(true_induction, fc_group)
  
  # Save per-cell summary table
  fwrite(true_induction,
         paste0("results/gene_induction_by_FC_groups/", c, ".tsv"),
         sep = '\t')
  
})))
