# ------------------------------------------------------------------------------
# Plot Interaction Frequency (IF) of Gene-Associated Loops by Expression Fold Change
# ------------------------------------------------------------------------------
# This script analyzes the relationship between **gene expression fold change (FC)**
# and **loop interaction frequency (IF)**. It performs the following steps:
#
# 1. Loads gene expression data and Hi-C loop interaction data.
# 2. Merges expression FC with loop IF FC for each gene.
# 3. Computes Pearson correlation between **expression FC** and **loop IF FC**.
# 4. Generates **boxplots** to visualize the relationship.
# 5. Saves correlation statistics and plots for each **cell line pair**.
# 6. Generates heatmaps summarizing correlation statistics across multiple comparisons.
#
# INPUT:
# - `cell_lines_pairs.txt`: List of cell line pairs to compare.
# - `gene_expression_CNVcorrected.tsv`: Corrected gene expression data.
# - `tssLoops_values_PP-loops.tsv`: Hi-C loop interaction frequency data.
#
# OUTPUT:
# - Boxplots showing **loop IF vs. gene expression FC** for each cell pair.
# - Heatmaps summarizing correlation p-values and rho values.
# - Summary statistics of correlation tests.
#
# AUTHOR: Gony Shanel
# DATE: March 2025
# ------------------------------------------------------------------------------

## Set the scripts directory as working directory.
## This script is designed with relative paths for the specific data structure
## with which it was uploaded to https://github.com/ElkonLab/micro-C.

# Setup ----

library(magrittr)
library(ggplot2)
library(dplyr)
library(data.table)
library(ggsci)
library(preprocessCore)
library(parallel)

# Arguments ----

pairs_path <- "../../data/cell_lines_pairs.txt"
exp_path <- "../../data//gene_expression_CNVcorrected.tsv"
loops_path <- "../../data/loops/tssLoops_values.tsv"

outdir <- "."
treatment <- "control"  # "control" | "nutlin"

# Commands ----

## Import and prepare data ----

loops_raw <- fread(loops_path, check.names = F)
loops_raw <- loops_raw[, !grepl("nutlin", colnames(loops_raw)), with = F]
colnames(loops_raw) <- sub("-control", "", colnames(loops_raw))

exp_raw <- fread(exp_path, check.names = F)
exp_raw <- exp_raw[, !grepl("nutlin", colnames(exp_raw)), with = F]
colnames(exp_raw) <- sub("-control", "", colnames(exp_raw))

pairs <- readLines(pairs_path)

# Initialize results storage
sum_stats <- data.frame()

# Pairwise Loop ----

pair <- pairs[1]  # for dev/test

sum_stats <- as.data.table(do.call(rbind, mclapply(pairs, function(pair) {
  dir.create(pair, recursive = T, showWarnings = F)
  
  cells <- strsplit(pair, "-")[[1]]
  cell1 <- cells[1]
  cell2 <- cells[2]
  
  ## Prepare data ----
  
  cols <- c(colnames(loops_raw)[1:7], "gene", "sym", cell1, cell2)
  loops <- loops_raw[, ..cols]
  loops$FC <- loops[[cell1]] / loops[[cell2]]
  loops$logFC <- log2(loops$FC)
  
  exp <- exp_raw[, c("gene", "sym", ..cell1, ..cell2)]
  exp$FC <- (exp[[cell1]] + 1) / (exp[[cell2]] + 1)
  exp$logFC <- log2(exp$FC)
  
  ## Merge and correlate ----
  
  exp_loops <- merge(
    exp[, c("gene", "logFC")],
    loops[, c(colnames(loops)[1:7], "gene", "logFC"), with = F],
    by = "gene",
    suffixes = c(".cpm", ".if")
  )
  
  exp_loops <- exp_loops[complete.cases(exp_loops), ]
  exp_loops$dis <- exp_loops$start2 - exp_loops$end1
  
  # Correlation: expression vs IF
  cor.IF <- cor.test(exp_loops$logFC.cpm, exp_loops$logFC.if, method = "pearson")
  r.IF <- format(cor.IF$estimate, digits = 2)
  pval.IF <- format(cor.IF$p.value, digits = 2)
  stats.IF <- paste0("r=", r.IF, "\np=", pval.IF)
  
  # Optional: correlation vs distance (disabled plot)
  cor.dis <- cor.test(exp_loops$logFC.cpm, exp_loops$dis, method = "pearson")
  r.dis <- format(cor.dis$estimate, digits = 2)
  pval.dis <- format(cor.dis$p.value, digits = 2)
  stats.dis <- paste0("r=", r.dis, "\np=", pval.dis)
  
  ## Bin by expression fold change ----
  
  exp_loops$round <- round(exp_loops$logFC.cpm)
  exp_loops$group <- sapply(exp_loops$round, function(x) {
    if (x < -6) {-6} else if (x > 6) {6} else {x}
  })
  
  data <- exp_loops %>%
    filter(group %in% c(-6:-2, 2:6))  # Keep |log2FC| >= 2
  data$round <- factor(data$group)
  
  ## Plot IF ~ Expression FC ----
  
  p.IF <- ggplot(data, aes(x = round, y = logFC.if, fill = round)) +
    geom_boxplot(outlier.shape = NA) +
    theme_bw() +
    scale_fill_brewer(palette = "RdBu", direction = -1) +
    theme(legend.position = "none", text = element_text(size = 8)) +
    xlab(expression('Gene expression [log'[2]*'FC]')) +
    ylab(expression('Loop intensity [log'[2]*'FC]')) +
    annotate("text", x = 7, y = -2.5, label = stats.IF, hjust = 0, vjust = 1, size = 2.5) +
    scale_x_discrete(labels = c(expression(phantom(x) <= -6), -5:-2, 2:5, expression(phantom(x) >= 6)))
  
  fn.IF <- paste0(outdir, "/", pair, "/gene-associated-loops_IF_pv-", pval.IF, "_r-", r.IF, ".png")
  ggsave(fn.IF, p.IF, width = 3, height = 2.5)
  
  # Optional plot for loop length vs gene FC group (disabled by default)
  # p.dis <- ggplot(data, aes(y = round, x = dis, fill = round)) +
  #   geom_boxplot(outlier.shape = NA) +
  #   ylim(0, 800000) +
  #   theme_bw() + scale_fill_brewer(palette = "RdBu") +
  #   theme(legend.position = "none",
  #         text = element_text(size = 8)) +
  #   ylab(expression('Gene expression log'[2]*FC)) +
  #   xlab("Loop length [bp]") +
  #   annotate(geom = "text", x = 8, y = 8, label = stats.dis,  size = 3) +
  #   scale_x_discrete(labels = c(expression(phantom(x) <=-6), -5:-2, 2:5, expression(phantom(x) >=6)))
  
  
  # Return stats
  c("cell1" = cell1, "cell2" = cell2, "N" = nrow(data), unlist(cor.IF))
}, mc.cores = 45)))

# Save Results ----

write.table(sum_stats, "summary_statistics.tsv",
            append = F, quote = F, sep = '\t', row.names = F, col.names = T)

# Summary Heatmaps ----

## Format stats ----
sum_stats$p.value <- as.numeric(sum_stats$p.value)
sum_stats$estimate.cor <- as.numeric(sum_stats$estimate.cor)
sum_stats$cell1 <- factor(sum_stats$cell1)
sum_stats$cell2 <- factor(sum_stats$cell2)

## Heatmap: p-values ----
ggplot(sum_stats, aes(cell1, cell2, fill = -log10(p.value))) +
  geom_tile() + theme_minimal() +
  scale_fill_distiller(palette = "Reds", direction = 1, na.value = "black") +
  labs(fill = expression("-log"[10]*P)) +
  geom_text(aes(label = round(-log10(p.value), 1)), size = 2) +
  theme(
    text = element_text(size = 8),
    legend.position = c(0.9, 0.7),
    panel.grid.major = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  ) +
  scale_y_discrete(limits = rev(levels(sum_stats$cell2))) +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 4))

ggsave("p.heatmap_pval_multiple_comparisons.png", width = 3, height = 3)

## Heatmap: correlation coefficients ----
ggplot(sum_stats, aes(cell1, cell2, fill = estimate.cor)) +
  geom_tile() + theme_minimal() +
  scale_fill_distiller(
    palette = "Reds",
    direction = 1,
    na.value = "black",
    limits = c(0, max(sum_stats$estimate.cor))
  ) +
  labs(fill = "r") +
  geom_text(aes(label = round(estimate.cor, 2)), size = 2) +
  theme(
    text = element_text(size = 8),
    legend.position = c(0.9, 0.7),
    panel.grid.major = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  ) +
  scale_y_discrete(limits = rev(levels(sum_stats$cell2))) +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 4))

ggsave("p.heatmap_r_multiple_comparisons.png", width = 3, height = 3)
