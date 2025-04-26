#!/usr/local/lib/R-3.5.3/lib/R/bin/Rscript

# ------------------------------------------------------------------------------
# Over-Representation Analysis of Gene Sets in CLICK Loop Clusters
#
# This script performs GO enrichment (Biological Process) on gene sets associated
# with Hi-C loop clusters defined by Expander's CLICK algorithm. It evaluates 
# functional enrichment using the clusterProfiler package and summarizes enriched
# terms across all clusters.
#
# Steps:
# 1. Load loop-gene mapping and CLICK cluster assignments.
# 2. Identify genes associated with each loop cluster.
# 3. Perform GO enrichment (Biological Process ontology) for each gene set.
# 4. Generate barplot and dotplot visualizations for enriched terms.
# 5. Save per-cluster enrichment tables and plots.
# 6. Summarize enrichment across all clusters and visualize the distribution.
#
# Inputs:
# - CLICK clusters: click_clusters.txt
# - TSS-loop interaction table: tssLoops_values.tsv
# - Output directory for enrichment results
#
# Outputs:
# - GO_BP_Cluster_*.tsv: Enrichment tables per cluster
# - GO_BP_barplot_*.png / GO_BP_dotplot_*.png: Visualization per cluster
# - GO_enrichment.rds: Combined enrichGO objects
# - all_enrichments.tsv: Summary table with fold change
# - p_barplot_sig_annots.png: Summary barplot of clusters with significant terms
#
# Author: Gony Shanel
# Date: March 2025
# ------------------------------------------------------------------------------

## Set the scripts directory as working directory.
## This script is designed with relative paths for the specific data structure
## with which it was uploaded to https://github.com/ElkonLab/micro-C.

# Functions ---------------------------------------------------------------

getRatio <- function(ratio){
  splt <- strsplit(ratio, "/")[[1]]
  r <- as.numeric(splt[1])/as.numeric(splt[2])
}

# Arguments ---------------------------------------------------------------

# ## Example
args <- c("../../data/loops/click_clusters.txt",
          "../../data/loops/tssLoops_values.tsv",
          ".") # For NAR publication example, set to script's directory

# # args <- commandArgs(trailingOnly = T)
# if (length(args) != 3){
#   stop("Usage: clusters_ORA.r {click_clusters} {tss_loops} {outDir}")
# }

click_path <- R.utils::getAbsolutePath(args[1])
tssLoops_path <- R.utils::getAbsolutePath(args[2])
outDir <- R.utils::getAbsolutePath(args[3])

# Commands ----------------------------------------------------------------
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(parallel)

dir.create("results", showWarnings = F)
dir.create("plots", showWarnings = F)

## Read files
tssLoops_raw = fread(tssLoops_path)
click_raw = fread(click_path)

## Shape data
tssLoops = tssLoops_raw[, !grepl("nutlin", colnames(tssLoops_raw)), with = F]
colnames(tssLoops) <- sub("-control", "", colnames(tssLoops))
click = click_raw[grepl("Cluster", click_raw$V2),]
cluster_names = unique(click$V2)
bg <- unlist(unique(tssLoops$gene))
Nbg <- length(bg) ## number of loops in background

## ORA analysis
clst = cluster_names[6]
## GO - Biological proccesses
BP_enrichments <- mclapply(cluster_names, function(clst){

  loopSet = unlist(click[click$V2 == clst, 1])  ## Cluster's loops
  geneSet = unique(unlist(tssLoops[tssLoops$peak_name %in% loopSet, "gene"]))  ## Loop cluster's associated gene set)

  ## ORA
  ego = enrichGO(gene = geneSet,
                 OrgDb = org.Hs.eg.db,
                 keyType = "ENSEMBL",
                 ont = "BP",
                 universe = bg)
  # BP_enrichments = c(BP_enrichments, ego)
  
  if (nrow(ego) != 0){
    ## Plot
    p = barplot(ego, 
                showCategory = 12,
                font.size = 7) + 
      theme(text = element_text(size = 8))
    p
    ggsave(paste0("plots/GO_BP_barplot_", clst, ".png"), width = 3.5, height = 3.5)
    
    dotplot(ego, 
            font.size = 7) +
      labs(x = "gene ratio") + 
      theme(text = element_text(size = 8))
    ggsave(paste0("plots/GO_BP_dotplot_", clst, ".png"), width = 3.5, height = 3.5)
    
    ## Save
    write.table(as.data.frame(ego),
                paste0("results/GO_BP_", clst, ".tsv"),
                append = F,
                quote = F,
                sep = '\t', 
                row.names = F,
                col.names = T)
  }
  ego
},
mc.cores = 20)
names(BP_enrichments) <- cluster_names
saveRDS(BP_enrichments, "GO_enrichment.rds")

## Summarize
all_enrichments <- as.data.frame(summary(merge_result(BP_enrichments)))

## Calculate fold change
geneRatio <- sapply(all_enrichments$GeneRatio, getRatio)
bgRatio   <- sapply(all_enrichments$BgRatio, getRatio)
all_enrichments$fold_change <- geneRatio/bgRatio 

## Save summarized results
write.table(all_enrichments,
            "all_enrichments.tsv",
            append = F, 
            quote = F,
            sep = '\t',
            row.names = F,
            col.names = T)

## 
culster_no <- sub("Cluster_", "", all_enrichments$Cluster)
all_enrichments$cluster_no <- factor(culster_no, levels = unique(culster_no), ordered = T)
p_col <- ggplot(all_enrichments, aes(x = cluster_no)) +
  geom_bar() +
  labs(title = "Significant enriched GO terms",
       subtitle = "(p-adjusted < 0.05)",
       x = "Cluster number",
       caption = paste0("Not showing clusters with zero results. Nbg = ", Nbg, ".")) +
  theme_bw() +
  theme(text = element_text(size = 8),
        axis.text.x = element_text(angle = 45),
        plot.caption = element_text(size = 6, hjust = 0))
  
ggsave("p_barplot_sig_annots.png", p_col, width = 3.5, height = 3)

