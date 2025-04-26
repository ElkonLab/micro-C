# ------------------------------------------------------------------------------
# Explore the Association Between TAD Organization and Gene Expression Patterns
# ------------------------------------------------------------------------------
# This script investigates how topologically associating domains (TADs) relate
# to gene expression regulation across multiple human cell lines. The analysis
# tests whether genes co-residing within the same TAD ("TAD neighbors") tend to:
# 
#   1. Be more differentially expressed across cell lines.
#   2. Show higher co-expression compared to non-neighboring genes at similar distances.
#   3. Be surrounded by other differentially expressed genes.
#
# Major steps:
# - Annotate TAD intervals with protein-coding genes (TSS mapping).
# - Identify genes co-residing within the same TAD.
# - Select the top differentially expressed genes between pairs of cell lines.
# - Match each TAD neighbor with a distance-matched non-neighbor control gene.
# - Compare gene expression fold change and correlation between neighbor vs. matched controls.
#
# INPUT:
# - TAD intervals and pooled TAD boundaries (from Arrowhead at 10kb resolution).
# - Gene expression data (raw + DESeq2 log2 fold change).
# - TSS annotations (Ensembl, protein-coding, hg38).
# - Gene-pair distances (precomputed genomic distances between genes).
#
# OUTPUT:
# - Boxplots and heatmaps of differential expression patterns.
# - Matched gene pair comparisons of TAD neighbors vs. non-neighbors.
# - Cumulative correlation plots across all cell lines.
#
# AUTHOR: Gony Shanel
# DATE: March 2025
# ------------------------------------------------------------------------------

## Set the scripts directory as working directory.
## This script is designed with relative paths for the specific data structure
## with which it was uploaded to https://github.com/ElkonLab/micro-C.

## This script requires bedtools to be installed (https://bedtools.readthedocs.io/en/latest/index.html)
path.bedtools <- "/path/to/bedtools"

# Functions ---------------------------------------------------------------

find_gene_neighbors <- function(GENE, TAD_DATA){
  tad_no <- TAD_DATA[gene == GENE, tad_id]
  unique(unlist(sapply(tad_no, function(n){
    TAD_DATA[tad_id == n & gene != GENE, gene]
  })))
}

pull_gene_pair_distance <- function(GENE1, GENE2, DT_DISTANCE){
  require(data.table)
  ROW <- DT_DISTANCE[gene1 == GENE1 & gene2 == GENE2]
  if (dim(ROW)[1] == 0){
    ROW <- DT_DISTANCE[gene1 == GENE2 & gene2 == GENE1]
  }
  if (dim(ROW)[1] == 0){
    return(NA)
  } else {
    return(ROW[[4]])
  }
}

match_control_genes <- function(GENE_NEIGHBORS, EXCLUDE=c()){

  as.data.table(do.call(rbind, mclapply(1:nrow(GENE_NEIGHBORS), function(i){
    
    ANCHOR_GENE <- unlist(GENE_NEIGHBORS[[i, 1]])
    NEIGH_GENES <- unique(unlist(GENE_NEIGHBORS[i,2])) # Anchor's TAD neighbors
    
    C_PAIRs <- GENE_PAIRS_DISTANCE[                   # Keep pairs of genes anchored to core gene that are not TAD neighbors
      distance > 10000 &
      (gene1 == ANCHOR_GENE | gene2 == ANCHOR_GENE) & # pair is anchored to ANCHOR gene
        (! gene1 %in% NEIGH_GENES) &        # Pair does not contain a neighbor
        (! gene2 %in% NEIGH_GENES) &         # Pair does not contain a neighbor
        (! gene1 %in% EXCLUDE | ! gene2 %in% EXCLUDE)]
    
    CONTROL_PAIRS <- data.table("anchor" = "dummy", "neighbor" = "dummy", "match" = "dummy")

    for (NEIGH in NEIGH_GENES){
      
      # Pull distance between anchor gene and the neighbor gene to be matched
      DIS <- pull_gene_pair_distance(ANCHOR_GENE, NEIGH, GENE_PAIRS_DISTANCE)
      # Look at anchored pairs by difference between their dostnace and the distance to match
      C_PAIRs$dis_diff <- C_PAIRs$distance-DIS
      C_PAIRs <- C_PAIRs[dis_diff < 1] %>% arrange(desc(dis_diff))
      if (nrow(C_PAIRs) == 0){next}
      # take closest match
      I=1
      C_PAIR <- C_PAIRs[I, c(gene1, gene2)]
      MATCH <- C_PAIR[!C_PAIR == ANCHOR_GENE]
      # Ensure pair was not picked for another neighbor
      
      while (MATCH %in% CONTROL_PAIRS$match){
        I=I+1
        if (I > nrow(C_PAIRs)){break}
        C_PAIR <- C_PAIRs[I, c(gene1, gene2)]
        MATCH <- C_PAIR[!C_PAIR == ANCHOR_GENE]
      }
      # Aggregate result
      if (!MATCH %in% CONTROL_PAIRS$match){
        CONTROL_PAIRS <- rbind(CONTROL_PAIRS, 
                               data.table("anchor" = ANCHOR_GENE,
                                          "neighbor" = NEIGH,
                                          "match" = MATCH))
      }
    }
    CONTROL_PAIRS[-1]
  }, 
  mc.cores = 45)))
}

pull_gene_FC <-  function(GENE, deseq_data){
  deseq_data[gene == GENE, log2FoldChange]
}

F_MATCHED_NEIGHBORS <- function(){
  lapply(CELL_PAIRS, function(PAIR){
    
    CELL1 <- SPLIT_CELLS[[PAIR]][1]
    CELL2 <- SPLIT_CELLS[[PAIR]][2]
    TOP_FC_GENES_1 <- TOP_FC_GENES[[PAIR]][[1]]
    TOP_FC_GENES_2 <- TOP_FC_GENES[[PAIR]][[2]]
    NEIGHBOUR_GENES_1 <- NEIGHBOUR_GENES[[PAIR]][cell == CELL1, gene]
    NEIGHBOUR_GENES_2 <- NEIGHBOUR_GENES[[PAIR]][cell == CELL2, gene]
    
    
    ## Find non-neighbor genes with distance similar to neighbors
    
    GENE_NEIGHBORS_1 <- data.table("gene" = TOP_FC_GENES_1$gene,
                                   "neighbor" = lapply(TOP_FC_GENES_1$gene,
                                                       find_gene_neighbors,
                                                       TAD_DATA = TAD_GENES[[CELL1]]))
    GENE_NEIGHBORS_2 <- data.table("gene" = TOP_FC_GENES_2$gene,
                                   "neighbor" = lapply(TOP_FC_GENES_2$gene,
                                                       find_gene_neighbors,
                                                       TAD_DATA = TAD_GENES[[CELL2]]))
    
    
    CONTROL_PAIRS_1 <- match_control_genes(GENE_NEIGHBORS = GENE_NEIGHBORS_1, 
                                           EXCLUDE = c(NEIGHBOUR_GENES_1, TOP_FC_GENES_1))
    CONTROL_PAIRS_2 <- match_control_genes(GENE_NEIGHBORS = GENE_NEIGHBORS_2, 
                                           EXCLUDE = c(NEIGHBOUR_GENES_2, TOP_FC_GENES_2))
    
    list(CONTROL_PAIRS_1, CONTROL_PAIRS_2)
  }
  )
}

find_genes_in_range <- function(GENE, CELL){
  GENE_PAIRS_IN_RANGE <-  GENE_PAIRS_DISTANCE[(gene1 == GENE |
                                                 gene2 == GENE) &
                                                distance > RANGE[1] &
                                                distance < RANGE[2]]
  if (nrow(GENE_PAIRS_IN_RANGE) > 0){
    NEIGHBORS <- find_gene_neighbors(GENE, TAD_GENES[[CELL]])
    GENES_IN_RANGE <- c(GENE_PAIRS_IN_RANGE$gene1, GENE_PAIRS_IN_RANGE$gene2)
    GENES_IN_RANGE <- GENES_IN_RANGE[GENES_IN_RANGE != GENE]
    GENES_IN_RANGE <- data.frame("cell" = CELL, 
                                 "gene" = GENES_IN_RANGE, 
                                 "group" = ifelse(GENES_IN_RANGE %in% NEIGHBORS,
                                                  "neighbor",
                                                  "non-neighbor"))
    GENES_IN_RANGE$log2FoldChange <- sapply(GENES_IN_RANGE$gene,
                                            pull_gene_FC,
                                            deseq_data = FC_DATA)
    GENES_IN_RANGE
  }
}

# for (GENE in TOP_FC_GENES_1[[1]]){
#   print(GENE)
#   print(find_genes_in_range(GENE, CELL1))
# }

# Arguments ---------------------------------------------------------------

EXP_PATH <- "../data/gene_expression_CNVcorrected.tsv"
TSS_PATH <- "../data/ensembl_mart_hg38_1st_TSS_protein_coding.txt"
TAD_FILES <- grep("control", list.files("../data/TAD_interval_files/", full.names = T), value = T)
TADS <- lapply(TAD_FILES, fread, col.names = c("chr", "start", "end"))
POOLED_TADS_PATH <- "../data/pooled_TAD_intervals.bed"
CELLS <- readLines("../data/cell_lines.txt")
CELL_PAIRS <- readLines("../data/cell_lines_pairs.txt")
DESEQ_DIR <- "../data/deseq2/paires/"
GENE_PAIR_DISTANCE <- "../data/distance_between_gene_pairs.tsv"
N_DIFF_GENES <- 200

# For testing
CELL <- CELLS[1]
PAIR <- CELL_PAIRS[[1]]

# Commands ----------------------------------------------------------------

## Set environment ----

library(data.table)
library(dplyr)
options(bedtools.path = path.bedtools)
library(bedtoolsr)
library(ggplot2)
library(ggsignif)
library(parallel)
library(ggpubr)

dir.create("plots")
dir.create("results")

## Prepare data ####
## Load gene expression data
DESEQ_FILES <- list.files(DESEQ_DIR,
                          "deseq2.res.Cells.txt",
                          full.names = T,
                          recursive = T)
DESEQ_DATA <- lapply(DESEQ_FILES, fread) ; names(DESEQ_DATA) <- CELL_PAIRS
EXP_RAW <- fread(EXP_PATH)
EXP_RAW <- EXP_RAW[, !grepl("nutlin", colnames(EXP_RAW)), with = F]
colnames(EXP_RAW) <- sub("-control", "", colnames(EXP_RAW))
#
TSS_RAW <- fread(TSS_PATH)

names(TADS) <- CELLS
TADS <- lapply(TADS, function(DT){
  DT <- DT[(end-start) < 1000000] ## Filter out TADs longer than 2 Mbp
  DT %>%
    arrange(as.numeric(sub("chr", "", chr)), start, end)
  cbind(DT, "tad_id" = 1:nrow(DT))
}) # add TAD IDs
POOLED_TADS <- fread(POOLED_TADS_PATH)

GENE_PAIRS_DISTANCE <- fread(GENE_PAIR_DISTANCE)

## Describe TADs ----
ALL_TADS <- as.data.table(do.call(rbind, mclapply(CELLS, function(CELL){
  cbind(TADS[[CELL]], cell = CELL)
})))

### Number of TADs per cell ----
median(table(ALL_TADS$cell)); sd(table(ALL_TADS$cell))
ggplot(ALL_TADS, aes(cell)) +
  geom_bar() + theme_bw() +
  coord_flip() +
  ylab("TADs [N]") +
  theme(text = element_text(size = 8),
        axis.title.y = element_blank())
ggsave("plots/p.barplot.TADs_N_per_cell.png",
       height = 2, width = 2)

### TAD length ----
ALL_TADS$length <- ALL_TADS[, end-start]
median(unique(ALL_TADS[, c(1:3, 6)])$length); sd(unique(ALL_TADS[, c(1:3, 6)])$length)

ggplot(ALL_TADS, aes(cell, length/1000)) +
  geom_boxplot(outlier.shape = NA) + theme_bw() +
  coord_flip(ylim = c(0, quantile(ALL_TADS$length, 0.99)/1000)) +
  ylab("TAD length [kbp]") +
  theme(text = element_text(size = 8),
        axis.title.y = element_blank())
ggsave("plots/p.boxplot.TAD_length.png",
       height = 2, width = 2)
  

## Associate TADs with genes ----

TAD_GENES <- lapply(CELLS, function(CELL){
  
  TAD <- TADS[[CELL]]
  TAD_GENES <- as.data.table(bt.intersect(a = TAD,
                                          b = TSS_RAW,
                                          loj = T))
  colnames(TAD_GENES) <- c(colnames(TAD), colnames(TSS_RAW))
  TAD_GENES[gene == ".", "gene"] = NA
  TAD_GENES[sym == ".", "sym"] = NA
  TAD_GENES[, -c(5:7)]
})
names(TAD_GENES) <- CELLS

## Plot the number of genes per TAD
GENES_PER_TAD <- lapply(CELLS, function(CELL){
  TAD <- TADS[[CELL]]
  GENES_PER_TAD <- as.data.table(bt.intersect(a = TAD,
                                          b = TSS_RAW,
                                          c = T))
  colnames(GENES_PER_TAD) <- c(colnames(TAD), "N_genes")
  GENES_PER_TAD
})
names(GENES_PER_TAD) <- CELLS

data_to_plot <- as.data.table(do.call(rbind, lapply(CELLS, function(CELL){
  cbind(GENES_PER_TAD[[CELL]], "cell" = CELL)
})))

mean(data_to_plot$N_genes)
sd(data_to_plot$N_genes)

data_to_plot$bin <- data_to_plot$N_genes
data_to_plot$bin[data_to_plot$bin>4] <- 5

ggplot(data_to_plot, aes(bin)) +
  geom_histogram(bins = 6) + theme_bw() +
  scale_x_continuous(breaks = 0:5,
                     labels=c(as.character(0:4), 
                              expression(phantom(x) >=5))) +
  labs(x = "genes per TAD [N]",
       y = "TAD count [N]") +
  facet_wrap(~ cell, nrow = 2) +
  theme(text = element_text(size = 8))

ggsave("plots/p.boxplot_promoters_per_TAD.png",
       height = 2, width = 5)

## Associate gene id with TADS that harbor multiple genes (neighborhood)
NEIGHBOURHOODS <- lapply(CELLS, function(CELL){
  right_join(TAD_GENES[[CELL]], GENES_PER_TAD[[CELL]][N_genes > 1])
})
names(NEIGHBOURHOODS) <- CELLS

## Genes that co-reside with highly deferential genes are also differential ####

## Preapre cell names for each pair
SPLIT_CELLS <- mclapply(CELL_PAIRS, function(PAIR){
  CELL1 <- strsplit(PAIR, "-")[[1]][[1]]
  CELL2 <- strsplit(PAIR, "-")[[1]][[2]]
  c(CELL1, CELL2)
})
names(SPLIT_CELLS) <- CELL_PAIRS

## Shape DESEQ (comparison of expression data) data for each pair
DESEQ <- mclapply(CELL_PAIRS, function(PAIR){
  
  CELL1 <- SPLIT_CELLS[[PAIR]][1]
  CELL2 <- SPLIT_CELLS[[PAIR]][2]
  
  DESEQ <- DESEQ_DATA[[PAIR]]
  names(DESEQ)[1] <- "gene"
  DESEQ$gene <- sub("\\..*", "", DESEQ$gene)
  DESEQ
  
})
names(DESEQ) <- CELL_PAIRS

## Find most differential genes that have TAD neighbors
TOP_FC_GENES <- mclapply(CELL_PAIRS, function(PAIR){
  
  CELL1 <- SPLIT_CELLS[[PAIR]][1]
  CELL2 <- SPLIT_CELLS[[PAIR]][2]
  DESEQ <- DESEQ[[PAIR]]
  
  ## find extreme deferential genes
  
  FILTERED_DESEQ_1 <- DESEQ[gene %in% NEIGHBOURHOODS[[CELL1]]$gene] %>% arrange(log2FoldChange)
  FILTERED_DESEQ_2 <- DESEQ[gene %in% NEIGHBOURHOODS[[CELL2]]$gene] %>% arrange(desc(log2FoldChange))
  
  TOP_FC_GENES_1 <- FILTERED_DESEQ_1[1:N_DIFF_GENES]
  TOP_FC_GENES_2 <- FILTERED_DESEQ_2[1:N_DIFF_GENES]
  
  list(TOP_FC_GENES_1, TOP_FC_GENES_2)
},
mc.cores = 15)
names(TOP_FC_GENES) <- CELL_PAIRS

## Find the neighbors of differential genes
# PAIR <- CELL_PAIRS[[2]]
# NEIGHBOUR_GENES <- lapply(CELL_PAIRS, function(PAIR){
# 
#   CELL1 <- SPLIT_CELLS[[PAIR]][1]
#   CELL2 <- SPLIT_CELLS[[PAIR]][2]
#   D <- DESEQ[[PAIR]]
#   TOP_FC_GENES_1 <- TOP_FC_GENES[[PAIR]][[1]]
#   TOP_FC_GENES_2 <- TOP_FC_GENES[[PAIR]][[2]]
# 
#   ## Find diff genes TAD neighbors
#   NEIGHBOUR_GENES <-
#     rbind(data.table("gene" = unlist(sapply(TOP_FC_GENES_1$gene,
#                                             find_gene_neighbors,
#                                             TAD_DATA = TAD_GENES[[CELL1]])),
#                      "cell" = CELL1),
#           data.table("gene" = unlist(sapply(TOP_FC_GENES_2$gene,
#                                             find_gene_neighbors,
#                                             TAD_DATA = TAD_GENES[[CELL2]])),
#                      "cell" = CELL2))
#   unique(left_join(NEIGHBOUR_GENES, D))
# 
# })
# names(NEIGHBOUR_GENES) <- CELL_PAIRS
# 
# saveRDS(NEIGHBOUR_GENES, "results/NEIGHBOUR_GENES.rds")
NEIGHBOUR_GENES <- readRDS("results/NEIGHBOUR_GENES.rds")

## Compare expression fold change between beighbors
COMPARISONS <- mclapply(CELL_PAIRS, function(PAIR){
  
  CELL1 <- SPLIT_CELLS[[PAIR]][1]
  CELL2 <- SPLIT_CELLS[[PAIR]][2]
  neighbors <- NEIGHBOUR_GENES[[PAIR]]
  
  C <- t.test(neighbors[cell == CELL1, log2FoldChange],
              neighbors[cell == CELL2, log2FoldChange], 
              paired = F, alternative = "less")
  
  c("cell1" = CELL1, "cell2" = CELL2, unlist(C))

})
names(COMPARISONS) <- CELL_PAIRS

## Plot comprisons in boxplots
dir.create("plots/pairwise_compare_diff_genes_neighbours_FC/")
DATA_TO_PLOT <- mclapply(CELL_PAIRS, function(PAIR){
  
  ## Shape
  CELL1 <- SPLIT_CELLS[[PAIR]][1]
  CELL2 <- SPLIT_CELLS[[PAIR]][2]
  TOP_FC_GENES_1 <- TOP_FC_GENES[[PAIR]][[1]]
  TOP_FC_GENES_2 <- TOP_FC_GENES[[PAIR]][[2]]  
  C <- COMPARISONS[[PAIR]]
  NEIGH_GENES <- NEIGHBOUR_GENES[[PAIR]]
  
  ## Add core genes
  NEIGH_GENES$type <- "neighbor"
  NEIGH_GENES <- rbind(NEIGH_GENES,
                           rbind(cbind(TOP_FC_GENES_1,
                                       "cell" = CELL1,
                                       "type" = "core"),
                                 cbind(TOP_FC_GENES_2, 
                                       "cell" = CELL2,
                                       "type" = "core")))
  LEVELS <- c(paste(c(CELL1, CELL2), "core"),
              paste(c(CELL1, CELL2), "neighbor"))
  NEIGH_GENES$group <- factor(paste(NEIGH_GENES$cell, NEIGH_GENES$type),
                              levels = LEVELS,
                              ordered = T)
  
  ## Plot
  xlabs <- paste0(c(paste0(CELL1, "\ncore"),
                    paste0(CELL2, "\ncore"),
                    paste0(CELL1, "\nneighbors"),
                    paste0(CELL2, "\nneighbors")),
                  "\n(N=", table(NEIGH_GENES$group), ")")
  
  ggplot(NEIGH_GENES, aes(x = group, y = log2FoldChange, fill = group)) +
    geom_boxplot(outlier.shape = NA) + theme_bw(base_size = 8) +
    scale_y_continuous(limits = quantile(NEIGH_GENES$log2FoldChange,
                                         c(0.1, 0.9), na.rm = T)) +
    geom_signif(comparisons = list(c(LEVELS[3], LEVELS[4])),
                annotation = formatC(as.numeric(C["p.value"]) ,3),
                textsize = 2.5, 
                extend_line = -0.02,
                tip_length = 0,
                y_position = quantile(NEIGH_GENES$log2FoldChange,
                                      0.9, na.rm = T)-3) +
    scale_fill_manual(values = c("grey70", "grey50", "coral", "cyan3")) +
    scale_x_discrete(labels = xlabs) +
    ylab(paste0("log2 fold change (",
                CELL2, "/", CELL1, ")")) +
    theme(axis.title.x = element_blank(),
          legend.position = "none")
  
  ggsave(paste0("plots/pairwise_compare_diff_genes_neighbours_FC/", PAIR, ".png"),
         height = 2.5, width = 2.5)
  
  NEIGH_GENES
},
mc.cores = 45)
names(DATA_TO_PLOT) <- CELL_PAIRS


## Summarize comparisons
SUM_COMPARISONS <- as.data.table(do.call(rbind, COMPARISONS))
SUM_COMPARISONS$p.value <- as.numeric(SUM_COMPARISONS$p.value)

ggplot(SUM_COMPARISONS, aes(cell1, cell2, fill = -log10(p.value))) +
  geom_tile() + theme_minimal() + 
  scale_fill_distiller(palette = "Reds", 
                       direction = 1, 
                       na.value = "black") +
  labs(fill = expression("-log"[10]*P)) +
  theme(text = element_text(size = 8),
        legend.position = c(0.9, 0.7),
        panel.grid.major = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45,
                                   vjust = 0.8)) +
  scale_y_discrete(limits = rev(unique(SUM_COMPARISONS$cell2))) +
  geom_text(aes(label = round(-log10(p.value))), size = 2.5)

ggsave("plots/pairwise_compare_diff_genes_neighbors_FC.png",
       height = 3, width = 3)


####
## Find for extreme diff genes a counterpart gene with a distance similar 
## to the distance of it's neighbor that is not in the same TAD and have 
## distance equal or smaller, and compare fold change. 


# MATCHED_NEIGHBORS <- F_MATCHED_NEIGHBORS()
# names(MATCHED_NEIGHBORS) <- CELL_PAIRS
# 
# saveRDS(MATCHED_NEIGHBORS, "results/MATCHED_NEIGHBORS.rds")
MATCHED_NEIGHBORS <- readRDS("results/MATCHED_NEIGHBORS.rds")

## Plot distances to neighbors compared to matched controls #### 
dir.create("plots/neighbours_compared_to_close_non-neighbour/distance_neighbor_match/",
           recursive = T, showWarnings = F)
mclapply(CELL_PAIRS, function(PAIR){
  
  CELL1 <- SPLIT_CELLS[[PAIR]][1]
  CELL2 <- SPLIT_CELLS[[PAIR]][2]
  MATCHED_NEIGHBORS_1 <- MATCHED_NEIGHBORS[[PAIR]][[1]]
  MATCHED_NEIGHBORS_2 <- MATCHED_NEIGHBORS[[PAIR]][[2]]
  
  DISTANCES_1 <- 
    rbind(
      data.table(
        "distance" = sapply(1:nrow(MATCHED_NEIGHBORS_1), function(i){
          ANCHOR <- MATCHED_NEIGHBORS_1[i, anchor]
          NEIGHBOR <- MATCHED_NEIGHBORS_1[i, neighbor]
          abs(TSS_RAW[gene == ANCHOR, start] - TSS_RAW[gene == NEIGHBOR, start])
        }),
        "group" = "neighbor"),
      data.table(
        "distance" = sapply(1:nrow(MATCHED_NEIGHBORS_1), function(i){
          ANCHOR <- MATCHED_NEIGHBORS_1[i, anchor]
          MATCH <- MATCHED_NEIGHBORS_1[i, match]
          abs(TSS_RAW[gene == ANCHOR, start] - TSS_RAW[gene == MATCH, start])
        }),
        "group" = "match")
      )
  
  DISTANCES_2 <- 
    rbind(
      data.table(
        "distance" = sapply(1:nrow(MATCHED_NEIGHBORS_2), function(i){
          ANCHOR <- MATCHED_NEIGHBORS_2[i, anchor]
          NEIGHBOR <- MATCHED_NEIGHBORS_2[i, neighbor]
          abs(TSS_RAW[gene == ANCHOR, start] - TSS_RAW[gene == NEIGHBOR, start])
        }),
        "group" = "neighbor"),
      data.table(
        "distance" = sapply(1:nrow(MATCHED_NEIGHBORS_2), function(i){
          ANCHOR <- MATCHED_NEIGHBORS_2[i, anchor]
          MATCH <- MATCHED_NEIGHBORS_2[i, match]
          abs(TSS_RAW[gene == ANCHOR, start] - TSS_RAW[gene == MATCH, start])
        }),
        "group" = "match")
      )
  
  # Plot
  DISTANCES_TO_PLOT <- rbind(  cbind(DISTANCES_1, "cell" = CELL1),
                               cbind(DISTANCES_2, "cell" = CELL2))
  DISTANCES_TO_PLOT$distance <- as.numeric(DISTANCES_TO_PLOT$distance)
  DISTANCES_TO_PLOT$group <- factor(DISTANCES_TO_PLOT$group,
                                    levels = c("neighbor", "match"),
                                    ordered = T)
  
  ggplot(DISTANCES_TO_PLOT, aes(distance, color = group)) +
    geom_density() + theme_bw() +
    scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
    scale_color_brewer(palette = "Set1") +
    theme(text = element_text(size = 8)) +
    facet_wrap(~ cell)
  ggsave(paste0("plots/neighbours_compared_to_close_non-neighbour/distance_neighbor_match/",
                PAIR, ".png"),
         height = 2, width = 5)
},
mc.cores = 45)

## Compare neighbors to matched controls #### 
dir.create("plots/neighbours_compared_to_close_non-neighbour/compare_neighbors_to_match/",
           recursive = T, showWarnings = F)
# NEIGHBOUR_COMPARED_TO_MATCHED <- as.data.table(do.call(rbind, lapply(CELL_PAIRS, function(PAIR){
# 
#   CELL1 <- SPLIT_CELLS[[PAIR]][1]
#   CELL2 <- SPLIT_CELLS[[PAIR]][2]
#   FC_DATA <- DESEQ[[PAIR]]
#   MATCHED_NEIGHBORS_1 <- MATCHED_NEIGHBORS[[PAIR]][[1]]
#   MATCHED_NEIGHBORS_2 <- MATCHED_NEIGHBORS[[PAIR]][[2]]
# 
#   ## Shape data
#   MATCHES_1 <- left_join(data.table("gene" = MATCHED_NEIGHBORS_1$match,
#                                     "cell" = CELL1,
#                                     "type" = "match",
#                                     "group" = paste(CELL1, "match")),
#                          FC_DATA)
#   MATCHES_2 <- left_join(data.table("gene" = MATCHED_NEIGHBORS_2$match,
#                                     "cell" = CELL2,
#                                     "type" = "match",
#                                     "group" = paste(CELL2, "match")),
#                          FC_DATA)
# 
#   DATA <- rbind(DATA_TO_PLOT[[PAIR]], # neighbors data
#                 MATCHES_1,
#                 MATCHES_2)
#   LEVELS <- levels(DATA$group)
# 
#   ## Statistical test
#   C_NEIGHBOURS <- compare_means(log2FoldChange ~ group,
#                                 DATA[group %in% c(LEVELS[3], LEVELS[4])])
#   C_MATCHES <- compare_means(log2FoldChange ~ group,
#                              DATA[group %in% c(LEVELS[5], LEVELS[6])])
# 
#   ## Plot
#   CC <- c(CELL1, CELL2)
#   XLABS <- paste0(c(paste(CC, "core", sep = '\n'),
#                     paste(CC, "neighbor", sep = '\n'),
#                           paste(CC, "match", sep = '\n')),
#                   "\n(N=", table(DATA$group), ")")
#   Q <- quantile(DATA$log2FoldChange,
#                 na.rm = T)
# 
#   ggplot(DATA, aes(group, log2FoldChange, fill = group)) +
#     geom_boxplot(outlier.shape = NA) + theme_bw() +
#     geom_signif(comparisons = list(c(LEVELS[3], LEVELS[4]),
#                                    c(LEVELS[5], LEVELS[6])),
#                 annotations = c(C_NEIGHBOURS$p.adj, C_MATCHES$p.adj),
#                 y_position = Q[5]*.75,
#                 tip_length = 0, test = "t.test",
#                 textsize = 2.5) +
#     scale_x_discrete(labels = XLABS) +
#     ylab(paste0("log2 mRNA fold\nchange [", CELL1, "/", CELL2, "]")) +
#     scale_fill_manual(values = c("grey70", "grey50", "pink4", "cyan4", "pink", "cyan")) +
#     # scale_fill_manual(values = c("grey70", "grey50", "red", "red", "blue", "blue")) +
#     theme(text = element_text(size = 8),
#           axis.title.x = element_blank(),
#           legend.position = "none")
#   ggsave(paste0("plots/neighbours_compared_to_close_non-neighbour/compare_neighbors_to_match/",
#                 PAIR, ".png"),
#          height = 2, width = 3)
# 
#   ## print statistics
#   rbind(c("pair" = PAIR, "type" = "neighbor", C_NEIGHBOURS),
#         c("pair" = PAIR, "type" = "match", C_MATCHES))
# 
# })))
# saveRDS(NEIGHBOUR_COMPARED_TO_MATCHED, 
#         "results/NEIGHBOUR_COMPARED_TO_MATCHED.rds")
NEIGHBOUR_COMPARED_TO_MATCHED <- readRDS("results/NEIGHBOUR_COMPARED_TO_MATCHED.rds")

## Within range 2 ----

## Find 75 precentile of TAD's length

POOLED_TADS$length <- POOLED_TADS[, V3-V2]
CUTOFF <- quantile(POOLED_TADS$length, 0.75)
MAX <- max(POOLED_TADS$length)
##
RANGE <- c(350000, 1000000)
OUTDIR <- paste0("plots/neighbours_compared_to_close_non-neighbour/compare_neighbors_to_match_range_", RANGE[1], "_", RANGE[2], "/")
dir.create(OUTDIR,
           recursive = T, showWarnings = F)
COMPARE_GENES_IN_RANGE <- as.data.table(do.call(rbind, lapply(CELL_PAIRS, function(PAIR){
  
  CELL1 <- SPLIT_CELLS[[PAIR]][1]
  CELL2 <- SPLIT_CELLS[[PAIR]][2]
  FC_DATA <- DESEQ[[PAIR]]
  MATCHED_NEIGHBORS_1 <- MATCHED_NEIGHBORS[[PAIR]][[1]]
  MATCHED_NEIGHBORS_2 <- MATCHED_NEIGHBORS[[PAIR]][[2]]
  
  ## Filter by range
  MATCHED_NEIGHBORS_1$dis <- mclapply(1:nrow(MATCHED_NEIGHBORS_1), function(i){
    pull_gene_pair_distance(MATCHED_NEIGHBORS_1[i,1],
                            MATCHED_NEIGHBORS_1[i,2],
                            GENE_PAIRS_DISTANCE)
  }, mc.cores = 45)
  MATCHED_NEIGHBORS_1 <- MATCHED_NEIGHBORS_1[dis > RANGE[1] & dis < RANGE[2]]
  
  MATCHED_NEIGHBORS_2$dis <- mclapply(1:nrow(MATCHED_NEIGHBORS_2), function(i){
    pull_gene_pair_distance(MATCHED_NEIGHBORS_1[i,1],
                            MATCHED_NEIGHBORS_1[i,2],
                            GENE_PAIRS_DISTANCE)
  }, mc.cores = 45)
  MATCHED_NEIGHBORS_2 <- MATCHED_NEIGHBORS_2[dis > RANGE[1] & dis < RANGE[2]]
  
  ## Shape data
  MATCHES_1 <- left_join(data.table("gene" = MATCHED_NEIGHBORS_1$match,
                                    "cell" = CELL1,
                                    "type" = "match",
                                    "group" = paste(CELL1, "match")),
                         FC_DATA)
  MATCHES_2 <- left_join(data.table("gene" = MATCHED_NEIGHBORS_2$match,
                                    "cell" = CELL2,
                                    "type" = "match",
                                    "group" = paste(CELL2, "match")),
                         FC_DATA)
  
  DATA <- rbind(DATA_TO_PLOT[[PAIR]], # neighbors data
                MATCHES_1,
                MATCHES_2)
  LEVELS <- levels(DATA$group)
  
  ## Statistical test
  C_NEIGHBORS <- compare_means(log2FoldChange ~ group,
                                DATA[group %in% c(LEVELS[3], LEVELS[4])])
  C_MATCHES <- compare_means(log2FoldChange ~ group,
                             DATA[group %in% c(LEVELS[5], LEVELS[6])])
  
  ## Plot
  CC <- c(CELL1, CELL2)
  XLABS <- paste0(c(paste(CC, "core", sep = '\n'),
                    paste(CC, "neighbor", sep = '\n'),
                    paste(CC, "match", sep = '\n')),
                  "\n(N=", table(DATA$group), ")")
  Q <- quantile(DATA$log2FoldChange,
                na.rm = T)
  
  ggplot(DATA, aes(group, log2FoldChange, fill = group)) +
    geom_boxplot(outlier.shape = NA) + theme_bw() +
    geom_signif(comparisons = list(c(LEVELS[3], LEVELS[4]),
                                   c(LEVELS[5], LEVELS[6])),
                annotations = c(C_NEIGHBORS$p.adj, C_MATCHES$p.adj),
                y_position = Q[5]*.75,
                tip_length = 0, test = "t.test",
                textsize = 2.5) +
    scale_x_discrete(labels = XLABS) +
    ylab(paste0("log2 mRNA fold\nchange [", CELL1, "/", CELL2, "]")) +
    scale_fill_manual(values = c("grey70", "grey50", "pink4", "cyan4", "pink", "cyan")) +
    # scale_fill_manual(values = c("grey70", "grey50", "red", "red", "blue", "blue")) +
    theme(text = element_text(size = 8),
          axis.title.x = element_blank(),
          legend.position = "none")
  ggsave(paste0(OUTDIR,
                PAIR, ".png"),
         height = 2, width = 3)
  
  ## print statistics
  data.table("pair" = PAIR, "cell1" = CELL1, "cell2" = CELL2,
             rbind(c("type" = "neighbors", unlist(C_NEIGHBORS)),
                   c("type" = "non-neighbors", unlist(C_MATCHES))))
  
})))

## Plot summary plots
COMPARE_GENES_IN_RANGE$p.adj <- as.numeric(unlist(COMPARE_GENES_IN_RANGE$p.adj))

## Summarize p-values
ggplot(COMPARE_GENES_IN_RANGE[COMPARE_GENES_IN_RANGE$type == "neighbors",],
       aes(cell1, cell2)) +
  geom_tile(aes(fill = -log10(p.adj))) + theme_minimal() + 
  labs(fill = expression("-log"[10]*P)) +
  scale_fill_gradientn(colours = c("white", "red"),
                       limits = c(0,
                                  max(-log10(COMPARE_GENES_IN_RANGE$p.adj)))) +
  theme(text = element_text(size = 8),
        legend.position = c(0.9, 0.7),
        panel.grid.major = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45,
                                   vjust = 0.8)) +
  scale_y_discrete(limits = rev(unique(SUM_COMPARISONS$cell2))) +
  geom_text(aes(label = round(-log10(p.adj))), size = 2.5)

ggsave(paste0(OUTDIR, "p.heatmap_pval_neighbors.png"),
       height = 3, width = 3)

ggplot(COMPARE_GENES_IN_RANGE[COMPARE_GENES_IN_RANGE$type == "non-neighbors",],
       aes(cell1, cell2)) +
  geom_tile(aes(fill = -log10(p.adj))) + theme_minimal() + 
  labs(fill = expression("-log"[10]*P)) +
  scale_fill_gradientn(colours = c("white", "red"),
                       limits = c(0,
                                  max(-log10(COMPARE_GENES_IN_RANGE$p.adj)))) +
  theme(text = element_text(size = 8),
        legend.position = c(0.9, 0.7),
        panel.grid.major = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45,
                                   vjust = 0.8)) +
  scale_y_discrete(limits = rev(unique(SUM_COMPARISONS$cell2))) +
  geom_text(aes(label = round(-log10(p.adj))), size = 2.5)

ggsave(paste0(OUTDIR, "p.heatmap_pval_non-neighbors.png"),
       height = 3, width = 3)

## Compare neighbors to non-neoghbors
compare_means(p.adj ~ type,
              COMPARE_GENES_IN_RANGE)

ggplot(COMPARE_GENES_IN_RANGE, aes(type, -log10(p.adj), fill = type)) +
  geom_boxplot(outlier.shape = NA) + theme_bw()+
  geom_hline(yintercept = -log10(0.05)) +
  coord_cartesian(ylim = c(0, quantile(-log10(COMPARE_GENES_IN_RANGE$p.adj), 0.99))) +
  scale_fill_brewer(palette = "Set1") +
  labs(y = "-log10 P-value") +
  theme(text = element_text(size = 8),
        axis.title.x = element_blank(),
        legend.position = "none") 

ggsave(paste0(OUTDIR, "p.boxplot_compare_neighbors_to_non-neighbor.png"),
       height = 2, width = 2)

##### test if neighbours have stronger effect ----

dir.create("plots/neighbours_compared_to_close_non-neighbour/compare_neighbors_to_match_2/",
           recursive = T, showWarnings = F)
NEIGHBOUR_COMPARED_TO_MATCHED_2 <- as.data.table(do.call(rbind, lapply(CELL_PAIRS, function(PAIR){
  
  CELL1 <- SPLIT_CELLS[[PAIR]][1]
  CELL2 <- SPLIT_CELLS[[PAIR]][2]
  FC_DATA <- DESEQ[[PAIR]]
  MATCHED_NEIGHBORS_1 <- MATCHED_NEIGHBORS[[PAIR]][[1]]
  MATCHED_NEIGHBORS_2 <- MATCHED_NEIGHBORS[[PAIR]][[2]]
  
  ## Shape data
  MATCHES_1 <- left_join(data.table("gene" = MATCHED_NEIGHBORS_1$match,
                                    "cell" = CELL1,
                                    "type" = "match",
                                    "group" = paste(CELL1, "match")),
                         FC_DATA)
  MATCHES_2 <- left_join(data.table("gene" = MATCHED_NEIGHBORS_2$match,
                                    "cell" = CELL2,
                                    "type" = "match",
                                    "group" = paste(CELL2, "match")),
                         FC_DATA)
  
  DATA <- rbind(DATA_TO_PLOT[[PAIR]], # neighbors data
                MATCHES_1,
                MATCHES_2)
  LEVELS <- c(paste(CELL1, c("core", "neighbor", "match")),
              paste(CELL2, c("core", "neighbor", "match")))
  DATA$group <- factor(DATA$group,
                       levels = LEVELS,
                       ordered = T)
  
  ## Statistical test
  C1 <- wilcox.test(DATA[group == LEVELS[[2]], log2FoldChange],
                    DATA[group == LEVELS[[3]], log2FoldChange],
                    alternative = "less")
  C2 <- wilcox.test(DATA[group == LEVELS[[5]], log2FoldChange],
                    DATA[group == LEVELS[[6]], log2FoldChange],
                    alternative = "greater")  
  ## Plot

  XLABS <- paste0(c(paste(CELL1, c("HCTGs", "intra-TAD\nmate", "extra-TAD\nmatched"), sep = "\n"),
                    paste(CELL2, c("HCTGs", "intra-TAD\nmate", "extra-TAD\nmatched"), sep = "\n")),
                 "\n(n=", table(DATA$group), ")")
  Q <- quantile(DATA$log2FoldChange,
                na.rm = T)
  
  ggplot(DATA, aes(group, log2FoldChange, fill = group)) +
    geom_boxplot(outlier.shape = NA) + theme_bw() +
    geom_signif(comparisons = list(c(LEVELS[2], LEVELS[3]),
                                   c(LEVELS[5], LEVELS[6])),
                annotations = c(ifelse(C1$p.value < 0.05, formatC(C1$p.value, digits = 2), "NS"),
                                ifelse(C2$p.value < 0.05, formatC(C2$p.value, digits = 2), "NS")),
                y_position = Q[5]*.5,
                tip_length = 0, test = "t.test",
                textsize = 2.5) +
    scale_x_discrete(labels = XLABS) +
    ylab(paste0("log2 mRNA fold\nchange [", CELL2, "/", CELL1, "]")) +
    scale_fill_manual(values = c("grey70", "pink", "pink4", "grey50", "cyan", "cyan4")) +
    # scale_fill_manual(values = c("grey70", "grey50", "red", "red", "blue", "blue")) +
    theme(text = element_text(size = 8),
          axis.title.x = element_blank(),
          legend.position = "none")
  ggsave(paste0("plots/neighbours_compared_to_close_non-neighbour/compare_neighbors_to_match_2/",
                PAIR, ".png"),
         height = 2, width = 3.2)
  
  ## print statistics
  rbind(c("pair" = PAIR, "cell" = CELL1, C1),
        c("pair" = PAIR, "cell" = CELL2, C2))
  
})))

###### Plot summary ----
dt <- as.data.table(t(apply(NEIGHBOUR_COMPARED_TO_MATCHED_2, 1, unlist)))
dt$p.value <- as.numeric(dt$p.value)

ggplot(dt, aes("p.value", p.value)) +
  geom_boxplot() + theme_bw() +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  geom_hline(yintercept = 0.05) +
  xlab("comparisons (N=90)") +
  theme(axis.text.x = element_blank())
  


# NEIGHBOUR_COMPARED_TO_MATCHED <- as.data.table(do.call(rbind, mclapply(CELL_PAIRS, function(PAIR){
#
#   CELL1 <- SPLIT_CELLS[[PAIR]][1]
#   CELL2 <- SPLIT_CELLS[[PAIR]][2]
#   FC_DATA <- DESEQ[[PAIR]]
#   MATCHED_NEIGHBORS_1 <- MATCHED_NEIGHBORS[[PAIR]][[1]]
#   MATCHED_NEIGHBORS_2 <- MATCHED_NEIGHBORS[[PAIR]][[2]]
#
#   ## Shape data
#   MATCHES_1 <- left_join(rbind(data.table("gene" = MATCHED_NEIGHBORS_1$match,
#                                           "cell" = CELL1,
#                                           "type" = "match",
#                                           "group" = paste(CELL1, "match")),
#                                data.table("gene" = MATCHED_NEIGHBORS_1$neighbor,
#                                           "cell" = CELL1,
#                                           "type" = "neighbor",
#                                           "group" = paste(CELL1, "neighbor"))),
#                          FC_DATA)
#   MATCHES_2 <- left_join(rbind(data.table("gene" = MATCHED_NEIGHBORS_2$match,
#                                           "cell" = CELL2,
#                                           "type" = "match",
#                                           "group" = paste(CELL2, "match")),
#                                data.table("gene" = MATCHED_NEIGHBORS_2$neighbor,
#                                           "cell" = CELL2,
#                                           "type" = "neighbor",
#                                           "group" = paste(CELL2, "neighbor"))),
#                          FC_DATA)
#   NEIGHBOR_DATA <- DATA_TO_PLOT[[PAIR]]
#
#   DATA <- rbind(NEIGHBOR_DATA[type == "core"], # neighbors data
#                 MATCHES_1,
#                 MATCHES_2)
#
#   LEVELS <- levels(DATA$group)
#   ## Plot
#   CC <- c(CELL1, CELL2)
#   XLABS <- paste0(c(paste(CC, "core", sep = '\n'),
#                     paste(CC, "neighbor", sep = '\n'),
#                     paste(CC, "match", sep = '\n')),
#                   "\n(N=", table(DATA$group), ")")
#   COLORS <- c("grey70", "grey50", "pink4", "cyan4", "pink", "cyan")
#   Q <- quantile(DATA$log2FoldChange,
#                 na.rm = T)
#
#   ggplot(DATA[type != "core"], aes(group, log2FoldChange, fill = group)) +
#     geom_boxplot(outlier.shape = NA) + theme_bw() +
#     geom_signif(comparisons = list(c(LEVELS[3], LEVELS[4]),
#                                    c(LEVELS[5], LEVELS[6])),
#                 y_position = Q[5]*.75,
#                 tip_length = 0, test = "t.test",
#                 textsize = 2.5) +
#     scale_x_discrete(labels = XLABS[-c(1:2)]) +
#     scale_fill_manual(values = COLORS[-c(1:2)]) +
#     theme(text = element_text(size = 8),
#           axis.title.x = element_blank(),
#           legend.position = "none")
#   ggsave(paste0("plots/neighbours_compared_to_close_non-neighbour/compare_neighbors_to_match/",
#                 PAIR, ".png"),
#          height = 2, width = 2.5)
#
#   ## Statistical test
#   C_NEIGHBOURS <- compare_means(log2FoldChange ~ group,
#                                 DATA[group %in% c(LEVELS[3], LEVELS[4])],
#                                 method = "t.test")
#   C_MATCHES <- compare_means(log2FoldChange ~ group,
#                              DATA[group %in% c(LEVELS[5], LEVELS[6])],
#                              method = "t.test")
#
#   rbind(c("pair" = PAIR, "type" = "neighbor", C_NEIGHBOURS),
#         c("pair" = PAIR, "type" = "match", C_MATCHES))
#
# },
# mc.cores = 45)))

## plot comparison
NEIGHBOUR_COMPARED_TO_MATCHED$p.adj <- as.numeric(NEIGHBOUR_COMPARED_TO_MATCHED$p.adj)
NEIGHBOUR_COMPARED_TO_MATCHED$type <- factor(NEIGHBOUR_COMPARED_TO_MATCHED$type,
                                             levels = c("neighbor", "match"),
                                             ordered = T)
c <- compare_means(p.adj ~ type,
                   NEIGHBOUR_COMPARED_TO_MATCHED)
ggplot(NEIGHBOUR_COMPARED_TO_MATCHED, aes(type, -log10(p.adj), fill = type)) +
  geom_boxplot() + theme_bw() +
  geom_hline(yintercept = -log10(0.05)) +
  # geom_signif(comparisons = list(c(levels(NEIGHBOUR_COMPARED_TO_MATCHED$type))),
  #             annotations = c$p.adj,
  #             test = "t.test",
  #             tip_length = 0,
  #             extend_line = -0.01,
  #             textsize = 2) +
  ylab(expression(-log[10]*"(q-value)")) +
  scale_fill_brewer(palette = "Set1") +
  theme(text = element_text(size = 8),
        legend.position = "none",
        axis.title.x = element_blank())
ggsave("plots/neighbours_compared_to_close_non-neighbour/p.boxplot_sum_comparisons.png",
       height = 2, width = 2)


##  Test if co-residing genes co-express ----

### Find TAD neighbor pairs and controls per cell line ----
# NEIGHBOR_PAIRS_AND_CONTROLS_PER_CELL <- as.data.table(do.call(rbind, lapply(CELLS, function(CELL){
#
#   ## Find pairs of TAD neighbor genes
#   SET_NEIGHBORHOODS <- TAD_GENES[[CELL]][gene %in% EXP_RAW$gene] ## Use only genes with expression data]
#   SET_NEIGHBORHOODS <- SET_NEIGHBORHOODS[tad_id %in% tad_id[duplicated(tad_id)]] ## Take only TADs with more than 1 gene
#   NEIBOUR_PAIRS <- as.data.table(do.call(rbind,
#                                          lapply(unique(SET_NEIGHBORHOODS$tad_id), function(TAD) {
#                                            data.frame(tad=TAD,
#                                                       gene = t(combn(SET_NEIGHBORHOODS$gene[SET_NEIGHBORHOODS$tad_id==TAD], 2)))
#                                            })))
#   ## Pull pair's distance
#   NEIBOUR_PAIRS$distance <- unlist(mclapply(1:nrow(NEIBOUR_PAIRS), function(i){
#     row <- NEIBOUR_PAIRS[i,]
#     pull_gene_pair_distance(row[[2]], row[[3]], GENE_PAIRS_DISTANCE)
#   }, mc.cores = 45))
#
#   ## Filter out paires closer than data resolution
#   NEIBOUR_PAIRS <- NEIBOUR_PAIRS[distance > 10000]
#   names(NEIBOUR_PAIRS)[2:3] <- c("gene1", "gene2")
#
#   ## Filter genes with expression data
#   GENE_PAIRS_DISTANCE_WITH_EXP_DATA <- GENE_PAIRS_DISTANCE[gene1 %in% EXP_RAW$gene &
#                                                              gene2 %in% EXP_RAW$gene]
#   ## Add id
#   NEIBOUR_PAIRS$neighbors.id <- paste0(CELL, "_", 1:nrow(NEIBOUR_PAIRS))
#
#   ## Find controls
#   CONTROL_PAIRS <- as.data.table(do.call(rbind, mclapply(1:nrow(NEIBOUR_PAIRS), function(i){
#
#     row <- NEIBOUR_PAIRS[i,]
#
#     ALL_TAD_GENES <- SET_NEIGHBORHOODS[tad_id == row[[1]], gene]
#     FILTERED_PAIRS_DISTANCES <- GENE_PAIRS_DISTANCE_WITH_EXP_DATA[(gene1 %in% row | gene2 %in% row)]  ## pair contain either genes
#     FILTERED_PAIRS_DISTANCES <- FILTERED_PAIRS_DISTANCES[ distance < row$distance] ## distance is equal or lower to the distance of the pair
#     FILTERED_PAIRS_DISTANCES <- FILTERED_PAIRS_DISTANCES[!(gene1 %in% ALL_TAD_GENES & gene2 %in% ALL_TAD_GENES)] %>% ## Gene is not in the same TAD
#       arrange(desc(distance))
#
#     MATCH1 <- FILTERED_PAIRS_DISTANCES[gene2 == row[[2]] | gene2 == row[[2]]][1,]
#     MATCH2 <- FILTERED_PAIRS_DISTANCES[gene2 == row[[3]] | gene2 == row[[3]]][1,]
#
#     RESULTS <- cbind(rbind(MATCH1, MATCH2),
#                      "neighbors.id" = row$neighbors.id)
#
#     RESULTS <- RESULTS[complete.cases(RESULTS)]
#
#   }, mc.cores = 45)))
#   ## Remove NAs
#   CONTROL_PAIRS <- CONTROL_PAIRS[!is.na(CONTROL_PAIRS[[1]])]
#
#   ## Plot
#   ALL <- cbind(rbind(cbind(NEIBOUR_PAIRS[,-1], group = "neighbor"),
#                      cbind(CONTROL_PAIRS[,-1], group = "control")),
#                cell = CELL)
#   data_to_plot <- ALL[ALL$neighbors.id %in% ALL$neighbors.id[duplicated(ALL$neighbors.id)]] ## Remove non-controlled pairs
#
#   ggplot(data_to_plot, aes(distance, color = group)) +
#            geom_density() + theme_bw() +
#     labs(x = "distance [Kbp]") +
#     scale_color_brewer(palette = "Set1", direction = -1, labels = c("match", "neighbor")) +
#     theme(text = element_text(size = 8),
#           legend.key.size = unit(0.5, "cm"),
#           legend.title = element_blank())
#   ggsave(paste0("plots/neighbours_compared_to_close_non-neighbour/distance_neighbor_match_per_cell/",
#                 CELL, ".png"),
#          height = 2, width = 3)
#   ## Aggregate
#   ALL
# })))
# write.table(NEIGHBOR_PAIRS_AND_CONTROLS_PER_CELL,
#             "all_neighbor_pairs_and_controls_per_cell.tsv",
#             append = F, quote = F, sep = '\t', row.names = F, col.names = T)
NEIGHBOR_PAIRS_AND_CONTROLS_PER_CELL <- fread("all_neighbor_pairs_and_controls_per_cell.tsv")

#### Summarize per cell line results ----

##### Neighbor pairs
PAIRS_PER_CELL <- as.data.table(table(NEIGHBOR_PAIRS_AND_CONTROLS_PER_CELL[group == "neighbor", cell]))
median(PAIRS_PER_CELL$N); sd(PAIRS_PER_CELL$N)
ggplot(NEIGHBOR_PAIRS_AND_CONTROLS_PER_CELL[group == "neighbor"], aes(cell)) +
  geom_bar() + theme_bw() +
  ylab("TAD neighbor gene pairs [N]") +
  coord_flip() +
  theme(text = element_text(size = 8),
        axis.title.y = element_blank())
ggsave("plots/p.barplot.TAD_residing_gene_pairs_N.png",
       height = 2, width = 2)

##### Control pairs
PAIRS_PER_CELL <- as.data.table(table(NEIGHBOR_PAIRS_AND_CONTROLS_PER_CELL[group == "control", cell]))
median(PAIRS_PER_CELL$N); sd(PAIRS_PER_CELL$N)
ggplot(NEIGHBOR_PAIRS_AND_CONTROLS_PER_CELL[group == "control"], aes(cell)) +
  geom_bar() + theme_bw() +
  ylab("control gene pairs [N]") +
  coord_flip() +
  theme(text = element_text(size = 8),
        axis.title.y = element_blank())
ggsave("plots/p.barplot.TAD_residing_gene_pairs_N_control.png",
       height = 2, width = 2)

### Pool pairs and correlate expression over 10 cell lines ----
NEIGHBOR_PAIRS_AND_CONTROLS_PER_CELL <- NEIGHBOR_PAIRS_AND_CONTROLS_PER_CELL[distance < 1000000]
ALL_NEIGHBOR_PAIRS_AND_CONTROLS <- NEIGHBOR_PAIRS_AND_CONTROLS_PER_CELL[NEIGHBOR_PAIRS_AND_CONTROLS_PER_CELL$neighbors.id %in%  
                                                          NEIGHBOR_PAIRS_AND_CONTROLS_PER_CELL$neighbors.id[duplicated(NEIGHBOR_PAIRS_AND_CONTROLS_PER_CELL$neighbors.id)]]
ALL_NEIGHBOR_PAIRS_AND_CONTROLS$group <- factor(ALL_NEIGHBOR_PAIRS_AND_CONTROLS$group,
                                                levels = c("neighbor", "control"),
                                                ordered = T)

ggplot(ALL_NEIGHBOR_PAIRS_AND_CONTROLS, aes(distance/1000000, color = group)) +
  geom_density() + theme_bw() +
  labs(x = "distance [Mbp]") +
  scale_color_brewer(palette = "Set1", 
                     labels = c("intra-TAD", "extra-TAD")) +
  theme(text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_blank(),
        legend.position = c(0.7, 0.8))

ggsave("plots/neighbours_compared_to_close_non-neighbour/p.density.matched_neighbors_pooled.png",
       height = 2, width = 2)

### Calculate correlations ----
CORRELATIONS_ALL_NEIGHBOR_PAIRS_AND_CONTROLS <- as.data.table(do.call(rbind, 
                                            mclapply(1:nrow(ALL_NEIGHBOR_PAIRS_AND_CONTROLS), function(i){
  row <- ALL_NEIGHBOR_PAIRS_AND_CONTROLS[i,]
  V1 <- EXP_RAW[gene == row[[1]], 3:12]
  V2 <- EXP_RAW[gene == row[[2]], 3:12]
  
  c(row, cor.test(unlist(V1), unlist(V2), method = "spearman"))
}, mc.cores = 45)))


### Plot results ----
CORRELATIONS_ALL_NEIGHBOR_PAIRS_AND_CONTROLS$estimate <- unlist(CORRELATIONS_ALL_NEIGHBOR_PAIRS_AND_CONTROLS$estimate)
CORRELATIONS_ALL_NEIGHBOR_PAIRS_AND_CONTROLS$group <- factor(unlist(CORRELATIONS_ALL_NEIGHBOR_PAIRS_AND_CONTROLS$group),
                                                levels = c("neighbor", "control"),
                                                ordered = T,
                                                labels = c("intra-TAD\nmate",
                                                           "extra-TAD\nmatched"))
c <- compare_means(estimate ~ group, CORRELATIONS_ALL_NEIGHBOR_PAIRS_AND_CONTROLS)
ggplot(CORRELATIONS_ALL_NEIGHBOR_PAIRS_AND_CONTROLS, aes(estimate, color = group)) +
  stat_ecdf() + theme_bw() +
  scale_color_brewer(palette = "Set1") +
  # annotate(geom = "text",
  #          x = -0.4, y = 0.75,
  #          label = paste0(c$method,
  #                         "\nP-value=", c$p.adj),
  #          size = 2.5) +
  labs(x = "coefficient of correlation",
       y = "cumulative density") +
  theme(text = element_text(size = 8),
        legend.position = c(0.3, 0.75),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black",
                                         size = 0.2))
ggsave("plots/p.cumulative.expression_correlation.png",
       height = 2, width = 2)