# ------------------------------------------------------------------------------
# Analysis of Promoter-Promoter (PP) Loops and Their Gene Expression Signatures
# ------------------------------------------------------------------------------
# This script explores the biological significance of promoter-promoter (PP) loops 
# identified in Hi-C data. Specifically, it evaluates whether genes connected by PP loops:
#
# 1. Exhibit higher co-expression across multiple cell lines.
# 2. Show elevated expression levels compared to distance-matched gene pairs.
# 3. Are over-represented compared to randomized loop structures.
#
# Additional goals include:
# - Annotating chromatin loops by their proximity to gene promoters (PP, PE, EE).
# - Quantifying the occurrence of PP loops across cell lines.
# - Comparing expression correlation and gene activity between looped and control gene pairs.
# - Estimating the statistical enrichment of PP loops using permutation tests.
#
# INPUT:
# - `tssLoops_values.tsv`: Loop anchors linked to genes.
# - `pooled_loops_binary.tsv`: Binary matrix of loop presence across cell lines.
# - `gene_expression_CNVcorrected.tsv`: Gene expression data across cell lines.
# - `distance_between_gene_pairs.tsv`: Genomic distances between gene pairs.
# - `ensembl_hg38_protein_coding_first_TSS_5000win.bed`: Promoter coordinates.
#
# OUTPUT:
# - Annotated loop files with promoter classifications (PP/PE/EE).
# - Distance and correlation comparisons between test and control gene pairs.
# - Statistical plots: barplots, boxplots, density plots.
# - Permutation-based enrichment analysis of PP loops.
#
# AUTHOR: Gony Shanel
# DATE: March 2025
# ------------------------------------------------------------------------------

## Set the scripts directory as working directory.
## This script is designed with relative paths for the specific data structure
## with which it was uploaded to https://github.com/ElkonLab/micro-C.

## This script requires bedtools to be installed (https://bedtools.readthedocs.io/en/latest/index.html)
path.bedtools <- "/path/to/bedtools"

# Function ----------------------------------------------------------------

annotate_loop <- function(bedpe,
                          tss.file="../../data/ensembl_hg38_protein_coding_first_TSS_5000win.bed"){
  annot1 <- intersect_with_tss(bedpe[,1:3], tss.file)
  annot1 <- ifelse(annot1 == TRUE, "P", "E")
  annot2 <- intersect_with_tss(bedpe[,4:6], tss.file)
  annot2 <- ifelse(annot2 == TRUE, "P", "E")
  annot <- paste0(annot1, annot2)
}

gene_pairs_correlation <- function(pairs){
  unlist(lapply(1:nrow(pairs),
                function(i){
                  GENE1 <- pairs[[1]][i] 
                  GENE2 <- pairs[[2]][i]
                  EXP1 <- sapply(CELLS, function(CELL) unlist(EXP_RAW[gene == GENE1, CELL, with = F]))
                  EXP2 <- sapply(CELLS, function(CELL) unlist(EXP_RAW[gene == GENE2, CELL, with = F]))
                  as.numeric(cor(EXP1, EXP2))
                }))
}

###############
# for (i in 1:nrow(pairs)){
#   print(i)
#   GENE1 <- pairs[[1]][i] 
#   GENE2 <- pairs[[2]][i]
#   EXP1 <- sapply(CELLS, function(CELL) unlist(EXP_RAW[gene == GENE1, CELL, with = F]))
#   EXP2 <- sapply(CELLS, function(CELL) unlist(EXP_RAW[gene == GENE2, CELL, with = F]))
#   as.numeric(cor(EXP1, EXP2))
# }
######################

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

find_control_pair <- function(gene1, gene2, n=1){
  
  PAIRs <- data.table("gene1" = "DUMMY", "gene2" = "DUMMY")
  dis <- pull_gene_pair_distance(gene1, gene2, GENE_PAIRS_DISTANCE) # The genomic distance between the genes
  
  while (nrow(PAIRs) != n+1){
    i=0
    I_ROW <- GENE_PAIRS_DISTANCE[distance == dis, which = T][[1]] # A row number containing pair with the same distance
    C_ROW <- GENE_PAIRS_DISTANCE[I_ROW+i]
    
    PAIR <- C_ROW[, c("gene1", "gene2")]
    while (
      ## assure control pair was not asigned already (for multiple output)
      nrow(PAIRs[gene1 == PAIR[[1]] & gene2 == PAIR[[2]]]) > 0 |
      # assure the control row does not contain the input genes
      sum(grepl(paste0(gene1, "|", gene2) ,C_ROW)) > 0 |
      # assure the control row is not a pp loop
      sum(sapply(1:nrow(ALL_PP_PAIRS), function(I){sum(grepl(paste0(C_ROW[[2]], "|", C_ROW[[3]]), ALL_PP_PAIRS[I])) == 2})) != 0 
    ){
      i = i+1
      C_ROW <- GENE_PAIRS_DISTANCE[I_ROW+i]
      PAIR <- C_ROW[, c("gene1", "gene2")]
    } 
    PAIRs <- rbind(PAIRs, PAIR)
  }
  PAIRs[-1]
}

# Arguments ---------------------------------------------------------------

GENE_PAIRS_DISTANCE_PATH <- "../../data/distance_between_gene_pairs.tsv"
TSS_LOOPS_PATH <- "../../data/loops/tssLoops_values.tsv"
LOOPS_BINARY_PATH <- "../../data/loops/pooled_loops_binary.tsv"
TSS_WITH_5K_WINDOWS_PATH <- "../../data/ensembl_hg38_protein_coding_first_TSS_5000win.bed"
CELL_LINES_PATH <- "../../data/cell_lines.txt"
CHROM_SIZE_PATH <- "../../data/hg38.chrom.sizes"
EXP_PATH <- "../../data/gene_expression_CNVcorrected.tsv"

LOOP_OCCURENCE_THRESHOLD <- 7

# Commands ----------------------------------------------------------------
library(ggplot2)
library(ggsignif)
library(RColorBrewer)
library(data.table)
library(parallel)
options(bedtools.path = path.bedtools)
library(bedtoolsr)
source("../../_functions_hic.r")
library(dplyr)

## Load data ----
TSS_LOOPS <- fread(TSS_LOOPS_PATH)

TSS_RAW <- fread(TSS_WITH_5K_WINDOWS_PATH)

LOOPS_BINARY <- fread(LOOPS_BINARY_PATH)
LOOPS_BINARY <- LOOPS_BINARY[, !grepl("nutlin", colnames(LOOPS_BINARY)), with = F]
colnames(LOOPS_BINARY) <- sub("-control", "", colnames(LOOPS_BINARY))

CELLS <- readLines(CELL_LINES_PATH)
CHROM_SIZE <- fread(CHROM_SIZE_PATH,
                    col.names = c("chr", "size"))

EXP_RAW <- fread(EXP_PATH)
EXP_RAW <- EXP_RAW[, !grepl("nutlin", colnames(EXP_RAW)), with = F]
colnames(EXP_RAW) <- sub("-control", "", colnames(EXP_RAW))
EXP_RAW$mean_exp <- rowSums(EXP_RAW[, 3:12])

GENE_PAIRS_DISTANCE <- fread(GENE_PAIRS_DISTANCE_PATH)
GENE_PAIRS_DISTANCE <- GENE_PAIRS_DISTANCE[gene1 %in% EXP_RAW$gene & gene2 %in% EXP_RAW$gene] # Filter out genes with no expression data
GENE_PAIRS_DISTANCE <- GENE_PAIRS_DISTANCE[order(distance)]

## Find promoter-promoter loops ----
TSS_LOOPS$annotation <- annotate_loop(bedpe = TSS_LOOPS)
PP_LOOPS <- aggregate(.~peak_name, TSS_LOOPS[annotation == "PP", c("peak_name", "gene")], paste)
ALL_PP_PAIRS <- unique(as.data.table(do.call(rbind, 
                                      lapply(PP_LOOPS$peak_name, 
                                             function(loop){
                                               genes <- PP_LOOPS[PP_LOOPS$peak_name == loop, "gene"]
                                               t(combn(x = unlist(genes), m = 2))
                                             }))))
## Find constitutive loops ----
PP_FREQ <- data.table(LOOPS_BINARY,
                      "frequency" = rowSums(LOOPS_BINARY[,8:17]))[peak_name %in% PP_LOOPS$peak_name]
COMMON_PP_LOOPS <- merge(PP_FREQ[frequency > LOOP_OCCURENCE_THRESHOLD], 
                         TSS_LOOPS[, c("peak_name", "gene", "sym")])
COMMON_PP_LOOPS <- COMMON_PP_LOOPS %>%
  arrange(as.numeric(sub("chr", "", chr1)), start1)
write.table(COMMON_PP_LOOPS,
            "constitutive_PP_loops.tsv",
            append = F, 
            quote = F, 
            sep = '\t',
            row.names = F,
            col.names = T)
## Break triplets into pairs 
PP_PAIRS <- unique(as.data.table(do.call(rbind, 
                                  lapply(COMMON_PP_LOOPS$peak_name, 
                                         function(loop){
                                           t(combn(x = unlist(PP_LOOPS[PP_LOOPS$peak_name == loop, "gene"]), m = 2))
                                         }))))
## Filter out pairs that reside in the same anchor
G1_TSS <- unlist(lapply(PP_PAIRS[[1]], function(GENE){TSS_RAW[gene == GENE, start]}))
G2_TSS <- unlist(lapply(PP_PAIRS[[2]], function(GENE){TSS_RAW[gene == GENE, start]}))
G1_G2_DISTANCE <- abs(G1_TSS-G2_TSS)
PP_PAIRS <- PP_PAIRS[G1_G2_DISTANCE > 5000]
## Filter out genes with no expression data
PP_PAIRS <- unique(PP_PAIRS[V1 %in% EXP_RAW$gene & V2 %in% EXP_RAW$gene])

# Generate a control group of gene pairs
# # Run once. Than use saved SDR for future runs to save time
# C_PP_PAIRS <- as.data.table(do.call(rbind,
#                                     mclapply(1:nrow(PP_PAIRS),
#                                              function(i){
#                                                find_control_pair(gene1 = PP_PAIRS[i, V1],
#                                                                  gene2 = PP_PAIRS[i, V2], n=10)
#                                              },
#                                              mc.cores = 40)))
# saveRDS(C_PP_PAIRS, "500_control_loops.rds")
C_PP_PAIRS <- readRDS("500_control_loops.rds")
names(C_PP_PAIRS) <- c("V1", "V2")

LABS <- c("test" = paste0("test\n(N=", nrow(PP_PAIRS),")"),
          "control" =  paste0("control\n(N=", nrow(C_PP_PAIRS),")"))
DATA_TO_PLOT <- rbind(cbind(PP_PAIRS, "group" = LABS["test"]),
                      cbind(C_PP_PAIRS, "group" = LABS["control"]))
DATA_TO_PLOT$distance <- sapply(1:nrow(DATA_TO_PLOT), function(i){
  gene1_coord <- TSS_RAW[gene == DATA_TO_PLOT[i,1], start]
  gene2_coord <- TSS_RAW[gene == DATA_TO_PLOT[i,2], start]
  dis <- abs(gene1_coord - gene2_coord)
})
DATA_TO_PLOT$group <- factor(DATA_TO_PLOT$group, levels = c(LABS["test"], LABS["control"]), ordered = T)
ggplot(DATA_TO_PLOT, aes(distance, color = group)) +
  geom_density() + theme_bw() +
  theme(text = element_text(size = 8)) +
  scale_color_brewer(palette = "Set1")
ggsave("plots/p.density_PPloop_length.png",
       height = 2, width = 3.5)

### For each cell line: how many loops are PP and how many are PE and EE? For each category plot distance distribution ####

LOOPS_BINARY$annotation <- annotate_loop(LOOPS_BINARY)
fwrite(LOOPS_BINARY, 
       "binary_loops_annotated.tsv",
       sep = '\t')

# Plot
LOOPS_BINARY$annotation <- sub("PE", "EP", LOOPS_BINARY$annotation)

ggplot(LOOPS_BINARY, aes(annotation, fill = annotation)) +
  geom_bar() + theme_bw() +
  geom_text(aes(label = ..count..), stat = "count", 
            vjust = -0.1, colour = "black", size = 2.5) +
  ylim(0, max(table(LOOPS_BINARY$annotation))*1.03) +
  labs(x = "loop type",
       y= "loops [N]") +
  scale_x_discrete(labels = c("DD", "PD", "PP")) +
  scale_fill_brewer(palette = "Set1") +
  theme(text = element_text(size = 8),
        legend.position = "none")
ggsave("plots/p.barplot_pooled_loops_annotation.png",
       height = 2, width = 1.5)

SUM_LOOP_TYPES <- as.data.table(do.call(rbind, lapply(CELLS, function(CELL){
  CELL_LOOPS <- as.data.table(table(LOOPS_BINARY[LOOPS_BINARY[[CELL]] == 1, annotation]))
  CELL_LOOPS <- cbind(CELL_LOOPS, "portion" = round(with(CELL_LOOPS, N/sum(N)), 2), "cell" = CELL)
})))

SUM_LOOP_TYPES$V1 <- factor(sub("EP", "PD", sub("EE", "DD", SUM_LOOP_TYPES$V1)), 
                            levels = c("DD", "PD", "PP"),
                            ordered = T)

### Number of loops for each type grouped by cell line 

ggplot(SUM_LOOP_TYPES, aes(x = cell, y = N, fill = V1)) +
  geom_bar(position = "dodge", stat = "identity") + theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Cell line",
       y = "genes [N]", 
       fill = "loop type") +
  theme(text = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, 
                                   vjust = 0.7),
        legend.position = "top",
        legend.key.size = unit(0.3, 'cm'))

ggsave("plots/p.barplot_loop_types_per_cell.png",
       height = 2, width = 3)

## Number of loops in each group across cell lines

ggplot(SUM_LOOP_TYPES, aes(V1, N, fill = V1)) +
  geom_boxplot() + theme_bw() +
  geom_text(data = aggregate(N ~ V1, SUM_LOOP_TYPES, mean),
            aes(label = round(N), y = N + 1200), size = 2) +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "loop type",
       y = "loops [N]", 
       fill = "loop type") +
  theme(text = element_text(size = 8),
        legend.position = "none")

ggsave("plots/p.boxplot_loop_type_numbers_across_cells.png",
       height = 2, width = 1.5)

## The fraction of the loop type out of all detected loop across cell lines

ggplot(SUM_LOOP_TYPES, aes(V1, portion, fill = V1)) +
  geom_boxplot() + theme_bw() +
  geom_text(data = aggregate(portion ~ V1, SUM_LOOP_TYPES, mean),
            aes(label = round(portion, 2), y = portion + 0.06), size = 2) +
  scale_fill_brewer(palette = "Set1") +
  ylim(0, 1) +
  labs(x = "loop type",
       y = "fraction", 
       fill = "loop type") +
  theme(text = element_text(size = 8),
        legend.position = "none")

ggsave("plots/p.boxplot_loop_type_portion_across_cells.png",
       height = 2, width = 1.5)

####  Does the expression of common loop-linked genes correlate between pairs of cells? ####

PP_PARIS <- PP_PAIRS[V1 %in% EXP_RAW$gene & V2 %in% EXP_RAW$gene]
DT_CORRELATION <- as.data.table(rbind(cbind("correlation" = gene_pairs_correlation(pairs = PP_PAIRS),
                                            "group" = "test"),
                                      cbind("correlation" = gene_pairs_correlation(C_PP_PAIRS),
                                            "group" = "control")))

DT_CORRELATION$correlation <- as.numeric(DT_CORRELATION$correlation)
DT_CORRELATION$group <- factor(DT_CORRELATION$group, levels = c("test", "control"), ordered = T)

ggplot(DT_CORRELATION, aes(group, correlation, fill = group)) +
  geom_boxplot() + theme_bw() +
  geom_signif(comparisons = list(c("test", "control")),
              textsize = 2) +
  scale_x_discrete(labels =paste0(levels(DT_CORRELATION$group), "\n(N=", table(DT_CORRELATION$group), ")")) + 
  scale_fill_brewer(palette = "Set1") +
  ylim(-1, 1.2) +
  ylab("correlation [r]") +
  theme(text = element_text(size = 8),
        axis.title.x = element_blank(),
        legend.position = "none")

ggsave("plots/p.boxplot_PPloop_associated_gene_expresssion_correlation.png",
       height = 2, width = 1.5)

#### Is the expression level of common PP-loop assocaited gene pairs different then random pairs with same distance? ####

DT <- as.data.table(rbind(cbind("mean_exp" = log(EXP_RAW[gene %in% union(PP_PAIRS$V1, PP_PAIRS$V2), mean_exp]), "group" = "test"),
                          cbind("mean_exp" = log(EXP_RAW[gene %in% union(C_PP_PAIRS[[1]], C_PP_PAIRS[[1]]), mean_exp]), "group" = "control")))
DT$mean_exp <- as.numeric(DT$mean_exp)

ggplot(DT, aes(group, log(mean_exp), fill = group)) +
  geom_boxplot(outlier.shape = NA) + theme_bw() +
  geom_signif(comparisons = list(c("test", "control")),
              textsize = 2) +
  scale_fill_brewer(palette = "Set1") +
  scale_x_discrete(labels = paste0(levels(DT_CORRELATION$group),
                                   "\n(N=", table(DT_CORRELATION$group), ")")) + 
  ylim(0, max(log(DT$mean_exp))+1) +
  ylab("mean expression [log CPM]") +
  theme(text = element_text(size = 8),
        axis.title.x = element_blank(),
        legend.position = "none") 

ggsave("plots/p.boxplot_PPloop_associated_gene_mean_expresssion.png",
       height = 2, width = 1.5)

#### Are there more promoter loops than expected? ####

N_P_LOOPS <- sum(LOOPS_BINARY$annotation != "EE")

## calculate the background distribution of the number of promoter loops by shifting all
# loops 10k-1M bp and anottation tssLoops 10,000 times  
# system("/specific/elkon/gonyshanel/2020-03_hsieh/hic/annotation/tads/enrichment_for_promoters/n_P-TADs_permutations.R")
N_PERM_P_DISTRIBUTION <- as.integer(readLines("../n_P-loops_permutations.txt"))

P <- sum(N_PERM_P_DISTRIBUTION > N_P_LOOPS)/length(N_PERM_P_DISTRIBUTION)

DATA <- data.table("N" = c(N_PERM_P_DISTRIBUTION, N_P_LOOPS))

ggplot(DATA, aes(N)) +
  geom_density() + theme_bw() +
  geom_segment(aes(xend = N_P_LOOPS, yend = 0.000015, x = N_P_LOOPS, y = 0.00025),
               arrow = arrow(length = unit(0.1, "cm"))) +
  labs(x = "promoter loops [N]",
       y = "density") +
  theme(text = element_text(size = 8))

ggsave("plots/p.density_promoter_loops_n.png",
       height = 1, width = 2.5)
#### Are there more PP loops then expected within th P-loops population

## calculate the distribution of the number of promoter-promoter loops over 10,000 permutations

P_LOOPS <- LOOPS_BINARY[annotation != "EE"]
PP_PORTION <- with(P_LOOPS, sum(annotation == "PP")/length(annotation))
CHROMs <- unique(LOOPS_BINARY$chr1)

## Generate permutations
# PERM_P_LOOPS <- mclapply(1:10000, function(i){
#   as.data.table(do.call(rbind, lapply(CHROMs, function(CHR){
# 
#     CHR_SIZE <- CHROM_SIZE[chr == CHR, size]
#     ## Shape data
#     LOOPS <- P_LOOPS[chr1 == CHR, -c(7:17)]
#     LOOPS$dis <- LOOPS[, start2-start1]
#     LOOPS$shuffle_dis <- sample(LOOPS$dis)
#     ## PP and PE loops are anchored to anchor1
#     LOOPS[annotation == "PP" | annotation == "PE", 5:6] <- LOOPS[annotation == "PP" | annotation == "PE", .(start1+shuffle_dis, end1+shuffle_dis)]
#     ## EP loops are anchored to anchor2
#     LOOPS[annotation == "EP", 2:3] <- LOOPS[annotation == "EP", .(start2-shuffle_dis, end2-shuffle_dis)]
#     ## Remove out of bound loops
#     LOOPS <- LOOPS[start1 > 0 & end2 < CHR_SIZE]
# 
#     LOOPS[, 1:6]
#   })))},
# mc.cores = 30)
# saveRDS(PERM_P_LOOPS, "data/permutated_P-loops.rds")
readRDS("data/permutated_P-loops.rds")

# PERM_P_LOOPS_ANNOTS <- lapply(PERM_P_LOOPS,
#                                 annotate_loop)
# saveRDS(PERM_P_LOOPS_ANNOTS, "data/permutated_P-loops_annotations.rds")
PERM_P_LOOPS_ANNOTS <-readRDS("data/permutated_P-loops_annotations.rds")

PERM_PP_PORTION <- sapply(PERM_P_LOOPS_ANNOTS, function(V){sum(V=="PP")/length(V)})

## calculate p value
P <- sum(PERM_PP_PORTION > PP_PORTION)/length(PERM_PP_PORTION)
## plot
DATA <- data.table("N" = c(PERM_PP_PORTION, PP_PORTION))

ggplot(DATA, aes(N)) +
  geom_density() + theme_bw() +
  geom_segment(aes(xend = PP_PORTION, yend = 0.000015, x = PP_PORTION, y = 90),
               arrow = arrow(length = unit(0.1, "cm"))) +
  labs(x = "PP-loop density [PP-loops/P-loops]",
       y = "density") +
  theme(text = element_text(size = 8))

ggsave("plots/p.density_PP_loop_density.png",
       height = 1, width = 2.5)

