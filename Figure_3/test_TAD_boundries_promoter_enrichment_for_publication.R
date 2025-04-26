#### Are there more promoter TADs than expected? ####

# ------------------------------------------------------------------------------
# This script evaluates whether **promoter-associated TADs (P-TADs)** are 
# more frequent than expected by chance. It does this by:
# 1. Annotating observed TADs to check how many contain promoters.
# 2. Performing a **permutation test** to compare observed vs. expected P-TAD counts.
# 3. Extending this to **promoter-promoter TADs (PP-TADs)** within the P-TAD population.
#
# INPUT:
# - `pooled_TAD_10kbp_boundaries.bedpe`: Defines the boundaries of pooled TADs.
# - `ensembl_mart_hg38_1st_TSS_protein_coding.txt`: Promoter locations.
# - `hg38.chrom.sizes`: Chromosome sizes for defining valid genome regions.
# - `data/permutated_P-tads.rds`: Randomly shuffled P-TADs for permutation testing.
# - `data/permutated_P-tads_annotations.rds`: Annotations of shuffled TADs.
#
# OUTPUT:
# - Density plots showing observed vs. expected P-TADs and PP-TADs.
# - P-values assessing significance of enrichment.
#
# AUTHOR: Gony Shanel
# DATE: March 2025
# ------------------------------------------------------------------------------

## Set the scripts directory as working directory.
## This script is designed with relative paths for the specific data structure
## with which it was uploaded to https://github.com/ElkonLab/micro-C.

## This script requires bedtools to be installed (https://bedtools.readthedocs.io/en/latest/index.html)
path.bedtools <- "/path/to/bedtools"

# Arguments ---------------------------------------------------------------

# List of chromosomes (autosomes + X)
CHROMs <- paste0("chr", c(1:22, "X"))

# Path to pooled TAD boundaries file
path.pooledTADs <- "../data/pooled_TAD_10kbp_boundaries.bedpe"

# Path to TSS file
path.tss <- "../data/ensembl_mart_hg38_1st_TSS_protein_coding.txt"

# Commands ----------------------------------------------------------------

## Set environment----

# Load required libraries
library(data.table)   # Efficient data handling
library(ggplot2)      # Visualization
source("../_functions_hic.r")  # Custom Hi-C functions
library(bedtoolsr)    # Wrapper for bedtools genomic operations
options(bedtools.path = path.bedtools)

# Set working directory and create output directories (if they don't exist)
dir.create("plots", showWarnings = F)
dir.create("data", showWarnings = F)

## Load data ----

# Load pooled TAD boundaries
pooledTADs <- fread(path.pooledTADs)
names(pooledTADs) <- c("chr1", "start1", "end1", "chr2", "start2", "end2", "tad")

# Load chromosome sizes for filtering valid genomic regions
CHROM_SIZE <- fread("../data/hg38.chrom.sizes",
                    col.names = c("chr", "size"))

## Portion of P-TADs ----

# Intersect TADs with promoter regions and compute the percentage of P-TADs
dt <- setDT(bt.pairtobed(pooledTADs,
                         fread(path.tss)[, 1:4]))

# Calculate the proportion of promoter-associated TADs
round(nrow(unique(dt[,1:6]))/nrow(pooledTADs)*100)

## Annotate TADs ----

# Assign functional annotations (PP, PE, EP, DD) to TADs
pooledTADs$annotation <- annotate_loops(pooledTADs)

# Count the number of **non-DD** (promoter-containing) TADs
N_P_TADS <- sum(pooledTADs$annotation != "DD")

## Permutation Test: Background Distribution ----

# Shift all TADs randomly (10k-1M bp) and recalculate P-TAD counts 10,000 times
# This system call runs an external script that generates these permutations:
# system("/specific/elkon/gonyshanel/2020-03_hsieh/hic/annotation/tads/enrichment_for_promoters/n_P-TADs_permutations.R")

# Load precomputed permutation results (expected P-TAD counts from shuffled data)
N_PERM_P_DISTRIBUTION <- as.integer(readLines("../data/n_P-boundaries_permutations.txt"))

# Compute p-value: Fraction of permutations where shuffled P-TAD count exceeds observed count
P <- sum(N_PERM_P_DISTRIBUTION > N_P_TADS)/length(N_PERM_P_DISTRIBUTION)

# Prepare data for visualization
DATA <- data.table("N" = c(N_PERM_P_DISTRIBUTION, N_P_TADS))

# Plot density of expected vs. observed P-TAD counts
ggplot(DATA, aes(N)) +
  geom_density() + theme_bw() +
  geom_segment(aes(xend = N_P_TADS, yend = 0.000015, x = N_P_TADS, y = 0.00025),
               arrow = arrow(length = unit(0.1, "cm"))) +
  labs(x = "Promoter TADs [N]",
       y = "Density") +
  theme(text = element_text(size = 8))

# Save plot
ggsave("plots/p.density_promoter_tads_n.png", height = 1, width = 2.5)

#### Are there more PP-TADs than expected within the P-TAD population? ####

## Compute observed **PP-TAD proportion** (within promoter-associated TADs)
P_TADS <- pooledTADs[annotation != "DD"]
PP_PORTION <- with(P_TADS, sum(annotation == "PP")/length(annotation))

## Permutation Test: Background Distribution of PP-TADs ----
# # Run once. Than use saved SDR for future runs to save time
# # Generate 10,000 shuffled versions of P-TADs by randomizing their positions
# PERM_P_TADS <- mclapply(1:10000, function(i){
#   as.data.table(do.call(rbind, lapply(CHROMs, function(CHR){
#     CHR_SIZE <- CHROM_SIZE[chr == CHR, size]
#     
#     ## Extract P-TADs for the current chromosome
#     TADS <- P_TADS[chr1 == CHR]
#     
#     ## Compute original TAD distances
#     TADS$dis <- TADS[, start2-start1]
#     
#     ## Shuffle TAD distances randomly
#     TADS$shuffle_dis <- sample(TADS$dis)
#     
#     ## PP and PE TADs are anchored to start1
#     TADS[annotation == "PP" | annotation == "PE", 5:6] <- TADS[annotation == "PP" | annotation == "PE", .(start1+shuffle_dis, end1+shuffle_dis)]
#     
#     ## EP TADs are anchored to start2
#     TADS[annotation == "EP", 2:3] <- TADS[annotation == "EP", .(start2-shuffle_dis, end2-shuffle_dis)]
#     
#     ## Remove out-of-bound shuffled TADs
#     TADS <- TADS[start1 > 0 & end2 < CHR_SIZE]
#     
#     TADS[, 1:6]  # Return only the first 6 columns
#   })))
# }, mc.cores = 40)
# 
# # Save shuffled P-TADs for future runs
# saveRDS(PERM_P_TADS, "data/permutated_P-tads.rds")

# Read shuffled P-TADs from previous run
PERM_P_TADS <- readRDS("data/permutated_P-tads.rds")

# # Annotate shuffled P-TADs
# PERM_P_TADS_ANNOTS <- lapply(PERM_P_TADS, annotate_loop)
# 
# # Save shuffled P-TAD annotations for future runs
# saveRDS(PERM_P_TADS_ANNOTS, "data/permutated_P-tads_annotations.rds") 

# Read shuffled P-TAD annotations from previous runs
PERM_P_TADS_ANNOTS <-readRDS("data/permutated_P-tads_annotations.rds")

# Compute PP-TAD proportions in shuffled data
PERM_PP_PORTION <- sapply(PERM_P_TADS_ANNOTS, function(V){sum(V=="PP")/length(V)})

## Compute p-value for PP-TAD enrichment
P <- sum(PERM_PP_PORTION > PP_PORTION)/length(PERM_PP_PORTION)

## Plot density of expected vs. observed PP-TAD proportions

DATA <- data.table("N" = c(PERM_PP_PORTION, PP_PORTION))

ggplot(DATA, aes(N)) +
  geom_density() + theme_bw() +
  geom_segment(aes(xend = PP_PORTION, yend = 0.000015, x = PP_PORTION, y = 90),
               arrow = arrow(length = unit(0.1, "cm"))) +
  labs(x = "Fruction of PP-TADs",
       y = "Density") +
  theme(text = element_text(size = 8))

# Save plot
ggsave("plots/p.density_PP_loop_density.png", height = 1, width = 2.5)
