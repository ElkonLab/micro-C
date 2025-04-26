# ------------------------------------------------------------------------------
# Descriptive Analysis of TADs Identified by Arrowhead at 10 kbp Resolution
# ------------------------------------------------------------------------------
# This script performs a **descriptive analysis of Topologically Associating Domains (TADs)**:
# 1. **Boundary Expansion:** Expands TAD boundaries to different sizes for comparison.
# 2. **Shared TAD Analysis:** Computes overlap between:
#    - Control samples (different cell lines)
#    - Control vs. Nutlin-treated samples (same cell line)
# 3. **Effect of Boundary Size:** Examines how varying boundary size affects results.
#
# INPUT:
# - `sample_list.txt`: List of all samples used.
# - `TAD_directories.txt`: Paths to TAD data for each sample.
# - `interval_10000_blocks.bed`: TAD boundary data per sample.
# - `ensembl_mart_hg38_1st_TSS_protein_coding.txt`: Promoter locations.
# - External tools: `bedtools`.
#
# OUTPUT:
# - Overlap statistics between TAD sets.
# - Heatmaps and bar plots showing shared TAD proportions.
# - Boxplots assessing the effect of boundary size.
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

control <- function(cell){paste0(cell, "-control")}
nutlin <- function(cell){paste0(cell, "-nutlin")}

# Arguments ---------------------------------------------------------------

samples <- readLines("../data/samples.txt")

tadDirs <- list.files("../data/TAD_interval_files/", "interval", full.names = T)
names(tadDirs) <- sub("_mer.*", "", sub(".*00/", "", tadDirs))

dir.duplicates_tads <- "../data/biological_duplicates_TADs/"
dup_samples <- list.dirs(dir.duplicates_tads, full.names = F)[-1]

boundary_sizes <- seq(10, 80, 5)

# Commands ----------------------------------------------------------------

## Set environment
library(data.table)
library(ggplot2)
library(dplyr)
library(parallel)
source("../_functions.r")

boundary_size <- boundary_sizes[1]
i <- 1
for (boundary_size in boundary_sizes){
  print(paste0("Running calculation for ", boundary_size, "kbp boundary size"))
# Generate boundary files with different boundary size

mclapply(1:length(tadDirs),
       function(i){
         inPath <- tadDirs[1]
         boundaries <- fread(inPath)
         boundaries <- cbind(boundaries, boundaries)
         
         addition <- boundary_size*1000/2
         
         boundaries[[3]] <- boundaries[[2]]+addition
         boundaries[[2]] <- boundaries[[2]]-addition
         boundaries[[5]] <- boundaries[[6]]-addition
         boundaries[[6]] <- boundaries[[6]]+addition
          
         names(boundaries) <- paste0("V", 1:6)
         boundaries <- boundaries %>%
           arrange(as.integer(sub("chr", "", V1)), V2)
         outPath <- paste0(dir, "/boundaries-", boundary_size, "kbp_10000_blocks.bedpe")
         write.tsv(boundaries, outPath)
       },
       mc.cores = 20)

#### Shared TADs between control samples ----
##### Find the proportion of overlapping TADs between pairs of samples

control_samples <- samples[grepl("control", samples)]
control_samples <- sub("-control", "", control_samples)
control_tadDirs <- tadDirs[grepl("control", tadDirs)]

TAD_interval_paths <- paste0(control_tadDirs, "/boundaries-", boundary_size, "kbp_10000_blocks.bedpe")
names(TAD_interval_paths) <- control_samples

pairs <- as.data.frame(t(combn(control_samples, 2)))

## Find the number of common TADs between pairs of samples 
overlap <- mclapply(1:nrow(pairs),
                    function(i){
                      pair <- unlist(pairs[i,])
                      fn_temp <- paste0(i, ".temp")
                      cmd <- paste0(path.bedtools,
                                    " pairtopair -a ",
                                    TAD_interval_paths[pair[1]], " -b ",
                                    TAD_interval_paths[pair[2]], " -f 0.05 > ",
                                    fn_temp)
                      system(cmd)
                      df <- fread(fn_temp)
                      unlink(fn_temp)
                      df[!duplicated(df)]
                    }, 
                    mc.cores = 45)
pairs$intersction <- sapply(overlap, nrow)

## Find the union of TADS between pairs of samples
# union = n1 + n2 - intersection
TAD_intervals <- mclapply(TAD_interval_paths, fread)
names(TAD_intervals) <- control_samples

unions <- mclapply(1:nrow(pairs),
                    function(i){
                      pair <- unlist(pairs[i,])
                      n1 <- nrow(TAD_intervals[[pair[1]]])        
                      n2 <- nrow(TAD_intervals[[pair[2]]])
                      inter <- pairs$intersction[i]
                      u <- n1 + n2 - inter
                    }, 
                    mc.cores = 45)
pairs$union <- unlist(unions)
## Calculate overlap portion
pairs$overlap <- with(pairs, intersction/union)
## Save
write.tsv(pairs,
          paste0("tables/compare_cell_lines_", boundary_size, "kbp_boundaries.tsv"), 
          header = F)

##### Plot heatmap

## Shape data
ggplot(pairs[,c(1:2,5)], aes(V1, V2, fill = overlap)) +
  geom_tile() + theme_minimal() + 
  scale_fill_distiller(palette = "Reds", 
                       direction = 1, 
                       na.value = "black", limits = c(0, 1)) +
  theme(text = element_text(size = 8),
        legend.position = c(0.9, 0.7),
        panel.grid.major = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.key.width=unit(0.1,"cm"),
        legend.key.height=unit(0.3,"cm"),
        axis.text.x = element_text(angle = 45,
                                   vjust = 0.8)) +
  scale_y_discrete(limits = rev(unique(pairs$V2))) +
  geom_text(aes(label = paste0(round(overlap*100), "%")), size = 2)

ggsave(paste0("plots/p.heatmap_overlapping_tads_", boundary_size, "k_boundaries.png"),
       height = 2.5, width = 2.5)

#### Shared TADs between control and treated samples ----

##### Find the proportion of overlapping TADs between pairs of control (c) - nutlin (n) pairs

pairs <- data.table("cell" = unique(sub("-.*", "", samples)))


## Find the number of common TADs 
overlap <- mclapply(1:nrow(pairs),
                          function(i){
                            cell <- pairs$cell[i]
                            fn_temp <- paste0(i, ".temp")
                            cmd <- paste0(path.bedtools,
                                          " pairtopair -a ",
                                          tadDirs[control(cell)], "/boundaries-", boundary_size, "kbp_10000_blocks.bedpe -b ",
                                          tadDirs[nutlin(cell)], "/boundaries-", boundary_size, "kbp_10000_blocks.bedpe -f 0.05 > ",
                                          fn_temp)
                            system(cmd)
                            df <- fread(fn_temp)
                            unlink(fn_temp)
                            df[!duplicated(df)]
                            }, 
                          mc.cores = 10)
pairs$intersction <- sapply(overlap, nrow)

## Find the union of TADS between pairs of samples
# union = Nn + Nc - intersection
TAD_intervals <- mclapply(paste0(tadDirs, "/boundaries-", boundary_size, "kbp_10000_blocks.bedpe"), fread)
names(TAD_intervals) <- names(tadDirs)

unions <- mclapply(1:nrow(pairs),
                    function(i){
                      cell <- pairs$cell[i]
                      Nn <- nrow(as.data.table(TAD_intervals[control(cell)]))
                      Nc <- nrow(as.data.table(TAD_intervals[nutlin(cell)]))
                      inter <- pairs$intersction[i]
                      u <- Nn + Nc - inter
                    }, 
                    mc.cores = 10)
pairs$union <- unlist(unions)
## Calculate overlap portion
pairs$overlap <- with(pairs, intersction/union)
## Save
write.tsv(pairs,
          paste0("tables/compare_control_nutlin_", boundary_size, "kbp_boundaries.tsv"), 
          header = F)

##### Plot barplot

ggplot(pairs, aes(cell, overlap, fill = overlap)) +
  geom_col() + theme_bw() +
  scale_fill_distiller(palette = "Reds", 
                       direction = 1, 
                       na.value = "black", limits = c(0, 1)) +
  labs(title = "Shared TADs in control-nutlin pairs",
       y = "overlap [intersection/union]") +
  theme(text = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.8),
        legend.position = "none") 
ggsave(paste0("plots/p.barplot_overlapping_TADs_c-t_", boundary_size, "_boundary.png"),
       height = 2, width = 2.5)


#### Shared TADs between biological duplicates

##### Find the proportion of overlapping TADs between pairs of biological duplicates

rep1 <- function(sample){paste0(sample, "_rep1")}
rep2 <- function(sample){paste0(sample, "_rep2")}

dup_samples <- data.table(sample = grep("IMR90|SKNSH", samples, invert = T, value = T))
dup_names <- grep("IMR90|SKNSH", list.dirs(dir.duplicates_tads, full.names = F), invert = T, value = T)[-1]

## Find the number of common TADs 
overlap <- mclapply(1:nrow(dup_samples),
                    function(i){
                      sample <- dup_samples[i, 1]
                      fn_temp <- paste0(i, ".temp")
                      cmd <- paste0(path.bedtools,
                                    " pairtopair -a ",
                                    dir.duplicates_tads, sample, "_rep1/boundaries-", boundary_size, "kbp_10000_blocks.bedpe -b ",
                                    dir.duplicates_tads, sample, "_rep2/boundaries-", boundary_size, "kbp_10000_blocks.bedpe -f 0.05 > ",
                                    fn_temp)
                      system(cmd)
                      df <- fread(fn_temp)
                      unlink(fn_temp)
                      df[!duplicated(df)]
                    }, 
                    mc.cores = 10)
dup_samples$intersction <- sapply(overlap, nrow)

## Find the union of TADS between pairs of samples
# union = N1 + N2 - intersection
dup_tad_dirs <- grep("IMR90|SKNSH", list.dirs(dir.duplicates_tads, full.names = T), invert = T, value = T)[-1]
TAD_intervals <- mclapply(paste0(dup_tad_dirs, "/boundaries-", boundary_size, "kbp_10000_blocks.bedpe"), fread)
names(TAD_intervals) <- dup_names

unions <- mclapply(1:nrow(dup_samples),
                   function(i){
                     s <- dup_samples[i, 1]
                     N1 <- nrow(as.data.table(TAD_intervals[paste0(s, "_rep1")]))
                     N2 <- nrow(as.data.table(TAD_intervals[paste0(s, "_rep2")]))
                     inter <- dup_samples$intersction[i]
                     u <- N1 + N2 - inter
                   }, 
                   mc.cores = 10)
dup_samples$union <- unlist(unions)
## Calculate overlap portion
dup_samples$overlap <- with(dup_samples, intersction/union)
## Save
write.tsv(dup_samples,
          paste0("tables/compare_biological_duplicates_", boundary_size, "kbp_boundaries.tsv"), 
          header = F)

##### Plot barplot

ggplot(dup_samples, aes(sample, overlap, fill = overlap)) +
  geom_col() + theme_bw() +
  scale_fill_distiller(palette = "Reds", 
                       direction = 1, 
                       na.value = "black", limits = c(0, 1)) +
  labs(title = "Shared TADs in control-nutlin pairs",
       y = "overlap [intersection/union]") +
  theme(text = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
        legend.position = "none") 
ggsave(paste0("plots/p.barplot_overlapping_TADs_between_biological_duplicates_", boundary_size, "_boundary.png"),
       height = 2.5, width = 3)

}

#### Examine the effect of boundary size on the results

## Aggregate data
overlap <- data.table()

for (i in seq(10, 80, 10)){
  print(i)
  # read cell-pair data
  inPath <- paste0("tables/compare_cell_lines_", i, "kbp_boundaries.tsv")
  in_df <- fread(inPath)
  out_df <- data.table("boundary_size" = as.character(i),
                       "overlap" = in_df$V5,
                       "comparison group" = "cell-pairs")
  overlap <- rbind(overlap, out_df)
  # read control-nutlin data
  inPath <- paste0("tables/compare_control_nutlin_", i, "kbp_boundaries.tsv")
  in_df <- fread(inPath)
  out_df <- data.table("boundary_size" = as.character(i),
                       "overlap" = in_df$V4,
                       "comparison group" = "control-nutlin")
  overlap <- rbind(overlap, out_df)
  # read biological duplicates data
  inPath <- paste0("tables/compare_biological_duplicates_", i, "kbp_boundaries.tsv")
  in_df <- fread(inPath)
  out_df <- data.table("boundary_size" = as.character(i),
                       "overlap" = in_df$V4,
                       "comparison group" = "biological duplicates")
  overlap <- rbind(overlap, out_df)
}
## Shape
overlap$boundary_size <- factor(overlap$boundary_size,
                                levels = unique(overlap$boundary_size),
                                ordered = T)
overlap$`comparison group` <- factor(overlap$`comparison group`,
                                     levels = c("cell-pairs", "control-nutlin", "biological duplicates"),
                                     labels = c("different\ncell-lines",
                                                "control vs.\nnutlin",
                                                "biological\nduplicates"),
                                     ordered = T)
## Plot
ggplot(overlap, aes(boundary_size, overlap, fill = `comparison group`)) +
  geom_boxplot(outlier.shape = NA) + theme_bw() +
  # geom_smooth(method = "loess", aes(group=`comparison group`), size = 0.5, color="black") +
  labs(x = "Boundary size [kbp]",
       y = "Overlap [intersection / union]",
       fill = "comparison") +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  theme(text = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.title = element_blank(),
        legend.direction="horizontal",
        legend.position = c(0.6, 0.1),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.background = element_rect(linetype="solid", 
                                         colour ="black", 
                                         size = 0.2))
## Save
ggsave("plots/p.boxplot_overlap_by_boundary_size.png",
       height = 2.5, width = 3.5)

## Using TAD interval overlap instead of TAD boundaries overlap ----
# Compare TAD overlap between cell lines to TAD overlap between control and nutlin-treated samples using TAD interval overlap instead of TAD boundaries overlap

## Create files with overlapping tads merged
pairs <- data.table("cell" = unique(sub("-.*", "", samples)))
mclapply(1:nrow(pairs),
         function(i){
           cell <- pairs$cell[i]
           cmd <- paste0(path.bedtools,
                         " sort -i ",
                         tadDirs[control(cell)], 
                         "/interval_10000_blocks.bed | ",
                         path.bedtools,
                         " merge -i stdin >",
                         tadDirs[control(cell)], "/interval_10000_blocks.merged.bed")
           system(cmd)
           
           cmd <- paste0(path.bedtools,
                         " sort -i ",
                         tadDirs[nutlin(cell)], 
                         "/interval_10000_blocks.bed | ",
                         path.bedtools,
                         " merge -i stdin >",
                         tadDirs[nutlin(cell)], "/interval_10000_blocks.merged.bed")
           system(cmd)
         }, 
         mc.cores = 40)

### Compare cell lines ----

control_samples <- samples[grepl("control", samples)]
control_samples <- sub("-control", "", control_samples)
control_tadDirs <- tadDirs[grepl("control", tadDirs)]

TAD_interval_paths <- paste0(control_tadDirs, "/interval_10000_blocks.merged.bed")
names(TAD_interval_paths) <- control_samples

pairs <- as.data.frame(t(combn(control_samples, 2)))

## Find the number of common TADs between pairs of samples 
i=1
overlap <- mclapply(1:nrow(pairs),
                    function(i){
                      pair <- unlist(pairs[i,])
                      fn_temp <- paste0(i, ".temp")
                      cmd <- paste0(path.bedtools,
                                    " intersect -a ",
                                    TAD_interval_paths[pair[1]], " -b ",
                                    TAD_interval_paths[pair[2]], " -wo > ",
                                    fn_temp)
                      system(cmd)
                      df <- fread(fn_temp)
                      unlink(fn_temp)
                      df[!duplicated(df)]
                    }, 
                    mc.cores = 45)
pairs$intersction <- sapply(overlap, function(df){sum(df$V7)})

## Find the union of TADS between pairs of samples
# union = n1 + n2 - intersection

TAD_intervals <- mclapply(TAD_interval_paths, fread)
names(TAD_intervals) <- control_samples

unions <- mclapply(1:nrow(pairs),
                   function(i){
                     pair <- unlist(pairs[i,])
                     N1 <- sum(as.data.table(TAD_intervals[[pair[1]]])[, .SD[[3]] - .SD[[2]]])
                     N2 <- sum(as.data.table(TAD_intervals[[pair[2]]])[, .SD[[3]] - .SD[[2]]])
                     inter <- pairs$intersction[i]
                     u <- N1 - inter + N2 
                   }, 
                   mc.cores = 45)
pairs$union <- unlist(unions)
## Calculate overlap portion
pairs$overlap <- with(pairs, intersction/union)
## Save
write.tsv(pairs,
          paste0("tables/compare_cell_lines_by_intervals_overlap.tsv"), 
          header = F)

### Control vs Nutlin ----

cell_lines <- data.table("cell" = unique(sub("-.*", "", samples)))

control <- function(cell){paste0(cell, "-control")}
nutlin <- function(cell){paste0(cell, "-nutlin")}

## Find the number of common TADs 
i=1
overlap <- mclapply(1:nrow(cell_lines),
                    function(i){
                      cell <- cell_lines$cell[i]
                      fn_temp <- paste0(i, ".temp")
                      cmd <- paste0(path.bedtools,
                                    " intersect -a ",
                                    tadDirs[control(cell)], "/interval_10000_blocks.merged.bed -b ",
                                    tadDirs[nutlin(cell)], "/interval_10000_blocks.merged.bed -wo > ",
                                    fn_temp)
                      system(cmd)
                      df <- fread(fn_temp)
                      unlink(fn_temp)
                      df[!duplicated(df)]
                    }, 
                    mc.cores = 40)
cell_lines$intersction <- sapply(overlap, function(df){sum(df$V7)})

## Find the union of TADS between pairs of samples
# union = Nn + Nc - intersection
TAD_intervals <- mclapply(paste0(tadDirs, "/interval_10000_blocks.merged.bed"), fread)
names(TAD_intervals) <- names(tadDirs)

unions <- mclapply(1:nrow(cell_lines),
                   function(i){
                     cell <- cell_lines$cell[i]
                     Nn <- sum(as.data.table(TAD_intervals[control(cell)])[, .SD[[3]] - .SD[[2]]])
                     Nc <- sum(as.data.table(TAD_intervals[nutlin(cell)])[, .SD[[3]] - .SD[[2]]])
                     inter <- cell_lines$intersction[i]
                     u <- Nn - inter + Nc 
                   }, 
                   mc.cores = 40)
cell_lines$union <- unlist(unions)
## Calculate overlap portion
cell_lines$overlap <- with(cell_lines, intersction/union)
## Save
write.tsv(cell_lines,
          paste0("tables/compare_control_nutlin_by_intervals_overlap.tsv"), 
          header = F)


### Biological duplicates ----

dup_samples <- data.table(sample = grep("IMR90|SKNSH", samples, invert = T, value = T))
dup_names <- grep("IMR90|SKNSH", list.dirs(dir.duplicates_tads, full.names = F), invert = T, value = T)[-1]

## Create files with overlapping tads merged
s=dup_names[1]
mclapply(dup_names,
         function(s){
           cell <- dup_names[i]
           cmd <- paste0(path.bedtools,
                         " sort -i /specific/elkon/gonyshanel/2020-03_hsieh/hic/annotation/juicer/duplicates/arrowhead/",
                         s, 
                         "/interval_10000_blocks.bed | ",
                         path.bedtools,
                         " merge -i stdin > /specific/elkon/gonyshanel/2020-03_hsieh/hic/annotation/juicer/duplicates/arrowhead/",
                         s, 
                         "/interval_10000_blocks.merged.bed")
           system(cmd)
         }, 
         mc.cores = 40)

## Find the number of common TADs 
overlap <- mclapply(1:nrow(dup_samples),
                    function(i){
                      sample <- dup_samples[i, 1]
                      fn_temp <- paste0(i, ".temp")
                      cmd <- paste0(path.bedtools,
                                    " intersect -a /specific/elkon/gonyshanel/2020-03_hsieh/hic/annotation/juicer/duplicates/arrowhead/",
                                    sample, "_rep1/interval_10000_blocks.merged.bed -b /specific/elkon/gonyshanel/2020-03_hsieh/hic/annotation/juicer/duplicates/arrowhead/",
                                    sample, "_rep2/interval_10000_blocks.merged.bed -wo > ",
                                    fn_temp)
                      system(cmd)
                      df <- fread(fn_temp)
                      unlink(fn_temp)
                      df[!duplicated(df)]
                    }, 
                    mc.cores = 40)

dup_samples$intersction <- sapply(overlap, function(df){sum(df$V7)})

## Find the union of TADS between pairs of samples
# union = N1 + N2 - intersection
dup_tad_dirs <- grep("IMR90|SKNSH", list.dirs(dir.duplicates_tads, full.names = T), invert = T, value = T)[-1]
TAD_intervals <- mclapply(paste0(dup_tad_dirs, "/interval_10000_blocks.merged.bed"), fread)
names(TAD_intervals) <- dup_names

unions <- mclapply(1:nrow(dup_samples),
                   function(i){
                     s <- dup_samples[i, 1]
                     N1 <- sum(as.data.table(TAD_intervals[paste0(s, "_rep1")])[, .SD[[3]] - .SD[[2]]])
                     N2 <- sum(as.data.table(TAD_intervals[paste0(s, "_rep2")])[, .SD[[3]] - .SD[[2]]])
                     inter <- dup_samples$intersction[i]
                     u <- N1 - inter + N2
                   }, 
                   mc.cores = 40)
dup_samples$union <- unlist(unions)
## Calculate overlap portion
dup_samples$overlap <- with(dup_samples, intersction/union)
## Save
write.tsv(dup_samples,
          paste0("tables/compare_biological_duplicates_by_intervals_overlap.tsv"), 
          header = F)


## Compare TAD conservation between bio duplicates to conservation between control vs. nitlin and between cell lines (by control sample) ----

data_to_plot <- rbind(data.table(type = "cell-lines", overlap = pairs$overlap),
                      data.table(type = "control-nutlin", overlap = cell_lines$overlap),
                      data.table(type = "biological-duplicates", overlap = dup_samples$overlap))
data_to_plot$type <- factor(data_to_plot$type, 
                            levels = c("biological-duplicates", "cell-lines", "control-nutlin"),
                            labels = c("biological\nduplicates", "cell\nlines", "control\nnutlin-3a"),
                            ordered = T)


ggplot(data_to_plot, aes(type, overlap, fill = type)) +
  geom_boxplot() + theme_bw() +
  labs(x = "Comparison",
       y = "Overlap [intersection/union]") +
  scale_fill_brewer(palette = "Set1") +
  theme(text = element_text(size = 8),
        # axis.title.x = element_blank(),
        legend.position = "none")
  
ggsave("plots/p.boxplot_overlap_by_TAD_intervals.png",
       width = 2.5, height = 2)

data_to_plot[, mean(overlap), by = type]
