
# pairToCells -------------------------------------------------------------
 
# Input:
#   A pair of cell-lines separated by a "-".
# 
# Output:
#   Define cell1 and cell2.
 
pairToCells <- function(pair){
  cells <- strsplit(pair, "-")[[1]]
  cell1 <<- cells[1]
  cell2 <<- cells[2]
}

# bedpe2juicebox ----------------------------------------------------------
 
# Input:
#   Coordinates of DNA-DNA interaction (2 anchors) in a vector format (eg. a row of a 6 col data.frame)
# 
# Output:
#   Coordinates in juicebox-hic format

bedpe2juicebox <- function(coordinates){
  c <- coordinates
  anchor1 <- paste0(c[1], ":", c[2], "-", c[3])
  anchor2 <- paste0(c[4], ":", c[5], "-", c[6])
  c(anchor1 = anchor1, anchor2 = anchor2)
}

# find_TSSloops -----------------------------------------------------------

# Input:
#   A data.frame containing paired interactions in it's first 6 columns
# 
# Output:
#   A data.frame containing input loops that had a TSS within 5000 bp from their anchors with gene symbol
 
find_TSSloops <- function(loops, silence=F){
  library(data.table)
  library(bedtoolsr)
  options(bedtools.path = "/specific/elkon/gonyshanel/tools/")
  # Annotate TSS associtaed loops (window = 5000).
  
  ## Import TSS file
  paths <- read.delim("/specific/elkon/gonyshanel/2020-03_hsieh/paths.tsv",
                      row.names = 1, stringsAsFactors = F)
  tss_path <- paths["TSS_file", "path"]
  
  ## import files
  tss <- fread(tss_path)
  tss[tss==""] <- NA
  ## Intersect, sort and remove duplicates
  TSSloops <- bt.sort(unique(bt.pairtobed(loops, tss)))
  names(TSSloops) <- c(names(loops), c("chr", "start", "end", "gene", "sym"))
  
  ## Write gene associated loops list
  data_to_write <- subset(TSSloops, select = -c(chr, start, end))
  
  if (silence == F){
    ## Summarize
    n_loops_in <- nrow(loops)
    n_loops_out <- nrow(unique(TSSloops[, 1:6]))
    n_genes <- length(unique(TSSloops$sym))
    
    cat(paste0("\nAnnotate gene-associated loops:\n\nInput:\n",
               n_loops_in," loops\n", "\nOuput:\n",
               n_loops_out, " loops\n", n_genes, " genes\n\n"))
  }
  data_to_write
}

# loops2Anchors -----------------------------------------------------------

# Input:
#   A data.frame containing paired interactions in it's first 6 columns
# 
# Output:
#   A data.frame containing the anchors of the input loops in .bed format

loops2Anchors <- function(loops, remove.duplicates=T, nameing = "prefix"){
  loops <- data.frame(loops, check.names = F)
  ## Split to anchors
  if (nameing == "prefix"){
    anchors_L <- cbind(loops[, -c(4:6)], "anchor_name" = paste("L", loops[,7], sep = "_"))
    anchors_R <- cbind(loops[, -c(1:3)], "anchor_name" = paste("R", loops[,7], sep = "_"))
  } else if (nameing == "suffix") {
    anchors_L <- cbind(loops[, -c(4:6)], "anchor_name" = paste(loops[,7], "L", sep = "_"))
    anchors_R <- cbind(loops[, -c(1:3)], "anchor_name" = paste(loops[,7], "R", sep = "_"))
  } else {
    stop("nameing must be either \"prefix\" or \"suffix\"")
  }
  ## Fix names
  names(anchors_R)[1:4] <- c("chr", "start", "end", "peak_name")
  names(anchors_L)[1:4] <- c("chr", "start", "end", "peak_name")
  ## Bind anchor lists
  anchors <- rbind(anchors_L, anchors_R)
  ## Remove duplicates
  if (remove.duplicates == T){
    anchors <- anchors[!duplicated(anchors[1:3]),]
  }
  anchors
}


# intersect_with_tss ------------------------------------------------------

# Input:
#   A 3 column data.frame in .bed format
#   Optional: a custum TSS file
# Output:
#   A logical vector containing information on intersection of 
#   interval with known TSS location (within 5kbp window by defult).

intersect_with_tss <- function(intervals,
                               tss.file="/specific/elkon/gonyshanel/data/public/ensembl_hg38_protein_coding_first_TSS_5000win.bed"){
  require(bedtoolsr)
  require(data.table)
  options(bedtools.path = "/specific/elkon/gonyshanel/tools/")
  
  df <- data.frame(intervals[,1:3])
  tss <- fread(tss.file)
  
  inter <- bt.intersect(a = df, b = tss, c = T)
  out <- inter[,4] != 0
}

# annotate_anchors --------------------------------------------------------

annotate_anchors <- function(anchors, mode="prom",
                             tss.file="/specific/elkon/gonyshanel/data/public/ensembl_hg38_protein_coding_first_TSS_5000win.bed"){
  # mode option: "prom", "non_prom", paths to a custom file in bed format
  df <- data.frame(anchors)
  if (mode == "prom"){
    # Promoters
    out <- intersect_with_tss(df, tss.file)
  } else if (mode == "non_prom"){
    # Non promoters
    out <- !intersect_with_tss(df, tss.file)
  } else {
    out <- intersect_with_tss(df, mode, tss.file)
  }
}
# intraTAD ----------------------------------------------------------------

# Input:
#   1. A data.frame containing paired interactions in it's first 6 columns
#   2. 
#
# Output:
#   A data.frame containing input loops that had a TSS within 5000 bp from their anchors with gene symbol

# intraTAD_loops <- function(loops, TAD.interval.path){
#   library(bedtoolsr)
#   interval <- read.delim(TAD.interval.path)
# }

# getLoopFDR -------------------------------------------------------------

# Input:
#   loop numbers
# 
# Output:
#   data frame containing loop FDR and source

getLoopFDR <- function(peaks, threshold = 0.1){
  library(data.table)
  library(parallel)
  paths = read.delim("/specific/elkon/gonyshanel/2020-03_hsieh/paths.tsv", 
                     header = T, 
                     row.names = 1,
                     stringsAsFactors = F)
  ## Load data
  loop_files <- fread(paths["mustache_5k_loops_paths", "path"], 
                      header = F,
                      stringsAsFactors = F)
  loop_index = fread(paths["mustache_5k_index", "path"], 
                     stringsAsFactors = F)
  samples <- apply(loop_files, 1, 
                   function(x) fread(x[2], stringsAsFactors = F))
  names <- loop_files$V1
  cols <- colnames(samples[[1]])[1:6]
  
  ## Container
  fdrs <- data.frame(matrix(nrow = 0, ncol = 3))
  colnames(fdrs) <- c("peak_name", "FDR", "sample")
  
  peak <- peaks[1]
  for (peak in peaks){
    if (!is.na(as.numeric(peak))){peak = paste0("peak_", peak)}
    print(peak)
    ## Shape data
    loop = loop_index[peak_name == peak,]
    colnames(loop)[1:6] <- cols
    ## 
    peak_fdrs <- mclapply(1:length(names), function(i){
      # Check if loops in sample
      loop_data <- merge(loop, samples[i])
      if (nrow(loop_data) > 0){
        row <- data.frame("peak_name" = peak,
                          "FDR" = loop_data$FDR,
                          "sample" = names[i])
        }
      }, 
      mc.cores = 20)
    fdrs <- rbind(fdrs,
                  Reduce(rbind, peak_fdrs[!sapply(fdrs, is.null)]))
    }
  
  fdrs <- fdrs[fdrs$FDR < threshold,]
}

# getGeneAssociatedLoopInterval -------------------------------------------------------------

# Input:
#   1. Gene symbol
#   2. Optional: a TSS loop list (defult: mustache_5k_0.1FDR)
# 
# Output:
  # An interval in UCSC format, containing all the loops associated with the input gene

getGeneAssociatedLoopInterval <- function(gene_symbol, margin=0, tssLoops=NA, format="UCSC"){
  # format options: "UCSC" | "bed"
  suppressMessages(require(data.table))
  if (is.character(tssLoops)){
    loops_raw <- fread(tssLoops)
  } else if (is.data.frame(tssLoops)){
    loops_raw <- tssLoops
  } else if (is.na(tssLoops)) {
    loops_raw <- fread("/specific/elkon/gonyshanel/2020-03_hsieh/hic/annotation/mustache/5000_0.1fdr/pooled_loops_values_normalized_TSSloops.tsv")
  } else {stop("Error: Unrecognzied tssLoops data format")}
  loops <- loops_raw[sym == gene_symbol]
  chr   <- min(loops$chr1)
  start <- min(loops$start1-margin)
  end   <- min(loops$end2+margin)
  if (format == "UCSC"){
    interval <- paste0(sub("chr", "", chr), ":", start, "-", end)
  } else if (format == "bed"){
    interval <- c(chr, start, end)
  } else {stop("Error: Unrecognized output format")}
}

# annotate_loops -------------------------------------------------------------

# Input:
# Data.frame with paired interactions
# 
# Output:
# A vector with the anchor annotations (promoter of distal)

annotate_loops <- function(bedpe,
                           tss.file="../../data/ensembl_hg38_protein_coding_first_TSS_5000win.bed"){
  annot1 <- intersect_with_tss(bedpe[,1:3], tss.file)
  annot1 <- ifelse(annot1 == TRUE, "P", "D")
  annot2 <- intersect_with_tss(bedpe[,4:6], tss.file)
  annot2 <- ifelse(annot2 == TRUE, "P", "D")
  paste0(annot1, annot2)
}
