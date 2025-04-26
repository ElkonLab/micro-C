

library(data.table)

setwd("/specific/elkon/gonyshanel/2020-03_hsieh/hic/annotation/ab_compartments/cool/gene_compartments/100k/")


paths <- read.delim("/specific/elkon/gonyshanel/2020-03_hsieh/paths.tsv",
                    row.names = 1, stringsAsFactors = F)

samples <- readLines(paths["all_samples", "path"])
inDir <- "/specific/elkon/gonyshanel/2020-03_hsieh/hic/annotation/ab_compartments/cool/gene_compartments/100k/"

files <- mclapply(list.files(inDir, full.names = T), fread)

## Values
compValues <- as.data.table(do.call(cbind, lapply(files, function(f)f[["eig"]])))
names(compValues) <- samples

compValues <- data.table(files[[1]][,1:3], compValues, files[[1]][,4:5])

write.table(compValues, "pooled_eig.tsv",
            append = F, quote = F, sep = '\t', row.names = F, col.names = T)

## A/B compartments
comps <- as.data.table(do.call(cbind, lapply(files, function(f)f[["compartment"]])))
names(comps) <- samples

comps <- data.table(files[[1]][,1:3], comps, files[[1]][,4:5])

write.table(comps, "pooled_compatmentalization.tsv",
            append = F, quote = F, sep = '\t', row.names = F, col.names = T)
