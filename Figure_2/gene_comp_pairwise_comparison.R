
# ------------------------------------------------------------------------------
# Gene Compartment and Expression Analysis
# ------------------------------------------------------------------------------
# This script compares gene expression changes (from DESeq2) with chromatin 
# compartment shifts (eigenvector differences) across multiple cell line pairs.
#
# INPUT FILES:
# - paths.tsv : Contains paths to various datasets (cell lines, cell pairs)
# - DESeq2 differential expression results
# - Chromatin compartment eigenvector data
#
# OUTPUTS:
# - Correlation plots between gene expression and compartment shifts
# - Boxplots of compartment switching effects on expression
# - Heatmaps summarizing correlation results across cell pairs
#
# AUTHOR: Gony Shanel
# DATE: March 2025
# ------------------------------------------------------------------------------

## Set the scripts directory as working directory.
## This script is designed with relative paths for the specific data structure
## with which it was uploaded to https://github.com/ElkonLab/micro-C. 

# Arguments ---------------------------------------------------------------

CELLS <- readLines("../data/cell_lines.txt") # Load list of cell lines
CELL_PAIRS <- readLines("../data/cell_lines_pairs.txt") # Load cell line pairs
TSS_PATH <- "../data/ensembl_mart_hg38_1st_TSS_protein_coding.txt"

# Define directories for input data
DESEQ_DIR <- "../data/deseq2/paires/"
COMP_DIR <- "../data/gene_compartments/"

COMP_RES <- "100k" # Resolution of compartment analysis


# Commands ----------------------------------------------------------------
library(data.table)
library(dplyr)
library(ggplot2)
library(ggsignif)

## Load cell type specific compartments data ----

COMP_FILES <- list.files(paste0(COMP_DIR, COMP_RES), 
                         "control",
                         full.names = T)
GENE_COMPs <- lapply(COMP_FILES, fread) ; names(GENE_COMPs) <- CELLS

## Load gene expression data ----
DESEQ_FILES <- list.files(DESEQ_DIR,
                          "deseq2.res.Cells.txt",
                          full.names = T,
                          recursive = T)
DESEQ_DATA <- lapply(DESEQ_FILES, fread) ; names(DESEQ_DATA) <- CELL_PAIRS

## Analysis ----

dEIG_mRNA_FC_CORRELATION <- data.table()
COMP_SWITCHING <- data.table()

PAIR <- CELL_PAIRS[[1]] # Define a PAIR for testing purposes

for (PAIR in CELL_PAIRS){
  
  CELL1 <- strsplit(PAIR, "-")[[1]][[1]]
  CELL2 <- strsplit(PAIR, "-")[[1]][[2]]
  
  DESEQ <- DESEQ_DATA[[PAIR]][, c("V1", "log2FoldChange", "pvalue", "padj")]
  names(DESEQ)[1] <- "gene"
  
  # Compare differential expression with compartment score
  
  ## Shape data
  
  ### DESeq2 data
  DESEQ$gene <- sub("\\..*", "", DESEQ$gene)
  
  ### Eignevector data
  EIG1 <- GENE_COMPs[[CELL1]]
  EIG2 <- GENE_COMPs[[CELL2]]
  EIG_DIFF <- data.table("gene" = EIG1$gene, 
                         "eig_diff" = (EIG2$eig - EIG1$eig))
  
  ### Merge differential expression with eigenvictor data
  DT <- inner_join(DESEQ,
                   EIG_DIFF)
  DT <- DT[complete.cases(DT)]

  C <- cor.test(DT$log2FoldChange, DT$eig_diff, method = "spearman")
  dEIG_mRNA_FC_CORRELATION <- rbind(dEIG_mRNA_FC_CORRELATION, 
                                     c("cell1" = CELL1, "cell2" = CELL2, C))
  
  ## Compare differential expression to difference in compartment score
  DT <- DT[padj > quantile(DT$padj, 0.01)] # remove outliers for visualization
  ggplot(DT, aes(log2FoldChange, -log10(padj), color = eig_diff)) +
    geom_point(size = 0.5) + theme_bw() +
    labs(x = paste0("mRNA log fold change\n(", CELL2, "/", CELL1,")"),
         y = expression(-log[10](q-value)),
         color = paste0("\u0394 compartment score", 
                        "\n(", CELL2, "-",
                        CELL1, ")")) + 
    scale_color_gradient2(low = "cornflowerblue",
                          mid = "white",
                          high = "red",
                          midpoint = 0) +
    geom_vline(xintercept = c(-5,5)) +
    geom_hline(yintercept = -log10(0.05)) +
    guides(color = guide_colorbar(title.position = "right")) +
    theme(text = element_text(size = 8),
          legend.key.width=unit(0.1,"cm"),
          legend.key.height=unit(0.5,"cm"),
          legend.box.spacing = margin(0.5),
          legend.title = element_text(angle = 90, vjust = 0, hjust = 0.5),
          legend.box.just = "center") 
  ggsave(paste0(WD, "/plots/gene_comps/expression_to_comp/", PAIR, ".png"),
         height = 2, width = 2.5)
  
  ## Compare compartment score with gene expression
  DT <- data.table("gene" = EIG1$gene,
                   "rank1" = rank(EIG1$eig, ties.method = "first"),
                   "rank2" = rank(EIG2$eig, ties.method = "first"))
  DT <- left_join(DT, DESEQ[padj < 0.05])
  DT <- DT[complete.cases(DT)]
  DT <- DT[padj < 0.05]
  
  ggplot(DT, aes(rank1, rank2, color = log2FoldChange)) +
    geom_point(size = 0.5) + theme_bw() + 
    scale_color_gradient2(low = "cornflowerblue",
                          mid = "white",
                          high = "red",
                          midpoint = 0) +
    labs(x = paste0("rank(", CELL1, ")"),
         y = paste0("rank(", CELL2, ")"),
         ) +
    theme(text = element_text(size = 8),
          legend.key.width=unit(0.1,"cm"),
          legend.key.height=unit(0.5,"cm"),
          legend.box.spacing = margin(0.5),
          legend.title = element_blank()) +
    coord_cartesian(ylim = c(0,20000), xlim = c(0,20000))
  ggsave(paste0(WD, "/plots/gene_comps/comp_rank_to_expression/", PAIR, ".png"),
         height = 2, width = 2.5)
  
  ## Compare gene compartmentalization 
  COMPARTMENTALIZATION <- left_join(GENE_COMPs[[CELL1]][, c("gene", "compartment")],
                                    GENE_COMPs[[CELL2]][, c("gene", "compartment")], 
                                    by = "gene")
  COMPARTMENTALIZATION$group <- factor(with(COMPARTMENTALIZATION, paste0(compartment.x, compartment.y)),
                                       levels = c("AA", "BB", "AB", "BA"),
                                       ordered = T)
  COMPARTMENTALIZATION <- left_join(COMPARTMENTALIZATION, DESEQ[,c("gene", "log2FoldChange")])
  COMPARTMENTALIZATION <- COMPARTMENTALIZATION[complete.cases(COMPARTMENTALIZATION)]
  COMPARTMENTALIZATION$log2FoldChange <- log2(1/(2^COMPARTMENTALIZATION$log2FoldChange))
  Q <- quantile(COMPARTMENTALIZATION$log2FoldChange)
  IQR <- Q[4]-Q[2]
  C <- wilcox.test(COMPARTMENTALIZATION[group == "AB",   log2FoldChange],
              COMPARTMENTALIZATION[group == "BA", log2FoldChange],
              paired = F,
              alternative = "greater")
  
  # Aggregate data
  v.frequencies <- c("AB(expA>expB)" = nrow(COMPARTMENTALIZATION[group == "AB" & log2FoldChange > 0]),
                     "BA(expA>expB)" = nrow(COMPARTMENTALIZATION[group == "BA" & log2FoldChange > 0]),
                     "AB(expA<expB)" = nrow(COMPARTMENTALIZATION[group == "AB" & log2FoldChange < 0]),
                     "BA(expA<expB)" = nrow(COMPARTMENTALIZATION[group == "BA" & log2FoldChange < 0]))
  COMP_SWITCHING <- rbind(COMP_SWITCHING, c("cell1" = CELL1, "cell2" = CELL2, C, v.frequencies))
  
  
  # Plot
  ggplot(COMPARTMENTALIZATION, aes(group, log2FoldChange, fill = group)) +
    geom_boxplot(outlier.shape = NA) + theme_bw() +
    # scale_fill_manual(values = c("grey70", "grey50", "coral", "cyan3")) +
    scale_fill_manual(values = c("grey70", "grey50", "red", "blue")) +
    coord_cartesian(ylim = c(quantile(COMPARTMENTALIZATION$log2FoldChange, 0.001),
                             quantile(COMPARTMENTALIZATION$log2FoldChange, 0.9999))) +
    geom_signif(comparisons = list(c("AB", "BA")),
                tip_length = 0,
                annotation = formatC(C$p.value,3),
                y_position = quantile(COMPARTMENTALIZATION$log2FoldChange, 0.99),
                textsize = 2.5, 
                extend_line = -0.02) +
    labs(x = paste0("gene compartmentalization\n(", CELL1, ",", CELL2, ")"),
         y = paste0("log2 expression fold change\n(",
                    CELL1, "/", CELL2, ")")) +
    scale_x_discrete(labels = paste0(levels(COMPARTMENTALIZATION$group), 
                                     "\n(N=", 
                                     table(COMPARTMENTALIZATION$group), ")"))+
    theme(text = element_text(size = 8),
          legend.position = "none")
  ggsave(paste0(WD, "plots/gene_comps/comp_switching/p.boxplot_", PAIR,"_expFC_comp_switching.png"),
         height = 2, width = 2.5)
}

## Plot correlation summary
dEIG_mRNA_FC_CORRELATION$cell1 <- factor(dEIG_mRNA_FC_CORRELATION$cell1, 
                                         levels = unique(dEIG_mRNA_FC_CORRELATION$cell1),
                                         ordered = T)
dEIG_mRNA_FC_CORRELATION$cell2 <- factor(dEIG_mRNA_FC_CORRELATION$cell2, 
                                         levels = unique(dEIG_mRNA_FC_CORRELATION$cell2),
                                         ordered = T)
# p-value
ggplot(dEIG_mRNA_FC_CORRELATION, aes(cell1, cell2, fill = -log10(p.value))) +
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
        legend.key.width=unit(0.1,"cm"),
        legend.key.height=unit(0.3,"cm"),
        axis.text.x = element_text(angle = 45,
                                   vjust = 0.8)) +
  scale_y_discrete(limits = rev(levels(dEIG_mRNA_FC_CORRELATION$cell2))) +
  geom_text(aes(label = round(-log10(p.value))), size = 2)

ggsave(paste0(WD, "/plots/gene_comps/p.heatmap_comp-exp_correlation_pvalue_summary.png"),
       height = 2, width = 2)

# correlation coefficient
ggplot(dEIG_mRNA_FC_CORRELATION, aes(cell1, cell2, fill = estimate)) +
  geom_tile() + theme_minimal() + 
  scale_fill_distiller(palette = "Reds", 
                       direction = 1, 
                       na.value = "black",
                       limits = c(0, 1), 
                       breaks = c(0, 0.5, 1)) +
  labs(fill = "r") +
  theme(text = element_text(size = 8),
        legend.position = c(0.9, 0.7),
        panel.grid.major = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.key.width=unit(0.1,"cm"),
        legend.key.height=unit(0.3,"cm"),
        axis.text.x = element_text(angle = 45,
                                   vjust = 0.8)) +
  scale_y_discrete(limits = rev(levels(dEIG_mRNA_FC_CORRELATION$cell2))) +
  geom_text(aes(label = round(estimate, 2)), size = 1.25)

ggsave(paste0(WD, "/plots/gene_comps/p.heatmap_comp-exp_correlation_coefficient_summary.png"),
       height = 2, width = 2)

## Plot compartment switching summary
DATA <- unique(COMP_SWITCHING[, .(cell1, cell2, p.value, 
                                  "effect_size" = (`AB(expA>expB)`+`BA(expA<expB)`)/(`AB(expA>expB)`+`BA(expA>expB)`+`AB(expA<expB)`+`BA(expA<expB)`))])

# Plot p-values
ggplot(DATA, aes(cell1, cell2, fill = -log10(p.value))) +
  geom_tile() + theme_minimal() + 
  scale_fill_distiller(palette = "Reds", 
                       direction = 1, 
                       na.value = "black") +
  labs(fill = expression("-log"[10]*P)) +
  theme(text = element_text(size = 8),
        legend.position = c(0.9, 0.7),
        panel.grid.major = element_blank(),
        axis.title.x = element_blank(),
        legend.key.width=unit(0.1,"cm"),
        legend.key.height=unit(0.3,"cm"),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45,
                                   vjust = 0.8)) +
  scale_y_discrete(limits = rev(unique(DATA$cell2))) +
  geom_text(aes(label = round(-log10(p.value))), size = 2)

ggsave(paste0(WD, "/plots/gene_comps/p.heatmap_comp_switching_wilcoxon_summary.png"),
       height = 2, width = 2)

# Plot effect size

ggplot(DATA, aes(cell1, cell2, fill = effect_size*100)) +
  geom_tile() + theme_minimal() + 
  scale_fill_distiller(palette = "Reds", 
                       direction = 1, 
                       na.value = "black",
                       limits = c(0,100),
                       breaks = c(0,100)) +
  labs(fill = "effect\nsize [%]") +
  theme(text = element_text(size = 8),
        legend.position = c(0.9, 0.7),
        panel.grid.major = element_blank(),
        axis.title.x = element_blank(),
        legend.key.width=unit(0.1,"cm"),
        legend.key.height=unit(0.3,"cm"),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45,
                                   vjust = 0.8)) +
  scale_y_discrete(limits = rev(unique(DATA$cell2))) +
  geom_text(aes(label = round(effect_size*100)), size = 2)

ggsave(paste0(WD, "/plots/gene_comps/p.heatmap_comp_switching_effect_size_summary.png"),
       height = 2, width = 2)

