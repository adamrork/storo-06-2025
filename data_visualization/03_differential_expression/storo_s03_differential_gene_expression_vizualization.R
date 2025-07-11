###############
### PURPOSE ###
###############

# This script generates a set of figures used to ensure #
# that the differential gene expression analysis was #
# properly conducted. Sets of figures to inspect the #
# outcomes of each pairwise contrast are also generated. #

######################
### INITIALIZATION ###
######################

# Load libraries #
library(tidyverse)
library(pheatmap)
library(ggpubr)
library(EnhancedVolcano)
library(RColorBrewer)
library(DESeq2)

# Set working directory #
setwd("~/Desktop/general/projects/clients/2025/002_steven_toro/")

# Load functions, contrast list, and differential gene expression analysis data #
source("storo-06-2025/functions/functions.R")
load("Data/storo_rdata/storo_s03_differential_gene_expression.RData")

# Set seed for reproducibility #
set.seed(123)

#########################
### FIGURE GENERATION ###
#########################

# Create a data frame of library size normalization factors #
allGroups <- gsub("_[1-3]$", "", names(dds@colData$sizeFactor))
normFactors <- data.frame(Names = names(dds@colData$sizeFactor),
                          Factors = dds@colData$sizeFactor,
                          Groups = allGroups)

# Create a dot plot of library size normalization factors #
colnames(normFactors)[1] <- "Samples"
nfPlot <- ggplot(normFactors, aes(x = Factors, y = Samples, fill = Groups)) +
            geom_dotplot() +
            coord_flip() +
            ggtitle("Normalization Factors") +
            theme(plot.title = element_text(hjust = 0.5, size = 15),
                  axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot to a PDF #
ggsave("Normalization Factors (DESeq2).pdf", plot = nfPlot, device = "pdf",
       path = "Figures/differential_expression/quality_assessment/",
       width = 8.5, height = 7, units = "in", dpi = 400)

# Create a gene-wise dispersion plot and save it to a PDF #
pdf("Figures/differential_expression/quality_assessment/Gene-wise Dispersion Estimates.pdf",
    width = 8.5, height = 7)

  plotDispEsts(dds)
    title("Gene-wise Dispersion Estimates")

dev.off()

#############
# PCA PLOTS #
#############

# Create a PCA plot annotated by timepoint and treatment (based on expression of all genes) #
pcaAll <- ggplot(pcaData, aes(x = PC1, y = PC2, color = Day, shape = Treatment)) +
            geom_point(size = 4) +
            scale_color_brewer(palette = "Dark2") +
            scale_shape_manual(values = c(18, 15, 17, 19)) +
            ggtitle("Treatments & Timepoints PCA (All Genes)") +
            xlab("PC1: 58% Variance") +
            ylab("PC1: 17% Variance") +
            theme_bw() +
            theme(plot.title = element_text(size = 17, hjust = 0.5),
                  axis.text = element_text(size = 8),
                  axis.title = element_text(size = 12))

# Create a PCA plot annotated by timepoint and treatment (based on expression of de genes) #
de.pcaAll <- ggplot(de.pcaData, aes(x = PC1, y = PC2, color = Day, shape = Treatment)) +
               geom_point(size = 4) +
               scale_color_brewer(palette = "Dark2") +
               scale_shape_manual(values = c(18, 15, 17, 19)) +
               ggtitle("Treatments & Timepoints PCA (DE Genes)") +
               xlab("PC1: 60% Variance") +
               ylab("PC1: 17% Variance") +
               theme_bw() +
               theme(plot.title = element_text(size = 17, hjust = 0.5),
                     axis.text = element_text(size = 8),
                     axis.title = element_text(size = 12))

# Save the "all genes" PCA plot to a PDF #
ggsave("PCA Plot (All Genes).pdf", plot = pcaAll, device = "pdf",
       path = "Figures/differential_expression/pca_plots/",
       width = 8.5, height = 7, units = "in", dpi = 400)

# Save the "de genes" PCA plot to a PDF #
ggsave("PCA Plot (DE Genes).pdf", plot = de.pcaAll, device = "pdf",
       path = "Figures/differential_expression/pca_plots/",
       width = 8.5, height = 7, units = "in", dpi = 400)

############
# MA PLOTS #
############

# Create a list to store our MA plots #
maList <- list()

# For each contrasts... #
for (ctr in seq_along(geneList)) {

  # Assign the data frame to a temporary variable #
  tmp <- geneList[[ctr]]

  # Create a cleaned up name from the original contrast using our custom rename function #
  tmp.name <- renameContrast(name = names(geneList)[[ctr]], addition = "MA Plot") %>% gsub("[ ]+", " ", .)

  # Dynamically generate y-min and y-max settings #
  ymin <- round(min(tmp$log2FoldChange) - 0.5)
  ymax <- round(max(tmp$log2FoldChange) + 0.5)

  # Create the MA plots #
  maList[[ctr]] <- ggmaplot(tmp,
                            top = 0,
                            fdr = 0.05, fc = 2.0,
                            size = 1) +
                     rremove("legend") +
                     ggtitle(tmp.name) +
                     xlab("log2(Mean Expression)") +
                     ylab("log2(Fold Change)") +
                     ylim(ymin, ymax) +
                     theme(plot.title = element_text(size = 20, hjust = 0.5),
                           axis.title = element_text(size = 16),
                           axis.text = element_text(size = 14))

  # Create a new filename after the plot title #
  filename <- paste0(tmp.name, ".pdf")

  # Save the plot to a PDF #
  ggsave(filename, plot = maList[[ctr]], device = "pdf",
         path = "Figures/differential_expression/ma_plots/",
         width = 8.5, height = 7, units = "in", dpi = 400)
}

##################
# VOLCANO PLOTS  #
##################

# Create a list to store our volcano plots #
voList <- list()

# For each contrasts... #
for (ctr in seq_along(geneList)) {

  # Assign the data frame to a temporary variable #
  tmp <- geneList[[ctr]]

  # Create a cleaned up name from the original contrast using our custom rename function #
  tmp.name <- renameContrast(name = names(geneList)[[ctr]], addition = "Volcano Plot") %>% gsub("[ ]+", " ", .)

  # Create the volcano plots #
  voList[[ctr]] <- EnhancedVolcano(tmp,
                                   lab = rownames(tmp),
                                   x = 'log2FoldChange',
                                   y = 'padj',
                                   pCutoff = 0.05,
                                   FCcutoff = 2,
                                   col = c("#999999", "#A6CEE3", "#FB9A99", "#E41A1C"),
                                   legendPosition = "right",
                                   legendLabSize = 10,
                                   legendIconSize = 3.0,
                                   legendLabels = c("| logFC | < 2 & FDR > 0.05",
                                                    "| logFC | > 2 & FDR > 0.05",
                                                    "| logFC | < 2 & FDR < 0.05",
                                                    "| logFC | > 2 & FDR < 0.05"),
                                   cutoffLineType = "dotted",
                                   xlab = "log2(Fold Change)",
                                   ylab = "-log10(FDR)",
                                   title = NULL,
                                   subtitle = NULL,
                                   caption = NULL,
                                   labSize = 0,
                                   borderWidth = 0.5,
                                   ) +
                     ggtitle(tmp.name) +
                     labs(color = "Significance") +
                     theme_bw() +
                     theme(plot.title = element_text(size = 16, face = "plain", hjust = 0.5),
                           axis.title = element_text(size = 15),
                           axis.ticks = element_line(linewidth = 0.5),
                           axis.text.x = element_text(size = 10),
                           axis.text.y = element_text(size = 10))

  # Create a new filename after the plot title #
  filename <- paste0(tmp.name, ".pdf")

  # Save the plot to a PDF #
  ggsave(filename, plot = voList[[ctr]], device = "pdf",
         path = "Figures/differential_expression/volcano_plots/",
         width = 8.5, height = 7,  units = "in", dpi = 400)
}

#############
# HEAT MAPS #
#############

# Gather N random genes and the N genes with the highest and lowest variance (N = size of DEG superset) #
sampleSize <- length(degVector)

random_genes <- sample(rownames(nrmDDS))[1:sampleSize]
lowvar_genes <- rownames(assay(nrmDDS)[order(apply(assay(nrmDDS), 1, var), decreasing = F)[1:sampleSize],])
hivar_genes <- rownames(assay(nrmDDS)[order(apply(assay(nrmDDS), 1, var), decreasing = T)[1:sampleSize],])

# Generate a heat map using normalized counts from N randomly selected genes #
randHM <- pheatmap(assay(nrmDDS)[random_genes,],
                   cluster_rows = TRUE, cluster_cols = TRUE,
                   show_rownames = FALSE, show_colnames = TRUE,
                   annotation_col = groupDF)

# Generate a heat map using normalized counts from the N genes with the highest variance #
hvarHM <- pheatmap(assay(nrmDDS)[lowvar_genes,],
                   cluster_rows = TRUE, cluster_cols = TRUE,
                   show_rownames = FALSE, show_colnames = TRUE,
                   annotation_col = groupDF)

# Generate a heat map using normalized counts from the N genes with the lowest variance #
lvarHM <- pheatmap(assay(nrmDDS)[hivar_genes,],
                   cluster_rows = TRUE, cluster_cols = TRUE,
                   show_rownames = FALSE, show_colnames = TRUE,
                   annotation_col = groupDF)

# Generate a heat map using normalized counts from the N genes which are DE in at least one contrast #
degsHM <- pheatmap(assay(nrmDDS)[degVector,],
                   cluster_rows = TRUE, cluster_cols = TRUE,
                   show_rownames = FALSE, show_colnames = TRUE,
                   annotation_col = groupDF)

# We only want to save this plot - the others were for demonstration and exploration purposes #
ggsave("Heatmap of DE Genes.pdf", plot = degsHM, device = "pdf",
       path = "Figures/differential_expression/heatmaps/",
       width = 8.5, height = 7, units = "in", dpi = 400)

# Clear the environment for the next project #
rm(list = ls())
