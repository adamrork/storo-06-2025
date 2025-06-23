###############
### PURPOSE ###
###############

# This script creates several combined ridgeline plots, #
# box plots and gene expression heatmaps from the results #
# of clusterProfiler functional enrichment and other analyses. #

######################
### INITIALIZATION ###
######################

# Load libraries #
library(tidyverse)
library(ggridges)
library(RColorBrewer)
library(DOSE)
library(ComplexHeatmap)
library(cowplot)

# Set working directory #
setwd("~/Desktop/general/projects/clients/2025/002_steven_toro/")

# Import some plotting themes and functions #
source("Scripts/functions/functions.R")
source("Data/storo_lists/lists_and_vectors.R")
load("Data/storo_rdata/storo_s05_clusterProfiler_enrichment.RData")

# Set seed for reproducibility #
set.seed(123)

#############################################
### DEFINE TERMS AND PATHWAYS OF INTEREST ###
#############################################

# Define all GO terms of interest (comma-delimited; spelling and case-sensitive) #
terms <- list("intermediate filament organization",
              "cell chemotaxis",
              "leukocyte migration")

# Uncomment the line your terms belong to. By default, we will use GO terms. #
type <- "GO"
#type <- "KEGG"

#################################
# RIDGELINE + BOXPLOT + HEATMAP #
#################################

# Create an empty data frame to append results to #
rdg.df <- data.frame()

# For each GO description... #
for (term in terms) {

  # For each contrast ... #
  for (ctr in deContrasts) {

    if ( type == "GO" ) {
      enrich.er <- goGSEnr[[ctr]]
    } else if ( type == "KEGG" ) {
      enrich.er <- kgGSEnr[[ctr]]
    }

    # If the description can be found in said contrasts' GSEA results object #
    if ( term %in% enrich.er@result$Description ) {

      # Obtain the relevant GO Category ID #
      term.ids <- filter(enrich.er, Description == term)$ID

      # Gather all genes associated with said ID #
      gene.ids <- geneInCategory(enrich.er)[term.ids][[1]]

      # Gather all genes and their logFC values based on the aforementnoned set of genes #
      genes.logfc <- enrich.er@geneList[names(enrich.er@geneList) %in% gene.ids]

      # Create a data frame from the above data and nullify the rows #
      tmp.df <- data.frame(ENTREZID = gene.ids, CONTRAST = ctr, CATEGORY = term.ids, log2FC = genes.logfc)
      rownames(tmp.df) <- NULL

      # Append to data frame rdg.df #
      rdg.df <- rbind(rdg.df, tmp.df)
    }
  }

  # Add some tidy axis labels to the data frame #
  rdg.df$Contrast <- gsub("_", " ", rdg.df$CONTRAST) %>%
    gsub("vs", "\nvs.", .) %>%
    gsub("PRG4.VEGF", "PRG4-VEGF", .) %>%
    gsub("D0 D00", "Day 0", .)

  # Create a data frame with expression data for genes associated with enriched terms #
  go2expr.df <- left_join(rdg.df, symb2expr, by = "ENTREZID")
  go2expr.df <- go2expr.df[!duplicated(go2expr.df[,"SYMBOL"]),]

  # Examine all genes associated with our particular term. #
  heatmap.genes <- go2expr.df[which(go2expr.df$ENTREZID %in% gene.ids),]

  # Subset the data frame to contain only gene symbols and expression data #
  heatmap.mt <- as.matrix(heatmap.genes[,c(8:28)])
  rownames(heatmap.mt) <- heatmap.genes$SYMBOL

  # Create a color palette for the heatmap #
  hm.colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)

  # Create a ridgeline plot #
  nes.rdg.plot <- ggplot(rdg.df, aes(x = log2FC, y = Contrast, fill = Contrast)) +
    geom_density_ridges(alpha = 0.7) +
    theme(plot.title = element_text(size = 17, hjust = 0.5),
          axis.title = element_text(size = 19),
          axis.text = element_text(size = 18),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.position = "bottom") +
    guides(fill = guide_legend(nrow = 3),
           fill = "none", color = "none")

  # Create a box or violin plot #
  nes.bov.plot <- ggplot(rdg.df, aes(x = log2FC, y = Contrast, fill = Contrast)) +
    geom_boxplot(width = 0.5, alpha = 1, aes(fill = Contrast)) + # switch between geom_boxplot and geom_violin
    theme_bw() +
    theme(plot.title = element_text(size = 17, hjust = 0.5),
          axis.title.x = element_text(size = 19),
          axis.text.x = element_text(size = 18),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    guides(fill = "none", color = "none")

  # Create a heat map #
  nes.hm.plot <- ComplexHeatmap::pheatmap(heatmap.mt, color = hm.colors,
                                          cluster_rows = FALSE, cluster_cols = FALSE,
                                          show_rownames = TRUE, show_colnames = TRUE,
                                          fontsize = 11.5, angle_col = "315",
                                          heatmap_legend_param = list(title = "log2(CPM)"))

  # Dynamically create some plot titles #
  tmp.title <- renameContrast(name = term,
                              addition = " - logFC Distributions and Expression of Associated Genes",
                              flip = FALSE) %>% gsub("[ ]+", " ", .)

  plot_title <- ggdraw() + draw_label(tmp.title, size = 30)

  # Display the ridgeline plot, box plot, and heatmap plot side-by-side #
  nes.plot <- plot_grid(plot_title,
                          plot_grid(nes.rdg.plot, nes.bov.plot, grid.grabExpr(draw(nes.hm.plot)), ncol = 3),
                        nrow = 2, rel_heights=c(0.1, 1))

  # Create file names #
  filename <- paste0(tmp.title, " (GSEA).pdf")

  # Save the plot to a PDF #
  ggsave(filename, plot = nes.plot, device = "pdf",
         path = "Figures/clusterprofiler/logfc_ridgeline_bar_heatmap_plots/",
         width = 28, height = 12, units = "in", dpi = 400)

  # Clear the "rdf.df" data frame
  rdg.df <- data.frame()
}
