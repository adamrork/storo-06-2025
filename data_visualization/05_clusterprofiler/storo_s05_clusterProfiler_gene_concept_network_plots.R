###############
### PURPOSE ###
###############

# This script creates several gene concept network plots from #
# the results of clusterProfiler functional enrichment analyses. #

######################
### INITIALIZATION ###
######################

# Load libraries #
library(tidyverse)
library(ggnewscale)
library(enrichplot)
library(DOSE)
library(org.Mm.eg.db)

# Set working directory #
setwd("~/Desktop/general/projects/clients/2025/002_steven_toro/")

# Import some plotting themes and functions #
source("Scripts/functions/functions.R")
source("Data/storo_lists/lists_and_vectors.R")
load("Data/storo_rdata/storo_s05_clusterProfiler_enrichment.RData")

# Set seed for reproducibility #
set.seed(123)

##############################
# GENE CONCEPT NETWORK PLOTS #
##############################

for ( ctr in deContrasts ) {

  for ( strategy in annotStrategies ) {

    for ( reg in regDirection ) {

      # Assign a data frame to a temporary variable #
      if ( strategy == "goEnr" ) {

        enrich.er <- goEnr[[reg]][[ctr]]
        method <- "Enr"
        type <- "GO"
        filepath <- "Figures/clusterprofiler/gene_concept_networks/enrichment/go/"

      } else if ( strategy == "goGSEnr" ) {

        enrich.er <- goGSEnr[[ctr]]
        method <- "GSEnr"
        type <- "GO"
        filepath <- "Figures/clusterprofiler/gene_concept_networks/gsea/go/"

      } else {
        next
      }

      # Only parse GSEA data once #
      if ( method == "GSEnr" && reg == "DOWN" ) {
        next
      }

      # Move on if there is no data to analyze #
      if ( !dim(enrich.er)[1] > 0 || is.list(enrich.er) ) {
        next
      }

      # Dynamically create some plot titles #
      if ( reg == "UP" ) {
        tmp.title <- renameContrast(name = ctr, addition = " - Top Enriched GO Terms", flip = FALSE) %>% gsub("[ ]+", " ", .)

      } else if ( reg == "DOWN" ) {
        tmp.title <- renameContrast(name = ctr, addition = " - Top Enriched GO Terms", flip = TRUE) %>% gsub("[ ]+", " ", .)
      }

      # Create GO Gene Set Enrichment Gene Concept Network Plots #
      gcn.plot <- setReadable(enrich.er, org.Mm.eg.db) %>%
        cnetplot(foldChange = geneVector.lst[[ctr]],
                showCategory = 5,
                node_label = "all", layout = "layout_nicely",
                color_category = "#4DAF4A", color_gene = "") +
        ggtitle(tmp.title) +
        theme(plot.title = element_text(hjust = 0.5, size = 25)) +
        labs(size = "No. Genes") +
        scale_color_continuous(name = "logFC", low = "#FFDBD3", high = "#D1001F")

      # Create file names #
      if ( method == "Enr" ) {
        filename <- paste0(tmp.title, " (Enr) Gene Concept Network.pdf")
      } else if ( method == "GSEnr" ) {
        filename <- paste0(tmp.title, " (GSEA) Gene Concept Network.pdf")
      }

      # Save the plot to a PDF #
      ggsave(filename, plot = gcn.plot, device = "pdf",
                 path = filepath,
                 width = 14, height = 12, units = "in", dpi = 400)
    }
  }
}
