###############
### PURPOSE ###
###############

# This script creates several running score vs. ranked list plots #
# from the results of clusterProfiler functional enrichment analyses. #

######################
### INITIALIZATION ###
######################

# Load libraries #
library(tidyverse)
library(RColorBrewer)
library(enrichplot)

# Set working directory #
setwd("~/Desktop/general/projects/clients/2025/002_steven_toro/")

# Import some plotting themes and functions #
source("Scripts/functions/functions.R")
source("Data/storo_lists/lists_and_vectors.R")
load("Data/storo_rdata/storo_s05_clusterProfiler_enrichment.RData")

# Set seed for reproducibility #
set.seed(123)

#######################
# RUNNING SCORE PLOTS #
#######################

for ( ctr in deContrasts ) {

  for ( strategy in annotStrategies ) {

    # Assign a data frame to a temporary variable #
    if ( strategy == "goGSEnr" ) {

      enrich.er <- goGSEnr[[ctr]]
      type <- "GO"
      filepath <- "Figures/clusterprofiler/running_score_plots//gsea/go/"

    } else if ( strategy == "kgGSEnr" ) {

      enrich.er <- kgGSEnr[[ctr]]
      type <- "KEGG"
      filepath <- "Figures/clusterprofiler/running_score_plots/gsea/kegg/"

    } else {
      next
    }

    # Move on there is no data to analyze #
    if ( !dim(enrich.er)[1] > 0 || is.list(enrich.er) ) {
      next
    }

    # Dynamically create some plot titles #
    if ( type == "GO" ) {
      tmp.title <- renameContrast(name = ctr, addition = " - Top Enriched GO Terms", flip = FALSE) %>% gsub("[ ]+", " ", .)

    } else if ( type == "KEGG" ) {
      tmp.title <- renameContrast(name = ctr, addition = " - Top Enriched KEGG Pathways", flip = FALSE) %>% gsub("[ ]+", " ", .)
    }

    # Create a color palette #
    colors <- brewer.pal(3, "Set2")

    # Create Running Score Plots #
    gsea.plot <- gseaplot2(enrich.er, geneSetID = 1:3, pvalue_table = FALSE, color = colors, title = tmp.title)

    # Create file names #
    filename <- paste0(tmp.title, " (GSEA) Running Score Plot.pdf")

    # Save the plot to a PDF #
    ggsave(filename, plot = gsea.plot, device = "pdf",
           path = filepath,
           width = 12, height = 10, units = "in", dpi = 400)
  }
}
