###############
### PURPOSE ###
###############

# This script creates several enrichment maps from the #
# results of clusterProfiler functional enrichment analyses. #

######################
### INITIALIZATION ###
######################

# Load libraries #
library(tidyverse)
library(ggnewscale)
library(enrichplot)

# Set working directory #
setwd("~/Desktop/general/projects/clients/2025/002_steven_toro/")

# Import some plotting themes and functions #
source("storo-06-2025/functions/functions.R")
source("storo-06-2025/other/lists_and_vectors.R")
load("Data/storo_rdata/storo_s05_clusterProfiler_enrichment.RData")

# Set seed for reproducibility #
set.seed(123)

########################
# ENRICHMENT MAP PLOTS #
########################

for ( ctr in deContrasts ) {

  for ( strategy in annotStrategies ) {

    for ( reg in regDirection ) {

      # Assign a data frame to a temporary variable #
      if ( strategy == "goEnr" ) {

        enrich.er <- goEnr[[reg]][[ctr]]
        method <- "Enr"
        type <- "GO"
        filepath <- "Figures/clusterprofiler/enrichment_maps/enrichment/go/"

      } else if ( strategy == "kgEnr" ) {

        enrich.er <- kgEnr[[reg]][[ctr]]
        method <- "Enr"
        type <- "KEGG"
        filepath <- "Figures/clusterprofiler/enrichment_maps/enrichment/kegg/"

      } else if ( strategy == "goGSEnr" ) {

        enrich.er <- goGSEnr[[ctr]]
        method <- "GSEnr"
        type <- "GO"
        filepath <- "Figures/clusterprofiler/enrichment_maps/gsea/go/"

      } else if ( strategy == "kgGSEnr" ) {

        enrich.er <- kgGSEnr[[ctr]]
        method <- "GSEnr"
        type <- "KEGG"
        filepath <- "Figures/clusterprofiler/enrichment_maps/gsea/kegg/"
      }

      # Only parse GSEA data once #
      if ( method == "GSEnr" && reg == "DOWN" ) {
        next
      }

      # Also break if there is no data to analyze #
      if ( dim(enrich.er)[1] > 0 && !is.list(enrich.er) ) {

        # Generate -log10(FDR) values and identify a midpoint for our -log10(FDR) scale #
        enrich.er@result$log.padjust <- -log10(enrich.er@result$p.adjust)
        med.logfdr <- median(enrich.er@result$log.padjust)

      } else {
        break
      }

      # Dynamically create some plot titles #
      if ( type == "GO") {

        tmp.addition <- paste0(" - Top Enriched GO Terms")

        if ( reg == "UP" ) {
          tmp.title <- renameContrast(name = ctr, addition = tmp.addition, flip = FALSE) %>% gsub("[ ]+", " ", .)

        } else if ( reg == "DOWN" ) {
          tmp.title <- renameContrast(name = ctr, addition = tmp.addition, flip = TRUE) %>% gsub("[ ]+", " ", .)
        }

      } else if ( type == "KEGG" ) {

        tmp.addition <- paste0(" - Top Enriched KEGG Pathways")

        if ( reg == "UP" ) {
          tmp.title <- renameContrast(name = ctr, addition = tmp.addition, flip = FALSE) %>% gsub("[ ]+", " ", .)

        } else if ( reg == "DOWN" ) {
          tmp.title <- renameContrast(name = ctr, addition = tmp.addition, flip = TRUE) %>% gsub("[ ]+", " ", .)
        }
      }

      emap.plot <- pairwise_termsim(enrich.er) %>%
        emapplot(layout = "star", showCategory = 20) +
        ggtitle(tmp.title) + labs(size = "No. Genes") +
        theme(plot.title = element_text(hjust = 0.5, size = 25),
              legend.title = element_text(size = 14),
              legend.text = element_text(size = 14)) +
        scale_color_continuous(name = "FDR", low = "#FFDBD3", high = "#D1001F")

      # Create file names #
      if ( method == "Enr" ) {

        filename <- paste0(tmp.title, " (Enr) Enrichment Map.pdf")

      } else if ( method == "GSEnr" ) {

        filename <- paste0(tmp.title, " (GSEA) Enrichment Map.pdf")
      }

      # Save the plot to a PDF #
      # If "Viewport has zero dimension(s)" error, plot is empty - just create it and move on #
      try(ggsave(filename, plot = emap.plot, device = "pdf", path = filepath,
             width = 14, height = 10, units = "in", dpi = 400))

    }
  }
}
