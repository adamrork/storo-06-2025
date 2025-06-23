###############
### PURPOSE ###
###############

# This script creates several ridgeline plots from the
# results of clusterProfiler functional enrichment analyses. #

######################
### INITIALIZATION ###
######################

# Load libraries #
library(tidyverse)
library(enrichplot)

# Set working directory #
setwd("~/Desktop/general/projects/clients/2025/002_steven_toro/")

# Import some plotting themes and functions #
source("Scripts/functions/functions.R")
source("Data/storo_lists/lists_and_vectors.R")
load("Data/storo_rdata/storo_s05_clusterProfiler_enrichment.RData")

# Set seed for reproducibility #
set.seed(123)

###################
# RIDGELINE PLOTS #
###################

for ( ctr in deContrasts ) {

  for ( strategy in annotStrategies ) {

    # Assign a data frame to a temporary variable #
    if ( strategy == "goGSEnr" ) {

      enrich.er <- goGSEnr[[ctr]]
      type <- "GO"
      filepath <- "Figures/clusterprofiler/logfc_ridgeline_plots/gsea/go/"

    } else if ( strategy == "kgGSEnr" ) {

      enrich.er <- kgGSEnr[[ctr]]
      type <- "KEGG"
      filepath <- "Figures/clusterprofiler/logfc_ridgeline_plots/gsea/kegg/"

    } else {
      next
    }

    # Move on there is no data to analyze #
    if ( dim(enrich.er)[1] > 0 && !is.list(enrich.er) ) {

      # Generate -log10(FDR) values #
      enrich.er@result$log.padjust <- -log10(enrich.er@result$p.adjust)

    } else {

      next
    }

    # Dynamically create some plot titles #
    if ( type == "GO" ) {
      tmp.title <- renameContrast(name = ctr, addition = " - Top Enriched GO Terms", flip = FALSE) %>% gsub("[ ]+", " ", .)

    } else if ( type == "KEGG" ) {
      tmp.title <- renameContrast(name = ctr, addition = " - Top Enriched KEGG Pathways", flip = FALSE) %>% gsub("[ ]+", " ", .)
    }

    # Create Ridgeline Plots #
    rdg.plot <- ridgeplot(enrich.er, showCategory = 10, fill = "log.padjust") +
      ggtitle(tmp.title) +
      xlab("log2(Fold Change)") + ylab("Enriched Term") +
      scale_fill_continuous(name =  "-log10(FDR)", low = "#FFDBD3", high = "#D1001F") +
      theme(plot.title = element_text(hjust = 0.5, size = 22),
            axis.title = element_text(size = 19),
            axis.text = element_text(size = 12),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 14))

    # Create file names #
    filename <- paste0(tmp.title, " (GSEA) Ridgeline Plot.pdf")

    # Save the plot to a PDF #
    ggsave(filename, plot = rdg.plot, device = "pdf", path = filepath,
           width = 14, height = 10, units = "in", dpi = 400)
  }
}
