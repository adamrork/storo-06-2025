###############
### PURPOSE ###
###############

# This script creates several KEGG pathview and GraphViz plots #
# from the results of clusterProfiler functional enrichment analyses. #

######################
### INITIALIZATION ###
######################

# Load libraries #
library(tidyverse)
library(pathview)

# Set working directory #
setwd("~/Desktop/general/projects/clients/2025/002_steven_toro/")
maindir <- getwd()

# Import some plotting themes and functions #
source("storo-06-2025/functions/functions.R")
source("storo-06-2025/other/lists_and_vectors.R")
load("Data/storo_rdata/storo_s05_clusterProfiler_enrichment.RData")

# Set seed for reproducibility #
set.seed(123)

#############################################
### DEFINE TERMS AND PATHWAYS OF INTEREST ###
#############################################

# Must be spelled exactly as pathview expects #
keggpaths <- list("mmu04310")

#######################
# KEGG PATHVIEW PLOTS #
#######################

for ( ctr in deContrasts ) {

  for ( strategy in annotStrategies ) {

    for ( reg in regDirection ) {

      # Assign a data frame to a temporary variable #
      if ( strategy == "kgEnr" ) {

        enrich.er <- kgEnr[[reg]][[ctr]]
        method <- "Enr"
        gv.filepath <- paste0(maindir, "/Figures/clusterprofiler/graphviz_pathview_plots/enrichment/graphviz/")
        pv.filepath <- paste0(maindir, "/Figures/clusterprofiler/graphviz_pathview_plots/enrichment/pathview/")

      } else if ( strategy == "kgGSEnr" ) {

        enrich.er <- kgGSEnr[[ctr]]
        method <- "GSEnr"
        gv.filepath <- paste0(maindir, "/Figures/clusterprofiler/graphviz_pathview_plots/gsea/graphviz/")
        pv.filepath <- paste0(maindir, "/Figures/clusterprofiler/graphviz_pathview_plots/gsea/pathview/")

      } else {

        next
      }

      # Move on if there is no data to analyze #
      if ( !dim(enrich.er)[1] > 0 || is.list(enrich.er) ) {
        next
      }

      # Dynamically create some file names #
      if ( method == "Enr" ) {

        geneVector <- de.geneVector.lst[[reg]][[ctr]]

        if ( reg == "UP") {
          tmp.title <- renameContrast(name = ctr, addition = " - Enr", flip = FALSE) %>% gsub("[ ]+", " ", .)
        } else if ( reg == "DOWN" ) {
          tmp.title <- renameContrast(name = ctr, addition = " - Enr", flip = TRUE) %>% gsub("[ ]+", " ", .)
        }

      } else if ( method == "GSEnr" ) {

        geneVector <- geneVector.lst[[ctr]]

        if ( reg == "UP") {
          tmp.title <- renameContrast(name = ctr, addition = " - GSEA", flip = FALSE) %>% gsub("[ ]+", " ", .)
        } else if ( reg == "DOWN" ) {
          next
        }
      }

      # For each pathway in the vector or list... #
      for ( pathway in keggpaths ) {

        setwd(gv.filepath)
        # Print the relevant pathways to image files #
        pathview(gene.data  = geneVector,
                 pathway.id = pathway,
                 species = "mmu", kegg.native = TRUE,
                 pdf.size = c(7, 7),
                 out.suffix = tmp.title)

        setwd(pv.filepath)
        # Print the relevant pathways to image files #
        pathview(gene.data  = geneVector,
                 pathway.id = pathway,
                 species = "mmu", kegg.native = FALSE,
                 pdf.size = c(7, 7),
                 out.suffix = tmp.title)
      }
    }
  }
}
