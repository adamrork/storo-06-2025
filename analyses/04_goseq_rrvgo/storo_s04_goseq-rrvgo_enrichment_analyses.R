###############
### PURPOSE ###
###############

# This script performs gene ontology enrichment #
# analyses via goseq and semantic similarity #
# clustering of enriched terms via rrvgo. All #
# GO categories are examined and are clustered #
# at three different levels for all contrasts #
# having a sufficient number of DEGs. #

######################
### INITIALIZATION ###
######################

# Load all necessary libraries #
library(tidyverse)
library(goseq)
library(rrvgo)
library(GenomicFeatures)

# Set seed for reproducibility #
set.seed(123)

# Set working directory and other important directories #
setwd("~/Desktop/general/projects/clients/2025/002_steven_toro/")

table.dir <- "~/Desktop/general/projects/clients/2025/002_steven_toro/Tables/"

# Load differential gene expression analysis data and list of contrasts #
source("Data/storo_lists/lists_and_vectors.R")
source("Scripts/functions/functions.R")
load("Data/storo_rdata/storo_s03_differential_gene_expression.RData")

# Define some important variables #
min.Genes <- 50
min.Terms <- 2
up.logFC <- 2
dn.logFC <- -2
de.FDR <- 0.05
go.FDR <- 0.05

# Create some empty lists to populate with enrichment results #
de.data <- list()
de.genes <- list()
pwf.lst <- list()
goData.lst <- list()
reducedTerms.lst <- list()
rrvgoData.lst <- list()

#################
# GO Enrichment #
#################

# For each direction of differential gene expression... #
for (reg in regDirection) {

  de.data[[reg]] <- list()
  de.genes[[reg]] <- list()
  pwf.lst[[reg]] <- list()

  filteredContrasts <- vector()

  # For each contrast for which differentially expressed genes were examined... #
  for ( ctr in seq_along(deContrasts) ) {

    # Convert data to a standard data frame #
    de.data[[reg]][[ctr]] <- data.frame(get(deContrasts[[ctr]]))

    # If we are examining genes downregulated in Y in an X vs. Y contrast #
    if (reg == "DOWN") {

      # Create a positional vector of downregulated (1) and non-downregulated (0) genes #
      de.genes[[reg]][[ctr]] <- as.integer(!is.na(de.data[[reg]][[ctr]]$padj)
                                           & de.data[[reg]][[ctr]]$padj < de.FDR
                                           & de.data[[reg]][[ctr]]$log2FoldChange < dn.logFC)

    # If we are examining genes upregulated in X in an X vs. Y contrast #
    } else if (reg == "UP") {

      # Create a positional vector of upregulated (1) and non-upregulated (0) genes #
      de.genes[[reg]][[ctr]] <- as.integer(!is.na(de.data[[reg]][[ctr]]$padj)
                                           & de.data[[reg]][[ctr]]$padj < de.FDR
                                           & de.data[[reg]][[ctr]]$log2FoldChange > up.logFC)
    }

    # Add gene names to vector of differential expression statuses #
    names(de.genes[[reg]][[ctr]]) <- row.names(de.data[[reg]][[ctr]])

    # In general, I would probably not try to perform GO enrichment analyses with fewer than 50 DE genes #
    if ( sum(de.genes[[reg]][[ctr]]) >= min.Genes ) {

      # Calculate a probability weighting function using mm10 gene lengths #
      # Don't fail over "plot.new() : figure margins to large" errors #
      try(pwf.lst[[reg]][[ctr]] <- nullp(de.genes[[reg]][[ctr]],
                                     "mm10", "ensGene",
                                     bias.data = NULL, plot.fit = TRUE))

    } else {

      cat("Contrast", deContrasts[[ctr]], "had too few DE Genes\n")
      filteredContrasts <- append(filteredContrasts, deContrasts[[ctr]])
      next
    }

    # Perform a GO enrichment analysis with GOseq #
    de.goResults <- goseq(pwf.lst[[reg]][[ctr]],
                         "mm10",
                         "ensGene",
                         test.cats=c("GO:BP", "GO:MF", "GO:CC"))

    # Correct over-represented p-values via Benjamini-Hochberg #
    de.goResults$FDR <- p.adjust(de.goResults$over_represented_pvalue, "BH")
    de.goResults <- de.goResults[which(de.goResults$FDR < go.FDR),]

    # Define additional lists to store enrichment and clustering data for each GO category #
    goData.lst[[reg]][[ctr]] <- list()
    reducedTerms.lst[[reg]][[ctr]] <- list()
    rrvgoData.lst[[reg]][[ctr]] <- list()

    # For each GO category... #
    for ( ctg in goCategories ) {

      # Perform semantic similarity clustering for enriched GO terms #
      goData.lst[[reg]][[ctr]][[ctg]] <- de.goResults[which(de.goResults$ontology == ctg),]

      if ( nrow(goData.lst[[reg]][[ctr]][[ctg]]) >= min.Terms ) {

        de.simData <- calculateSimMatrix(goData.lst[[reg]][[ctr]][[ctg]]$category,
                                         orgdb = "org.Mm.eg.db",
                                         ont = ctg,
                                         method = "Rel")

       # Use negative log-transformed FDRs as scores #
       de.scores <- setNames(-log10(goData.lst[[reg]][[ctr]][[ctg]]$FDR), goData.lst[[reg]][[ctr]][[ctg]]$category)

       # Define additional lists to store clustering data at various thresholds #
       reducedTerms.lst[[reg]][[ctr]][[ctg]] <- list()
       rrvgoData.lst[[reg]][[ctr]][[ctg]] <- list()

       # Use a few clustering thresholds when reducing the similarity matrix to enable exploratory analyses... #
       for ( lvl in clusterLevels ) {

         # Group terms based on our clustering thresholds #
         reducedTerms.lst[[reg]][[ctr]][[ctg]][[lvl]] <- reduceSimMatrix(de.simData,
                                                                  de.scores,
                                                                  threshold = (as.numeric(lvl) / 100),
                                                                  orgdb = "org.Mm.eg.db")

         # Make a results data frame with terms, parent terms, gene counts, and enrichment FDRs #
         rrvgoData.lst[[reg]][[ctr]][[ctg]][[lvl]] <- left_join(goData.lst[[reg]][[ctr]][[ctg]],
                                                                reducedTerms.lst[[reg]][[ctr]][[ctg]][[lvl]],
                                                                by = "term")[c("category", "term",
                                                                               "parent", "parentTerm",
                                                                               "numDEInCat", "numInCat", "FDR")]

         # Sort said data frame by FDR (lowest to highest) #
         rrvgoData.lst[[reg]][[ctr]][[ctg]][[lvl]] <- rrvgoData.lst[[reg]][[ctr]][[ctg]][[lvl]][order(rrvgoData.lst[[reg]][[ctr]][[ctg]][[lvl]]$FDR,
                                                                                                      decreasing = FALSE),]

         # Convert terms to factors in case we need them in that format for cleaner plotting #
         rrvgoData.lst[[reg]][[ctr]][[ctg]][[lvl]]$term <- factor(rrvgoData.lst[[reg]][[ctr]][[ctg]][[lvl]]$term,
                                                    levels = rrvgoData.lst[[reg]][[ctr]][[ctg]][[lvl]]$term)
        }

      # If there are too few terms to perform clustering, move to the next category #
      } else {

        de.simData < "NA"
      }
    }
  }

  # Filter out null list elements and name all level-2 list elements (i.e. contrasts) #
  rrvgoData.lst[[reg]] <- compact(rrvgoData.lst[[reg]])
  names(rrvgoData.lst[[reg]]) <- deContrasts[!(deContrasts %in% filteredContrasts)]
}

# Write results to tables #
for (reg in regDirection) {

  for ( ctr in seq_along(deContrasts) ) {

    i <- ctr
    ctr <- deContrasts[[ctr]]

    go.tmp.table <- data.frame()
    rv.tmp.table <- data.frame()

    for ( ctg in goCategories ) {
        go.tmp.table <- rbind(go.tmp.table, goData.lst[[reg]][[i]][[ctg]])
        rv.tmp.table <- rbind(rv.tmp.table, rrvgoData.lst[[reg]][[ctr]][[ctg]][["90"]])
    }

      if ( reg == "UP") {
        tmp.name <- renameContrast(name = ctr, flip = FALSE) %>% gsub("[ ]+", " ", .)
        go.filename <- paste0("Tables/goseq_rrvgo_enrichment/", tmp.name, " - GOSEQ Results.tsv")
        rv.filename <- paste0("Tables/goseq_rrvgo_enrichment/", tmp.name, " - RRVGO Results.tsv")

      } else if ( reg == "DOWN" ) {
        tmp.name <- renameContrast(name = ctr, flip = TRUE) %>% gsub("[ ]+", " ", .)
        go.filename <- paste0("Tables/goseq_rrvgo_enrichment/", tmp.name, " - GOSEQ Results.tsv")
        rv.filename <- paste0("Tables/goseq_rrvgo_enrichment/", tmp.name, " - RRVGO Results.tsv")
      }

      write.table(x = go.tmp.table, file = go.filename, sep = "\t", row.names = FALSE)
      write.table(x = rv.tmp.table, file = rv.filename, sep = "\t", row.names = FALSE)
  }
}

# Save all RData #
save.image(file = "Data/storo_rdata/storo_s04_goseq-rrvgo_enrichment.RData")

# Print any warnings and sessionInfo #
writeLines(capture.output(warnings()), "Logs/storo_s04_goseq-rrvgo_enrichment_warnings.txt")
writeLines(capture.output(sessionInfo()), "Logs/storo_s04_goseq-rrvgo_enrichment_sessionInfo.txt")

# Clear the environment for the next project #
rm(list = ls())
