###############
### PURPOSE ###
###############

# This script performs several types of functional #
# enrichment analyses via clusterProfiler including #
# gene ontology, KEGG pathway, and gene set enrichment #
# for all contrasts having a sufficient number of DEGs. #

######################
### INITIALIZATION ###
######################

# Load all necessary libraries #
library(tidyverse)
library(DOSE)
library(clusterProfiler)
library(DESeq2)
library(org.Mm.eg.db)

# Set seed for reproducibility #
set.seed(123)

# Set working directory #
setwd("~/Desktop/general/projects/clients/2025/002_steven_toro/")

# Load differential gene expression analysis data #
source("Data/storo_lists/lists_and_vectors.R")
source("Scripts/functions/functions.R")
load("Data/storo_rdata/storo_s03_differential_gene_expression.RData")

# Define some important variables #
min.Genes <- 25
min.Terms <- 2
up.logFC <- 2
dn.logFC <- -2
de.FDR <- 0.05
go.FDR <- 0.05

# Create some empty lists to populate with enrichment results #
enterezAll <- list()
geneVector.lst <- list()
goGSEnr <- list()
kgGSEnr <- list()

enterezDiffExpr <- list()
de.geneVector.lst <- list()
goEnr <- list()
kgEnr <- list()

################################
# GENE SET ENRICHMENT ANALYSIS #
################################

# For each contrast... #
for (ctr in deContrasts) {

  # Print a message to our warnings log to describe which gene set we are analyzing #
  warning(paste0("Now analyzing contrast ", ctr))

  # Convert our DESeq object to a data frame from easier handling #
  tmp <- data.frame(get(ctr))
  tmp[,1] <- rownames(tmp)
  colnames(tmp)[1] <- "ENSEMBL"

  # Convert ENSEMBL IDs to ENTEREZ IDs #
  tmpEntrez <- bitr(tmp$ENSEMBL,
                    fromType = "ENSEMBL",
                    toType = c("ENTREZID", "SYMBOL"),
                    OrgDb = org.Mm.eg.db)

  # If we have not stored the IDs which cannot be converted
  if ( exists("unmappable_ids") == FALSE ) {

    # Do so, for future reference #
    unmappable_ids <- tmp$ENSEMBL[!tmp$ENSEMBL %in% tmpEntrez$ENSEMBL]
  }

  # If we have not yet created an expression matrix w/ gene symbols... #
  if ( exists("symb2expr") == FALSE ) {

    # Create a normalized gene expression matrix with gene symbols #
    geneExpr <- data.frame(ENSEMBL = rownames(assay(nrmDDS)), assay(nrmDDS))

    # We only need to do this once; the contrast doesn't matter because all contrasts report all genes #
    symb2expr <- left_join(tmpEntrez, geneExpr, by = "ENSEMBL")
  }

  # Create a new data frame with Enterez IDs included #
  enterezAll[[ctr]] <- drop_na(left_join(tmpEntrez, tmp, by = "ENSEMBL"))

  # Create named vectors of logFC values
  geneVector.lst[[ctr]] <- enterezAll[[ctr]]$log2FoldChange
  names(geneVector.lst[[ctr]]) <- enterezAll[[ctr]]$ENTREZID

  # A small number of Ensembl IDs map to a yet smaller set of Enterez IDs. Deduplication is necessary. #
  geneVector.lst[[ctr]] <- geneVector.lst[[ctr]][!duplicated(names(geneVector.lst[[ctr]]))]

  # Sort the vector by logFC just in case it's not already sorted. #
  geneVector.lst[[ctr]] <- sort(geneVector.lst[[ctr]], decreasing = T)

  # Perform GO and KEGG Gene Set Enrichment analyses #
  goGSEnr[[ctr]] <- gseGO(geneVector.lst[[ctr]], ont = "ALL", org.Mm.eg.db, eps = 0)
  kgGSEnr[[ctr]] <- gseKEGG(geneVector.lst[[ctr]], organism = "mmu", eps = 0)
}

###############################
# GO/KEGG ENRICHMENT ANALYSIS #
###############################

# For each direction of regulation... #
for (reg in regDirection) {

  # For each contrast... #
  for (ctr in deContrasts) {

    # Print a message to our warnings log to describe which gene set we are analyzing #
    warning(paste0("Now analyzing the ", reg, "REGULATED gene set of contrast ", ctr))

    # Create nested lists for our results #
    enterezDiffExpr[[reg]][[ctr]] <- list()
    de.geneVector.lst[[reg]][[ctr]] <- list()
    goEnr[[reg]][[ctr]] <- list()
    kgEnr[[reg]][[ctr]] <- list()

    # Gather all upregulated and downregulated genes #
    if (reg == "DOWN") {

      enterezDiffExpr[[reg]][[ctr]] <- enterezAll[[ctr]][enterezAll[[ctr]]$padj < 0.05
                                                         & enterezAll[[ctr]]$log2FoldChange < -2,]
    } else if (reg == "UP") {

      enterezDiffExpr[[reg]][[ctr]] <- enterezAll[[ctr]][enterezAll[[ctr]]$padj < 0.05
                                                         & enterezAll[[ctr]]$log2FoldChange > 2,]
    }

    # Create a vector of differentially expressed genes and a background list #
    de.geneVector.lst[[reg]][[ctr]] <- enterezDiffExpr[[reg]][[ctr]]$log2FoldChange
    names(de.geneVector.lst[[reg]][[ctr]]) <- enterezDiffExpr[[reg]][[ctr]]$ENTREZID

    # Implement a manual filter for gene sets that pass the gene set threshold #
    # but nevertheless display a low-quality PWF or model convergence warnings #
    gene.set <- paste0(ctr, "_", reg)

    # Here I would not perform GO enrichment analyses with fewer than 25 DE genes or when pwf plot quality is low #
    if ( length(de.geneVector.lst[[reg]][[ctr]]) >= min.Genes && ! gene.set %in% low.qual.pwfs ) {

      # Perform GO and KEGG Enrichment analyses #
      goEnr[[reg]][[ctr]] <- enrichGO(names(de.geneVector.lst[[reg]][[ctr]]), ont = "ALL", org.Mm.eg.db)
      kgEnr[[reg]][[ctr]] <- enrichKEGG(names(de.geneVector.lst[[reg]][[ctr]]), organism = "mmu")
    } else {

      warning(paste0("Contrast ", ctr, " has too few DE Genes or a low-quality PWF plot"))
      next
    }
  }
}

for ( ctr in deContrasts ) {

  if ( !is.list(goGSEnr[[ctr]]) && !is.list(kgGSEnr[[ctr]]) ) {

  tmp.name <- renameContrast(name = ctr, flip = FALSE) %>% gsub("[ ]+", " ", .)
  goEnr.filename <- paste0("Tables/clusterprofiler/gsea/", tmp.name, " - clusterProfiler GO-GSEA Results.tsv")
  kgEnr.filename <- paste0("Tables/clusterprofiler/gsea/", tmp.name, " - clusterProfiler KEGG-GSEA Results.tsv")

  write.table(x = as.data.frame(goGSEnr[[ctr]]@result), file = goEnr.filename, sep = "\t")
  write.table(x = as.data.frame(kgGSEnr[[ctr]]@result), file = kgEnr.filename, sep = "\t")

  }

  for (reg in regDirection) {

    if ( !is.list(goEnr[[reg]][[ctr]]) && !is.list(kgEnr[[reg]][[ctr]]) ) {

      if ( reg == "UP") {
        tmp.name <- renameContrast(name = ctr, flip = FALSE) %>% gsub("[ ]+", " ", .)
        goEnr.filename <- paste0("Tables/clusterprofiler/enrichment/", tmp.name, " - clusterprofiler GO-Enr Results.tsv")
        kgEnr.filename <- paste0("Tables/clusterprofiler/enrichment/", tmp.name, " - clusterprofiler KEGG-Enr Results.tsv")

      } else if ( reg == "DOWN" ) {
        tmp.name <- renameContrast(name = ctr, flip = TRUE) %>% gsub("[ ]+", " ", .)
        goEnr.filename <- paste0("Tables/clusterprofiler/enrichment/", tmp.name, " - clusterprofiler GO-Enr Results.tsv")
        kgEnr.filename <- paste0("Tables/clusterprofiler/enrichment/", tmp.name, " - clusterprofiler KEGG-Enr Results.tsv")
      }

      write.table(x = as.data.frame(goEnr[[reg]][[ctr]]), file = goEnr.filename, sep = "\t")
      write.table(x = as.data.frame(kgEnr[[reg]][[ctr]]), file = kgEnr.filename, sep = "\t")
    }
  }
}

# Save all RData #
save.image(file = "Data/storo_rdata/storo_s05_clusterProfiler_enrichment.RData")

# Print any warnings and sessionInfo #
writeLines(capture.output(warnings()), "Logs/storo_s05_clusterProfiler_enrichment_warnings.txt")
writeLines(capture.output(sessionInfo()), "Logs/storo_s05_clusterProfiler_enrichment_sessionInfo.txt")

# Clear the environment for the next project #
rm(list = ls())
