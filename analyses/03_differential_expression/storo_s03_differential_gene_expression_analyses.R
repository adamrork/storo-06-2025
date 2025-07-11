###############
### PURPOSE ###
###############

# This script performs differential gene expression #
# analyses (DGEA) via DESeq2. It also generates some #
# summary statistics for post-DGEA quality assessment. #

######################
### INITIALIZATION ###
######################

# Load libraries #
library(tidyverse)
library(DESeq2)

# Set seed for reproducibility #
set.seed(123)

# Set working directory and other important directories #
setwd("~/Desktop/general/projects/clients/2025/002_steven_toro/Data/storo_counts/counts/")

# Load contrasts #
source("../../../storo-06-2025/functions/functions.R")
source("../../../storo-06-2025/other/lists_and_vectors.R")

#################
# PREPROCESSING #
#################

# Populate a data frame with counts data #
datafiles <- list.files(".", full.names = FALSE)

for (i in seq_along(datafiles)) {
  dfile <- datafiles[[i]]
  if (i == 1) {
    counts.df <- read.delim(dfile, skip = 1)
    lengths.df <- counts.df[,c(1, 6)]
    counts.df <- counts.df[,c(1, 7)]
  } else {
    tmp <- read.delim(dfile, skip = 1)
    tmp <- tmp[,c(1, 7)]
    counts.df <- left_join(counts.df, tmp, "Geneid")
  }
  colnames(counts.df)[i + 1] <- gsub(".counts.tsv", "", dfile)
}

# Clean up the data frames #
colnames(counts.df)[1] <- "GeneID"
counts.df[,1] <- gsub('^gene:', "", counts.df[,1])
counts.df <- data.frame(counts.df[,-1], row.names = counts.df[,1])

colnames(lengths.df)[1] <- "GeneID"
lengths.df[,1] <- gsub('^gene:', "", lengths.df[,1])

# Create a vector of gene lengths for downstream analyses #
geneLengths <- lengths.df$Length
names(geneLengths) <- lengths.df$GeneID

# Generate a data frame with all relevant sample substrings #
l.sstr <- colnames(counts.df) # Library string #
g.sstr <- gsub("_[1-3]$", "", l.sstr) # Group substring #
t.sstr <- gsub("_.*", "", l.sstr) # Treatment substring #
d.sstr <- gsub("_[1-3]$", "", l.sstr) %>% gsub(".*_", "", .) # Day (timepoint) substring #
r.sstr <- (gsub(".*_", "",  l.sstr)) # Replicate substring #
samples.df <- data.frame(cbind(l.sstr, g.sstr, t.sstr, d.sstr, r.sstr)) # Group data frame #
colnames(samples.df) <- c("Library", "Group", "Treatment", "Day", "Replicate")

################################
# DIFFERENTIAL GENE EXPRESSION #
################################

# Reset working directory #
setwd("~/Desktop/general/projects/clients/2025/002_steven_toro/")

# Create a DESeq dataset w/ a cell-means design to enable flexible pairwise contrasts #
unfiltered_dds <- DESeqDataSetFromMatrix(countData = counts.df,
                                         colData = samples.df,
                                         design = ~0 + Group)

# Filter out genes with extremely low expression (helpful for QC - otherwise, not strictly necessary) #
expressed <- rowSums(counts(unfiltered_dds) >= 10) >= 3 # Set to ">= 0 & >= 0" for no filtering.
filtered_dds <- unfiltered_dds[expressed,]

# Run DESeq #
dds <- DESeq(filtered_dds)

# Create a model matrix and coefficient vectors for each sample #
model_matrix <- model.matrix(design(dds), colData(dds))
D0_D00 <- colMeans(model_matrix[dds$Group == "DAY0_D00", ])
NT_D10 <- colMeans(model_matrix[dds$Group == "NT_D10", ])
NT_D24 <- colMeans(model_matrix[dds$Group == "NT_D24", ])
PRG4_D10 <- colMeans(model_matrix[dds$Group == "PRG4_D10", ])
PRG4_D24 <- colMeans(model_matrix[dds$Group == "PRG4_D24", ])
PRG4.VEGF_D10 <- colMeans(model_matrix[dds$Group == "PRG4.VEGF_D10", ])
PRG4.VEGF_D24 <- colMeans(model_matrix[dds$Group == "PRG4.VEGF_D24", ])

# Extract results for all contrasts of interest #

# NT_D10 vs. DAY0_D00 #
NT_D10_vs_D0_D00 <- results(dds, contrast = NT_D10 - D0_D00, alpha = 0.05)

NT_D10_vs_D0_D00 <- NT_D10_vs_D0_D00[order(NT_D10_vs_D0_D00$log2FoldChange,
                                           decreasing = TRUE),]
summary(NT_D10_vs_D0_D00)

# NT_D24 vs. DAY0_D00 #
NT_D24_vs_D0_D00 <- results(dds, contrast = NT_D24 - D0_D00, alpha = 0.05)

NT_D24_vs_D0_D00 <- NT_D24_vs_D0_D00[order(NT_D24_vs_D0_D00$log2FoldChange,
                                           decreasing = TRUE),]
summary(NT_D24_vs_D0_D00)

# PRG4_D10 vs. NT_D10 #
PRG4_D10_vs_NT_D10 <- results(dds, contrast = PRG4_D10 - NT_D10, alpha = 0.05)

PRG4_D10_vs_NT_D10 <- PRG4_D10_vs_NT_D10[order(PRG4_D10_vs_NT_D10$log2FoldChange,
                                               decreasing = TRUE),]
summary(PRG4_D10_vs_NT_D10)

# PRG4_D24 vs. NT_D24 #
PRG4_D24_vs_NT_D24 <- results(dds, contrast = PRG4_D24 - NT_D24, alpha = 0.05)

PRG4_D24_vs_NT_D24 <- PRG4_D24_vs_NT_D24[order(PRG4_D24_vs_NT_D24$log2FoldChange,
                                               decreasing = TRUE),]
summary(PRG4_D24_vs_NT_D24)

# PRG4.VEGF_D10 vs. NT_D10 #
PRG4.VEGF_D10_vs_NT_D10 <- results(dds, contrast = PRG4.VEGF_D10 - NT_D10, alpha = 0.05)

PRG4.VEGF_D10_vs_NT_D10 <- PRG4.VEGF_D10_vs_NT_D10[order(PRG4.VEGF_D10_vs_NT_D10$log2FoldChange,
                                                         decreasing = TRUE),]
summary(PRG4.VEGF_D10_vs_NT_D10)

# PRG4.VEGF_D24 vs. NT_D24 #
PRG4.VEGF_D24_vs_NT_D24 <- results(dds, contrast = PRG4.VEGF_D24 - NT_D24, alpha = 0.05)

PRG4.VEGF_D24_vs_NT_D24 <- PRG4.VEGF_D24_vs_NT_D24[order(PRG4.VEGF_D24_vs_NT_D24$log2FoldChange,
                                                         decreasing = TRUE),]
summary(PRG4.VEGF_D24_vs_NT_D24)

# PRG4.VEGF_D10 vs. PRG4_D10 #
PRG4.VEGF_D10_vs_PRG4_D10 <- results(dds, contrast = PRG4.VEGF_D10 - PRG4_D10, alpha = 0.05)

PRG4.VEGF_D10_vs_PRG4_D10 <- PRG4.VEGF_D10_vs_PRG4_D10[order(PRG4.VEGF_D10_vs_PRG4_D10$log2FoldChange,
                                                             decreasing = TRUE),]
summary(PRG4.VEGF_D10_vs_PRG4_D10)

# PRG4.VEGF_D24 vs. PRG4_D24 #
PRG4.VEGF_D24_vs_PRG4_D24 <- results(dds, contrast = PRG4.VEGF_D24 - PRG4_D24, alpha = 0.05)

PRG4.VEGF_D24_vs_PRG4_D24 <- PRG4.VEGF_D24_vs_PRG4_D24[order(PRG4.VEGF_D24_vs_PRG4_D24$log2FoldChange,
                                                             decreasing = TRUE),]
summary(PRG4.VEGF_D24_vs_PRG4_D24)

# NT_D24 vs. NT_D10 #
NT_D24_vs_NT_D10 <- results(dds, contrast = NT_D24 - NT_D10, alpha = 0.05)

NT_D24_vs_NT_D10 <- NT_D24_vs_NT_D10[order(NT_D24_vs_NT_D10$log2FoldChange,
                                           decreasing = TRUE),]
summary(NT_D24_vs_NT_D10)

# PRG4_D24 vs. PRG4_D10 #
PRG4_D24_vs_PRG4_D10 <- results(dds, contrast = PRG4_D24 - PRG4_D10, alpha = 0.05)

PRG4_D24_vs_PRG4_D10 <- PRG4_D24_vs_PRG4_D10[order(PRG4_D24_vs_PRG4_D10$log2FoldChange,
                                                   decreasing = TRUE),]
summary(PRG4_D24_vs_PRG4_D10)

# PRG4.VEGF_D24 vs. PRG4.VEGF_D10 #
PRG4.VEGF_D24_vs_PRG4.VEGF_D10 <- results(dds, contrast = PRG4.VEGF_D24 - PRG4.VEGF_D10, alpha = 0.05)

PRG4.VEGF_D24_vs_PRG4.VEGF_D10 <- PRG4.VEGF_D24_vs_PRG4.VEGF_D10[order(PRG4.VEGF_D24_vs_PRG4.VEGF_D10$log2FoldChange,
                                                                       decreasing = TRUE),]
summary(PRG4.VEGF_D24_vs_PRG4.VEGF_D10)

######################
# SUMMARY STATISTICS #
######################

# Create a superset of all differentially expressed genes and a list of DEG sets #
degVector <- vector()
geneList <- list()
de.geneList<- list()

for (ctr in deContrasts) {

  tmp <- get(ctr)

  degs <- rownames(tmp[!is.na(tmp$padj)
                   & abs(tmp$log2FoldChange) > 2.0
                   & tmp$padj < 0.05,])


  geneList[[ctr]] <- tmp[!is.na(tmp$padj),]

  de.geneList[[ctr]] <- tmp[!is.na(tmp$padj)
                        & abs(tmp$log2FoldChange) > 2.0
                        & tmp$padj < 0.05,]

  degVector <- append(degs, degVector)

  tmp.name <- renameContrast(name = ctr) %>% gsub("[ ]+", " ", .)

  # Write results to tables #
  filename <- paste0("Tables/differential_expression/", tmp.name, " - Differenial Expression Results.tsv")
  write.table(x = data.frame(get(ctr)), file = filename, sep = "\t")
}

# Generate some summary statistics for quality assessment #
groupDF <- as.data.frame(colData(dds)[,c("Treatment", "Day")])

# Deduplicate the degVector #
degVector <- unique(degVector)

# Generate transformed counts #
nrmDDS <- normTransform(dds)
vsdDDS <- vst(dds, blind = FALSE)
rlgDDS <- rlog(dds, blind = FALSE)

de.nrmDDS <- nrmDDS[which(rownames(nrmDDS) %in% degVector)]
de.vsdDDS <- vsdDDS[which(rownames(vsdDDS) %in% degVector)]
de.rlgDDS <- rlgDDS[which(rownames(rlgDDS) %in% degVector)]

# Generate distance matrices #
nrmDDS.dist <- dist(t(assay(nrmDDS)))
vsdDDS.dist <- dist(t(assay(vsdDDS)))
rlgDDS.dist <- dist(t(assay(rlgDDS)))

de.nrmDDS.dist <- dist(t(assay(de.nrmDDS)))
de.vsdDDS.dist <- dist(t(assay(de.vsdDDS)))
de.rlgDDS.dist <- dist(t(assay(de.rlgDDS)))

# Generate PCA data #
pcaData <- plotPCA(vsdDDS, intgroup = c("Treatment", "Day"), returnData = TRUE)
de.pcaData <- plotPCA(de.vsdDDS, intgroup = c("Treatment", "Day"), returnData = TRUE)

pcaData$Treatment <- gsub("DAY0", "Day 0", pcaData$Treatment) %>% gsub("PRG4.VEGF", "PRG4-VEGF", .)
pcaData$Day <- gsub("D", "", pcaData$Day) %>% gsub("00", "  0", .)

de.pcaData$Treatment <- gsub("DAY0", "Day 0", de.pcaData$Treatment) %>% gsub("PRG4.VEGF", "PRG4-VEGF", .)
de.pcaData$Day <- gsub("D", "", de.pcaData$Day) %>% gsub("00", "  0", .)

pcaVar <- round(100 * attr(pcaData, "percentVar"))
de.pcaVar <- round(100 * attr(de.pcaData, "percentVar"))

#####################
# CONCLUDE ANALYSIS #
#####################

# Save all RData #
save.image(file = "Data/storo_rdata/storo_s03_differential_gene_expression.RData")

# Print any warnings and sessionInfo #
writeLines(capture.output(warnings()), "Logs/storo_s03_differential_gene_expression_warnings.txt")
writeLines(capture.output(sessionInfo()), "Logs/storo_s03_differential_gene_expression_sessionInfo.txt")

# Clear the environment for the next project #
rm(list = ls())
