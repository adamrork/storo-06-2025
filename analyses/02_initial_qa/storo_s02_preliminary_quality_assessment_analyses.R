###############
### PURPOSE ###
###############

# This script calculates basic summary statistics on #
# all counts data before performing differential gene #
# expression analyses to ensure there are no obvious #
# quality issues with any of the libraries. #

######################
### INITIALIZATION ###
######################

# Load libraries #
library(tidyverse)
library(edgeR)
library(vsn)

# Set seed for reproducibility #
set.seed(123)

# Set working directory #
setwd("~/Desktop/general/projects/clients/2025/002_steven_toro/Data/storo_counts/counts/")

#################
# PREPROCESSING #
#################

# Populate an edgeR list with counts data #
datafiles <- list.files(".", full.names = FALSE)
countsList <- readDGE(datafiles, columns = c(1, 7), skip = 1)

# Clean up column names #
colnames(countsList) <- gsub(".counts", "", colnames(countsList))

# Add treatment and timepoint (day) information to the samples data frame #
countsList$samples$treatment <- as.factor(rep(c("Day0", "NT", "PRG4", "PRG4-VEGF"), times = c(3, 6, 6, 6)))
countsList$samples$day <- as.factor(c(rep("D00", times = 3), rep(c("D10", "D24"), times = 3, each = 3)))

######################
# SUMMARY STATISTICS #
######################

# Calculate raw and log-transformed CPM values for all counts data #
rawCPM <- cpm(countsList, log = FALSE)
logCPM <- cpm(countsList, log = TRUE)

# Gather some basic library size statistics #
avg.libSize <- mean(countsList$samples$lib.size) # 19.9 million reads
med.libSize <- median(countsList$samples$lib.size) # 19.6 million reads
rng.libSize <- range(countsList$samples$lib.size) # 16.4 million reads to 27.2 million reads

# Filter out genes which have low expression across all samples #
allGenes <- rownames(countsList$counts)
exprGenes <- rowSums(countsList$counts >= 10) >= 3
filteredCounts <- countsList[exprGenes, keep.lib.sizes = FALSE]

# Calculate raw and log-transformed CPM values for filtered counts data #
filtered.rawCPM <- cpm(filteredCounts, log = FALSE)
filtered.logCPM <- cpm(filteredCounts, log = TRUE)

# Cross-normalize samples #
filteredCounts <- calcNormFactors(filteredCounts, method = "TMM")
normFactors <- data.frame(Names = rownames(filteredCounts$samples),
                          Factors = filteredCounts$samples$norm.factors,
                          Groups = paste0(filteredCounts$samples$treatment, "_", filteredCounts$samples$day))

# Inspect the spread of normalization factors (all should be around 1) #
normRanges <- range(normFactors$Factors)

# Convert all CPM tables to formats suitable for plotting #
rawCPM.df <- data.frame(genes = rownames(rawCPM), rawCPM)
logCPM.df <- data.frame(genes = rownames(logCPM), logCPM)
filtered.rawCPM.df <- data.frame(genes = rownames(filtered.rawCPM), filtered.rawCPM)
filtered.logCPM.df <- data.frame(genes = rownames(filtered.logCPM), filtered.logCPM)

rawCPM.df <- pivot_longer(rawCPM.df, cols = 2:22, names_to = "Sample")
logCPM.df <- pivot_longer(logCPM.df, cols = 2:22, names_to = "Sample")
filtered.rawCPM.df <- pivot_longer(filtered.rawCPM.df, cols = 2:22, names_to = "Sample")
filtered.logCPM.df <- pivot_longer(filtered.logCPM.df, cols = 2:22, names_to = "Sample")

rawCPM.df$Group <- gsub("_[1-3]$", "", rawCPM.df$Sample)
logCPM.df$Group <- gsub("_[1-3]$", "", logCPM.df$Sample)
filtered.rawCPM.df$Group <- gsub("_[1-3]$", "", filtered.rawCPM.df$Sample)
filtered.logCPM.df$Group <- gsub("_[1-3]$", "", filtered.logCPM.df$Sample)

#####################
# CONCLUDE ANALYSIS #
#####################

# Save all RData #
setwd("~/Desktop/general/projects/clients/2025/002_steven_toro/")
save.image(file = "Data/storo_rdata/storo_s02_preliminary_quality_assessment.RData")

# Print any warnings and sessionInfo #
writeLines(capture.output(warnings()), "Logs/storo_s02_preliminary_quality_assessment_warnings.txt")
writeLines(capture.output(sessionInfo()), "Logs/storo_s02_preliminary_quality_assessment_sessionInfo.txt")

# Clear the environment for the next project #
rm(list = ls())
