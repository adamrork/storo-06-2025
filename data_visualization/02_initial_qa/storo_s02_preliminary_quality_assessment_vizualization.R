###############
### PURPOSE ###
###############

# This script generates a small set of figures used to ensure #
# the counts data are of reasonable quality before proceeding #
# with differential gene expression analyses. #

######################
### INITIALIZATION ###
######################

# Load libraries #
library(tidyverse)
library(RColorBrewer)
library(edgeR)
library(cowplot)

# Set working directory #
setwd("~/Desktop/general/projects/clients/2025/002_steven_toro/")

# Load differential gene expression analysis data #
load("Data/storo_rdata/storo_s02_preliminary_quality_assessment.RData")

# Set seed for reproducibility #
set.seed(123)

#########################
### FIGURE GENERATION ###
#########################

# Define color palettes #
dayColors <- rep(brewer.pal(3, "Dark2"), each = 3, times = 3)[-c(10, 11, 12, 19, 20, 21)]
treatmentColors <- rep(brewer.pal(4, "Dark2"), times = c(3, 6, 6, 6))
groupColorsAll <- rep(brewer.pal(7, "Dark2"), each = 3)
groupColorsUniq <- brewer.pal(7, "Dark2")

# Create treatment-timepoint group names #
groupNames <- unique(paste0(countsList$samples$treatment, "_", countsList$samples$day))

#################################
# GENE EXPRESSION DISTRIBUTIONS #
#################################

# Create box plots of filtered and unfiltered raw and normalized CPM values #
rawCPM.plot <- ggplot(rawCPM.df, aes(x = Sample, y = value, fill = Group)) +
  geom_boxplot() +
  ggtitle("Unfiltered CPM") + ylab("CPM") +
  ylim(0, 75) +
  theme(plot.title= element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 340))

logCPM.plot <- ggplot(logCPM.df, aes(x = Sample, y = value, fill = Group)) +
  geom_boxplot() +
  ggtitle("Unfiltered log2(CPM)") + ylab("log2(CPM)") +
  ylim(0, 20) +
  theme(plot.title= element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 340))

flt.rawCPM.plot <- ggplot(filtered.rawCPM.df, aes(x = Sample, y = value, fill = Group)) +
  geom_boxplot() +
  ggtitle("Filtered CPM") + ylab("CPM") +
  ylim(0, 75) +
  theme(plot.title= element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 340))

flt.logCPM.plot <- ggplot(filtered.logCPM.df, aes(x = Sample, y = value, fill = Group)) +
  geom_boxplot() +
  ggtitle("Filtered log2(CPM)") + ylab("log2(CPM)") +
  ylim(0, 20) +
  theme(plot.title= element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 340))

unfiltered.CPM.plot <- plot_grid(rawCPM.plot, logCPM.plot, nrow = 2)
filtered.CPM.plot <- plot_grid(flt.rawCPM.plot, flt.logCPM.plot, nrow = 2)

# Save the plots to PDFs #
ggsave("Unfiltered CPM Boxplot.pdf", plot = unfiltered.CPM.plot, device = "pdf",
       path = "Figures/initial_quality_assessment/",
       width = 14, height = 8, units = "in", dpi = 400)

ggsave("Filtered CPM Boxplot.pdf", plot = filtered.CPM.plot, device = "pdf",
       path = "Figures/initial_quality_assessment/",
       width = 14, height = 8, units = "in", dpi = 400)

# Create density plots of unfiltered log2-normalized CPM values from each sample and save plot to a PDF #
pdf("Figures/initial_quality_assessment/log2CPM Density Plots.pdf",
    width = 8, height = 8)

  plot(density(logCPM[,1]), col = groupColorsAll, main = "Unfiltered log2(CPM)", xlab = "log2(CPM)")

  for (i in 2:21) {
    den <- density(logCPM[,i])
    lines(den$x, den$y, col = groupColorsAll[i])
  }

  legend("topright", groupNames, text.col = groupColorsUniq)

  plot(density(filtered.logCPM[,1]), col = groupColorsAll, main = "Filtered log2(CPM)", xlab = "log2(CPM)")

  for (i in 2:21) {
    den <- density(filtered.logCPM[,i])
    lines(den$x, den$y, col = groupColorsAll[i])
  }

  legend("topright", groupNames, text.col = groupColorsUniq)

dev.off()

#########################
# NORMALIZATION FACTORS #
#########################
# Create a dot plot of library size normalization factors #
colnames(normFactors)[1] <- "Samples"
nfPlot <- ggplot(normFactors, aes(x = Factors, y = Samples, fill = Groups)) +
            geom_dotplot() +
            coord_flip() +
            ggtitle("Normalization Factors") +
            theme(plot.title = element_text(hjust = 0.5, size = 15),
                  axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot to a PDF #
ggsave("Normalization Factors (edgeR).pdf", plot = nfPlot, device = "pdf",
       path = "Figures/initial_quality_assessment/",
       width = 8, height = 8, units = "in", dpi = 400)

#############
# MDS PLOTS #
#############

# Create an MDS plot annotated by treatment and timepoint (based on expression of all genes and a filtered subset) #
par(mfrow = c(1, 2))

pdf("Figures/initial_quality_assessment/MDS Plot by Treatment.pdf",
    width = 8, height = 8)

  plotMDS(logCPM, labels = countsList$samples$treatment, col = treatmentColors, dim = c(1, 2))
    title("Treatments (Unfiltered) MDS")
  plotMDS(filtered.logCPM, labels = countsList$samples$treatment, col = treatmentColors, dim = c(1, 2))
    title("Treatments (Filtered) MDS")

dev.off()

pdf("Figures/initial_quality_assessment/MDS Plot by Timepoint.pdf",
    width = 8, height = 8)

  plotMDS(logCPM, labels = countsList$samples$day, col = dayColors, dim = c(1, 2))
    title("Timepoints (Unfiltered) MDS")
  plotMDS(filtered.logCPM, labels = countsList$samples$day, col = dayColors, dim = c(1, 2))
    title("Timepoints (Filtered) MDS")

dev.off()

# Clear the environment for the next project #
rm(list = ls())
