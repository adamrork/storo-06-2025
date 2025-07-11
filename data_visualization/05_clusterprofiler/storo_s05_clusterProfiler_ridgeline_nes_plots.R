###############
### PURPOSE ###
###############

# This script creates several NES ridgeline and dot plots from #
# the results of clusterProfiler functional enrichment analyses. #

# Design adapted from Figure 5 of "Saarensilta et al. 2025"
# https://doi.org/10.1038/s41598-025-91511-0 #

######################
### INITIALIZATION ###
######################

# Load libraries #
library(tidyverse)
library(ggridges)
library(cowplot)

# Set working directory #
setwd("~/Desktop/general/projects/clients/2025/002_steven_toro/")

# Import some plotting themes and functions #
source("storo-06-2025/functions/functions.R")
source("storo-06-2025/other/lists_and_vectors.R")
load("Data/storo_rdata/storo_s05_clusterProfiler_enrichment.RData")

# Set seed for reproducibility #
set.seed(123)

#############################
# NES RIDGELINE + DOT PLOTS #
#############################

tmp.go.nes <- data.frame()
tmp.kg.nes <- data.frame()

# For each contrast... #
for (ctr in deContrasts) {

  if ( dim(data.frame(goGSEnr[[ctr]]))[[1]] > 0 ) {
    tmp.goGSEnr <- data.frame(goGSEnr[[ctr]])
  } else {
    next
  }

  if ( dim(data.frame(kgGSEnr[[ctr]]))[[1]] > 0 ) {
    tmp.kgGSEnr <- data.frame(kgGSEnr[[ctr]])
  } else {
    next
  }

  tmp.goGSEnr$Contrast <- ctr
  tmp.kgGSEnr$Contrast <- ctr

  tmp.goGSEnr$log.padjust <- -log10(tmp.goGSEnr$p.adjust)
  tmp.kgGSEnr$log.padjust <- -log10(tmp.kgGSEnr$p.adjust)

  tmp.go.nes <- rbind(tmp.go.nes, tmp.goGSEnr)
  tmp.kg.nes <- rbind(tmp.kg.nes, tmp.kgGSEnr)
}

sub.go.nes <- rbind(head(tmp.go.nes[order(tmp.go.nes$NES, decreasing = TRUE),], 20),
                    head(tmp.go.nes[order(tmp.go.nes$NES, decreasing = FALSE),], 20))

sub.kg.nes <- rbind(head(tmp.kg.nes[order(tmp.kg.nes$NES, decreasing = TRUE),], 20),
                    head(tmp.kg.nes[order(tmp.kg.nes$NES, decreasing = FALSE),], 20))

go.nes.rdg.plot <- ggplot(tmp.go.nes, aes(x = NES, y = Contrast, fill = Contrast)) +
  geom_density_ridges(alpha = 0.7) +
  theme(plot.title = element_text(size = 19, hjust = 0.5),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 13),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14))

kg.nes.rdg.plot <- ggplot(tmp.kg.nes, aes(x = NES, y = Contrast, fill = Contrast)) +
  geom_density_ridges(alpha = 0.7) +
  theme(plot.title = element_text(size = 19, hjust = 0.5),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 13),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14))

go.nes.dot.plot <- ggplot(sub.go.nes,
  aes(x = NES, y = reorder(Description, NES), color = Contrast, size = setSize, alpha = log.padjust)) +
  geom_point() +
  ylab("GO Term") +
  labs(alpha = "-log10(FDR)") +
  scale_alpha() +
  theme_bw() +
  theme(plot.title = element_text(size = 19, hjust = 0.5),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 13),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14))

kg.nes.dot.plot <- ggplot(sub.kg.nes,
  aes(x = NES, y = reorder(Description, NES), color = Contrast, size = setSize, alpha = log.padjust)) +
  geom_point() +
  ylab("KEGG Pathway") +
  labs(alpha = "-log10(FDR)") +
  scale_alpha() +
  theme_bw() +
  theme(plot.title = element_text(size = 19, hjust = 0.5),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 13),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14))

go.plot_title <- ggdraw() + draw_label("Normalized Enrichment Score Distributions - GO (GSEA)", size = 20)
kg.plot_title <- ggdraw() + draw_label("Normalized Enrichment Score Distributions - KEGG (GSEA)", size = 20)

go.nes.plot <- plot_grid(go.plot_title,
                           plot_grid(go.nes.rdg.plot, go.nes.dot.plot, nrow = 2),
                         nrow = 2, rel_heights=c(0.1, 1))

kg.nes.plot <- plot_grid(kg.plot_title,
                           plot_grid(kg.nes.rdg.plot, kg.nes.dot.plot, nrow = 2),
                         nrow = 2, rel_heights=c(0.1, 1))

# Save the plots to PDFs #
ggsave("Normalized Enrichment Score Distributions - GO (GSEA).pdf", plot = go.nes.plot, device = "pdf",
       path = "Figures/clusterprofiler/nes_ridgeline_dotplots/",
       width = 16, height = 12, units = "in", dpi = 400)

ggsave("Normalized Enrichment Score Distributions - KEGG (GSEA).pdf", plot = kg.nes.plot, device = "pdf",
       path = "Figures/clusterprofiler/nes_ridgeline_dotplots/",
       width = 16, height = 12, units = "in", dpi = 400)
