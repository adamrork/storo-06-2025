###############
### PURPOSE ###
###############

# This script generates a set of GO enrichment dot #
# plots based on the results of the  goseq + rrvgo #
# analyses. Here, we only plot the terms clustered #
# at threshold = 0.9, as it seems most informative. #

######################
### INITIALIZATION ###
######################

# Load libraries #
library(tidyverse)
library(RColorBrewer)
library(ggplot2)

# Set working directory #
setwd("~/Desktop/general/projects/clients/2025/002_steven_toro/")

# Import some plotting themes and functions #
source("Data/storo_lists/lists_and_vectors.R")
source("Scripts/functions/functions.R")
load("Data/storo_rdata/storo_s04_goseq-rrvgo_enrichment.RData")

# Set seed for reproducibility #
set.seed(123)

#########################
### FIGURE GENERATION ###
#########################

# Create a color palette (Set1, but with a darker yellow) #
colors <- brewer.pal(n = 9, "Set1")
colors[6] <- "#E9D92E"

# For each direction of differential expression... #
for (i in seq_along(regDirection)) {

  # Create a list of contrasts having enriched GO Terms #
  ctrVector <- vector()
  dir <- regDirection[[i]]
  ctrVector <- names(rrvgoData.lst[[dir]])

  # For each such contrast... #
  for (ctr in ctrVector) {

    # For each GO category... #
    for (ont in goCategories) {

      # Assign the data frame to a temporary variable #
      tmp = rrvgoData.lst[[dir]][[ctr]][[ont]]$`90`

      # Create an add-on to the plot title #
      tmp.addition <- paste0(" - Top Enriched ", ont, " Terms")

      # Dynamically generate a title for each plot #
      if ( dir == "UP" ) {

        tmp.name <- renameContrast(name = ctr, addition = tmp.addition, flip = FALSE) %>%
                      gsub("[ ]+", " ", .)

      } else {

        tmp.name <- renameContrast(name = ctr, addition = tmp.addition, flip = TRUE) %>%
          gsub("[ ]+", " ", .)

      }

      # Get our subset of enriched terms #
      filteredTerms <- head(tmp[which(tmp$parentTerm != "NANA"),], n = 25)

      # Discretely adjust several plot size parameters based on how many terms we have. #
      # Most parameters should remain unchanged but are included in case they are needed. #
      # NOTE: Do not use continuous scaling ; it performs poorly relative to discrete scaling. #

      if ( dim(filteredTerms)[1] == 25 ) {
        ptitle = 18 ; ltitle = 12; xtitle = 16 ; ytitle = 16
        xtext = 10 ; ytext = 10 ; ltext = 10.5 ; csize = 3
        pwidth = 10 ; pheight = 10
      } else if (dim(filteredTerms)[1] > 20) {
        ptitle = 16 ; ltitle = 11 ; xtitle = 15 ; ytitle = 15
        xtext = 9 ; ytext = 9 ; ltext = 10 ; csize = 2.5
        pwidth = 10 ; pheight = 9
      } else if (dim(filteredTerms)[1] > 15) {
        ptitle = 16 ; ltitle = 11 ; xtitle = 15 ; ytitle = 15
        xtext = 9 ; ytext = 9 ; ltext = 10 ; csize = 2.5
        pwidth = 10 ; pheight = 8
      } else if (dim(filteredTerms)[1] > 10) {
        ptitle = 16 ; ltitle = 11 ; xtitle = 15 ; ytitle = 15
        xtext = 9 ; ytext = 9 ; ltext = 10 ; csize = 2.5
        pwidth = 10 ; pheight = 7
      } else if (dim(filteredTerms)[1] > 5) {
        ptitle = 16 ; ltitle = 11 ; xtitle = 15 ; ytitle = 15
        xtext = 9 ; ytext = 9 ; ltext = 10 ; csize = 2.5
        pwidth = 10 ; pheight = 6
      } else {
        ptitle = 16 ; ltitle = 11 ; xtitle = 15 ; ytitle = 15
        xtext = 9 ; ytext = 9 ; ltext = 10 ; csize = 2.5
        pwidth = 10 ; pheight = 5
      }

      # Create a GO Enrichment dot plot #
      goseq.rrvgo.plot <- filteredTerms %>%
        arrange(across(.cols=c("parentTerm", "FDR"))) %>%
        rowid_to_column() %>%
        ggplot(aes(x = -log10(FDR), y = reorder(str_wrap(term, width = 50), rowid),
                   colour = str_wrap(parentTerm, width = 50), size = numDEInCat)) +
          geom_point() +
          scale_color_manual(values = colors) +
          ggtitle(tmp.name) +
          labs(x = "-log10(FDR)", y = "GO Term",
               colour = "Representative GO Term",
               size = "Upregulated Genes \nin Category") +
          theme_light() +
          theme(plot.title = element_text(hjust = 0.5, size = ptitle),
                axis.title.x = element_text(size = xtitle),
                axis.title.y = element_text(size = ytitle),
                axis.text.x = element_text(size = xtext),
                axis.text.y = element_text(size = ytext),
                legend.title = element_text(size = ltitle),
                legend.text = element_text(size = ltext)) +
          guides(colour = guide_legend(override.aes = list(size = csize)))

      # Create a new filename after the plot title #
      filename <- paste0(tmp.name, ".pdf")

      # Save the plot to a PDF #
      ggsave(filename, plot = goseq.rrvgo.plot, device = "pdf",
             path = paste0("Figures/goseq_and_rrvgo/clustered_dotplots/", tolower(ont), "/"),
             width = pwidth, height = pheight, units = "in", dpi = 400)
    }
  }
}

# Clear the environment for the next project #
rm(list = ls())
