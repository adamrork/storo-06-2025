###############
### PURPOSE ###
###############

# This script creates several dotplots from the results #
# of clusterProfiler functional enrichment analyses. #

######################
### INITIALIZATION ###
######################

# Load libraries #
library(tidyverse)

# Set working directory #
setwd("~/Desktop/general/projects/clients/2025/002_steven_toro/")

# Import some plotting themes and functions #
source("Scripts/functions/functions.R")
source("Data/storo_lists/lists_and_vectors.R")
load("Data/storo_rdata/storo_s05_clusterProfiler_enrichment.RData")

# Set seed for reproducibility #
set.seed(123)

#############
# DOT PLOTS #
#############

for ( ctr in deContrasts ) {

  for ( strategy in annotStrategies ) {

    for ( reg in regDirection ) {

      # Assign a data frame to a temporary variable #
      if ( strategy == "goEnr" ) {

        enrich.df <- data.frame(goEnr[[reg]][[ctr]])
        method <- "Enr"
        type <- "GO"
        filepath <- "Figures/clusterprofiler/dotplots/enrichment/go/"

      } else if ( strategy == "kgEnr" ) {

        enrich.df <- data.frame(kgEnr[[reg]][[ctr]])
        method <- "Enr"
        type <- "KEGG"
        filepath <- "Figures/clusterprofiler/dotplots/enrichment/kegg/"

      } else if ( strategy == "goGSEnr" ) {

        enrich.df <- data.frame(goGSEnr[[ctr]])
        method <- "GSEnr"
        type <- "GO"
        filepath <- "Figures/clusterprofiler/dotplots/gsea/go/"

      } else if ( strategy == "kgGSEnr" ) {

        enrich.df <- data.frame(kgGSEnr[[ctr]])
        method <- "GSEnr"
        type <- "KEGG"
        filepath <- "Figures/clusterprofiler/dotplots/gsea/kegg/"

      }

      # Only parse GSEA data once #
      if ( method == "GSEnr" && reg == "DOWN") {
        break
      }

      # Also break if there is no data to analyze #
      if ( dim(enrich.df)[1] > 0 ) {

        # Generate -log10(FDR) values and identify a midpoint for our -log10(FDR) scale #
        enrich.df$log.padjust <- -log10(enrich.df$p.adjust)
        med.logfdr <- median(enrich.df$log.padjust)

      } else {
        break
      }

      # Dynamically create some plot titles #
      if ( type == "GO") {

        tmp.addition <- paste0(" - Top Enriched GO Terms")
        bp.tmp.addition <- paste0(" - Top Enriched BP Terms")
        mf.tmp.addition <- paste0(" - Top Enriched MF Terms")
        cc.tmp.addition <- paste0(" - Top Enriched CC Terms")

        if ( reg == "UP" ) {
          tmp.title <- renameContrast(name = ctr, addition = tmp.addition, flip = FALSE) %>% gsub("[ ]+", " ", .)
          bp.tmp.title <- renameContrast(name = ctr, addition = bp.tmp.addition, flip = FALSE) %>% gsub("[ ]+", " ", .)
          mf.tmp.title <- renameContrast(name = ctr, addition = mf.tmp.addition, flip = FALSE) %>% gsub("[ ]+", " ", .)
          cc.tmp.title <- renameContrast(name = ctr, addition = cc.tmp.addition, flip = FALSE) %>% gsub("[ ]+", " ", .)

        } else if ( reg == "DOWN" ) {
          tmp.title <- renameContrast(name = ctr, addition = tmp.addition, flip = TRUE) %>% gsub("[ ]+", " ", .)
          bp.tmp.title <- renameContrast(name = ctr, addition = bp.tmp.addition, flip = TRUE) %>% gsub("[ ]+", " ", .)
          mf.tmp.title <- renameContrast(name = ctr, addition = mf.tmp.addition, flip = TRUE) %>% gsub("[ ]+", " ", .)
          cc.tmp.title <- renameContrast(name = ctr, addition = cc.tmp.addition, flip = TRUE) %>% gsub("[ ]+", " ", .)
        }

      } else if ( type == "KEGG" ) {

        tmp.addition <- paste0(" - Top Enriched KEGG Pathways")

        if ( reg == "UP" ) {
          tmp.title <- renameContrast(name = ctr, addition = tmp.addition, flip = FALSE) %>% gsub("[ ]+", " ", .)

        } else if ( reg == "DOWN" ) {
          tmp.title <- renameContrast(name = ctr, addition = tmp.addition, flip = TRUE) %>% gsub("[ ]+", " ", .)
        }
      }

      # If we are plotting GO enrichment data... #
      if ( strategy == "goEnr" ) {

        # Create a base dot plot showing the top 25 terms sorted by RichFactor #
        enrich.df <- slice_max(enrich.df, order_by = RichFactor, n = 25, with_ties = FALSE)
        enr.dotplot <- ggplot(enrich.df, aes(x = RichFactor, y = reorder(Description, RichFactor),
          color = log.padjust, size = Count)) + xlab("RichFactor") + ggtitle(tmp.title)

        # Create a base dot plot showing the Top terms per category sorted by RichFactor #
        grp.enrich.df <- slice_max(group_by(enrich.df, ONTOLOGY), order_by = RichFactor, n = 25, with_ties = FALSE)
        grp.enrich.df <- drop_na(grp.enrich.df)

        grp.enr.dotplot <- ggplot(grp.enrich.df, aes(x = RichFactor, y = reorder(Description, RichFactor),
          color = log.padjust, size = Count)) + xlab("RichFactor") + ggtitle(tmp.title)

        # Create base dot plots for each category showing the top 25 terms per category, sorted by RichFactor #
        bp.enrich.df <- slice_max(enrich.df, order_by = RichFactor, n = 25, with_ties = FALSE) %>%
          filter(ONTOLOGY == "BP")
        mf.enrich.df <- slice_max(enrich.df, order_by = RichFactor, n = 25, with_ties = FALSE) %>%
          filter(ONTOLOGY == "MF")
        cc.enrich.df <- slice_max(enrich.df, order_by = RichFactor, n = 25, with_ties = FALSE) %>%
          filter(ONTOLOGY == "CC")

        bp.enr.dotplot <- ggplot(bp.enrich.df, aes(x = RichFactor, y = reorder(Description, RichFactor),
          color = log.padjust, size = Count)) + xlab("RichFactor") + ggtitle(bp.tmp.title)
        mf.enr.dotplot <- ggplot(mf.enrich.df, aes(x = RichFactor, y = reorder(Description, RichFactor),
          color = log.padjust, size = Count)) + xlab("RichFactor") + ggtitle(mf.tmp.title)
        cc.enr.dotplot <- ggplot(cc.enrich.df, aes(x = RichFactor, y = reorder(Description, RichFactor),
          color = log.padjust, size = Count)) + xlab("RichFactor") + ggtitle(cc.tmp.title)

      # If we are plotting KEGG enrichment data... #
      } else if ( strategy == "kgEnr" ) {

        # Create a base dot plot showing the top 25 pathways sorted by RichFactor #
        enrich.df <- slice_max(enrich.df, order_by = RichFactor, n = 25, with_ties = FALSE)

        enr.dotplot <- ggplot(enrich.df, aes(x = RichFactor, y = reorder(Description, RichFactor),
          color = log.padjust, size = Count)) + xlab("RichFactor") + ggtitle(tmp.title)

        # Create a base dot plot showing the Top pathways per category sorted by RichFactor #
        grp.enrich.df <- slice_max(group_by(enrich.df, category), order_by = RichFactor, n = 25, with_ties = FALSE)
        grp.enrich.df <- grp.enrich.df[!is.na(grp.enrich.df$category),]

        grp.enr.dotplot <- ggplot(grp.enrich.df, aes(x = RichFactor, y = reorder(Description, RichFactor),
                                             color = log.padjust, size = Count)) + xlab("RichFactor") + ggtitle(tmp.title)

      # If we are plotting GO gene set enrichment data... #
      } else if ( strategy == "goGSEnr" ) {

        # Create a base dot plot showing the top 25 terms sorted by enrichmentScore #
        enrich.df <- slice_max(enrich.df, order_by = enrichmentScore, n = 25, with_ties = FALSE)

        enr.dotplot <- ggplot(enrich.df, aes(x = enrichmentScore, y = reorder(Description, enrichmentScore),
                                             color = log.padjust, size = setSize)) + xlab("enrichmentScore") + ggtitle(tmp.title)

        # Create a base dot plot showing the Top terms per category sorted by enrichmentScore #
        grp.enrich.df <- slice_max(group_by(enrich.df, ONTOLOGY), order_by = enrichmentScore, n = 25, with_ties = FALSE)
        grp.enrich.df <- drop_na(grp.enrich.df)

        grp.enr.dotplot <- ggplot(grp.enrich.df, aes(x = enrichmentScore, y = reorder(Description, enrichmentScore),
                                                     color = log.padjust, size = setSize)) + xlab("enrichmentScore") + ggtitle(tmp.title)

        # Create base dot plots for each category showing the top 25 terms per category, sorted by enrichmentScore #
        bp.enrich.df <- slice_max(enrich.df, order_by = enrichmentScore, n = 25, with_ties = FALSE) %>%
          filter(ONTOLOGY == "BP")
        mf.enrich.df <- slice_max(enrich.df, order_by = enrichmentScore, n = 25, with_ties = FALSE) %>%
          filter(ONTOLOGY == "MF")
        cc.enrich.df <- slice_max(enrich.df, order_by = enrichmentScore, n = 25, with_ties = FALSE) %>%
          filter(ONTOLOGY == "CC")

        bp.enr.dotplot <- ggplot(bp.enrich.df, aes(x = enrichmentScore, y = reorder(Description, enrichmentScore),
          color = log.padjust, size = setSize)) + xlab("enrichmentScore") + ggtitle(bp.tmp.title)
        mf.enr.dotplot <- ggplot(mf.enrich.df, aes(x = enrichmentScore, y = reorder(Description, enrichmentScore),
          color = log.padjust, size = setSize)) + xlab("enrichmentScore") + ggtitle(mf.tmp.title)
        cc.enr.dotplot <- ggplot(cc.enrich.df, aes(x = enrichmentScore, y = reorder(Description, enrichmentScore),
          color = log.padjust, size = setSize)) + xlab("enrichmentScore") + ggtitle(cc.tmp.title)

      # If we are plotting KEGG gene set enrichment data... #
      } else if ( strategy == "kgGSEnr" ) {

        # Create a base dot plot showing the top 25 pathways sorted by enrichmentScore #
        enrich.df <- slice_max(enrich.df, order_by = enrichmentScore, n = 25, with_ties = FALSE)

        enr.dotplot <- ggplot(enrich.df, aes(x = enrichmentScore, y = reorder(Description, enrichmentScore),
                                             color = log.padjust, size = setSize)) + xlab("enrichmentScore") + ggtitle(tmp.title)
      }

      # If we are plotting GO data... #
      if ( type == "GO" ) {

        # Create fully annotated dot plots, one per category #
        bp.enr.dotplot <- bp.enr.dotplot +
          geom_point() +
          ylab("Enriched Term") + labs(color = "-log10(FDR)") +
          scale_colour_gradient2(midpoint = med.logfdr, low = "#FFDBD3", mid = "#E86E79", high = "#D1001F") +
          scale_y_discrete(labels = function(Description) str_wrap(Description, width = 50)) +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5, size = 22),
                axis.title = element_text(size = 19),
                axis.text = element_text(size = 12),
                legend.title = element_text(size = 14),
                legend.text = element_text(size = 14))

        mf.enr.dotplot <- mf.enr.dotplot +
          geom_point() +
          ylab("Enriched Term") + labs(color = "-log10(FDR)") +
          scale_colour_gradient2(midpoint = med.logfdr, low = "#FFDBD3", mid = "#E86E79", high = "#D1001F") +
          scale_y_discrete(labels = function(Description) str_wrap(Description, width = 50)) +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5, size = 22),
                axis.title = element_text(size = 19),
                axis.text = element_text(size = 12),
                legend.title = element_text(size = 14),
                legend.text = element_text(size = 14))

        cc.enr.dotplot <- cc.enr.dotplot +
          geom_point() +
          ylab("Enriched Term") + labs(color = "-log10(FDR)") +
          scale_colour_gradient2(midpoint = med.logfdr, low = "#FFDBD3", mid = "#E86E79", high = "#D1001F") +
          scale_y_discrete(labels = function(Description) str_wrap(Description, width = 50)) +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5, size = 22),
                axis.title = element_text(size = 19),
                axis.text = element_text(size = 12),
                legend.title = element_text(size = 14),
                legend.text = element_text(size = 14))

        # Create file names #
        if ( method == "Enr" ) {

          filename <- paste0(tmp.title, " (Enr) Dotplots.pdf")
          bp.filename <- paste0(bp.tmp.title, " (Enr) Dotplots.pdf")
          mf.filename <- paste0(mf.tmp.title, " (Enr) Dotplots.pdf")
          cc.filename <- paste0(cc.tmp.title, " (Enr) Dotplots.pdf")

        } else if ( method == "GSEnr" ) {

          filename <- paste0(tmp.title, " (GSEA) Dotplots.pdf")
          bp.filename <- paste0(bp.tmp.title, " (GSEA) Dotplots.pdf")
          mf.filename <- paste0(mf.tmp.title, " (GSEA) Dotplots.pdf")
          cc.filename <- paste0(cc.tmp.title, " (GSEA) Dotplots.pdf")
        }

        # Save each dot plot to a PDF #
        ggsave(bp.filename, plot = bp.enr.dotplot, device = "pdf",
               path = paste0(filepath, "bp/"),
               width = 14, height = 12, units = "in", dpi = 400)

        ggsave(mf.filename, plot = mf.enr.dotplot, device = "pdf",
               path = paste0(filepath, "mf/"),
               width = 14, height = 12, units = "in", dpi = 400)

        ggsave(cc.filename, plot = cc.enr.dotplot, device = "pdf",
               path = paste0(filepath, "cc/"),
               width = 14, height = 12, units = "in", dpi = 400)

        # Also, create a facet dot plot #
        grp.enr.dotplot <- grp.enr.dotplot +
          geom_point() +
          ylab("Enriched Term") +
          labs(color = "-log10(FDR)") +
          scale_colour_gradient2(midpoint = med.logfdr, low = "#FFDBD3", mid = "#E86E79", high = "#D1001F") +
          scale_y_discrete(labels = function(Description) str_wrap(Description, width = 50)) +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5, size = 22),
                axis.title = element_text(size = 19),
                axis.text = element_text(size = 12),
                legend.title = element_text(size = 14),
                legend.text = element_text(size = 14))

        go.facet.dotplot <- grp.enr.dotplot +
          facet_grid(vars(ONTOLOGY), scales = "free", space = "free_y") +
          theme(strip.text = element_text(size = 12))

        # Save the facet dot plot to a PDF #
        ggsave(filename, plot = go.facet.dotplot, device = "pdf",
               path = paste0(filepath, "all/"),,
               width = 14, height = 12, units = "in", dpi = 400)

      # If we are plotting KEGG data... #
      } else if ( type == "KEGG" ) {

        # Create fully annotated dot plots #
        enr.dotplot <- enr.dotplot +
          geom_point() +
          ylab("Enriched Term") +
          labs(color = "-log10(FDR)") +
          scale_colour_gradient2(midpoint = med.logfdr, low = "#FFDBD3", mid = "#E86E79", high = "#D1001F") +
          scale_y_discrete(labels = function(Description) str_wrap(Description, width = 50)) +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5, size = 22),
                axis.title = element_text(size = 19),
                axis.text = element_text(size = 12),
                legend.title = element_text(size = 14),
                legend.text = element_text(size = 14))

        # Create file names #
        if ( method == "Enr" ) {

          filename <- paste0(tmp.title, " (Enr) Dotplots.pdf")

        } else if ( method == "GSEnr" ) {

          filename <- paste0(tmp.title, " (GSEA) Dotplots.pdf")
        }

        # Save the plot to a PDF #
        ggsave(filename, plot = enr.dotplot, device = "pdf",
              path = paste0(filepath, "all/"),
              width = 14, height = 12, units = "in", dpi = 400)
      }
    }
  }
}
