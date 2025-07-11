###############
### PURPOSE ###
###############

# This script creates several barplots from the results #
# of clusterProfiler functional enrichment analyses. #

######################
### INITIALIZATION ###
######################

# Load libraries #
library(tidyverse)
library(ggnewscale)

# Set working directory #
setwd("~/Desktop/general/projects/clients/2025/002_steven_toro/")

# Import some plotting themes and functions #
source("storo-06-2025/functions/functions.R")
source("storo-06-2025/other/lists_and_vectors.R")
load("Data/storo_rdata/storo_s05_clusterProfiler_enrichment.RData")

# Set seed for reproducibility #
set.seed(123)

#############
# BAR PLOTS #
#############

for ( ctr in deContrasts ) {

  for ( strategy in annotStrategies ) {

    for ( reg in regDirection ) {

      # Assign a data frame to a temporary variable #
      if ( strategy == "goEnr" ) {

        enrich.df <- data.frame(goEnr[[reg]][[ctr]])
        method <- "Enr"
        type <- "GO"
        filepath <- "Figures/clusterprofiler/barplots/enrichment/go/"

      } else if ( strategy == "kgEnr" ) {

        enrich.df <- data.frame(kgEnr[[reg]][[ctr]])
        method <- "Enr"
        type <- "KEGG"
        filepath <- "Figures/clusterprofiler/barplots/enrichment/kegg/"

      } else if ( strategy == "goGSEnr" ) {

        enrich.df <- data.frame(goGSEnr[[ctr]])
        method <- "GSEnr"
        type <- "GO"
        filepath <- "Figures/clusterprofiler/barplots/gsea/go/"

      } else if ( strategy == "kgGSEnr" ) {

        enrich.df <- data.frame(kgGSEnr[[ctr]])
        method <- "GSEnr"
        type <- "KEGG"
        filepath <- "Figures/clusterprofiler/barplots/gsea/kegg/"
      }

      # Only parse GSEA data once #
      if ( method == "GSEnr" && reg == "DOWN" ) {
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

      # If we are plotting GO data... #
      if ( strategy == "goEnr" ) {

        # Create a base dot plot showing the top 25 terms sorted by RichFactor #
        enrich.df <- slice_max(enrich.df, order_by = RichFactor, n = 25, with_ties = FALSE)
        enr.barplot <- ggplot(enrich.df, aes(x = RichFactor, y = reorder(Description, RichFactor), fill = log.padjust)) +
          xlab("RichFactor") + ggtitle(tmp.title)

        # Create a base bar plot showing the top terms per category sorted by RichFactor #
        grp.enrich.df <- slice_max(group_by(enrich.df, ONTOLOGY), order_by = RichFactor, n = 25, with_ties = FALSE)
        grp.enrich.df <- drop_na(grp.enrich.df)

        grp.enr.barplot <- ggplot(grp.enrich.df, aes(x = RichFactor, y = reorder(Description, RichFactor), fill = log.padjust)) +
          xlab("RichFactor") + ggtitle(tmp.title)

        # Create a bar plot with a rich set of annotations #
        annot.enrich.df <- slice_max(group_by(enrich.df, ONTOLOGY), order_by = RichFactor, n = 25, with_ties = FALSE) %>%
          mutate(Description = factor(Description, Description))
        annot.enrich.df <- drop_na(annot.enrich.df)
        annot.enrich.df$ProportionDE <- sapply(grp.enrich.df$GeneRatio, function(x) eval(parse(text = x)))

        annot.enr.barplot <- ggplot(annot.enrich.df, aes(x = RichFactor, y = Description, fill = ONTOLOGY)) +
          xlab("RichFactor") + ggtitle(tmp.title) + labs(fill = "Ontology")

        xmax <- max(annot.enrich.df$RichFactor)

        # Create base bar plots for each category showing the top 25 terms per category, sorted by RichFactor #
        bp.enrich.df <- slice_max(enrich.df, order_by = RichFactor, n = 25, with_ties = FALSE) %>%
          filter(ONTOLOGY == "BP")
        mf.enrich.df <- slice_max(enrich.df, order_by = RichFactor, n = 25, with_ties = FALSE) %>%
          filter(ONTOLOGY == "MF")
        cc.enrich.df <- slice_max(enrich.df, order_by = RichFactor, n = 25, with_ties = FALSE) %>%
          filter(ONTOLOGY == "CC")

        bp.enr.barplot <- ggplot(bp.enrich.df, aes(x = RichFactor, y = reorder(Description, RichFactor), fill = log.padjust)) +
          xlab("RichFactor") + ggtitle(bp.tmp.title)
        mf.enr.barplot <- ggplot(mf.enrich.df, aes(x = RichFactor, y = reorder(Description, RichFactor), fill = log.padjust)) +
          xlab("RichFactor") + ggtitle(mf.tmp.title)
        cc.enr.barplot <- ggplot(cc.enrich.df, aes(x = RichFactor, y = reorder(Description, RichFactor), fill = log.padjust)) +
          xlab("RichFactor") + ggtitle(cc.tmp.title)

      # If we are plotting KEGG enrichment data... #
      } else if ( strategy == "kgEnr" ) {

        # Create a base dot plot showing the top 25 terms sorted by RichFactor #
        enrich.df <- slice_max(enrich.df, order_by = RichFactor, n = 25, with_ties = FALSE)
        enr.barplot <- ggplot(enrich.df, aes(x = RichFactor, y = reorder(Description, RichFactor), fill = log.padjust)) +
          xlab("RichFactor") + ggtitle(tmp.title)

        # Create a base bar plot showing the top terms per category sorted by RichFactor #
        grp.enrich.df <- slice_max(group_by(enrich.df, category), order_by = RichFactor, n = 25, with_ties = FALSE)
        grp.enrich.df <- grp.enrich.df[!is.na(grp.enrich.df$category),]

        grp.enr.barplot <- ggplot(grp.enrich.df, aes(x = RichFactor, y = reorder(Description, RichFactor), fill = log.padjust)) +
          xlab("RichFactor") + ggtitle(tmp.title)

        # Create a bar plot with a rich set of annotations #
        annot.enrich.df <- slice_max(group_by(enrich.df, category), order_by = RichFactor, n = 25, with_ties = FALSE) %>%
          mutate(Description = factor(Description, Description))
        annot.enrich.df <- annot.enrich.df[!is.na(annot.enrich.df$category),]
        annot.enrich.df$ProportionDE <- sapply(grp.enrich.df$GeneRatio, function(x) eval(parse(text = x)))

        annot.enr.barplot <- ggplot(annot.enrich.df, aes(x = RichFactor, y = Description, fill = category)) +
          xlab("RichFactor") + ggtitle(tmp.title) + labs(fill = "Category")

        xmax <- max(grp.enrich.df$RichFactor)

      # If we are plotting GO gene set enrichment data... #
      } else if ( strategy == "goGSEnr" ) {

        # Create a base dot plot showing the top 25 terms sorted by enrichmentScore #
        enrich.df <- slice_max(enrich.df, order_by = enrichmentScore, n = 25, with_ties = FALSE)
        enr.barplot <- ggplot(enrich.df, aes(x = enrichmentScore, y = reorder(Description, enrichmentScore), fill = log.padjust)) +
          xlab("RichFactor") + ggtitle(tmp.title)

        # Create a base bar plot showing the top terms per category sorted by enrichmentScore #
        grp.enrich.df <- slice_max(group_by(enrich.df, ONTOLOGY), order_by = enrichmentScore, n = 25, with_ties = FALSE)
        grp.enrich.df <- drop_na(grp.enrich.df)

        grp.enr.barplot <- ggplot(grp.enrich.df, aes(x = enrichmentScore, y = reorder(Description, enrichmentScore), fill = log.padjust)) +
          xlab("enrichmentScore") + ggtitle(tmp.title)

        # Create base bar plots for each category showing the top 25 terms per category, sorted by enrichmentScore #
        bp.enrich.df <- slice_max(enrich.df, order_by = enrichmentScore, n = 25, with_ties = FALSE) %>%
          filter(ONTOLOGY == "BP")
        mf.enrich.df <- slice_max(enrich.df, order_by = enrichmentScore, n = 25, with_ties = FALSE) %>%
          filter(ONTOLOGY == "MF")
        cc.enrich.df <- slice_max(enrich.df, order_by = enrichmentScore, n = 25, with_ties = FALSE) %>%
          filter(ONTOLOGY == "CC")

        bp.enr.barplot <- ggplot(bp.enrich.df, aes(x = enrichmentScore, y = reorder(Description, enrichmentScore), fill = log.padjust)) +
          xlab("enrichmentScore") + ggtitle(bp.tmp.title)
        mf.enr.barplot <- ggplot(mf.enrich.df, aes(x = enrichmentScore, y = reorder(Description, enrichmentScore), fill = log.padjust)) +
          xlab("enrichmentScore") + ggtitle(mf.tmp.title)
        cc.enr.barplot <- ggplot(cc.enrich.df, aes(x = enrichmentScore, y = reorder(Description, enrichmentScore), fill = log.padjust)) +
          xlab("enrichmentScore") + ggtitle(cc.tmp.title)

      # If we are plotting KEGG gene set enrichment data... #
      } else if ( strategy == "kgGSEnr" ) {

        # Create a base bar plot showing the top 25 pathways sorted by enrichmentScore #
        enrich.df <- slice_max(enrich.df, order_by = enrichmentScore, n = 25, with_ties = FALSE)

        enr.barplot <- ggplot(enrich.df, aes(x = enrichmentScore, y = reorder(Description, enrichmentScore), fill = log.padjust)) +
          xlab("enrichmentScore") + ggtitle(tmp.title)
      }

      # If we are plotting GO data... #
      if ( type == "GO" ) {

        # Create fully annotated bar plots #
        bp.enr.barplot <- bp.enr.barplot +
          geom_bar(stat = "identity") +
          ylab("Enriched Term") +
          labs(fill = "-log10(FDR)") +
          scale_fill_gradient2(midpoint = med.logfdr, low = "#FFDBD3", mid = "#E86E79", high = "#D1001F") +
          scale_y_discrete(labels = function(Description) str_wrap(Description, width = 50)) +
          scale_x_continuous(expand = c(0.005, 0)) +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5, size = 22),
                axis.title = element_text(size = 19),
                axis.text = element_text(size = 12),
                legend.title = element_text(size = 14),
                legend.text = element_text(size = 14))

        mf.enr.barplot <- mf.enr.barplot +
          geom_bar(stat = "identity") +
          ylab("Enriched Term") +
          labs(fill = "-log10(FDR)") +
          scale_fill_gradient2(midpoint = med.logfdr, low = "#FFDBD3", mid = "#E86E79", high = "#D1001F") +
          scale_y_discrete(labels = function(Description) str_wrap(Description, width = 50)) +
          scale_x_continuous(expand = c(0.005, 0)) +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5, size = 22),
                axis.title = element_text(size = 19),
                axis.text = element_text(size = 12),
                legend.title = element_text(size = 14),
                legend.text = element_text(size = 14))

        cc.enr.barplot <- cc.enr.barplot +
          geom_bar(stat = "identity") +
          ylab("Enriched Term") +
          labs(fill = "-log10(FDR)") +
          scale_fill_gradient2(midpoint = med.logfdr, low = "#FFDBD3", mid = "#E86E79", high = "#D1001F") +
          scale_y_discrete(labels = function(Description) str_wrap(Description, width = 50)) +
          scale_x_continuous(expand = c(0.005, 0)) +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5, size = 22),
                axis.title = element_text(size = 19),
                axis.text = element_text(size = 12),
                legend.title = element_text(size = 14),
                legend.text = element_text(size = 14))

        # Create file names #
        if ( method == "Enr" ) {

          filename <- paste0(tmp.title, " (Enr) Barplots.pdf")
          bp.filename <- paste0(bp.tmp.title, " (Enr) Barplots.pdf")
          mf.filename <- paste0(mf.tmp.title, " (Enr) Barplots.pdf")
          cc.filename <- paste0(cc.tmp.title, " (Enr) Barplots.pdf")

        } else if ( method == "GSEnr" ) {
          filename <- paste0(tmp.title, " (GSEA) Barplots.pdf")
          bp.filename <- paste0(bp.tmp.title, " (GSEA) Barplots.pdf")
          mf.filename <- paste0(mf.tmp.title, " (GSEA) Barplots.pdf")
          cc.filename <- paste0(cc.tmp.title, " (GSEA) Barplots.pdf")
        }

        # Save each bar plot to a PDF #
        ggsave(bp.filename, plot = bp.enr.barplot, device = "pdf",
               path = paste0(filepath, "bp/"),
               width = 14, height = 12, units = "in", dpi = 400)

        ggsave(mf.filename, plot = mf.enr.barplot, device = "pdf",
               path = paste0(filepath, "mf/"),
               width = 14, height = 12, units = "in", dpi = 400)

        ggsave(cc.filename, plot = cc.enr.barplot, device = "pdf",
               path = paste0(filepath, "cc/"),
               width = 14, height = 12, units = "in", dpi = 400)

        # Also, create a facet bar plot #
        grp.enr.barplot <- grp.enr.barplot +
          geom_bar(stat = "identity") +
          ylab("Enriched Term") +
          labs(fill = "-log10(FDR)") +
          scale_fill_gradient2(midpoint = med.logfdr, low = "#FFDBD3", mid = "#E86E79", high = "#D1001F") +
          scale_y_discrete(labels = function(Description) str_wrap(Description, width = 50)) +
          scale_x_continuous(expand = c(0.005, 0)) +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5, size = 22),
                axis.title = element_text(size = 19),
                axis.text = element_text(size = 12),
                legend.title = element_text(size = 14),
                legend.text = element_text(size = 14))

        go.facet.barplot <- grp.enr.barplot +
          facet_grid(vars(ONTOLOGY), scales = "free", space = "free_y") +
          theme(strip.text = element_text(size = 12))

        # Save the facet bar plot to a PDF #
        ggsave(filename, plot = go.facet.barplot, device = "pdf",
              path = paste0(filepath, "all/"),
              width = 14, height = 12, units = "in", dpi = 400)

        # If we are plotting KEGG data... #
      } else if ( type == "KEGG" ) {

        # Create fully annotated bar plots #
        enr.barplot <- enr.barplot +
          geom_bar(stat = "identity") +
          ylab("Enriched Term") +
          labs(fill = "-log10(FDR)") +
          scale_fill_gradient2(midpoint = med.logfdr, low = "#FFDBD3", mid = "#E86E79", high = "#D1001F") +
          scale_y_discrete(labels = function(Description) str_wrap(Description, width = 50)) +
          scale_x_continuous(expand = c(0.005, 0)) +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5, size = 22),
                axis.title = element_text(size = 19),
                axis.text = element_text(size = 12),
                legend.title = element_text(size = 14),
                legend.text = element_text(size = 14))

        # Create file names #
        if ( method == "Enr" ) {

          filename <- paste0(tmp.title, " (Enr) Barplots.pdf")

        } else if ( method == "GSEnr" ) {

          filename <- paste0(tmp.title, " (GSEA) Barplots.pdf")
        }

        # Save the plot to a PDF #
        ggsave(filename, plot = enr.barplot, device = "pdf",
               path = paste0(filepath, "all/"),,
               width = 14, height = 12, units = "in", dpi = 400)
      }

      if ( method == "Enr" ) {
        # Annotate our existing enrichment bar plot #
        annot.enr.barplot <- annot.enr.barplot +
          geom_bar(stat = "identity") +
          scale_y_discrete(labels = function(Description) str_wrap(Description, width = 50)) +
          scale_x_continuous(expand = c(0.005, 0), limits = c(0, xmax + 0.1)) +
          scale_fill_brewer(palette = "Dark2") +
          # Start a new scale for -log10(FDR) tiles #
          new_scale_fill() +
          geom_tile(aes(x = xmax + 0.02, y = Description, fill = log.padjust, width = 0.03)) +
          scale_fill_gradient(low = "#FFDBD3", high = "#D1001F", na.value = NA) +
          geom_text(aes(x = xmax + 0.02, y = Description, label = round(log.padjust, 2)),
                    size = 3, vjust = 0.5, angle = 0, fontface = "bold") +
          labs(fill = "-log10(FDR)") +
          # Start a new scale for DE gene proportion tiles #
          new_scale_fill() +
          geom_tile(aes(x = xmax + 0.05, y = Description, fill = ProportionDE), width = 0.03) +
          scale_fill_gradient(low = "#D3F3FF", high = "#0095D1", na.value = NA) +
          geom_text(aes(x = xmax + 0.05, y = Description, label = round(ProportionDE, 2)),
                    size = 3, vjust = 0.5, angle = 0, fontface = "bold") +
          labs(fill = "Proportion DE") +
          # Start a new scale for fold enrichment proportion tiles #
          new_scale_fill() +
          geom_tile(aes(x = xmax + 0.08, y = Description, fill = FoldEnrichment), width = 0.03) +
          scale_fill_gradient(low = "#FFFFAC", high = "#E7E732", na.value = NA) +
          geom_text(aes(x = xmax + 0.08, y = Description, label = round(FoldEnrichment, 2)),
                    size = 3, vjust = 0.5, angle = 0, fontface = "bold") +
          labs(fill = "Fold Enrichment") +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5, size = 22),
                axis.title = element_text(size = 19),
                axis.text = element_text(size = 12),
                legend.title = element_text(size = 14),
                legend.text = element_text(size = 14))

        # Save the plot to a PDF #
        ggsave(filename, plot = annot.enr.barplot, device = "pdf",
               path = paste0(filepath, "rich/"),
               width = 14, height = 12, units = "in", dpi = 400)
      }
    }
  }
}
