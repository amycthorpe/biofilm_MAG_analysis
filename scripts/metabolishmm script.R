# metabolishmm
# figure 4 - carbon, nitrogen and sulfur gene heatmaps
# extended data figure 4 - hydrology and oxygene gene heatmaps


# load libraries
library(ggplot2)
library(dplyr)
library(tidyverse)
library(reshape2)
library(viridisLite)
library(RColorBrewer)
library(patchwork)

set.seed(123)


# import data
metabolishmm <- read.csv("data/metabolishmm_results.csv", row.names=1, check.names=F, sep=",")

gtdb <- read.csv("data/checkm_gtdb.csv", row.names=1, check.names=F, sep=",") %>% 
  select(14:20)

coverage <- read.csv("data/finalbins_coverage.csv", row.names=1, check.names=F, sep=",")

# convert to relative abundance
coverage_rel <- coverage %>%
  mutate(across(everything(), ~ . / sum(.))) %>% 
  t(.) # samples as rows

# calculate mean across all samples
coverage_mean <- as.data.frame(colMeans(coverage_rel))
coverage_mean <- coverage_mean %>% rownames_to_column(., var = "label")
colnames(coverage_mean)[2] <- "mean_rel_abund"


# merge
metabolishmm_tax <- merge(metabolishmm, gtdb, by.y="row.names", by.x="bin")
metabolishmm_merged <- merge(metabolishmm_tax, coverage_mean, by.x="bin", by.y="label")


# calculate weighted phylum counts
metabolishmm_phylum <- metabolishmm_merged %>% 
  mutate(weighted_rel_abund=count*mean_rel_abund) %>% # count weighted by mean relative abundance
  group_by(phylum, gene_short, functional_group) %>% 
  summarise(phylum_count = sum(count), # total count by phylum
            phylum_abund = sum(weighted_rel_abund)) %>% # total phylum count weighted by mean relative abundance
  filter(phylum_abund != 0) # remove zeros for better visualisation

metabolishmm_phylum$phylum <- factor(metabolishmm_phylum$phylum,
                                      levels = rev(unique(metabolishmm_phylum$phylum))) # reverse order for plotting


# set colours
color_gradients <- list(
  carbon = brewer.pal(5, "Blues"), 
  hydrogen = brewer.pal(5, "Purples"), 
  nitrogen = brewer.pal(5, "Greens"),
  sulfur = brewer.pal(5, "Oranges"),
  oxygen = brewer.pal(5, "Greys")
)

ordered_groups <- c("carbon", "nitrogen", "hydrogen", "oxygen", "sulfur")


# plot heatmaps
plots <- lapply(ordered_groups, function(group) {
  group_data <- metabolishmm_phylum %>% filter(functional_group == group)
  
  ggplot(group_data, aes(x = gene_short, y = phylum, fill = phylum_abund)) +
    geom_tile() +
    scale_fill_gradientn(colors = color_gradients[[group]], name = "Weighted Count",
                         limits=c(0,0.4),
                         na.value="white") +
    labs(title = tools::toTitleCase(group)) + 
    theme_bw() +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          legend.title = element_blank(),
          legend.text = element_text(size = 5),
          axis.title = element_blank(),
          axis.text.y = element_text(size = 5, colour = "black"),
          axis.text.x = element_text(size = 5, colour = "black", angle = 45, hjust = 1),
          plot.title = element_blank()) + 
    scale_y_discrete(limits = levels(metabolishmm_phylum$phylum)) 
})

main_plot <- wrap_plots(plots[c(1,2,5)], ncol = 3) +
  plot_layout(axes="collect") &
  theme(legend.position="none")

main_plot

ggsave("output/metabolishmm_CNS.png", main_plot, width = 250, height = 120, units = "mm", dpi = 600)
ggsave("output/metabolishmm_CNS.pdf", main_plot, width = 180, height = 70, units = "mm", dpi = 600)

supp_plot <- wrap_plots(plots[c(3,4)], ncol = 2) +
  plot_layout(axes="collect") &
  theme(legend.position="bottom")

supp_plot

ggsave("output/metabolishmm_HO.png", supp_plot, width = 180, height = 120, units = "mm", dpi = 600)
ggsave("output/metabolishmm_HO.pdf", supp_plot, width = 140, height = 100, units = "mm", dpi = 600)
