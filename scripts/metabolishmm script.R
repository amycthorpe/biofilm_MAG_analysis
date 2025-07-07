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

coverage <- read.csv("data/finalbins_coverage.csv", row.names=1, check.names=F, sep=",")


# filter low quality MAGs
checkm <- read.csv("data/checkm_gtdb.csv", row.names=1, check.names=F, sep=",")

checkm <- mutate(checkm, Completeness_quality = case_when(
  Completeness >= 90 & Contamination <= 5 ~ "Near-complete MAGs",
  Completeness >= 70 & Completeness < 90 & Contamination <= 10 ~ "Medium-quality MAGs",
  TRUE ~ "Low-quality MAGs"
))

gtdb <- checkm %>% 
  subset(., Completeness_quality != "Low-quality MAGs") %>% 
  select(14:20)

coverage <- coverage[rownames(coverage) %in% rownames(gtdb), ]
metabolishmm <- metabolishmm[metabolishmm$bin %in% rownames(gtdb), ]


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


# genus level
# deal with unassigned genera
fill_unassigned <- function(row) {
  last_known <- NA
  for (i in seq_along(row)) {
    if (is.na(row[i]) || row[i] == "") {
      if (!is.na(last_known)) {
        row[i] <- paste0("unassigned ", last_known)
      } else {
        row[i] <- "unassigned root"
      }
    } else {
      last_known <- row[i]
    }
  }
  return(row)
}

gtdb2 <- as.data.frame(t(apply(gtdb, 1, fill_unassigned)), stringsAsFactors = FALSE)

# merge
metabolishmm_tax <- merge(metabolishmm, gtdb2, by.y="row.names", by.x="bin")
metabolishmm_merged <- merge(metabolishmm_tax, coverage_mean, by.x="bin", by.y="label")

# calculate weighted genera counts
metabolishmm_genus <- metabolishmm_merged %>% 
  mutate(weighted_rel_abund=count*mean_rel_abund) %>% # count weighted by mean relative abundance
  group_by(phylum, genus, gene_short, functional_group) %>% 
  summarise(genus_count = sum(count), # total count by genus
            genus_abund = sum(weighted_rel_abund)) %>% # total genus count weighted by mean relative abundance
  filter(genus_abund != 0) # remove zeros for better visualisation

metabolishmm_genus$genus <- factor(metabolishmm_genus$genus,
                                     levels = rev(unique(metabolishmm_genus$genus))) # reverse order for plotting

# subset to top genera
coverage_mean_genera <- coverage_rel %>% 
  t(.) %>% 
  merge(., gtdb2, by="row.names") %>% 
  group_by(genus) %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  column_to_rownames(., var="genus") %>% 
  mutate(mean_coverage = rowMeans(.)) %>%
  arrange(desc(mean_coverage)) %>%
  slice_head(n = 30)

metabolishmm_genus_sub <- metabolishmm_genus[metabolishmm_genus$genus %in% rownames(coverage_mean_genera), ]
metabolishmm_genus_sub$genus <- droplevels(metabolishmm_genus_sub$genus)

# plot heatmaps
plots <- lapply(ordered_groups, function(group) {
  group_data <- metabolishmm_genus_sub %>% filter(functional_group == group)
  
  ggplot(group_data, aes(x = gene_short, y = genus, fill = genus_abund)) +
    geom_tile() +
    scale_fill_gradientn(colors = color_gradients[[group]], name = "Weighted Count",
                         limits=c(0,0.04),
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
    scale_y_discrete(limits = levels(metabolishmm_genus_sub$genus)) 
})

# phylum colour plot
phylum_colour_map <- data.frame(
  Phylum = c("Pseudomonadota", "Myxococcota", "Acidobacteriota", "Eisenbacteria", "Gemmatimonadota", 
             "Nitrospirota", "Deinococcota", "Actinomycetota", "Chloroflexota", "Armatimonadota", 
             "Bacillota_A", "Bacillota_C", "Cyanobacteriota", "Bdellovibrionota", "Planctomycetota", 
             "Verrucomicrobiota", "Spirochaetota", "Chlamydiota", "Campylobacterota", "Bacillota", 
             "Patescibacteria", "Bacteroidota", "Unknown"),
  Colour = c("#1F77B4FF", "#FF7F0EFF", "#2CA02CFF", "#D3D3D3", "#D62728FF", "#87CEEB", "#4B0085", 
             "#52BD90", "#A3623A", "#FF69B4", "#FFF700", "#F7B6D2FF", "#9467BDFF", "#787276", 
             "#EFC050", "#98DF8AFF", "#AA336A", "#C5B0D5FF", "#808000", "#FF6F61", "#018571", 
             "#17BECFFF", "#C7C7C7FF"))

phylum_colour_map <- setNames(phylum_colour_map$Colour, phylum_colour_map$Phylum)

phyla <- ggplot(metabolishmm_genus_sub, aes(x = 1, y = genus, fill = phylum)) +
  geom_tile() +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="bottom",
        legend.text = element_text(size = 5),
        strip.background = element_blank(),
        strip.text = element_text(size = 6, colour = "black"),
        plot.title = element_text(size=6, colour="black", hjust=0.5),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 5, colour = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = phylum_colour_map, name = "Phylum")

phyla

main_plot <- wrap_plots(c(list(phyla), plots[c(1, 2, 5)], list(phyla), plots[c(3, 4)]), ncol = 4, widths = c(0.1, 1, 1, 1, 0.1, 1, 1)) +
  plot_layout(axes="collect") &
  theme(legend.position="none")

main_plot

ggsave("output/metabolishmm_genus.png", main_plot, width = 180, height = 180, units = "mm", dpi = 600)
ggsave("output/metabolishmm_genus.pdf", main_plot, width = 180, height = 180, units = "mm", dpi = 600)

main_plot <- wrap_plots(c(list(phyla), plots[c(1, 2, 5)], list(phyla), plots[c(3, 4)]), ncol = 4, widths = c(0.1, 1, 1, 1, 0.1, 1, 1)) +
  plot_layout(axes="collect") &
  theme(legend.position="right")

main_plot

ggsave("output/metabolishmm_genus_legend.png", main_plot, width = 180, height = 180, units = "mm", dpi = 600)
ggsave("output/metabolishmm_genus_legend.pdf", main_plot, width = 180, height = 180, units = "mm", dpi = 600)


# supp_plot <- wrap_plots(plots[c(3,4)], ncol = 2) +
#   plot_layout(axes="collect") &
#   theme(legend.position="bottom")
# 
# supp_plot
# 
# ggsave("output/metabolishmm_HO.png", supp_plot, width = 180, height = 120, units = "mm", dpi = 600)
# ggsave("output/metabolishmm_HO.pdf", supp_plot, width = 140, height = 100, units = "mm", dpi = 600)
