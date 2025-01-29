# metabolic
# extended data figure 3 - metabolic pathways alluvial plot


# load libraries
library(dplyr)
library(tidyverse)
library(ggalluvial)
library(ggthemes)

set.seed(123)


# import data
metabolic <- read.csv("data/metabolic_results.csv", row.names=1, check.names=F, sep=",")

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
metabolic_tax <- merge(metabolic, gtdb, by.x="bin", by.y="row.names")
metabolic_merged <- merge(metabolic_tax, coverage_mean, by.x="bin", by.y="label")


# calculate weighted phylum counts
metabolic_phylum <- metabolic_merged %>% 
  mutate(weighted_rel_abund=presence*mean_rel_abund) %>% # count weighted by mean relative abundance
  group_by(phylum, reaction) %>% 
  summarise(phylum_count = sum(presence), # total count by phylum
            phylum_abund = sum(weighted_rel_abund)) # total phylum count weighted by mean relative abundance

energy.flow <- metabolic_phylum

energy.flow$category <- ifelse(grepl("C-S", energy.flow$reaction), "Carbon", # add categories
                               ifelse(grepl("N-S", energy.flow$reaction), "Nitrogen",
                                      ifelse(grepl("S-S", energy.flow$reaction), "Sulfur",
                                             ifelse(grepl("O-S", energy.flow$reaction), "Others",""))))
energy.flow$category <- as.factor(energy.flow$category)

energy.flow$reaction <- sub(".*:", "", energy.flow$reaction) # simplify reaction names


# filter for better visualisation
energy.flow2 <- energy.flow %>%
  group_by(phylum) %>%
  mutate(total_count = sum(phylum_abund)) %>%
  ungroup() %>% 
  filter(total_count >= 0.1) %>% # remove phyla with total count < 0.1
  select(-total_count)

energy.flow3 <- energy.flow2 %>%
  group_by(reaction) %>%
  mutate(total_count = sum(phylum_abund)) %>%
  ungroup() %>% 
  filter(total_count >= 0.1) %>% # remove reactions with total count < 0.1
  select(-total_count)


# set colours
category_colour_map <- c(
  "Carbon" = "#6BAED6",
  "Nitrogen" = "#74C476",
  "Sulfur" = "#FD8D3C",
  "Others" = "grey80")

reaction_colour_map <- c(
  "Organic carbon oxidation" = "#08306B",
  "Carbon fixation" = "#08519C",
  "Ethanol oxidation" = "#2171B5",
  "Acetate oxidation" = "#4292C6",
  "Hydrogen generation" = "#6BAED6",
  "Fermentation" = "#9ECAE1",
  "Methanogenesis" = "#C6DBEF",
  "Methanotrophy" = "#DEEBF7",
  "Hydrogen oxidation" = "#F7FBFF",
  "Nitrogen fixation" = "#00441B",
  "Ammonia oxidation" = "#006D2C",
  "Nitrite oxidation" = "#238B45",
  "Nitrate reduction" = "#41AB5D",
  "Nitrite reduction" = "#74C476",
  "Nitric oxide reduction" = "#A1D99B",
  "Nitrous oxide reduction" = "#C7E9C0",
  "Nitrite ammonification" = "#E5F5E0",
  "Anammox" = "#F7FCF5",
  "Iron reduction" = "#252525",
  "Iron oxidation" = "#636363",
  "Arsenate reduction" = "#969696",
  "Arsenite oxidation" = "#CCCCCC",
  "Selenate reduction" = "#F7F7F7",
  "Sulfide oxidation" = "#7F2704",
  "Sulfur reduction" = "#A63603",
  "Sulfur oxidation" = "#D94801",
  "Sulfite oxidation" = "#F16913",
  "Sulfate reduction" = "#FD8D3C",
  "Sulfite reduction" = "#FDAE6B",
  "Thiosulfate oxidation" = "#FDD0A2",
  "Thiosulfate disproportionation 1" = "#FEE6CE",
  "Thiosulfate disproportionation 2" = "#FFF5EB")

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


# plot
plot <- ggplot(as.data.frame(energy.flow3),
               aes(y = phylum_abund, axis1= phylum, axis2 = reaction, axis3 = category)) +
  theme_bw () +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 6, colour = "black"),
        axis.text.x = element_text(size = 6, colour = "black"),
        axis.text.y = element_text(size = 5, colour = "black"),
        axis.ticks.x = element_blank(),
        legend.position = "none") + 
  geom_alluvium(aes(fill = phylum),
                width = 1/12, curve_type = "quintic") +
  geom_stratum(aes(fill = after_stat(stratum)),
               width = 1/2, color = "black") +
  geom_text(stat = "stratum", infer.label = TRUE, size = 1.55) +
  scale_x_discrete(limits = c("Phylum", "Reaction", "Category"), expand=c(0,0)) +
  #  scale_y_continuous(limits = c(0, 6000), breaks = seq(0, 5000, by = 1000), expand = c(0,0)) +
  scale_fill_manual(values = c(phylum_colour_map, category_colour_map, reaction_colour_map),
                    na.value = "grey80", limits = \(x) x[!is.na(x)]) +
  labs(y = "Weighted count")

plot

ggsave("output/alluvial_plot_weighted.png", plot, width = 90, height = 100, units = "mm", dpi = 600)
ggsave("output/alluvial_plot_weighted.pdf", plot, width = 90, height = 100, units = "mm", dpi = 600)
