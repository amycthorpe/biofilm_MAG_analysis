# levins index and occupancy
# figure 3 - generalists vs specialists and occupancy vs relative abundance


# load libraries
library(tidyverse)
library(ggplot2)
library(ggsci)
library(patchwork)
library(microeco)
library(MicroNiche)
library(tibble)
library(ggridges)

set.seed(123)


# import data
bins <- read.csv("data/finalbins_coverage.csv", check.names=F, sep=",") %>% 
  dplyr::rename("Taxon"="bin")

checkm <- read.csv("data/checkm_gtdb.csv", row.names=1, check.names=F, sep=",")

checkm <- mutate(checkm, Completeness_quality = case_when(
  Completeness >= 90 & Contamination <= 5 ~ "Near-complete MAGs",
  Completeness >= 70 & Completeness < 90 & Contamination <= 10 ~ "Medium-quality MAGs",
  TRUE ~ "Low-quality MAGs"
))

checkm <- subset(checkm, Completeness_quality != "Low-quality MAGs") # filter low quality MAGs
gtdb <- checkm[, 14:20]

coverage <- read.csv("data/finalbins_coverage.csv", row.names=1, check.names=F, sep=",")
coverage <- coverage[rownames(coverage) %in% rownames(checkm), ] # filter low quality MAGs
bins <- bins[bins$Taxon %in% rownames(checkm), ] # filter low quality MAGs


# calculate levins index
sampleinfo <- colnames(bins) %>% as.data.frame() %>% 
  rename("samples"=".") %>% 
  filter(samples != "Taxon")

levins_out <- levins.Bn(bins, length(sampleinfo$samples), sampleinfo$samples, q = 1.65) # run levins

levins_median <- median(levins_out$Bn) # get median Bn

levins_out <- levins_out %>% 
  mutate(category = ifelse(Bn > levins_median, "generalist", "specialist")) # categorise based on median Bn

write.csv(levins_out, "output/levins_median.csv")


# percent generalist vs specialist by taxonomy
levins_tax <- levins_out %>% 
  merge(., gtdb, by = "row.names")

levins_percent <- levins_tax %>% 
  group_by(phylum, category) %>% 
  summarise(count=n()) %>% # sum generalists vs specialists by phylum
  ungroup() %>% 
  group_by(phylum) %>% 
  mutate(percentage=(count/sum(count))*100) %>% # get percentage
  select(phylum, category, count, percentage)

levins_percent$phylum <- factor(levins_percent$phylum,
                                levels = rev(unique(levins_percent$phylum))) # reverse order for plotting

# plot levins proportions
levins_percent2 <- levins_percent %>%
  mutate(category = case_when(
    category == "generalist" ~ "b_generalist", # rename to control order in plot
    category == "specialist" ~ "a_specialist"
  ))

levins_plot <- ggplot(levins_percent2, aes(x = percentage, y = phylum, fill = category)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_text(size = 5, colour = "black")) +
  scale_fill_manual(values=c("b_generalist" = "seagreen3", "a_specialist" = "steelblue3")) +
  scale_x_continuous(labels = scales::percent_format(scale = 100)) 

levins_plot

ggsave("output/levins_barplot.png", levins_plot, width = 6, height = 4, units = "in", dpi = 600)
ggsave("output/levins_barplot.pdf", levins_plot, width = 90, height = 87, units = "mm", dpi = 600)


# occupancy
# convert to relative abundance
coverage_rel <- coverage %>%
  mutate(across(everything(), ~ . / sum(.))) %>% 
  t(.) # samples as rows

# calculate mean across all samples
coverage_mean <- as.data.frame(colMeans(coverage_rel))
coverage_mean <- coverage_mean %>% rownames_to_column(., var = "label")
colnames(coverage_mean)[2] <- "mean_rel_abund"

# merge
gtdb_coverage <- merge(coverage_mean, gtdb, by.y="row.names", by.x = "label")

# calculate presence absence
presence_absence <- (coverage > 0) * 1 # convert coverage to binary presence/absence

# sum number of samples each MAG had abundance > 0 in
occupancy <- as.data.frame(rowSums(presence_absence))

occupancy <- occupancy %>% rownames_to_column(., var = "label")
colnames(occupancy)[2] <- "occupancy"

# calculate number of samples
total_samples <- ncol(coverage)

# convert to percentage of total samples
occupancy$percentage <- (occupancy$occupancy / total_samples) * 100

# merge
occupancy_percentage_coverage <- merge(occupancy, coverage_mean, by="label")
occupancy_gtdb <- merge(occupancy_percentage_coverage, gtdb, by.x="label", by.y="row.names")


# set colours
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


# plot occupancy vs abundance log 
scientific_10_labels <- function(x) {
  parse(text = paste0("10^", floor(log10(abs(x))))) # convert labels to scientific notation
}

abundance_occupancy_log <- ggplot(occupancy_gtdb, aes(x = occupancy, y = mean_rel_abund, colour = phylum)) +
  geom_point() +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        axis.title = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black")) +
  labs(x = "Occupancy", y = "Mean relative abundance") +
  guides(colour = guide_legend(ncol=6)) +
  scale_color_manual(values = phylum_colour_map) +
  scale_x_continuous(limits = c(0, 450), breaks = seq(0, 450, by = 50)) +
  scale_y_continuous(limits = c(0, 0.02), breaks = seq(0, 0.02, by = 0.005), labels = scientific_10_labels) +
  scale_y_log10()

abundance_occupancy_log

ggsave("output/abundance_occupancy_log.png", abundance_occupancy_log, width = 9, height = 7, units = "in", dpi = 600)
ggsave("output/abundance_occupancy_log.pdf", abundance_occupancy_log, width = 90, height = 90, units = "mm", dpi = 600)
