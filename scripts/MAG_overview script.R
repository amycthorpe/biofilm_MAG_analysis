# MAG overview plots
# figure 1 - contamination vs completeness, GC content vs genome size, taxonomic novelty
# extended data figure 1 - genomic trait distributions


# load libraries
library(ggplot2)
library(dplyr)
library(ggpubr)
library(patchwork)
library(tidyverse)
library(reshape2)
library(pals)
library(ggsci)

set.seed(123)


# import data
checkm <- read.csv("data/checkm_gtdb.csv", row.names=1, check.names=F, sep=",")
checkm$Genome_Size_Mbp <- checkm$Genome_Size / 1000000 # convert units

coverage <- read.csv("data/finalbins_coverage.csv", row.names=1, check.names=F, sep=",")


# contamination and completeness thresholds
checkm <- mutate(checkm, Completeness_quality = case_when(
  Completeness >= 90 & Contamination <= 5 ~ "Near-complete MAGs",
  Completeness >= 70 & Completeness < 90 & Contamination <= 10 ~ "Medium-quality MAGs",
  TRUE ~ "Low-quality MAGs"
))


# set colours
quality_colors <- c("Near-complete MAGs" = "forestgreen", "Medium-quality MAGs" = "steelblue3", "Low-quality MAGs" = "violetred")

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


# contamination vs completeness plot
# scatter plot
scatter_plot <- ggplot(checkm, aes(x = Completeness, y = Contamination, color = Completeness_quality)) +
  geom_point(size=0.5) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black")) +
  labs(x = "Completeness (%)", y = "Contamination (%)") +
  guides(color = guide_legend(title = "")) +
  scale_color_manual(values = quality_colors) +
  scale_x_continuous(limits = c(50, 100), breaks = seq(0, 100, by = 10)) +
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 100, by = 5))

scatter_plot

# contamination histogram
checkm <- checkm %>% # set categories
  mutate(Contamination_category = case_when(
    Contamination >= 0 & Contamination < 2.5 ~ "0-2.5",
    Contamination >= 2.5 & Contamination < 5 ~ "2.5-5",
    Contamination >= 5 & Contamination < 7.5 ~ "5-7.5",
    Contamination >= 7.5 & Contamination < 10 ~ "7.5-10",
    Contamination >= 10 & Contamination < 12.5 ~ "10-12.5",
    Contamination >= 12.5 & Contamination < 15 ~ "12.5-15",
    Contamination >= 15 & Contamination < 17.5 ~ "15-17.5",
    Contamination >= 17.5 & Contamination < 20 ~ "17.5-20",
    Contamination >= 20 & Contamination <= 22.5 ~ "20-22.5",
    Contamination >= 22.5 & Contamination <= 25 ~ "22.5-25",
    Contamination >= 25 & Contamination <= 27.5 ~ "25-27.5",
    Contamination >= 27.5 & Contamination <= 30 ~ "27.5-30",
    TRUE ~ NA_character_
  ))

contamination_counts <- checkm %>% # sum number MAGs per category
  filter(!is.na(Contamination_category) & !is.na(Completeness_quality)) %>%
  group_by(Contamination_category, Completeness_quality) %>%
  summarise(count = n(), .groups = 'drop')

contam_bar_plot <- ggplot(contamination_counts, aes(x = Contamination_category, y = count, fill = Completeness_quality)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 6, colour = "black"),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 5, colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  guides(fill = guide_legend(title = "")) +
  labs(x = "Contamination level (%)", y = "Number of MAGs") +
  scale_fill_manual(values = quality_colors) +
  scale_y_continuous(limits = c(0, 600), breaks = seq(0, 600, by = 300)) +
  scale_x_discrete(limits = c("0-2.5", "2.5-5", "5-7.5", "7.5-10", "10-12.5", 
                              "12.5-15", "15-17.5", "17.5-20", "20-22.5", "22.5-25",
                              "25-27.5", "27.5-30")) +
  coord_flip()

contam_bar_plot

# completeness histogram
checkm <- checkm %>% # set categories
  mutate(Completeness_category = case_when(
    Completeness >= 55 & Completeness < 57.5 ~ "55-57.5",
    Completeness >= 57.5 & Completeness < 60 ~ "57.5-60",
    Completeness >= 60 & Completeness < 62.5 ~ "60-62.5",
    Completeness >= 62.5 & Completeness < 65 ~ "62.5-65",
    Completeness >= 65 & Completeness < 67.5 ~ "65-67.5",
    Completeness >= 67.5 & Completeness < 70 ~ "67.5-70",
    Completeness >= 70 & Completeness < 72.5 ~ "70-72.5",
    Completeness >= 72.5 & Completeness < 75 ~ "72.5-75",
    Completeness >= 75 & Completeness < 77.5 ~ "75-77.5",
    Completeness >= 77.5 & Completeness < 80 ~ "77.5-80",
    Completeness >= 80 & Completeness < 82.5 ~ "80-82.5",
    Completeness >= 82.5 & Completeness < 85 ~ "82.5-85",
    Completeness >= 85 & Completeness < 87.5 ~ "85-87.5",
    Completeness >= 87.5 & Completeness < 90 ~ "87.5-90",
    Completeness >= 90 & Completeness < 92.5 ~ "90-92.5",
    Completeness >= 92.5 & Completeness < 95 ~ "92.5-95",
    Completeness >= 95 & Completeness < 97.5 ~ "95-97.5",
    Completeness >= 97.5 & Completeness <= 100 ~ "97.5-100",
    TRUE ~ NA_character_
  ))

completeness_counts <- checkm %>% # sum number MAGs per category
  filter(!is.na(Completeness_category) & !is.na(Completeness_quality)) %>%
  group_by(Completeness_category, Completeness_quality) %>%
  summarise(count = n(), .groups = 'drop')

complete_bar_plot <- ggplot(completeness_counts, aes(x = Completeness_category, y = count, fill = Completeness_quality)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 6, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 5, colour = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(title = "")) +
  labs(x = "Completeness level (%)", y = "Number of MAGs") +
  scale_fill_manual(values = quality_colors) +
  scale_x_discrete(limits = c("55-57.5", "57.5-60", "60-62.5", "62.5-65", "65-67.5",
                              "67.5-70", "70-72.5", "72.5-75", "75-77.5", "77.5-80",
                              "80-82.5", "82.5-85", "85-87.5", "87.5-90", "90-92.5", 
                              "92.5-95", "95-97.5", "97.5-100"))

complete_bar_plot

# arrange plots
design <-
  "
AB           
CD           
"

combined_bar_plot <- list(
  complete_bar_plot,
  plot_spacer(),
  scatter_plot,
  contam_bar_plot,
  guide_area()
) |>
  wrap_plots() +
  plot_layout(axes = "collect", guides = "collect", widths = c(6, 1.5), heights = c(1.5, 6), design = design) &
  theme(legend.direction = "vertical", legend.box = "horizontal", legend.position = "none")

combined_bar_plot

ggsave("output/contam_compl_bar_plot.png", combined_bar_plot, width = 7, height = 7, units = "in", dpi = 600)
ggsave("output/contam_compl_bar_plot.pdf", combined_bar_plot, width = 70, height = 70, units = "mm", dpi = 600)


# checkm histograms 
# genome size
checkm_filter <- checkm %>% 
  filter(., Genome_Size < 11400000) # remove large MAG outlier due to very high proportion protist genes (for plotting only)

genome_size <- ggplot(checkm_filter, aes(x = Genome_Size_Mbp)) +
  geom_histogram(binwidth=0.5, colour="black", fill="lightblue") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black")) +
  scale_y_continuous(limits = c(0, 200), breaks = seq(0, 200, by = 50)) +
  labs(x = "Genome size (Mbp)", y = "Number of MAGs")

genome_size

ggsave("output/genome_size_histogram.png", genome_size, width = 90, height = 45, units = "mm", dpi = 600)
ggsave("output/genome_size_histogram.pdf", genome_size, width = 90, height = 45, units = "mm", dpi = 600)

# total coding seqs histogram
total_coding <- ggplot(checkm_filter, aes(x = Total_Coding_Sequences)) +
  geom_histogram(binwidth=500, colour="black", fill="lightblue") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black")) +
  scale_y_continuous(limits = c(0, 200), breaks = seq(0, 200, by = 50)) +
  labs(x = "Total coding sequences", y = "Number of MAGs")

total_coding

ggsave("output/total_coding_histogram.png", total_coding, width = 90, height = 45, units = "mm", dpi = 600)
ggsave("output/total_coding_histogram.pdf", total_coding, width = 90, height = 45, units = "mm", dpi = 600)

# total contigs histogram
contigs <- ggplot(checkm_filter, aes(x = Total_Contigs)) +
  geom_histogram(binwidth=100, colour="black", fill="lightblue") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black")) +
  scale_y_continuous(limits = c(0, 200), breaks = seq(0, 200, by = 50)) +
  labs(x = "Total contigs", y = "Number of MAGs")

contigs

ggsave("output/total_contigs_histogram.png", contigs, width = 90, height = 45, units = "mm", dpi = 600)
ggsave("output/total_contigs_histogram.pdf", contigs, width = 90, height = 45, units = "mm", dpi = 600)

# N50 histogram
N50 <- ggplot(checkm_filter, aes(x = Contig_N50)) +
  geom_histogram(binwidth=20000, colour="black", fill="lightblue") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black")) +
  scale_y_continuous(limits = c(0, 500), breaks = seq(0, 500, by = 100)) +
  labs(x = "Contig N50", y = "Number of MAGs")

N50

ggsave("output/N50_histogram.png", N50, width = 90, height = 45, units = "mm", dpi = 600)
ggsave("output/N50_histogram.pdf", N50, width = 90, height = 45, units = "mm", dpi = 600)

# coding density histogram
coding_density <- ggplot(checkm_filter, aes(x = Coding_Density)) +
  geom_histogram(binwidth=0.02, colour="black", fill="lightblue") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black")) +
  scale_y_continuous(limits = c(0, 400), breaks = seq(0, 400, by = 100)) +
  labs(x = "Coding density (%)", y = "Number of MAGs")

coding_density

ggsave("output/coding_density_histogram.png", coding_density, width = 90, height = 45, units = "mm", dpi = 600)
ggsave("output/coding_density_histogram.pdf", coding_density, width = 90, height = 45, units = "mm", dpi = 600)

# GC content histogram
GC <- ggplot(checkm_filter, aes(x = GC_Content)) +
  geom_histogram(binwidth=0.02, colour="black", fill="lightblue") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black")) +
  labs(x = "GC content (%)", y = "Number of MAGs") +
  scale_y_continuous(limits = c(0, 80), breaks = seq(0, 800, by = 20))

GC

ggsave("output/GC_histogram.png", GC, width = 90, height = 45, units = "mm", dpi = 600)
ggsave("output/GC_histogram.pdf", GC, width = 90, height = 45, units = "mm", dpi = 600)

# gene length histogram
gene_length <- ggplot(checkm_filter, aes(x = Average_Gene_Length)) +
  geom_histogram(binwidth=10, colour="black", fill="lightblue") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black")) +
  labs(x = "Average gene length (bp)", y = "Number of MAGs") +
  scale_y_continuous(limits = c(0, 200), breaks = seq(0, 200, by = 50))

gene_length

ggsave("output/gene_length_histogram.png", gene_length, width = 90, height = 45, units = "mm", dpi = 600)
ggsave("output/gene_length_histogram.pdf", gene_length, width = 90, height = 45, units = "mm", dpi = 600)

# contig length histogram
contig_length <- ggplot(checkm_filter, aes(x = Max_Contig_Length)) +
  geom_histogram(binwidth=40000, colour="black", fill="lightblue") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black")) +
  labs(x = "Maximum contig length (bp)", y = "Number of MAGs") +
  scale_y_continuous(limits = c(0, 500), breaks = seq(0, 500, by = 100))

contig_length

ggsave("output/contig_length_histogram.png", contig_length, width = 90, height = 45, units = "mm", dpi = 600)
ggsave("output/contig_length_histogram.pdf", contig_length, width = 90, height = 45, units = "mm", dpi = 600)


# GC content vs genome size vs abundance vs taxonomy
# convert to relative abundance
coverage_rel <- coverage %>%
  mutate(across(everything(), ~ . / sum(.))) %>%
  t(.) # samples as rows

# calculate mean across all samples
coverage_mean <- as.data.frame(colMeans(coverage_rel))
coverage_mean <- coverage_mean %>% rownames_to_column(., var = "label")
colnames(coverage_mean)[2] <- "mean_rel_abund"

# merge
gtdb_checkm_coverage <- merge(coverage_mean, checkm, by.y="row.names", by.x = "label")

# sum by phylum
gtdb_checkm_coverage_phylum <- gtdb_checkm_coverage %>% 
  select(phylum, where(is.numeric)) %>%
  group_by(phylum) %>% 
  summarise(across(everything(), sum))

# plot scatter
GC_size <- ggplot(gtdb_checkm_coverage, aes(x = Genome_Size_Mbp, y = GC_Content, colour = phylum, size = mean_rel_abund)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black")) +
  labs(x = "Genome size (Mbp)", y = "GC content (%)") +
  guides(colour = guide_legend(title = "Phylum"),
         size = guide_legend(title = "Mean relative abundance")) +
  scale_color_manual(values = phylum_colour_map) +
  scale_size_continuous(range = c(0.2, 2)) +
  scale_y_continuous(labels = function(x) x*100) +
  scale_x_continuous(limits=c(0,12), breaks=seq(0,12,by=2)) # remove limits to include biggest MAG outlier

GC_size

ggsave("output/GC_size.png", GC_size, width = 10, height = 6, units = "in", dpi = 600)
ggsave("output/GC_size.pdf", GC_size, width = 70, height = 50, units = "mm", dpi = 600)


# genome size and environmental metadata
coverage_melt <- melt(coverage_rel)

genome_size_coverage <- gtdb_checkm_coverage %>% 
  select(., c(label, Genome_Size_Mbp)) %>% 
  merge(., coverage_melt, by.x="label", by.y="Var2")

colnames(genome_size_coverage)[3] <- "metag_ID"
colnames(genome_size_coverage)[4] <- "rel_coverage"

weighted_gsize_coverage <- genome_size_coverage %>% 
  mutate(weighted_genome_size = rel_coverage * Genome_Size_Mbp) %>% 
  group_by(metag_ID) %>% 
  summarise(mean_weighted_size = sum(weighted_genome_size))

water_chem <- read.csv("data/water_chem.csv", check.names=F, row.names=1, sep=",") %>% 
  select(ends_with("_mean")) # only mean columns

new_names <- c(
  "Alkalinity",
  "Ammoniacal_N",
  "DOC",
  "Conductivity",
  "Nitrate_N",
  "Nitrite_N",
  "Orthophosphate",
  "DO",
  "SiO2",
  "Water_temp",
  "pH")

colnames(water_chem) <- new_names # tidy colnames

merged <- merge(weighted_gsize_coverage, water_chem, by.y="row.names", by.x="metag_ID")

plot <- ggplot(merged, aes(x = Nitrate_N, y = mean_weighted_size)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 10, colour = "black"),
        axis.text = element_text(size = 8, colour = "black")) +
  labs(x = "Nitrate-N", y = "Weighted genome size")

plot

model <- lm(Nitrate_N ~ mean_weighted_size, data = merged)
summary(model)

plot <- ggplot(merged, aes(x = Orthophosphate, y = mean_weighted_size)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 10, colour = "black"),
        axis.text = element_text(size = 8, colour = "black")) +
  labs(x = "Orthophosphate", y = "Weighted genome size")

plot

model <- lm(Orthophosphate ~ mean_weighted_size, data = merged)
summary(model)

plot <- ggplot(merged, aes(x = Nitrite_N, y = mean_weighted_size)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 10, colour = "black"),
        axis.text = element_text(size = 8, colour = "black")) +
  labs(x = "Nitrite-N", y = "Weighted genome size")

plot

model <- lm(Nitrite-N ~ mean_weighted_size, data = merged)
summary(model)


# taxonomic novelty
# subset to only near-complete (comment out to run on all MAGs)
#checkm <- subset(checkm, Completeness_quality == "Near-complete MAGs")

# subset to taxonomy
taxonomy <- checkm[, 14:20]

# label empty taxonomy
taxonomy_labelled <- taxonomy %>%
  mutate(across(c(kingdom, phylum, class, order, family, genus, species),
                ~ if_else(. == "", "unknown", .))) %>%
  mutate(across(c(order, family, genus, species),
                ~ if_else(. != "unknown", "known", .)))

# known vs unknown at each rank
taxonomy_transformed <- taxonomy_labelled %>%
  pivot_longer(cols = c(order, family, genus, species),
               names_to = "rank",
               values_to = "taxonomy")

# set order
rank_order <- c("order", "family", "genus", "species")
rank_order_rev <- rev(rank_order)

# plot proportion
counts <- taxonomy_transformed %>% # calculate known and unknown MAG counts for annotations
  group_by(rank, taxonomy) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(rank) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

known_v_unknown <- ggplot(taxonomy_transformed, aes(x = rank, fill = taxonomy)) +
  geom_bar(stat = "count", position = "fill") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title.x = element_text(size = 6, colour = "black"),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 5, colour = "black")) +
  labs(x = "Rank", y = "Proportion of MAGs (%)") +
  guides(fill = guide_legend(title = "")) +
  scale_x_discrete(limits = rank_order_rev,
                   labels = function(x) str_to_title(x)) +
  scale_fill_manual(values = c("lightblue", "grey"),
                    labels = function(x) str_to_title(x)) +
  scale_y_continuous(labels = function(x) paste0(x*100)) +
  coord_flip() +
  geom_text(data = counts, aes(x = rank, y = count / sum(count), 
                               label = paste0(count, "\n(", sprintf("%.1f", percentage), "%)")), 
            position = position_fill(vjust = 0.5), size = 1.5, color = "black")

known_v_unknown

ggsave("output/known_v_unknown.png", known_v_unknown, width = 6, height = 4, units = "in", dpi = 600)
ggsave("output/known_v_unknown.pdf", known_v_unknown, width = 70, height = 40, units = "mm", dpi = 600)
