# microtrait
# figure 5 - trait heatmap


# load libraries
library(microtrait)
library(readr)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(reshape2)

set.seed(123)


# import data
microtrait <- read.csv("data/microtrait_results.csv", row.names=1, check.names=F, sep=",") # microtrait output is normalised by genome size

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
microtrait_tax <- merge(microtrait, gtdb, by.x="id", by.y="row.names")
microtrait_merged <- merge(microtrait_tax, coverage_mean, by.x="id", by.y="label")


# calculate weighted phylum counts
trait_list <- colnames(microtrait_merged[2:101]) # get trait list

traits <- microtrait_merged %>% 
  mutate(across(all_of(trait_list), ~ . * mean_rel_abund)) %>% # count weighted by mean relative abundance
  group_by(phylum) %>% 
  summarise(across(all_of(trait_list), sum),) # sum trait presence by phylum


# prepare table for plotting
traits2 <- traits %>% 
  column_to_rownames(., var="phylum") %>% 
  select(where(~ sum(.) != 0)) %>% # remove traits that are all zero
  rownames_to_column(., var = "phylum")

traits_melt <- traits2 %>% 
  melt(.) %>% 
  filter(variable %in% trait_list) %>% 
  mutate(value=na_if(value,0)) %>% # zeros to NA for better visualisation
  separate(variable, into = c("level1", "level2", "level3","level4","level5", "level6"), # separate traits into their levels
           sep=":",
           extra="merge",
           fill="right") %>% 
  mutate(phylum=factor(phylum, levels=sort(unique(phylum), decreasing=T))) %>%  # reverse phylum order for plotting
  arrange(phylum)

# change NA to last non-NA value within the level columns
replace_na <- function(row) {
  last_non_na <- NA  # store last non-NA value
  for (i in 1:length(row)) {  # iterate from left to right across rows (take last non-NA value to the left)
    if (!is.na(row[i])) {  # if the value is not NA, update last_non_na
      last_non_na <- row[i]
    }
    if (is.na(row[i])) {  # if the value is NA, replace it with the last_non_na
      row[i] <- last_non_na
    }
  }
  return(row)
}

level_columns <- grep("^level", names(traits_melt))  # level columns

for (i in 1:nrow(traits_melt)) {
  traits_melt[i, level_columns] <- replace_na(traits_melt[i, level_columns]) # apply function
}

# harmonise level labels:
# resource use can all be shown at level4 and acquisition can be shown at level5
# what level to show stress tolerance at is more complicated:
# scavenging calculated for general and oxidative stress (general is more broad, just keep general)
# protection calculated for general, pH stress and temperature (sums are the same, just keep general)
# redox calculated for oxidative stress, oxygen limitation and pH stress (sums are the same for oxidative stress and oxygen limitation, keep oxidative stress and pH stress but label this way @ level 5)
# membrane stabilisation calculated for pH stress and temperature stress (label this way at level5)
# all others good at level 5

traits_melt2 <- traits_melt %>% 
  filter(!(level2 == "Specific" & level4 == "scavenging of reactive oxygen species")) %>% 
  filter(!(level2 == "Specific" & level5 == "protection, repair, degradation of denatured/misfolded proteins")) %>% 
  filter(!(level4 == "redox sensing" & level3 == "oxygen limitation")) %>% 
  mutate(level5 = ifelse(level4 == "redox sensing", paste(level3, level5, sep= " "), level5)) %>% 
  mutate(level5 = ifelse(level6 == "membrane stabilization", paste(level3, level5, sep = " "), level5)) %>% 
  mutate(level2 = ifelse(level2 %in% c("General", "Specific"), level1, level2))


# plot
ordered_groups <- c("Substrate assimilation", "Substrate degradation", "Substrate uptake", "Chemotrophy", "Phototrophy", "Stress Tolerance")

level_4_plots <- lapply(ordered_groups, function(group) {
  group_data <- traits_melt2 %>% filter(level2 == group)
  
  ggplot(group_data, aes(x = level4, y = phylum, fill = value)) +
    geom_tile() +
    scale_fill_gradientn(colors = brewer.pal(9, "Blues"),
                         limits = c(0,0.5),
                         name = "Weighted count",
                         na.value="white") +
    labs(title = group) + 
    theme_bw() +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position="bottom",
          strip.background = element_blank(),
          strip.text = element_text(size = 6, colour = "black"),
          plot.title = element_text(size=6, colour="black", hjust=0.5),
          axis.title = element_blank(),
          axis.text.y = element_text(size = 5, colour = "black"),
          axis.text.x = element_text(size = 5, colour = "black", angle = 45, hjust = 1)) 
})

level_5_plots <- lapply(ordered_groups, function(group) {
  group_data <- traits_melt2 %>% filter(level2 == group)
  
  ggplot(group_data, aes(x = level5, y = phylum, fill = value)) +
    geom_tile() +
    scale_fill_gradientn(colors = brewer.pal(9, "Blues"),
                         limits = c(0,0.5),
                         name = "Weighted count",
                         na.value="white") +
    labs(title = group) + 
    theme_bw() +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position="bottom",
          strip.background = element_blank(),
          strip.text = element_text(size = 6, colour = "black"),
          plot.title = element_text(size=6, colour="black", hjust=0.5),
          axis.title = element_blank(),
          axis.text.y = element_text(size = 5, colour = "black"),
          axis.text.x = element_text(size = 5, colour = "black", angle = 45, hjust = 1)) 
})


patchwork <- (level_4_plots[[1]] +
                level_4_plots[[2]] +
                level_4_plots[[3]] +
                level_4_plots[[4]] +
                level_5_plots[[5]] +
                level_5_plots[[6]]) +
  plot_layout(ncol=6, widths=c(7, 11, 37, 6, 5, 17), guides="collect") &
  theme(legend.position = "bottom")

patchwork

patchwork[[2]] = patchwork[[2]] + theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank())

patchwork[[3]] = patchwork[[3]] + theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank())

patchwork[[4]] = patchwork[[4]] + theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank())

patchwork[[5]] = patchwork[[5]] + theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank())

patchwork[[6]] = patchwork[[6]] + theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank())

patchwork

ggsave("output/microtrait_phylum_combined2.png", patchwork, width = 18, height = 8, units = "in", dpi = 600)
ggsave("output/microtrait_phylum_combined2.pdf", patchwork, width = 200, height = 120, units = "mm", dpi = 600)
