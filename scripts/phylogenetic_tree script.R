# MAG phylogenetic tree from fasta file
# figure 1 - phylogenetic tree


# load libraries
library(Biostrings)
library(R.utils)
library(phangorn)
library(tidytree)
library(treeio)
library(ggtree)
library(tidyr)
library(stringr)
library(dplyr)
library(seqinr)
library(ape)
library(msa)
library(ggtreeExtra)
library(ggsci)
library(ggnewscale)
library(reshape2)
library(pals)
library(scales)
library(tidyverse)

set.seed(123)


# import data
checkm <- read.csv("data/checkm_gtdb.csv", row.names=1, check.names=F, sep=",")

checkm <- mutate(checkm, Completeness_quality = case_when(
  Completeness >= 90 & Contamination <= 5 ~ "Near-complete MAGs",
  Completeness >= 70 & Completeness < 90 & Contamination <= 10 ~ "Medium-quality MAGs",
  TRUE ~ "Low-quality MAGs"
))

coverage <- read.csv("data/finalbins_coverage.csv", row.names=1, check.names=F, sep=",")

checkm <- subset(checkm, Completeness_quality != "Low-quality MAGs") # filter low quality MAGs
coverage <- coverage[rownames(coverage) %in% rownames(checkm_quality), ]

gtdb <- checkm[, 14:20] %>%
  rownames_to_column(., var="label")

# convert to relative abundance
coverage_rel <- coverage %>%
  mutate(across(everything(), ~ . / sum(.))) %>% 
  t(.) # samples as rows

# calculate mean across all samples
coverage_mean <- as.data.frame(colMeans(coverage_rel))
coverage_mean <- coverage_mean %>% rownames_to_column(., var = "label")
colnames(coverage_mean)[2] <- "mean_rel_abund"


# import fasta file
fasta_file <- "data/gtdbtk.bac120.user_msa.fasta.gz"
decompressed_file <- gunzip(fasta_file, remove = FALSE, temporary = TRUE)
fasta_data <- readAAStringSet(decompressed_file)
head(fasta_data)

fasta_data <- fasta_data[names(fasta_data) %in% rownames(checkm)]


# prepare tree data
# convert the AAStringSet to AAbin format (ape's format for AA sequences)
fasta_data_bin <- as.AAbin(fasta_data)

# calculate pairwise distances using the JTT (Jones-Taylor-Thornton) model
distance_matrix <- dist.ml(fasta_data_bin, model = "JTT")

# build the tree using the Neighbor-Joining method
tree_nj <- nj(distance_matrix)

# maximum-likelihood tree
aa_phyDat <- phyDat(fasta_data_bin, type = "AA") # create a phyDat object for amino acids

# estimate an ML tree using the JTT model
fit <- pml(tree_nj, aa_phyDat)
fit_opt <- optim.pml(fit, model = "JTT")  # optimize the ML tree

y <- fit_opt$tree
y_data <- y %>% as.treedata %>% as_tibble
head(y_data)


# join with taxonomy
y_data_tax <- left_join(y_data, gtdb, by="label") %>% 
  mutate(phylum = if_else(is.na(phylum), "Unknown", phylum))
y_tax_tree <- as.treedata(y_data_tax)

# join with coverage (mean relative abundance)
y_data_tax_cov <- left_join(y_data_tax, coverage_mean, by = "label")
y_tax_cov_tree <- as.treedata(y_data_tax_cov)


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


# plot
plot <- ggtree(y_tax_cov_tree, branch.length="none", layout="fan", open.angle = 8, aes(colour=phylum), size=0.25) +
  geom_tippoint(aes(colour=phylum), size=0.1) +
  scale_colour_manual(values=phylum_colour_map) +
  geom_fruit(geom=geom_bar,
             mapping=aes(y=label, x=mean_rel_abund, size=2),
             fill="darkblue",
             offset=0.05, 
             size=0.5,
             stat="identity",
             axis.params=list(axis="x",
                              text.angle=-90,
                              text.size=2.8,
                              hjust=-0.3,
                              vjust=0.3),
             grid.params=list()) +
  guides(colour=guide_legend(override.aes = list(size=3), ncol=3, title="Phylum"),
         fill=guide_legend(override.aes = list(size=3), ncol=3)) +
  theme(legend.title = element_text(size=14),
        legend.position="none",
        legend.text = element_text(size=12))

plot

ggsave("output/tree_plot.pdf", plot, width = 180, height = 170, units = "mm")
ggsave("output/tree_plot.pdf", plot, width = 150, height = 150, units = "mm")
