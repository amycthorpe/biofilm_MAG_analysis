# microbeannotator results


# load libraries
library(httr)
library(dplyr)
library(purrr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(reshape2)
library(pals)
library(ggsci)

set.seed(123)


# import database list
bbsdb <- read.csv("data/biofilm_KOs/protein.csv", row.names=1, check.names=F, sep=",")


# function to match gene symbols in KEGG
get_symbol_ko_numbers <- function(gene_name) {
  # KEGG API URL for symbol search
  url <- paste0("https://rest.kegg.jp/find/ko/", gene_name)
  
  # API request
  response <- tryCatch({
    GET(url)
  }, error = function(e) {
    message(paste("Error fetching data for", gene_name, ":", e$message))
    return(NULL)
  })
  
  # check if request failed
  if (is.null(response) || status_code(response) != 200) {
    return("no KO assigned")
  }
  
  # parse the response
  content <- content(response, "text", encoding = "UTF-8")
  
  # if empty response
  if (nchar(content) == 0) {
    return("no KO assigned")
  }
  
  # process each line to find exact symbol matches
  lines <- strsplit(content, "\n")[[1]]
  symbol_matches <- character(0)
  
  for (line in lines) {
    # split into KO and description parts
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) >= 2) {
      # check if gene_name appears as a standalone symbol in the description, looks for:
      # 1. start of string or space, followed by gene_name
      # 2. followed by end of string, space, comma, or semicolon
      if (grepl(paste0("(^| )", gene_name, "($| |,|;)"), parts[2])) {
        ko <- regmatches(parts[1], regexpr("K\\d{5}", parts[1]))
        if (length(ko) > 0) {
          symbol_matches <- c(symbol_matches, ko)
        }
      }
    }
  }
  
  # process results
  if (length(symbol_matches) > 0) {
    return(paste(unique(symbol_matches), collapse = "; "))
  } else {
    return("no KO assigned")
  }
}

# add KO numbers to dataframe using symbol matching
bbsdb <- bbsdb %>%
  mutate(KO_numbers = map_chr(gene_name, function(gene) {
    message(paste("Processing gene:", gene))
    ko_result <- get_symbol_ko_numbers(gene)
    Sys.sleep(0.5)
    return(ko_result)
  }))

# get KO list
assigned <- subset(bbsdb, KO_numbers != "no KO assigned")
ko_list <- unlist(strsplit(assigned$KO_numbers, ";"))
ko_list <- trimws(ko_list)
unique_ko_list <- unique(ko_list)
KOs <- as.data.frame(ko_list)
colnames(KOs) <- "BBSdb"  

write.csv(KOs, "data/biofilm_KOs/BBSdb_KO_list.csv", row.names = FALSE)

# save results
write.csv(bbsdb, "data/biofilm_KOs/symbol_ko_mappings.csv", row.names = FALSE)


# import KO lists
bbsdb <- read.csv("data/biofilm_KOs/symbol_ko_mappings.csv")
KO_lists <- read.csv("data/biofilm_KOs/KO_list.csv") # from KEGG pathways

KO_lists_long <- KO_lists %>% 
  pivot_longer(cols = everything(), 
               names_to = "pathway", 
               values_to = "KO") %>%
  filter(KO != "") %>% 
  distinct(KO, pathway, .keep_all = TRUE) %>% 
  mutate(pathway = case_when(
    grepl("biofilm_formation_", pathway) ~ "biofilm_formation",
    TRUE ~ pathway
  ))

final_list <- KO_lists_long %>%
  group_by(KO) %>%
  summarize(pathway = str_c(pathway, collapse = "; ")) %>%
  ungroup()

clean_ko <- function(x) {
  x %>%
    gsub("\u00A0", "", ., fixed = TRUE) %>% # remove hidden spaces
    trimws() %>%
    toupper()
}

final_list2 <- final_list %>%
  mutate(KO = clean_ko(KO)) 

write_excel_csv(final_list2, "data/biofilm_KOs/KO_list_pathways.csv")


# match with microbeannotator annotation results
files <- list.files(path = "data/biofilm_KOs/annotation_results", pattern = "\\.faa\\.annot$", full.names = TRUE)

annotation_list <- lapply(files, read.delim, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
names(annotation_list) <- tools::file_path_sans_ext(basename(files))

annotation_list <- annotation_list[sapply(annotation_list, nrow) > 0]

annotation_df <- bind_rows(annotation_list, .id = "bin")

merged_data <- annotation_df %>%
  left_join(final_list2, by = c("ko_number" = "KO")) %>% 
  mutate(biofilm_KOs = ifelse(!is.na(pathway), "biofilm_associated",  # if true, mark as biofilm associated
      NA_character_)) %>% # if false, keep as NA
  select(ko_number, bin, ko_product, pathway, biofilm_KOs)

merged_data$bin <- gsub("\\.faa$", "", merged_data$bin) # remove .faa from bin names

biofilm_ko_counts <- merged_data %>%
  group_by(bin) %>%
  summarise(total_ko_count = n_distinct(ko_number),
            biofilm_ko_count = n_distinct(ko_number[biofilm_KOs == "biofilm_associated"]),
            BBSdb_ko_count = n_distinct(ko_number[grepl("BBSdb", pathway)]),
            quorum_ko_count = n_distinct(ko_number[grepl("quorum_sensing", pathway)]),
            chemotaxis_ko_count = n_distinct(ko_number[grepl("chemotaxis", pathway)]),
            flagellar_ko_count = n_distinct(ko_number[grepl("flagellar_assembly", pathway)]),
            TCS_ko_count = n_distinct(ko_number[grepl("two_component_system", pathway)]),
            ABC_ko_count = n_distinct(ko_number[grepl("abc_transporters", pathway)]),
            formation_ko_count = n_distinct(ko_number[grepl("biofilm_formation", pathway)]),
            EPS_ko_count = n_distinct(ko_number[grepl("polysac_eps", pathway)]))

biofilm_ko_percents <- biofilm_ko_counts %>%
  mutate(biofilm_ko_percent = 100 * biofilm_ko_count / total_ko_count,
         BBSdb_ko_percent = 100 * BBSdb_ko_count / total_ko_count,
         quorum_ko_percent = 100 * quorum_ko_count / total_ko_count,
         chemotaxis_ko_percent = 100 * chemotaxis_ko_count / total_ko_count,
         flagellar_ko_percent = 100 * flagellar_ko_count / total_ko_count,
         TCS_ko_percent = 100 * TCS_ko_count / total_ko_count,
         ABC_ko_percent = 100 * ABC_ko_count / total_ko_count,
         formation_ko_percent = 100 * formation_ko_count / total_ko_count,
         EPS_ko_percent = 100 * EPS_ko_count / total_ko_count) %>% 
  select(bin, ends_with("_percent"))

write.csv(biofilm_ko_counts, "data/biofilm_KOs/biofilm_KO_counts.csv")
write.csv(biofilm_ko_percents, "data/biofilm_KOs/biofilm_KO_percents.csv")


# match with taxonomy and genome data
# import data
checkm <- read.csv("data/checkm_gtdb.csv", row.names=1, check.names=F, sep=",")
checkm$Genome_Size_Mbp <- checkm$Genome_Size / 1000000 # convert units

coverage <- read.csv("data/finalbins_coverage.csv", row.names=1, check.names=F, sep=",")

# convert to relative abundance
coverage_rel <- coverage %>%
  mutate(across(everything(), ~ . / sum(.))) %>%
  t(.) # samples as rows

# calculate mean across all samples
coverage_mean <- as.data.frame(colMeans(coverage_rel))
coverage_mean <- coverage_mean %>% rownames_to_column(., var = "bin")
colnames(coverage_mean)[2] <- "mean_rel_abund"

# filter low quality MAGs
checkm <- mutate(checkm, Completeness_quality = case_when(
  Completeness >= 90 & Contamination <= 5 ~ "Near-complete MAGs",
  Completeness >= 70 & Completeness < 90 & Contamination <= 10 ~ "Medium-quality MAGs",
  TRUE ~ "Low-quality MAGs"))

checkm <- subset(checkm, Completeness_quality != "Low-quality MAGs")
biofilm_ko_counts <- biofilm_ko_counts[biofilm_ko_counts$bin %in% rownames(checkm), ]

# merge
checkm <- rownames_to_column(checkm, var="bin")
merged_df <- merge(checkm, biofilm_ko_counts, by = "bin")
merged_df <- merge(merged_df, coverage_mean, by = "bin")

# normalise KO counts to total coding sequences
merged_df2 <- merged_df

ko_columns <- grep("_ko_count$", names(merged_df2), value = TRUE)
merged_df2[ko_columns] <- lapply(ko_columns, function(col) merged_df2[[col]] / merged_df2$Total_Coding_Sequences)


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


# plot scatter - biofilm associated
biofilm_plot <- ggplot(merged_df2, aes(x = Genome_Size_Mbp, y = biofilm_ko_count, colour = phylum, size = mean_rel_abund)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black")) +
  labs(x = "Genome size (Mbp)", y = "Biofilm associated KOs") +
  guides(colour = guide_legend(title = "Phylum"),
         size = guide_legend(title = "Mean relative abundance")) +
  scale_color_manual(values = phylum_colour_map) +
  scale_size_continuous(range = c(0.2, 2)) +
  scale_x_continuous(limits=c(0,12), breaks=seq(0,12,by=2)) # remove limits to include biggest MAG outlier

biofilm_plot

# ggsave("output/biofilm_KOs.png", biofilm_plot, width = 70, height = 50, units = "mm", dpi = 600)
# ggsave("output/biofilm_KOs.pdf", biofilm_plot, width = 70, height = 50, units = "mm", dpi = 600)


# plot scatter - quorum sensing
QS_plot <- ggplot(merged_df2, aes(x = Genome_Size_Mbp, y = quorum_ko_count, colour = phylum, size = mean_rel_abund)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black")) +
  labs(x = "Genome size (Mbp)", y = "KOs associated with\n quorum sensing") +
  guides(colour = guide_legend(title = "Phylum"),
         size = guide_legend(title = "Mean relative abundance")) +
  scale_color_manual(values = phylum_colour_map) +
  scale_size_continuous(range = c(0.2, 2)) +
  scale_x_continuous(limits=c(0,12), breaks=seq(0,12,by=2)) # remove limits to include biggest MAG outlier

QS_plot

# ggsave("output/quorum_KOs.png", QS_plot, width = 70, height = 50, units = "mm", dpi = 600)
# ggsave("output/quorum_KOs.pdf", QS_plot, width = 70, height = 50, units = "mm", dpi = 600)


# plot scatter - chemotaxis
chemotaxis_plot <- ggplot(merged_df2, aes(x = Genome_Size_Mbp, y = chemotaxis_ko_count, colour = phylum, size = mean_rel_abund)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black")) +
  labs(x = "Genome size (Mbp)", y = "KOs associated with\n chemotaxis") +
  guides(colour = guide_legend(title = "Phylum"),
         size = guide_legend(title = "Mean relative abundance")) +
  scale_color_manual(values = phylum_colour_map) +
  scale_size_continuous(range = c(0.2, 2)) +
  scale_x_continuous(limits=c(0,12), breaks=seq(0,12,by=2)) # remove limits to include biggest MAG outlier

chemotaxis_plot

# ggsave("output/chemotaxis_KOs.png", chemotaxis_plot, width = 70, height = 50, units = "mm", dpi = 600)
# ggsave("output/chemotaxis_KOs.pdf", chemotaxis_plot, width = 70, height = 50, units = "mm", dpi = 600)


# plot scatter - flagellar assembly
flagellar_plot <- ggplot(merged_df2, aes(x = Genome_Size_Mbp, y = flagellar_ko_count, colour = phylum, size = mean_rel_abund)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black")) +
  labs(x = "Genome size (Mbp)", y = "KOs associated with\n flagellar assembly") +
  guides(colour = guide_legend(title = "Phylum"),
         size = guide_legend(title = "Mean relative abundance")) +
  scale_color_manual(values = phylum_colour_map) +
  scale_size_continuous(range = c(0.2, 2)) +
  scale_x_continuous(limits=c(0,12), breaks=seq(0,12,by=2)) # remove limits to include biggest MAG outlier

flagellar_plot

# ggsave("output/flagellar_KOs.png", flagellar_plot, width = 70, height = 50, units = "mm", dpi = 600)
# ggsave("output/flagellar_KOs.pdf", flagellar_plot, width = 70, height = 50, units = "mm", dpi = 600)


# plot scatter - two component system
TCS_plot <- ggplot(merged_df2, aes(x = Genome_Size_Mbp, y = TCS_ko_count, colour = phylum, size = mean_rel_abund)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black")) +
  labs(x = "Genome size (Mbp)", y = "KOs associated with\n the two-component system") +
  guides(colour = guide_legend(title = "Phylum"),
         size = guide_legend(title = "Mean relative abundance")) +
  scale_color_manual(values = phylum_colour_map) +
  scale_size_continuous(range = c(0.2, 2)) +
  scale_x_continuous(limits=c(0,12), breaks=seq(0,12,by=2)) # remove limits to include biggest MAG outlier

TCS_plot

# ggsave("output/TCS_KOs.png", TCS_plot, width = 70, height = 50, units = "mm", dpi = 600)
# ggsave("output/TCS_KOs.pdf", TCS_plot, width = 70, height = 50, units = "mm", dpi = 600)


# plot scatter - ABC transporters
ABC_plot <- ggplot(merged_df2, aes(x = Genome_Size_Mbp, y = ABC_ko_count, colour = phylum, size = mean_rel_abund)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black")) +
  labs(x = "Genome size (Mbp)", y = "KOs associated with\n ABC transporters") +
  guides(colour = guide_legend(title = "Phylum"),
         size = guide_legend(title = "Mean relative abundance")) +
  scale_color_manual(values = phylum_colour_map) +
  scale_size_continuous(range = c(0.2, 2)) +
  scale_x_continuous(limits=c(0,12), breaks=seq(0,12,by=2)) # remove limits to include biggest MAG outlier

ABC_plot

# ggsave("output/ABC_KOs.png", ABC_plot, width = 70, height = 50, units = "mm", dpi = 600)
# ggsave("output/ABC_KOs.pdf", ABC_plot, width = 70, height = 50, units = "mm", dpi = 600)


# plot scatter - biofilm formation
formation_plot <- ggplot(merged_df2, aes(x = Genome_Size_Mbp, y = formation_ko_count, colour = phylum, size = mean_rel_abund)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black")) +
  labs(x = "Genome size (Mbp)", y = "KOs associated with\n biofilm formation") +
  guides(colour = guide_legend(title = "Phylum"),
         size = guide_legend(title = "Mean relative abundance")) +
  scale_color_manual(values = phylum_colour_map) +
  scale_size_continuous(range = c(0.2, 2)) +
  scale_x_continuous(limits=c(0,12), breaks=seq(0,12,by=2)) # remove limits to include biggest MAG outlier

formation_plot

# ggsave("output/formation_KOs.png", formation_plot, width = 70, height = 50, units = "mm", dpi = 600)
# ggsave("output/formation_KOs.pdf", formation_plot, width = 70, height = 50, units = "mm", dpi = 600)


# plot scatter - EPS
EPS_plot <- ggplot(merged_df2, aes(x = Genome_Size_Mbp, y = EPS_ko_count, colour = phylum, size = mean_rel_abund)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black")) +
  labs(x = "Genome size (Mbp)", y = "KOs associated with\n EPS biosynthesis") +
  guides(colour = guide_legend(title = "Phylum"),
         size = guide_legend(title = "Mean relative abundance")) +
  scale_color_manual(values = phylum_colour_map) +
  scale_size_continuous(range = c(0.2, 2)) +
  scale_x_continuous(limits=c(0,12), breaks=seq(0,12,by=2)) # remove limits to include biggest MAG outlier

EPS_plot

# ggsave("output/EPS_KOs.png", EPS_plot, width = 70, height = 50, units = "mm", dpi = 600)
# ggsave("output/EPS_KOs.pdf", EPS_plot, width = 70, height = 50, units = "mm", dpi = 600)


# coding seqs vs genome size
total_coding_size <- ggplot(merged_df2, aes(x = Genome_Size_Mbp, y = Total_Coding_Sequences, colour = phylum, size = mean_rel_abund)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black")) +
  labs(x = "Genome size (Mbp)", y = "Total coding sequences") +
  guides(colour = guide_legend(title = "Phylum"),
         size = guide_legend(title = "Mean relative abundance")) +
  scale_color_manual(values = phylum_colour_map) +
  scale_size_continuous(range = c(0.2, 2)) +
  scale_y_continuous(limits=c(0,10000), breaks=seq(0,10000,by=3000)) +
  scale_x_continuous(limits=c(0,12), breaks=seq(0,12,by=2)) # remove limits to include biggest MAG outlier

total_coding_size


# combined plot
combined <- wrap_plots(biofilm_plot, 
                       QS_plot, 
                       chemotaxis_plot,
                       flagellar_plot,
                       TCS_plot,
                       ABC_plot,
                       formation_plot,
                       EPS_plot,
                       total_coding_size,
                       ncol = 3) +
  plot_layout(guides = "collect", axes = "collect") &
  theme(legend.position="none")

combined

ggsave("output/biofilm_KOs_combined_totalcodingseqnorm.pdf", combined, width = 180, height = 120, units = "mm", dpi = 600)


# boxplots - facet by phylum
long_df <- merged_df2 %>%
  pivot_longer(cols = ends_with("_ko_count"),
               names_to = "KO_type",
               values_to = "KO_count") %>% 
  filter(!KO_type %in% c("BBSdb_ko_count", "total_ko_count")) %>%
  mutate(KO_type = recode(KO_type,
                          "ABC_ko_count" = "KOs associated with ABC transporters",
                          "BBSdb_ko_count" = "BBSdb KOs",
                          "biofilm_ko_count" = "Biofilm associated KOs",
                          "chemotaxis_ko_count" = "KOs associated with chemotaxis",
                          "EPS_ko_count" = "KOs associated with EPS biosynthesis",
                          "flagellar_ko_count" = "KOs associated with flagellar assembly",
                          "formation_ko_count" = "KOs associated with biofilm formation",
                          "quorum_ko_count" = "KOs associated with quorum sensing",
                          "TCS_ko_count" = "KOs associated with the two-component system"))

phylum_order <- long_df %>%
  group_by(phylum) %>%
  summarise(total_abund = sum(mean_rel_abund)) %>%
  arrange(desc(total_abund)) %>%
  filter(total_abund >= 0.01) %>%
  pull(phylum)

long_df <- long_df %>%
  filter(phylum %in% phylum_order) %>%
  mutate(phylum = factor(phylum, levels = phylum_order))

boxplots <- ggplot(long_df, aes(x = phylum, y = KO_count, fill = phylum)) +
  geom_boxplot(outlier.size = 0.5, size = 0.2) +
  theme_bw() +
  facet_wrap(~ KO_type, scales = "free_y", ncol=4) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 6, colour = "black"),
        axis.title.y = element_text(size = 6, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 5, colour = "black"),
        axis.text.x = element_text(size = 5, colour = "black", angle = 45, hjust = 1)) +
  labs(y = "KO count / total coding sequences") +
  scale_fill_manual(values = phylum_colour_map)

boxplots

ggsave("output/biofilm_KOs_boxplot_totalcodingseqnorm.pdf", boxplots, width = 180, height = 120, units = "mm", dpi = 600)


# normalise KO counts to coding density
merged_df3 <- merged_df

ko_columns <- grep("_ko_count$", names(merged_df3), value = TRUE)
merged_df3[ko_columns] <- lapply(ko_columns, function(col) merged_df3[[col]] / merged_df3$Coding_Density)


# plot scatter - biofilm associated
model <- lm(biofilm_ko_count ~ Genome_Size_Mbp, data = merged_df3)
summary(model)

eq_text <- paste0(
  "R² = ", round(summary(model)$r.squared, 3), "; ",
  ifelse(summary(model)$coefficients[2, 4] < 0.001, "p < 0.001",
         ifelse(summary(model)$coefficients[2, 4] < 0.01, "p < 0.01",
                ifelse(summary(model)$coefficients[2, 4] < 0.05, "p < 0.05",
                       paste0("p = ", signif(summary(model)$coefficients[2, 4], 2))))))

biofilm_plot <- ggplot(merged_df3, aes(x = Genome_Size_Mbp, y = biofilm_ko_count, colour = phylum, size = mean_rel_abund)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  geom_smooth(aes(group = 1), colour = "black", method = "lm", size = 0.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black")) +
  annotate("text", x = 0.5, y = max(merged_df3$biofilm_ko_count, na.rm = TRUE) * 0.95, 
           label = eq_text, hjust = 0, size = 2) + 
  labs(x = "Genome size (Mbp)", y = "Biofilm associated KOs") +
  guides(colour = guide_legend(title = "Phylum"),
         size = guide_legend(title = "Mean relative abundance")) +
  scale_color_manual(values = phylum_colour_map) +
  scale_size_continuous(range = c(0.2, 2)) +
  scale_x_continuous(limits=c(0,12), breaks=seq(0,12,by=2)) # remove limits to include biggest MAG outlier

biofilm_plot


# plot scatter - quorum sensing
model <- lm(quorum_ko_count ~ Genome_Size_Mbp, data = merged_df3)
summary(model)

eq_text <- paste0(
  "R² = ", round(summary(model)$r.squared, 3), "; ",
  ifelse(summary(model)$coefficients[2, 4] < 0.001, "p < 0.001",
         ifelse(summary(model)$coefficients[2, 4] < 0.01, "p < 0.01",
                ifelse(summary(model)$coefficients[2, 4] < 0.05, "p < 0.05",
                       paste0("p = ", signif(summary(model)$coefficients[2, 4], 2))))))

QS_plot <- ggplot(merged_df3, aes(x = Genome_Size_Mbp, y = quorum_ko_count, colour = phylum, size = mean_rel_abund)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  geom_smooth(aes(group = 1), colour = "black", method = "lm", size = 0.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black")) +
  annotate("text", x = 0.5, y = max(merged_df3$quorum_ko_count, na.rm = TRUE) * 0.95, 
           label = eq_text, hjust = 0, size = 2) + 
  labs(x = "Genome size (Mbp)", y = "KOs associated with\n quorum sensing") +
  guides(colour = guide_legend(title = "Phylum"),
         size = guide_legend(title = "Mean relative abundance")) +
  scale_color_manual(values = phylum_colour_map) +
  scale_size_continuous(range = c(0.2, 2)) +
  scale_x_continuous(limits=c(0,12), breaks=seq(0,12,by=2)) # remove limits to include biggest MAG outlier

QS_plot


# plot scatter - chemotaxis
model <- lm(chemotaxis_ko_count ~ Genome_Size_Mbp, data = merged_df3)
summary(model)

eq_text <- paste0(
  "R² = ", round(summary(model)$r.squared, 3), "; ",
  ifelse(summary(model)$coefficients[2, 4] < 0.001, "p < 0.001",
         ifelse(summary(model)$coefficients[2, 4] < 0.01, "p < 0.01",
                ifelse(summary(model)$coefficients[2, 4] < 0.05, "p < 0.05",
                       paste0("p = ", signif(summary(model)$coefficients[2, 4], 2))))))

chemotaxis_plot <- ggplot(merged_df3, aes(x = Genome_Size_Mbp, y = chemotaxis_ko_count, colour = phylum, size = mean_rel_abund)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  geom_smooth(aes(group = 1), colour = "black", method = "lm", size = 0.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black")) +
  annotate("text", x = 0.5, y = max(merged_df3$chemotaxis_ko_count, na.rm = TRUE) * 0.95, 
           label = eq_text, hjust = 0, size = 2) + 
  labs(x = "Genome size (Mbp)", y = "KOs associated with\n chemotaxis") +
  guides(colour = guide_legend(title = "Phylum"),
         size = guide_legend(title = "Mean relative abundance")) +
  scale_color_manual(values = phylum_colour_map) +
  scale_size_continuous(range = c(0.2, 2)) +
  scale_x_continuous(limits=c(0,12), breaks=seq(0,12,by=2)) # remove limits to include biggest MAG outlier

chemotaxis_plot


# plot scatter - flagellar assembly
model <- lm(flagellar_ko_count ~ Genome_Size_Mbp, data = merged_df3)
summary(model)

eq_text <- paste0(
  "R² = ", round(summary(model)$r.squared, 3), "; ",
  ifelse(summary(model)$coefficients[2, 4] < 0.001, "p < 0.001",
         ifelse(summary(model)$coefficients[2, 4] < 0.01, "p < 0.01",
                ifelse(summary(model)$coefficients[2, 4] < 0.05, "p < 0.05",
                       paste0("p = ", signif(summary(model)$coefficients[2, 4], 2))))))

flagellar_plot <- ggplot(merged_df3, aes(x = Genome_Size_Mbp, y = flagellar_ko_count, colour = phylum, size = mean_rel_abund)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  geom_smooth(aes(group = 1), colour = "black", method = "lm", size = 0.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black")) +
  annotate("text", x = 0.5, y = max(merged_df3$flagellar_ko_count, na.rm = TRUE) * 0.95, 
           label = eq_text, hjust = 0, size = 2) + 
  labs(x = "Genome size (Mbp)", y = "KOs associated with\n flagellar assembly") +
  guides(colour = guide_legend(title = "Phylum"),
         size = guide_legend(title = "Mean relative abundance")) +
  scale_color_manual(values = phylum_colour_map) +
  scale_size_continuous(range = c(0.2, 2)) +
  scale_x_continuous(limits=c(0,12), breaks=seq(0,12,by=2)) # remove limits to include biggest MAG outlier

flagellar_plot


# plot scatter - two component system
model <- lm(TCS_ko_count ~ Genome_Size_Mbp, data = merged_df3)
summary(model)

eq_text <- paste0(
  "R² = ", round(summary(model)$r.squared, 3), "; ",
  ifelse(summary(model)$coefficients[2, 4] < 0.001, "p < 0.001",
         ifelse(summary(model)$coefficients[2, 4] < 0.01, "p < 0.01",
                ifelse(summary(model)$coefficients[2, 4] < 0.05, "p < 0.05",
                       paste0("p = ", signif(summary(model)$coefficients[2, 4], 2))))))

TCS_plot <- ggplot(merged_df3, aes(x = Genome_Size_Mbp, y = TCS_ko_count, colour = phylum, size = mean_rel_abund)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  geom_smooth(aes(group = 1), colour = "black", method = "lm", size = 0.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black")) +
  annotate("text", x = 0.5, y = max(merged_df3$TCS_ko_count, na.rm = TRUE) * 0.95, 
           label = eq_text, hjust = 0, size = 2) + 
  labs(x = "Genome size (Mbp)", y = "KOs associated with\n the two-component system") +
  guides(colour = guide_legend(title = "Phylum"),
         size = guide_legend(title = "Mean relative abundance")) +
  scale_color_manual(values = phylum_colour_map) +
  scale_size_continuous(range = c(0.2, 2)) +
  scale_x_continuous(limits=c(0,12), breaks=seq(0,12,by=2)) # remove limits to include biggest MAG outlier

TCS_plot


# plot scatter - ABC transporters
model <- lm(ABC_ko_count ~ Genome_Size_Mbp, data = merged_df3)
summary(model)

eq_text <- paste0(
  "R² = ", round(summary(model)$r.squared, 3), "; ",
  ifelse(summary(model)$coefficients[2, 4] < 0.001, "p < 0.001",
         ifelse(summary(model)$coefficients[2, 4] < 0.01, "p < 0.01",
                ifelse(summary(model)$coefficients[2, 4] < 0.05, "p < 0.05",
                       paste0("p = ", signif(summary(model)$coefficients[2, 4], 2))))))

ABC_plot <- ggplot(merged_df3, aes(x = Genome_Size_Mbp, y = ABC_ko_count, colour = phylum, size = mean_rel_abund)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  geom_smooth(aes(group = 1), colour = "black", method = "lm", size = 0.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black")) +
  annotate("text", x = 0.5, y = max(merged_df3$ABC_ko_count, na.rm = TRUE) * 0.95, 
           label = eq_text, hjust = 0, size = 2) + 
  labs(x = "Genome size (Mbp)", y = "KOs associated with\n ABC transporters") +
  guides(colour = guide_legend(title = "Phylum"),
         size = guide_legend(title = "Mean relative abundance")) +
  scale_color_manual(values = phylum_colour_map) +
  scale_size_continuous(range = c(0.2, 2)) +
  scale_x_continuous(limits=c(0,12), breaks=seq(0,12,by=2)) # remove limits to include biggest MAG outlier

ABC_plot


# plot scatter - biofilm formation
model <- lm(formation_ko_count ~ Genome_Size_Mbp, data = merged_df3)
summary(model)

eq_text <- paste0(
  "R² = ", round(summary(model)$r.squared, 3), "; ",
  ifelse(summary(model)$coefficients[2, 4] < 0.001, "p < 0.001",
         ifelse(summary(model)$coefficients[2, 4] < 0.01, "p < 0.01",
                ifelse(summary(model)$coefficients[2, 4] < 0.05, "p < 0.05",
                       paste0("p = ", signif(summary(model)$coefficients[2, 4], 2))))))

formation_plot <- ggplot(merged_df3, aes(x = Genome_Size_Mbp, y = formation_ko_count, colour = phylum, size = mean_rel_abund)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  geom_smooth(aes(group = 1), colour = "black", method = "lm", size = 0.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black")) +
  annotate("text", x = 0.5, y = max(merged_df3$formation_ko_count, na.rm = TRUE) * 0.95, 
           label = eq_text, hjust = 0, size = 2) + 
  labs(x = "Genome size (Mbp)", y = "KOs associated with\n biofilm formation") +
  guides(colour = guide_legend(title = "Phylum"),
         size = guide_legend(title = "Mean relative abundance")) +
  scale_color_manual(values = phylum_colour_map) +
  scale_size_continuous(range = c(0.2, 2)) +
  scale_x_continuous(limits=c(0,12), breaks=seq(0,12,by=2)) # remove limits to include biggest MAG outlier

formation_plot


# plot scatter - EPS
model <- lm(EPS_ko_count ~ Genome_Size_Mbp, data = merged_df3)
summary(model)

eq_text <- paste0(
  "R² = ", round(summary(model)$r.squared, 3), "; ",
  ifelse(summary(model)$coefficients[2, 4] < 0.001, "p < 0.001",
         ifelse(summary(model)$coefficients[2, 4] < 0.01, "p < 0.01",
                ifelse(summary(model)$coefficients[2, 4] < 0.05, "p < 0.05",
                       paste0("p = ", signif(summary(model)$coefficients[2, 4], 2))))))

EPS_plot <- ggplot(merged_df3, aes(x = Genome_Size_Mbp, y = EPS_ko_count, colour = phylum, size = mean_rel_abund)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  geom_smooth(aes(group = 1), colour = "black", method = "lm", size = 0.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black")) +
  annotate("text", x = 0.5, y = max(merged_df3$EPS_ko_count, na.rm = TRUE) * 0.95, 
           label = eq_text, hjust = 0, size = 2) + 
  labs(x = "Genome size (Mbp)", y = "KOs associated with\n EPS biosynthesis") +
  guides(colour = guide_legend(title = "Phylum"),
         size = guide_legend(title = "Mean relative abundance")) +
  scale_color_manual(values = phylum_colour_map) +
  scale_size_continuous(range = c(0.2, 2)) +
  scale_x_continuous(limits=c(0,12), breaks=seq(0,12,by=2)) # remove limits to include biggest MAG outlier

EPS_plot


# coding density vs genome size
model <- lm(Coding_Density ~ Genome_Size_Mbp, data = merged_df3)
summary(model)

eq_text <- paste0(
  "R² = ", round(summary(model)$r.squared, 3), "; ",
  ifelse(summary(model)$coefficients[2, 4] < 0.001, "p < 0.001",
         ifelse(summary(model)$coefficients[2, 4] < 0.01, "p < 0.01",
                ifelse(summary(model)$coefficients[2, 4] < 0.05, "p < 0.05",
                       paste0("p = ", signif(summary(model)$coefficients[2, 4], 2))))))

coding_size <- ggplot(merged_df3, aes(x = Genome_Size_Mbp, y = Coding_Density, colour = phylum, size = mean_rel_abund)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  geom_smooth(aes(group = 1), colour = "black", method = "lm", size = 0.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black")) +
  annotate("text", x = 0.5, y = max(merged_df3$Coding_Density, na.rm = TRUE), 
           label = eq_text, hjust = 0, size = 2) + 
  labs(x = "Genome size (Mbp)", y = "Coding density (%)") +
  guides(colour = guide_legend(title = "Phylum"),
         size = guide_legend(title = "Mean relative abundance")) +
  scale_color_manual(values = phylum_colour_map) +
  scale_size_continuous(range = c(0.2, 2)) +
  scale_y_continuous(labels = function(x) x*100, limits=c(0.6,1), breaks=seq(0.6,1,by=0.1)) +
  scale_x_continuous(limits=c(0,12), breaks=seq(0,12,by=2)) # remove limits to include biggest MAG outlier

coding_size


# combined plot
combined <- wrap_plots(QS_plot, 
                       chemotaxis_plot,
                       flagellar_plot,
                       TCS_plot,
                       EPS_plot,
                       formation_plot,
                       ABC_plot,
                       coding_size,
                       ncol = 3) +
  plot_layout(guides = "collect", axes = "collect") &
  theme(legend.position="none")

combined

ggsave("output/biofilm_KOs_combined_codingdensitynorm.pdf", combined, width = 180, height = 120, units = "mm", dpi = 600)


# boxplots - facet by phylum
long_df <- merged_df3 %>%
  pivot_longer(cols = ends_with("_ko_count"),
               names_to = "KO_type",
               values_to = "KO_count") %>% 
  filter(!KO_type %in% c("BBSdb_ko_count", "total_ko_count")) %>%
  mutate(KO_type = recode(KO_type,
                          "ABC_ko_count" = "KOs associated with ABC transporters",
                          "BBSdb_ko_count" = "BBSdb KOs",
                          "biofilm_ko_count" = "Biofilm associated KOs",
                          "chemotaxis_ko_count" = "KOs associated with chemotaxis",
                          "EPS_ko_count" = "KOs associated with EPS biosynthesis",
                          "flagellar_ko_count" = "KOs associated with flagellar assembly",
                          "formation_ko_count" = "KOs associated with biofilm formation",
                          "quorum_ko_count" = "KOs associated with quorum sensing",
                          "TCS_ko_count" = "KOs associated with the two-component system"))

phylum_order <- long_df %>%
  group_by(phylum) %>%
  summarise(total_abund = sum(mean_rel_abund)) %>%
  arrange(desc(total_abund)) %>%
  filter(total_abund >= 0.01) %>%
  pull(phylum)

long_df <- long_df %>%
  filter(phylum %in% phylum_order) %>%
  mutate(phylum = factor(phylum, levels = phylum_order))

boxplots <- ggplot(long_df, aes(x = phylum, y = KO_count, fill = phylum)) +
  geom_boxplot(outlier.size = 0.5, size = 0.2) +
  theme_bw() +
  facet_wrap(~ KO_type, scales = "free_y", ncol=4) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 6, colour = "black"),
        axis.title.y = element_text(size = 6, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 5, colour = "black"),
        axis.text.x = element_text(size = 5, colour = "black", angle = 45, hjust = 1)) +
  labs(y = "KO count / coding density") +
  scale_fill_manual(values = phylum_colour_map)

boxplots

ggsave("output/biofilm_KOs_boxplot_codingdensitynorm.pdf", boxplots, width = 180, height = 120, units = "mm", dpi = 600)
