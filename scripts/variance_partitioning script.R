# variance partitioning and correlations
# figure 6 - variance partitioning
# extended data figure 5 - correlation matrix


# load libraries
library(ggplot2)
library(dplyr)
library(vegan)
library(tidyverse)
library(reshape2)
library(variancePartition)
library(pals)
library(ggsci)
library(corrplot)
library(RColorBrewer)
library(patchwork)
library(Hmisc)

set.seed(123)


# import data
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


# convert to relative abundance
coverage_rel <- coverage %>%
  mutate(across(everything(), ~ . / sum(.)))


# metadata
water_chem <- read.csv("data/supplementary/water_chem.csv", check.names=F, row.names=1, sep=",") %>% 
  select(ends_with("_mean")) # only mean columns

new_names <- c(
  "Alkalinity",
  "Ammonia_N_unionised",
  "Ammoniacal_N",
  "DOC",
  "Chloride",
  "Conductivity",
  "Nitrate_N",
  "Nitrite_N",
  "TN_oxidised",
  "TN",
  "Orthophosphate",
  "DO",
  "DO_sat",
  "TP",
  "SiO2",
  "Water_temp",
  "pH") 

colnames(water_chem) <- new_names # tidy colnames

land_cover <- read.csv("data/catchment_land_cover.csv", check.names=F, row.names=1, sep=",")
geology <- read.csv("data/catchment_geology.csv", check.names=F, row.names=1, sep=",")
chars <- read.csv("data/catchment_chars.csv", check.names=F, row.names=1, sep=",")

location <- read.csv("data/latitude.csv", row.names=1, check.names=F, sep=",")

# merge metadata
meta_merged1 <- merge(land_cover, geology, by="row.names") %>% 
  column_to_rownames(., var="Row.names")

colnames(meta_merged1) <- gsub("[ -]", "_", colnames(meta_merged1)) # replace spaces in colnames
colnames(meta_merged1) <- gsub(",", "", colnames(meta_merged1)) # remove commas

meta_merged2 <- merge(meta_merged1, chars, by="row.names") %>% 
  column_to_rownames(., var="Row.names")
meta_merged3 <- merge(meta_merged2, water_chem, by="row.names") %>% 
  column_to_rownames(., var="Row.names")
meta_merged3 <- merge(meta_merged3, location, by="row.names") %>% 
  column_to_rownames(., var="Row.names")


# correlation matrix
# include phylum coverage (mean relative summed by phylum)
phylum_coverage <- coverage_rel %>% 
  merge(., gtdb, by="row.names") %>% 
  group_by(phylum) %>% # group by phylum
  summarise(across(where(is.numeric), sum)) %>% # sum coverage by phylum
  column_to_rownames(., var="phylum") %>% 
  t(.)

# include metadata
meta_coverage <- meta_merged3 %>% 
  merge(., phylum_coverage, by="row.names") %>% 
  column_to_rownames(., var="Row.names")

# compute correlations
correlation_matrix <- cor(meta_coverage, use = "pairwise.complete.obs", method="pearson")
corrplot(correlation_matrix, method = "circle", tl.col = "black", tl.srt = 45)

# get r values
rcors <- rcorr(as.matrix(meta_coverage), type="pearson")
rvalues <- as.matrix(rcors[["r"]])
rvalues <- melt(rvalues) %>% 
  filter(value != 1)

# get p values
pvalues <- as.matrix(rcors[["P"]])
pvalues <- melt(pvalues) %>% 
  filter(across(everything(), ~ !is.na(.)))

cors <- merge(rvalues, pvalues, by=c("Var1", "Var2"))
colnames(cors)[3] <- "r_value"
colnames(cors)[4] <- "p_value"
#cors <- cors[, -4]

# add significance levels
cors$significance <- ifelse(cors$p_value <= 0.001, "***",
                            ifelse(cors$p_value <= 0.01, "**",
                                   ifelse(cors$p_value <= 0.05, "*",
                                          "")))

write.csv(cors, "output/correlations.csv")

# loop through p values to remove non-significant correlations
for (i in 1:nrow(pvalues)) {
  var1 <- pvalues$Var1[i]
  var2 <- pvalues$Var2[i]
  pval <- pvalues$value[i]
  
  # if p-value is >= 0.05, set corresponding correlation values to NA for plotting
  if (pval >= 0.05) {
    if (var1 %in% rownames(correlation_matrix) && var2 %in% colnames(correlation_matrix)) {
      correlation_matrix[var1, var2] <- NA
      correlation_matrix[var2, var1] <- NA # symmetric matrix
    }
  }
}


# plot correlation matrix
# set order of variables
water <- colnames(water_chem)
land_geo <- colnames(meta_merged1)
characteristics <- colnames(chars)
latitude <- colnames(location)
phyla <- colnames(phylum_coverage)

ordered_vars <- c(water, land_geo, characteristics, latitude, phyla)
ordered_vars <- intersect(ordered_vars, colnames(correlation_matrix))
correlation_matrix2 <- correlation_matrix[ordered_vars, ordered_vars]

# "Saltwater" , # only 3 non-zero values
# "Supra_littoral_rock" , # only 3 non-zero values
# "Littoral_rock" , # only 7 non-zero values
# "Littoral_sediment" , # only 2 non-zero values
# "Suburban" , # fails if included with urban due to co-correlation r=0.70 ***
# "Salt_pct" , # only 2 non-zero values

remove <- c("Saltwater", "Supra_littoral_rock", "Littoral_rock", "Littoral_sediment", "Salt_pct") # fewer than 7 non-zero values

correlation_matrix2 <- correlation_matrix2[!(rownames(correlation_matrix2) %in% remove), 
                                         !(colnames(correlation_matrix2) %in% remove)]
output_path <- "output/cor_matrix5.pdf"
pdf(file = output_path,
    width = 220 / 25.4,   # Convert mm to inches
    height = 220 / 25.4) 

corrplot(correlation_matrix2, 
         method = "circle", 
         tl.col = "black", 
         tl.srt = 90, 
         tl.cex = 0.4,
         na.label = " ",       # Use an empty string for NA (NA's are non-sig)
         na.label.col = "white",
         type = 'upper', diag = FALSE)

dev.off()


# variance partitioning
# subset to variables for analysis - remove any strong co-correlating variables (r>0.80 & ***) or those with few non-zero measurements otherwise model fails
variables <- c(
  "Alkalinity",
  # "Ammonia_N_unionised", # co-correlates with chloride (r=0.89 ***) and ammoniacal-N (r=0.86 ***)
  "Ammoniacal_N",
  "DOC",
  # "Chloride", # co-correlates with ammoniacal-N (r=0.86 ***) and conductivity (r=0.84 ***)
  "Conductivity",
  "Nitrate_N",
  "Nitrite_N",
  # "TN_oxidised", # co-correlates with nitrate-N (0.9999 ***) and TN (r=0.98 ***)
  # "TN", # co-correlates with nitrate-N (0.98 ***), also only 214 non-NA measurements
  "Orthophosphate",
  "DO",
  # "DO_sat", # co-correlates with DO (r=0.84 ***)
  # "TP", # co-correlates with orthophosphate (r=0.98 ***)
  "SiO2",
  "Water_temp",
  "pH",
  "Broadleaved_woodland" , 
  "Coniferous_woodland" , 
  "Arable_and_horticulture" , 
  "Improved_grassland" , 
  "Neutral_grassland" , 
  "Calcareous_grassland" , 
  "Acid_grassland" , 
  "Fen_marsh_and_swamp" , 
  "Heather" , 
  "Heather_grassland" , 
  "Bog" , 
  "Inland_rock" ,
  # "Saltwater" , # only 3 non-zero values
  "Freshwater" , 
  # "Supra_littoral_rock" , # only 3 non-zero values
  "Supra_littoral_sediment" , 
  # "Littoral_rock" , # only 7 non-zero values
  # "Littoral_sediment" , # only 2 non-zero values
  "Saltmarsh" , 
  "Urban" ,
  # "Suburban" , # fails if included with urban due to co-correlation r=0.70 ***
  "Siliceous_pct" , 
  "Calcareous_pct" , 
  # "Salt_pct" , # only 2 non-zero values
  "Peat_pct" ,
  "Chalk_pct" ,
  "altitude_m" , 
  "DSM_Shade_Average" , 
  "catchment_area_sqkm" , 
  "river_width" ,
  "river_depth" , 
  "strahler_number" , 
  "sum_WWTP_load_kg_sqkm" , 
  "lat")

meta_merged3_sub <- meta_merged3[colnames(meta_merged3) %in% variables]

# deal with NAs
non_na_counts <- as.data.frame(colSums(!is.na(meta_merged3_sub))) # check how many complete observations per variable
meta_merged3_sub <- na.omit(meta_merged3_sub) # keep only complete observations - 401 samples complete
coverage_rel_subset <- coverage_rel[colnames(coverage_rel) %in% rownames(meta_merged3_sub)] # subset to same samples

# set up formula
names <- colnames(meta_merged3_sub)
formatted_names <- paste(names, " +", sep = "")
cat(formatted_names, sep = "\n") # copy and paste output to below

form <- ~ 
  Broadleaved_woodland +
  Coniferous_woodland +
  Arable_and_horticulture +
  Improved_grassland +
  Neutral_grassland +
  Calcareous_grassland +
  Acid_grassland +
  Fen_marsh_and_swamp +
  Heather +
  Heather_grassland +
  Bog +
  Inland_rock +
  Freshwater +
  Supra_littoral_sediment +
  Saltmarsh +
  Urban +
  Siliceous_pct +
  Calcareous_pct +
  Peat_pct +
  Chalk_pct +
  altitude_m +
  DSM_Shade_Average +
  catchment_area_sqkm +
  river_width +
  river_depth +
  strahler_number +
  sum_WWTP_load_kg_sqkm +
  Alkalinity +
  Ammoniacal_N +
  DOC +
  Conductivity +
  Nitrate_N +
  Nitrite_N +
  Orthophosphate +
  DO +
  SiO2 +
  Water_temp +
  pH +
  lat

# align data
meta_merged3_sub <- meta_merged3_sub[order(rownames(meta_merged3_sub)), ] # sort samples alphabetically
coverage_rel_subset <- coverage_rel_subset[, order(colnames(coverage_rel_subset))] # sort samples alphabetically
all(rownames(meta_merged3_sub) == colnames(coverage_rel_subset))  # confirm they are in the same order

# run variance partitioning
varPart <- fitExtractVarPartModel(coverage_rel_subset, form, meta_merged3_sub)
write.csv(varPart, "output/varPart.csv")

# variance partitioning MAG stats
varPart <- subset(varPart, select = -Residuals)
MAG_sums <- as.data.frame(rowSums(varPart)) * 100 # % variance explained per MAG across all metadata
MAG_mean <- mean(MAG_sums$`rowSums(varPart)`) # 70.63% mean variance explained across all MAGs and all metadata
MAG_SD <- sd(MAG_sums$`rowSums(varPart)`) # 19.5%

# add taxonomy
varPart_tax <- as.data.frame(varPart) %>% 
  merge(., gtdb, by="row.names") %>% 
  group_by(phylum) %>% # group by phylum
  summarise(across(where(is.numeric), mean)) # calculate mean variance by phylum

varPart2 <- melt(varPart_tax) %>% 
  mutate(value = value * 100) %>% # convert to percentages
  filter(variable != "Residuals")

# set variable categories
varPart2$category <- ifelse(varPart2$variable %in% colnames(water_chem), "Water chemistry",
                            ifelse(varPart2$variable %in% colnames(geology), "Geology",
                                   ifelse(varPart2$variable %in% colnames(meta_merged1), "Land cover",
                                          "Watershed characteristics"))) 

varPart2$phylum <- factor(varPart2$phylum,
                          levels = rev(unique(varPart2$phylum))) # reverse order for plotting

varPart3 <- varPart2 %>% 
  group_by(phylum, category) %>% 
  summarise(value=sum(value)) # sum by category and phylum

varPart3 <- varPart3 %>%
  mutate(category = factor(category, levels = c("Watershed characteristics", "Geology", "Land cover", "Water chemistry"))) # set order for plotting


# plot all categores
plot <- ggplot(varPart3, aes(x = phylum, y = value, fill = category)) +
  geom_bar(stat = "identity", position="stack", color=NA, size=0) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size=5, colour="black"),
        legend.key.width = unit(2, "mm"),
        legend.key.height = unit(2, "mm"),
        axis.title.x = element_text(size = 6, colour = "black"),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 5, colour = "black")) +
  labs(x = "Phylum", y = "Variance explained (%)") +
  guides(fill = guide_legend(title = "")) +
  scale_y_continuous(limits=c(0,100), expand=c(0.02,0.02)) +
  scale_fill_manual(values=c("#FD8D3C", "#74C476", "#8682BC", "#6BAED6")) + 
  coord_flip()

plot

ggsave("output/variance_partitioning_categories.png", plot, width = 10, height = 6, units = "in", dpi = 600)
ggsave("output/variance_partitioning_categories.pdf", plot, width = 180, height = 80, units = "mm", dpi = 600)


# stats
sums <- varPart2 %>%
  group_by(phylum) %>%
  summarise(total_value = sum(value)) # total variance explained for each phylum

sums <- varPart2 %>%
  group_by(category, phylum) %>%
  summarise(total_value = sum(value)) # total variance explained by each category for each phylum

sums <- varPart2 %>%
  group_by(category, variable, phylum) %>%
  summarise(total_value = sum(value)) # total variance explained by each variable for each phylum


# plot by category - catchment geology
colourCount = length(unique(varPart2[varPart2$category == "Geology", "variable"]))
getPalette = colorRampPalette(brewer.pal(4, "Spectral"))

green <- brewer.pal(4, "Greens")

geology <- ggplot(varPart2 %>% filter(category=="Geology"), aes(x = phylum, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size=5, colour="black"),
        legend.key.width = unit(2, "mm"),
        legend.key.height = unit(2, "mm"),
        axis.title.x = element_text(size = 6, colour = "black"),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 5, colour = "black")) +
  labs(x = "Phylum", y = "Variance explained (%)") +
  guides(fill = guide_legend(title = "")) +
  scale_y_continuous(limits=c(0,65), breaks = seq(0, 65, by = 20), expand=c(0.02,0.02)) +
  scale_fill_manual(values = green) +
  coord_flip()

geology


# plot by category - catchment land cover
colourCount = length(unique(varPart2[varPart2$category == "Land cover", "variable"]))
getPalette = colorRampPalette(brewer.pal(16, "Spectral"))

purple <- colorRampPalette(brewer.pal(9, "Purples"))(16)

land <- ggplot(varPart2 %>% filter(category=="Land cover"), aes(x = phylum, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size=5, colour="black"),
        legend.key.width = unit(2, "mm"),
        legend.key.height = unit(2, "mm"),
        axis.title.x = element_text(size = 6, colour = "black"),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 5, colour = "black")) +
  labs(x = "Phylum", y = "Variance explained (%)") +
  guides(fill = guide_legend(title = "")) +
  scale_y_continuous(limits=c(0,30), breaks = seq(0, 30, by = 10), expand=c(0.02,0.02)) +
  scale_fill_manual(values = purple) +
  coord_flip()

land


# plot by category - catchment chars
colourCount = length(unique(varPart2[varPart2$category == "Watershed characteristics", "variable"]))
getPalette = colorRampPalette(brewer.pal(8, "Spectral"))

orange <- colorRampPalette(brewer.pal(9, "Oranges"))(8)

char <- ggplot(varPart2 %>% filter(category=="Watershed characteristics"), aes(x = phylum, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size=5, colour="black"),
        legend.key.width = unit(2, "mm"),
        legend.key.height = unit(2, "mm"),
        axis.title.x = element_text(size = 6, colour = "black"),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 5, colour = "black")) +
  labs(x = "Phylum", y = "Variance explained (%)") +
  guides(fill = guide_legend(title = "")) +
  scale_y_continuous(limits=c(0,20), breaks = seq(0, 20, by = 5), expand=c(0.02,0.02)) +
  scale_fill_manual(values = orange) +
  coord_flip()

char


# plot by category - water chemistry
colourCount = length(unique(varPart2[varPart2$category == "Water chemistry", "variable"]))
getPalette = colorRampPalette(brewer.pal(7, "Spectral"))

blue <- colorRampPalette(brewer.pal(9, "Blues"))(11)

chem <- ggplot(varPart2 %>% filter(category=="Water chemistry"), aes(x = phylum, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size=5, colour="black"),
        legend.key.width = unit(2, "mm"),
        legend.key.height = unit(2, "mm"),
        axis.title.x = element_text(size = 6, colour = "black"),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 5, colour = "black")) +
  labs(x = "Phylum", y = "Variance explained (%)") +
  guides(fill = guide_legend(title = "")) +
  scale_y_continuous(limits=c(0,10), breaks = seq(0, 10, by = 2), expand=c(0.02,0.02)) +
  scale_fill_manual(values = blue) +
  coord_flip()

chem


# combined plot
combined <- wrap_plots(geology,land,chem,char, ncol = 4) +
  plot_layout(axes="collect") &
  theme(legend.position="right")

combined

ggsave("output/variance_partitioning_combined_legend.pdf", combined, width = 384, height = 50, units = "mm", dpi = 600)

combined <- wrap_plots(geology,land,chem,char, ncol = 4) +
  plot_layout(axes="collect") &
  theme(legend.position="none")

ggsave("output/variance_partitioning_combined.pdf", combined, width = 184, height = 50, units = "mm", dpi = 600)
