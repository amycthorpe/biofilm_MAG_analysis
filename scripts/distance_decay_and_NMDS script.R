# distance decay and ordination
# extended data figure 2 - distance decay and beta diversity NMDS


# load libaries
library(ggplot2)
library(dplyr)
library(vegan)
library(geosphere)
library(tidyverse)
library(reshape2)
library(ggrepel)

set.seed(123)


# import data
checkm <- read.csv("data/checkm_gtdb.csv", row.names=1, check.names=F, sep=",")

checkm <- mutate(checkm, Completeness_quality = case_when(
  Completeness >= 90 & Contamination <= 5 ~ "Near-complete MAGs",
  Completeness >= 70 & Completeness < 90 & Contamination <= 10 ~ "Medium-quality MAGs",
  TRUE ~ "Low-quality MAGs"
))

checkm <- subset(checkm, Completeness_quality != "Low-quality MAGs") # filter low quality MAGs

coverage <- read.csv("data/finalbins_coverage.csv", row.names=1, check.names=F, sep=",")
coverage <- coverage[rownames(coverage) %in% rownames(checkm), ] # filter low quality MAGs


# convert to relative abundance
coverage_rel <- coverage %>%
  mutate(across(everything(), ~ . / sum(.))) %>% 
  t(.) # samples as rows


# distance matrix
bray_dist <- vegdist(coverage_rel, method = "bray")
bray_dist <- as.matrix(bray_dist)

# geographic matrix
# geo_dist <- distm(lat_long, fun=distGeo) # geographic distances in m calculated between pairwise sites using latitude and longitude
geo_dist <- read.csv("output/geo_dist_matrix.csv", row.names=1, check.names=F, sep=",")

# align matrices
geo_dist_ordered <- geo_dist[row.names(bray_dist), colnames(bray_dist)]

# melt and merge
geo_dist_melt <- melt(as.matrix(geo_dist_ordered))
bray_dist_melt <- melt(as.matrix(bray_dist))

distances <- merge(bray_dist_melt, geo_dist_melt, by = c("Var1", "Var2"))

distances <- distances %>%
  rename(
    bray_distance = value.x,
    geo_distance = value.y
  ) %>%
  mutate(geo_distance_km = geo_distance / 1000) # convert m to km

# plot scatter
distance_decay <- ggplot(distances, aes(x = geo_distance_km, y = bray_distance)) +
  geom_point(alpha = 0.1, size=0.5) +
  geom_smooth(method = "lm") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black")) +
  labs(x = "Geographic distance (km)", y = "Bray-Curtis dissimilarity") +
  scale_x_continuous(limits = c(0, 650), breaks = seq(0, 650, by = 100))

distance_decay

ggsave("output/distance_decay.png", distance_decay, width = 10, height = 6, units = "in", dpi = 600)
ggsave("output/distance_decay.pdf", distance_decay, width = 90, height = 90, units = "mm", dpi = 600)

# regression
model <- lm(geo_distance_km ~ bray_distance, data = distances)
summary(model)


# NMDS
nmds <- metaMDS(coverage_rel)

# extract scores
nmds_scores<-as.data.frame(nmds[["points"]])

# plot with latitude contour lines
location <- read.csv("data/latitude.csv", row.names=1, check.names=F, sep=",")

location <- location[match(rownames(coverage_rel), rownames(location)),]

surf <- ordisurf(nmds ~ location)

surf_grid <- expand.grid(x=surf$grid$x, y=surf$grid$y)
surf_grid$z <- as.vector(surf$grid$z)

plot<-ggplot(nmds_scores,aes(x=MDS1,y=MDS2))+
  theme_bw()+
  geom_point(size=0.1,alpha=0.7,stroke=1)+
  geom_contour(data=surf_grid, aes(x=x, y=y, z=z, colour=..level..)) +
  scale_color_gradient(low = "blue", high = "red", name="Latitude") +
  theme(axis.title = element_text(size=6,colour="black"),
        axis.text.x = element_text(size=5,colour="black"),
        axis.text.y = element_text(size=5,colour="black"),
        panel.border=element_rect(colour="black",fill=NA,size=0.6),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.text = element_text(size=5,colour="black"),
        legend.title = element_text(size=6,colour="black"))+
  geom_vline(xintercept=c(0,0), linetype="dashed", colour="black",size=0.6)+
  geom_hline(yintercept=c(0,0), linetype="dashed",colour="black",size=0.6)+
  labs(x="NMDS1",y="NMDS2")

plot

ggsave("output/NMDS_contour.png", plot, width = 6, height = 5, units = "in", dpi = 600)
ggsave("output/NMDS_contour.pdf", plot, width = 90, height = 90, units = "mm", dpi = 600)


# NMDS and envfit
# import metadata
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

water_chem <- water_chem %>% select(-TN) # TN has lots of NAs so remove

land_cover <- read.csv("data/catchment_land_cover.csv", check.names=F, row.names=1, sep=",")
geology <- read.csv("data/catchment_geology.csv", check.names=F, row.names=1, sep=",")
chars <- read.csv("data/catchment_chars.csv", check.names=F, row.names=1, sep=",")

# merge metadata
meta_merged1 <- merge(land_cover, geology, by="row.names") %>% 
  column_to_rownames(., var="Row.names")

colnames(meta_merged1) <- gsub("[ -]", "_", colnames(meta_merged1)) # replace spaces in colnames
colnames(meta_merged1) <- gsub(",", "", colnames(meta_merged1)) # remove commas

meta_merged2 <- merge(meta_merged1, chars, by="row.names") %>% 
  column_to_rownames(., var="Row.names")
meta_merged3 <- merge(meta_merged2, water_chem, by="row.names") %>% 
  column_to_rownames(., var="Row.names")

# calculate log ratio calcareous:silleous
constant <- 0.01 # avoid zero

meta_merged3$log_ratio_calc_sil <- log((meta_merged3$Calcareous_pct + constant) / 
                                         (meta_merged3$Siliceous_pct + constant))

columns_to_include <- colnames(meta_merged3)

formula <- as.formula(paste("nmds ~", paste(columns_to_include, collapse="+"))) # set envfit formula

meta_merged3 <- meta_merged3[match(rownames(coverage_rel), rownames(meta_merged3)),] # match order of rownames (samples)

nmds.fit <- envfit(formula,
                   data=meta_merged3,
                   perm=999,
                   na.rm=T)

nmds.fit

df<-nmds.fit$vectors # get datafrane of vector r and p values
df_r<-as.data.frame(df[["r"]])
df_pvals<-as.data.frame(df[["pvals"]])
df_arrows<-as.data.frame(df[["arrows"]])
nmds.fit_df<-cbind(df_arrows,df_r,df_pvals)
colnames(nmds.fit_df)[3] <- "rvalue"
colnames(nmds.fit_df)[4] <- "pvalue"

nmds.fit_df <- nmds.fit_df %>% # add significance levels
  mutate(significance = case_when(
    pvalue > 0.05 ~ "ns",
    pvalue <= 0.05 & pvalue > 0.01 ~ "*",
    pvalue <= 0.01 & pvalue > 0.001 ~ "**",
    pvalue <= 0.001 ~ "***"
  ))

write.csv(nmds.fit_df, "output/envfit.csv")

nmds_scores<-as.data.frame(nmds[["points"]]) # extract scores
nmds_scores_meta<-merge(nmds_scores,meta_merged3,by="row.names") %>% # combine with meta
  column_to_rownames(.,"Row.names")

vectors<-as.data.frame(scores(nmds.fit, "vectors")) * ordiArrowMul(nmds.fit) # set up vector arrows
arrowhead = arrow(length = unit(0.02, "npc"))

vectors <- rownames_to_column(vectors, var="variable") # vector groups
vectors$category <- ifelse(vectors$variable %in% colnames(water_chem), "Water chemistry",
                           ifelse(vectors$variable %in% c(colnames(geology), "log_ratio_calc_sil"), "Geology",
                                  ifelse(vectors$variable %in% colnames(land_cover), "Land cover",
                                         "Watershed characteristics")))
vectors <- column_to_rownames(vectors, var="variable")

# filter vectors
# nmds.fit_sub <- nmds.fit_df %>% 
#   filter(., pvalue <= 0.05, rvalue >= 0.05) # only significant

# vectors_sub <- vectors[rownames(vectors) %in% rownames(nmds.fit_sub), ] 
vectors_sub <- vectors[rownames(vectors) %in% colnames(water_chem), ] # focus on water chem

# plot
plot<-ggplot(nmds_scores_meta,aes(x=MDS1,y=MDS2))+
  theme_bw()+
  geom_point(size=0.1,alpha=0.7,stroke=1, aes(colour=log_ratio_calc_sil))+
  theme(axis.title = element_text(size=6,colour="black"),
        axis.text.x = element_text(size=5,colour="black"),
        axis.text.y = element_text(size=5,colour="black"),
        panel.border=element_rect(colour="black",fill=NA,size=0.6),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.text = element_text(size=5,colour="black"),
        legend.title = element_text(size=6,colour="black"))+
  geom_vline(xintercept=c(0,0), linetype="dashed", colour="black",size=0.6)+
  geom_hline(yintercept=c(0,0), linetype="dashed",colour="black",size=0.6)+
  geom_segment(aes(x=0,y=0,xend=NMDS1,yend=NMDS2),data=vectors_sub,size=0.4,arrow=arrowhead, colour="black")+
  # geom_label(data=vectors_sub,aes(x=1.1*NMDS1,y=1.1*NMDS2),label=row.names(vectors_sub),size=2,color="black")+
  geom_text_repel(data = vectors_sub,aes(x = 1 * NMDS1, y = 1 * NMDS2, label = row.names(vectors_sub)),size=2,colour = "black",
    force = 5, segment.color=NA) +
  labs(x="NMDS1",y="NMDS2") +
  scale_colour_viridis_c(option="D", na.value="grey40", name="Calcareous:Siliceous")

plot

ggsave("output/NMDS_geology_water.png", plot, width = 6, height = 5, units = "in", dpi = 600)
ggsave("output/NMDS_geology_water.pdf", plot, width = 120, height = 90, units = "mm", dpi = 600)
