# distance decay and ordination
# extended data figure 2 - distance decay and beta diversity NMDS


# load libaries
library(ggplot2)
library(dplyr)
library(vegan)
library(geosphere)
library(tidyverse)
library(reshape2)

set.seed(123)


# import data
coverage <- read.csv("data/finalbins_coverage.csv", row.names=1, check.names=F, sep=",")


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
