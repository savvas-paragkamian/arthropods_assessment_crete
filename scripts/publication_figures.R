#!/usr/bin/Rscript

# load packages and functions
library(tidyverse)
library(ggnewscale)
library(sf)
library(raster)
source("functions.R")

g_base <- g_base()

# load data
locations_shp <- sf::st_read("../data/arthropods_occurrences/arthropods_occurrences.shp")
locations_spatial <- sf::st_read("../results/locations_spatial/locations_spatial.shp")
locations_grid <- sf::st_read("../results/locations_grid/locations_grid.shp") 
crete_shp <- sf::st_read("../data/crete/crete.shp")
crete_peaks <- read_delim("../data/crete_mountain_peaks.csv", delim=";", col_names=T) %>% 
    st_as_sf(coords=c("X", "Y"),
             remove=F,
             crs="WGS84")
endemic_species <- read_delim("../results/endemic_species_assessment.tsv", delim="\t")
clc_crete_shp <- st_read("../data/clc_crete_shp/clc_crete_shp.shp")
natura_crete <- sf::st_read("../data/natura2000/natura2000_crete.shp")
wdpa_crete <- sf::st_read("../data/wdpa_crete/wdpa_crete.shp")
natura_crete_land <- st_intersection(natura_crete, crete_shp)

# raster DEM hangling
dem_crete <- raster("../data/dem_crete/dem_crete.tif")
dem_crete_pixel <- as(dem_crete, "SpatialPixelsDataFrame")
dem_crete_df <- as.data.frame(dem_crete_pixel) %>% filter(dem_crete>0)


# split the SPA SCI

natura_crete_land_sci <- natura_crete_land %>% filter(SITETYPE=="B")
## Hotspots and threatspots
endemic_hotspots <- st_read("../results/endemic_hotspots/endemic_hotspots.shp")
threatspots <- st_read("../results/threatspots/threatspots.shp")
threatspots_lt <- threatspots %>% 
    filter(pc_thrt>= quantile(pc_thrt,0.90))

locations_inland <- st_join(locations_shp, crete_shp, left=F)


# Colorblind palette
palette.colors(palette = "Okabe-Ito") 
# Crete figures
## Fig1a
crete_base <- ggplot() +
    geom_sf(crete_shp, mapping=aes()) +
    geom_sf(natura_crete_land_sci,
            mapping=aes(fill="Natura2000 SAC"),
            alpha=1,
            colour="transparent",
            show.legend=T) +
    geom_sf(locations_inland,
            mapping=aes(colour="Sampling sites"),
            size=0.1,
            alpha=0.5,
            show.legend=T) +
    geom_sf(crete_peaks,
            mapping=aes(color = "Places"),
            size=1,
            alpha=1,
            show.legend=T) +
    geom_label(data = crete_peaks, 
               mapping=aes(x = X, y = Y, label = name),
               size = 1.5,
               nudge_x = 0.05,
               nudge_y=0.05, label.padding = unit(0.1, "lines"))+ 
    scale_fill_manual(values = c("Natura2000 SAC" = "#56B4E9"),
                      guide = "legend") +
    scale_colour_manual(values = c("Sampling sites" = "#CC79A7",
                                   "Places"="#D55E00"),
                        guide = "legend") +
    guides(fill = guide_legend(override.aes = list(color = "#56B4E9", alpha=1)),
           colour = guide_legend(override.aes = list(size = c(1,1),
                                                     alpha=c(1,1),
                                                     fill="transparent")))+
    coord_sf(crs="WGS84") +
    theme_bw()+
    theme(axis.title=element_blank(),
          axis.text=element_text(colour="black"),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.box.background = element_blank())


ggsave("../figures/Fig1a.tiff", 
       plot=crete_base, 
       height = 10, 
       width = 20,
       dpi = 600, 
       units="cm",
       device="tiff")

ggsave("../figures/Fig1a.png", 
       plot=crete_base, 
       height = 10, 
       width = 20,
       dpi = 600, 
       units="cm",
       device="png")

## fig1b
crete_hotspot <- ggplot() +
    geom_sf(crete_shp, mapping=aes()) +
    geom_sf(natura_crete_land_sci,
            mapping=aes(fill="Natura2000 SAC"),
            alpha=1,
            colour="transparent",
            show.legend=T) +
    scale_fill_manual(values = c("Natura2000 SAC" = "#56B4E9"))+
    new_scale_fill() +
    geom_sf(endemic_hotspots, mapping=aes(fill=n_species),
            alpha=0.6,
            colour="transparent",
            na.rm = FALSE,
            show.legend=T) +
    scale_fill_gradient(low="#F0E442",
                        high="#D55E00")+
    geom_sf(crete_peaks,
            mapping=aes(),
            colour="#D55E00",
            size=1,
            alpha=1,
            show.legend=F) +
    geom_label(data = crete_peaks, 
               mapping=aes(x = X, y = Y, label = name),
               size = 1.5,
               nudge_x = 0.05,
               nudge_y=0.05, label.padding = unit(0.1, "lines"))+ 
    coord_sf(crs="WGS84") +
    theme_bw()+
    theme(axis.title=element_blank(),
          axis.text=element_text(colour="black"),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.box.background = element_blank())

ggsave("../figures/Fig1b.tiff", 
       plot=crete_hotspot, 
       height = 10, 
       width = 20,
       dpi = 600, 
       units="cm",
       device="tiff")

ggsave("../figures/Fig1b.png", 
       plot=crete_hotspot, 
       height = 10, 
       width = 20,
       dpi = 600, 
       units="cm",
       device="png")

## fig1c
crete_threat <- ggplot() +
    geom_sf(crete_shp, mapping=aes()) +
    geom_sf(natura_crete_land_sci,
            mapping=aes(fill="Natura2000 SAC"),
            alpha=1,
            colour="transparent",
            show.legend=T) +
    scale_fill_manual(values = c("Natura2000 SAC" = "#56B4E9"))+
    new_scale_fill() +
    geom_sf(threatspots_lt, mapping=aes(fill=pc_thrt),
            alpha=0.6,
            colour="transparent",
            size=0.1,
            na.rm = FALSE,
            show.legend=T) +
    scale_fill_gradient(low="#E69F00",
                        high="#CC79A7",
                        breaks = c(26, 36,45,50),
                        labels = c(25,35,45, 51))+
    geom_sf(crete_peaks,
            mapping=aes(),
            colour="#D55E00",
            size=1,
            alpha=1,
            show.legend=F) +
    geom_label(data = crete_peaks, 
               mapping=aes(x = X, y = Y, label = name),
               size = 1.5,
               nudge_x = 0.05,
               nudge_y=0.05, label.padding = unit(0.1, "lines"))+ 
    coord_sf(crs="WGS84") +
    theme_bw()+
    theme(axis.title=element_blank(),
          axis.text=element_text(colour="black"),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.box.background = element_blank())

ggsave("../figures/Fig1c.tiff", 
       plot=crete_threat, 
       height = 10, 
       width = 20,
       dpi = 600, 
       units="cm",
       device="tiff")

ggsave("../figures/Fig1c.png", 
       plot=crete_threat, 
       height = 10, 
       width = 20,
       dpi = 600, 
       units="cm",
       device="png")

## fig1d

crete_corine <- ggplot() +
    geom_sf(crete_shp, mapping=aes()) +
    geom_sf(clc_crete_shp,
            mapping=aes(fill=LABEL1),
            alpha=1,
            colour="transparent",
            show.legend=T) +
    geom_sf(natura_crete_land_sci,
            mapping=aes(color="Natura2000 SAC"),
            linewidth=0.6,
            fill=NA,
            alpha=1,
            show.legend=T) +
    geom_sf(crete_peaks,
            mapping=aes(),
            color = "#D55E00",
            size=1,
            alpha=1,
            show.legend=F) +
    geom_label(data = crete_peaks, 
               mapping=aes(x = X, y = Y, label = name),
               size = 1.5,
               nudge_x = 0.05,
               nudge_y=0.05, label.padding = unit(0.1, "lines"))+ 
    scale_fill_manual(values = c("Artificial surfaces"="#000000",
                                 "Agricultural areas"="#E69F00",
                                 "Forest and semi natural areas" = "#009E73",
                                 "Water bodies" = "#0072B2",
                                 "Natura2000 SAC"=NA),
                      guide = "legend") +
    scale_colour_manual(values = c("Natura2000 SAC" = "#999999"),
                        guide = "legend") +
    guides(fill = guide_legend(override.aes = list(color = "transparent", alpha=1) ),
           colour = guide_legend(override.aes = list(alpha=1, fill="transparent") ) )+
    coord_sf(crs="WGS84") +
    theme_bw()+
    theme(axis.title=element_blank(),
          axis.text=element_text(colour="black"),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.box.background = element_blank(),
          legend.key.size = unit(6, "mm"), 
          legend.text=element_text(size=8))


ggsave("../figures/Fig1d.tiff", 
       plot=crete_corine, 
       height = 10, 
       width = 20,
       dpi = 600, 
       units="cm",
       device="tiff")

ggsave("../figures/Fig1d.png", 
       plot=crete_corine, 
       height = 10, 
       width = 20,
       dpi = 600, 
       units="cm",
       device="png")
