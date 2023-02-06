#!/usr/bin/Rscript

# load packages and functions
library(tidyverse)
library(ggnewscale)
library(sf)
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

# split the SPA SCI

natura_crete_land_sci <- natura_crete_land %>% filter(SITETYPE=="B")
## Hotspots and threatspots
endemic_hotspots <- st_read("../results/endemic_hotspots/endemic_hotspots.shp")
threatspots <- st_read("../results/threatspots/threatspots.shp")

locations_inland <- st_join(locations_shp, crete_shp, left=F)


# Colorblind palette
palette.colors(palette = "Okabe-Ito") 
# Crete figures

crete_base <- ggplot() +
    geom_sf(crete_shp, mapping=aes()) +
    geom_sf(natura_crete_land_sci,
            mapping=aes(fill="Natura2000 SAC"),
#            alpha=0.2,
#            size=0.1,
            colour="transparent",
            show.legend=T) +
    geom_sf(locations_inland,
            mapping=aes(colour="Occurrences"),
            size=0.1,
            alpha=0.5,
            show.legend=T) +
#    geom_sf(crete_peaks,
#            mapping=aes(color = "Peaks"),
#            size=0.5,
#            alpha=0.5,
#            show.legend=T) +
    geom_label(data = crete_peaks, 
               mapping=aes(x = X, y = Y, label = name),
               size = 1, nudge_x = 0.05, nudge_y=0.05)+ 
    scale_fill_manual(values = c("Natura2000 SAC" = "#E69F00"),
                      guide = "legend") +
    scale_colour_manual(values = c("Occurrences" = "#CC79A7"),
                        guide = "legend") +
    guides(fill = guide_legend(override.aes = list(color = "#E69F00", alpha=1) ),
           colour = guide_legend(override.aes = list(size = 1,alpha=1, fill="transparent") ) )+
    coord_sf(crs="WGS84") +
    theme_bw()+
    theme(axis.title=element_blank(),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.box.background = element_blank())


ggsave("../figures/Fig1.tiff", 
       plot=crete_base, 
       height = 10, 
       width = 20,
       dpi = 600, 
       units="cm",
       device="tiff")
