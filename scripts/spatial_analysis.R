#!/usr/bin/Rscript

## Script name: spacial_analysis.R
##
## Purpose of script: Spatial analysis using data from the Copernicus system
## for elevation, habitats and corine land change
##
## How to run:
## Rscript spacial_analysis.R
##
## Execution time: 
##
## Author: Savvas Paragkamian
##
## Date Created: 2022-12-22

library(tidyverse)
library(sf)
library(raster)
source("functions.R")

# Load data
g_base <- g_base()
locations_shp <- sf::st_read("../data/arthropods_occurrences/arthropods_occurrences.shp")
crete_shp <- sf::st_read("../data/crete/crete.shp")
endemic_species <- read_delim("../results/endemic_species_paca.tsv", delim="\t")
natura_crete <- sf::st_read("../data/natura2000/natura2000_crete.shp")
wdpa_crete <- sf::st_read("../data/wdpa_crete/wdpa_crete.shp")

natura_crete_land <- st_intersection(natura_crete, crete_shp)

# split the SPA SCI

natura_crete_land_sci <- natura_crete_land %>% filter(SITETYPE=="B")

# Spatial data

## Habitats
habitats_crete <- raster("../data/habitats_crete/habitats_crete.tif")
locations_shp$habitats <- raster::extract(habitats_crete, locations_shp, cellnumbers=F)

habitats_meta <- read_delim("../data/habitats_crete/CLC2018_CLC2018_V2018_20_QGIS.txt", delim=",", col_names=F)

colnames(habitats_meta) <- c("code", "r","g","b","255", "description")

habitats_meta$id <- seq(1:nrow(habitats_meta))

habitats_meta$hex <- rgb(habitats_meta$r,
                         habitats_meta$g,
                         habitats_meta$b,
                         maxColorValue=255)


habitats_crete_pixel <- as(habitats_crete, "SpatialPixelsDataFrame")
habitats_crete_df <- as.data.frame(habitats_crete_pixel) %>%
    left_join(habitats_meta, by=c("habitats_crete"="id"))

habitats_crete_df$description <- factor(habitats_crete_df$description, 
                                        levels=unique(habitats_crete_df$description))

habitats_crete_df$hex <- factor(habitats_crete_df$hex, 
                                        levels=unique(habitats_crete_df$hex))


colors_all <- setNames(habitats_meta$hex,habitats_meta$description)
colors_crete <- colors_all[colors_all %in% unique(habitats_crete_df$hex)]
g_hab <- ggplot() +
    geom_raster(habitats_crete_df, mapping=aes(x=x, y=y, fill=description)) +
    scale_fill_manual(values=colors_crete)+
    theme(legend.position="bottom", legend.margin=margin()) +
    guides(fill=guide_legend(nrow=8,byrow=TRUE, title="")) +
    coord_equal() 

#    geom_sf(locations_shp, mapping=aes(),color="blue", size=0.1, alpha=0.2)

ggsave("../plots/crete_habitats.png",
       plot=g_hab,
       width = 50,
       height = 30,
       units='cm', 
       device = "png",
       dpi = 300)

dim_x <- res(habitats_crete)[1]
dim_y <- res(habitats_crete)[2]
habitats_summary <- as.data.frame(habitats_crete) %>% 
    group_by(habitats_crete) %>% 
    tally() %>% 
    mutate(area=n * dim_x * dim_y) %>%
    left_join(habitats_meta, by=c("habitats_crete"="id"))
## Dem

dem_crete <- raster("../data/dem_crete/dem_crete.tif")

locations_shp$elevation <- raster::extract(dem_crete, locations_shp, cellnumbers=F)

dem_crete_pixel <- as(dem_crete, "SpatialPixelsDataFrame")
dem_crete_df <- as.data.frame(dem_crete_pixel)


g_dem <- g_base +
    geom_raster(dem_crete_df, mapping=aes(x=x, y=y, fill=dem_crete))+
    geom_sf(locations_shp, mapping=aes(),color="blue", size=0.1, alpha=0.2)

ggsave("../plots/crete_dem.png", plot=g_dem, device="png")


g_ele <- g_base + 
    geom_sf(locations_shp, mapping=aes(color=elevation), size=0.1, alpha=0.2) +
    scale_color_gradientn(colours = terrain.colors(10)) 

ggsave("../plots/crete_occurrences_dem.png", plot=g_ele, device="png")

# Export of locations shapefile with all the spatial metadata
st_write(locations_shp,"../results/locations_spatial/locations_spatial.shp", append=TRUE)
