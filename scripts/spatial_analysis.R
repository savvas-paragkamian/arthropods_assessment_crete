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

## Spatial data

habitats_crete <- raster("../data/habitats_crete/habitats_crete.tif")
habitats_meta <- read_delim("../data/habitats_crete/CLC2018_CLC2018_V2018_20_QGIS.txt", delim=",", col_names=F)
habs = raster::extract(nlcd, zion, df = TRUE, factors = TRUE)
habitats <- raster("~/Downloads/u2018_clc2018_v2020_20u1_raster100m/DATA/U2018_CLC2018_V2020_20u1.tif")


crete_epsg <- st_transform(crete_shp, crs="EPSG:3035")

habitats_crete <- crop(habitats, crete_epsg)
wgs84 = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
habitats_crete <- projectRaster(habitats_crete, crs = wgs84, method = "ngb")
writeRaster(habitats_crete, filename="../data/habitats_crete/habitats_crete.tif")

dem_crete <- raster("../data/dem_crete/dem_crete.tif")

dem_crete_pixel <- as(dem_crete, "SpatialPixelsDataFrame")
dem_crete_df <- as.data.frame(test_spdf)


g_dem <- g_base +
    geom_raster(dem_crete_df, mapping=aes(x=x, y=y, fill=dem_crete))+
    geom_sf(locations_shp, mapping=aes(),color="blue", size=0.1, alpha=0.2)

ggsave("../plots/crete_dem.png", plot=g_dem, device="png")


locations_shp$elevation <- raster::extract(dem_crete, locations_shp, cellnumbers=F)


g_ele <- g_base + 
    geom_sf(locations_shp, mapping=aes(color=elevation), size=0.1, alpha=0.2) +
    scale_color_gradientn(colours = terrain.colors(10)) 

ggsave("../plots/crete_occurrences_dem.png", plot=g_ele, device="png")


st_write(locations_shp,"../results/locations_spatial/locations_spatial.shp")
