#!/usr/bin/Rscript

## Script name: spacial_analysis.R
##
## Purpose of script: Spatial analysis using data from the Copernicus system
## for elevation, habitats and corine land change. Main output is an
## enriched locations file which is exported for further analysis.
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
library(units)
source("functions.R")

# Load data
g_base <- g_base()
locations_shp <- sf::st_read("../data/arthropods_occurrences/arthropods_occurrences.shp")
crete_shp <- sf::st_read("../data/crete/crete.shp")
endemic_species <- read_delim("../results/endemic_species_assessment.tsv", delim="\t")
endemic_hotspots <- st_read("../results/endemic_hotspots/endemic_hotspots.shp")
threatspots <- st_read("../results/threatspots/threatspots.shp")
threatspots_lt <- threatspots %>% 
    filter(pc_thrt>= quantile(pc_thrt,0.90))

natura_crete <- sf::st_read("../data/natura2000/natura2000_crete.shp")
wdpa_crete <- sf::st_read("../data/wdpa_crete/wdpa_crete.shp")

wdpa_crete_wildlife <- wdpa_crete %>% filter(DESIG_ENG=="Wildlife Refugee")
natura_crete_land <- st_intersection(natura_crete, crete_shp)

# split the SPA SCI

natura_crete_land_sci <- natura_crete_land %>% filter(SITETYPE=="B")

# Spatial data

locations_shp <- st_join(locations_shp, natura_crete_land_sci, left=T)
## CORINE Land Cover nomenclature
clc_crete_shp <- st_read("../data/clc_crete_shp/clc_crete_shp.shp")

clc_crete_colors <- clc_crete_shp %>% 
    st_drop_geometry() %>%
    distinct(LABEL3,hex)
## create the named vector in order to plot with correct colors
##
colors_all <- setNames(clc_crete_colors$hex,clc_crete_colors$LABEL3)

g_clc_s <- g_base + 
    geom_sf(clc_crete_shp, mapping=aes(fill=LABEL3,color=LABEL3))+
    scale_fill_manual(values=colors_all)+
    scale_color_manual(values=colors_all)+
    theme(legend.position="bottom", legend.margin=margin()) +
    guides(fill=guide_legend(nrow=8,byrow=TRUE, title=""), color="none") 

ggsave("../plots/clc_crete_shp.png", 
       plot=g_clc_s,
       width = 50,
       height = 30,
       units='cm',
       device = "png",
       dpi = 300)


locations_shp <- st_join(locations_shp, clc_crete_shp, left=T)

### Summary of overlaps of CLC with areas

list_sf <- list(natura2000 = natura_crete_land_sci,
                wildlife = wdpa_crete_wildlife,
                hotspots = endemic_hotspots,
                threatspots = threatspots_lt)

clabels <- c("LABEL1", "LABEL2", "LABEL3")

### Here we calculate the overlap of areas with CLC
clc_area_ovelaps <- area_overlap_combination(clc_crete_shp, list_sf, clabels)

### data transformation
clc_area_ovelaps_df <- convert_nested_l_df(clc_area_ovelaps)

### Here we calculate the area of each LABEL for Crete
clc_crete_summary <- lapply(clabels, function(x) spatial_area_summary(clc_crete_shp, x))

### Here we merge the total area with the overlaps and 
### print the output

for (i in seq_along(clc_crete_summary)){
    
    merged <- clc_crete_summary[[i]] %>%
        left_join(clc_area_ovelaps_df[[i]])

    write_delim(merged,
                paste0("../results/clc_crete_", 
                       names(clc_area_ovelaps_df)[i],
                       ".tsv",sep=""),
                delim="\t")
}

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
st_write(locations_shp,
         "../results/locations_spatial/locations_spatial.shp",
         layer_options = "ENCODING=UTF-8", 
         delete_layer=T, 
         delete_dsn = TRUE)
         append=TRUE)

st_write(locations_shp,"../results/locations_spatial/locations_spatial.csv")
