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
library(terra)
library(tidyterra)
library(units)
source("functions.R")

# Load data
g_base <- g_base()
locations_shp <- sf::st_read("../data/arthropods_occurrences/arthropods_occurrences.shp")
crete_shp <- sf::st_read("../data/crete/crete.shp")
endemic_species <- read_delim("../results/endemic_species_assessment.tsv", delim="\t")
endemic_hotspots <- st_read("../results/endemic_hotspots/endemic_hotspots.shp")
threatspots <- st_read("../results/threatspots/threatspots.shp")
threatspots_lt <- threatspots |> 
    filter(pc_thrt>= quantile(pc_thrt,0.90))

natura_crete <- sf::st_read("../data/natura2000/natura2000_crete.shp")
wdpa_crete <- sf::st_read("../data/wdpa_crete/wdpa_crete.shp")

wdpa_crete_wildlife <- wdpa_crete |> filter(DESIG_ENG=="Wildlife Refugee")
natura_crete_land <- st_intersection(natura_crete, crete_shp)

# split the SPA SCI

natura_crete_land_sci <- natura_crete_land |> filter(SITETYPE=="B")

# Spatial data

locations_shp <- st_join(locations_shp, natura_crete_land_sci, left=T)
## CORINE Land Cover nomenclature
clc_crete_shp <- st_read("../data/clc_crete_shp/clc_crete_shp.shp")

clc_crete_colors <- clc_crete_shp |> 
    st_drop_geometry() |>
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
    
    merged <- clc_crete_summary[[i]] |>
        left_join(clc_area_ovelaps_df[[i]])

    write_delim(merged,
                paste0("../results/clc_crete_", 
                       names(clc_area_ovelaps_df)[i],
                       ".tsv",sep=""),
                delim="\t")
}

## Dem

dem_crete <- rast("../data/dem_crete/dem_crete.tif")

locations_shp$elevation <- terra::extract(dem_crete, locations_shp, cellnumbers=F)

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


## Crete land use change
colors_clc_label2_v <- c("Artificial, non-agricultural vegetated areas"="#000000",
"Open spaces with little or no vegetation"="#A77300",
"Pastures"="#98BA6A",
"Mine, dump and construction sites"="#46B0D3",
"Forests"="#07A07D",
"Scrub and/or herbaceous vegetation associations"="#98CA53",
"Urban fabric" ="#A4A869",
"Inland waters"="#1370A1",
"Permanent crops"="#AE6120",
"Industrial, commercial and transport units"="#D06C5B",
"Heterogeneous agricultural areas"="#BE81A3",
"Arable land"="#999999")

hilda_cat <- data.frame(hilda_id = c("11","22","33","44","55","66","77"),
                        hilda_name=c("urban","cropland","pasture/rangeland",
                                     "forest", "unmanaged grass/shrubland","sparse/no vegetation", "water"),
                        hilda_hex=c("#000000","#AE6120","#98BA6A","#07A07D","#BE81A3","#999999", "#1370A1"))

hilda_cat_v <- c("urban"="#000000",
                 "cropland"="#AE6120",
                 "pasture/rangeland"="#98BA6A",
                 "forest"="#07A07D",
                 "unmanaged grass/shrubland"="#BE81A3",
                 "sparse/no vegetation"="#999999",
                 "water"="#1370A1")

hilda_1994 <- rast("../data/hildap_GLOB-v1.0_lulc-states_crete/crete_hilda_plus_1994_states_GLOB-v1-0_wgs84-nn.tif")
hilda_1994_df <- terra::as.data.frame(hilda_1994, xy=TRUE, cells=TRUE) |>
    filter(`hilda_plus_1994_states_GLOB-v1-0_wgs84-nn`>0) |>
    mutate(hilda_1994=as.character(`hilda_plus_1994_states_GLOB-v1-0_wgs84-nn`)) |>
    left_join(hilda_cat, by=c("hilda_1994"="hilda_id"))

hilda_1994_df$hilda_name <- factor(hilda_1994_df$hilda_name, levels=as.character(unique(sort(hilda_1994_df$hilda_name))))

g_hilda_1994 <- g_base +
    geom_raster(hilda_1994_df,
                mapping=aes(x=x, y=y, fill=hilda_name)) +
    scale_fill_manual(values=hilda_cat_v)

ggsave("../plots/crete_hilda_1994.png", plot=g_hilda_1994, device="png")

#### Hilda analysis
hilda_path <- "../data/hildap_GLOB-v1.0_lulc-states_crete/"
hilda_files <- list.files(hilda_path)

for (i in hilda_files){

    i=100
    filename <- hilda_files[i]
    
    raster_file <- rast(paste0(hilda_path,hilda_files[i],sep=""))
    
    raster_name <- paste0("hilda_",gsub(".*([0-9]{4}).*", "\\1", filename),sep="")
    
    raster_df <- terra::as.data.frame(raster_file, xy=TRUE, cells=TRUE)

    raster_df <- raster_df |>
        mutate(hilda_id=as.character(raster_df[,4])) |>
        filter(raster_df[,4]>0) |>
        left_join(hilda_cat, by=c("hilda_id"="hilda_id"))
    
    raster_df$hilda_name <- factor(raster_df$hilda_name, levels=as.character(unique(sort(raster_df$hilda_name))))
    
    g_hilda_map <- g_base +
        geom_raster(raster_df,
                    mapping=aes(x=x, y=y, fill=hilda_name)) +
        scale_fill_manual(values=hilda_cat_v) +
        guides(fill = guide_legend(nrow = 1)) +
        ggtitle(raster_name)+
        theme(axis.title=element_blank(),
              legend.position="bottom",
              legend.key.size = unit(4, "mm"), 
              legend.text=element_text(size=8),
              legend.title=element_blank())
    
    ggsave(paste0("../plots/crete_",raster_name,"_map.png",sep=""),
           plot=g_hilda_map,
           height = 10, 
           width = 20,
           dpi = 300, 
           units="cm",
           device="png")
    
    hilda_sum <- zonal(cellSize(raster_file), raster_file, "sum") |> 
        mutate(area_m2=units::set_units(area,m^2)) |>
        mutate(area=units::set_units(area/10^6, km^2)) 
    
    hilda_sum <- hilda_sum |>
        mutate(hilda_id=as.character(hilda_sum[,1])) |>
        filter(hilda_sum[,1]>0) |>
        left_join(hilda_cat) 
    
    hilda_sum_g <- ggplot()+
        geom_col(hilda_sum,
                 mapping= aes(x=area,
                              y="",
                              fill = hilda_name),
                 position = position_stack()) +
        scale_fill_manual(values=hilda_cat_v) +
        theme_bw()+
        theme(legend.position='none',
            panel.border = element_blank(),
            panel.grid.major = element_blank(), #remove major gridlines
            panel.grid.minor = element_blank()) #remove minor gridlines
    
    ggsave(paste0("../plots/crete_",raster_name,"_bar.png",sep=""),
           plot=hilda_sum_g,
           height = 10, 
           width = 20,
           dpi = 300, 
           units="cm",
           device="png")
}
# Export of locations shapefile with all the spatial metadata
st_write(locations_shp,
         "../results/locations_spatial/locations_spatial.shp",
         layer_options = "ENCODING=UTF-8", 
         delete_layer=T, 
         delete_dsn = TRUE)
         append=TRUE)

st_write(locations_shp,"../results/locations_spatial/locations_spatial.csv")



