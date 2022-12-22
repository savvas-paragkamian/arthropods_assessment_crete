#!/usr/bin/Rscript

## Script name: code_snippets.R
##
## Purpose of script: here are some parts of code to crop spatial 
## objects around the area of Crete and other test code.
##
## Author: Savvas Paragkamian
##
## Date Created: 2022-12-22

library(tidyverse)
library(readxl)
library(rredlist)
library(taxize)
library(sf)
library(raster)
library(units)
#library(rgdal)
library(ConR)
library(vegan)



## Load shapefiles and crop only the Crete overlap

periphereies_shp <- sf::st_read("~/Documents/spatial_data/periphereies/periphereies.shp")

#### crete is the 12th region in this shapefile
#### https://geodata.gov.gr/en/dataset/28121643-d977-48eb-a8ca-a6fac6b4af6d/resource/7c80a2c1-93b7-4814-9fc4-245e775acaa6/download/periphereies.zip
#
crete_shp <- periphereies_shp[12,] %>% st_transform(., "WGS84")

st_write(crete_shp,"../data/crete/crete.shp",
         layer_options = "ENCODING=UTF-8", delete_layer = TRUE)

# natura2000 shapefile downloaded from here https://www.eea.europa.eu/data-and-maps/data/natura-13 
natura2000 <- sf::st_read("~/Downloads/Natura2000_end2021_Shapefile/Natura2000_end2021_epsg3035.shp") 
                   %>% st_transform(., "WGS84")

natura_crete <- st_crop(natura2000,
                        y=st_bbox(c(xmin=23,ymin=34,xmax=27,ymax=36),
                                  crs="WGS84"))

st_write(natura_crete,"../data/natura2000/natura2000_crete.shp",
         layer_options = "ENCODING=UTF-8", delete_layer = TRUE)

# World Database of Protected Areas
wdpa_gr_0 <- sf::st_read("~/Documents/spatial_data/WDPA_WDOECM_Dec2022_Public_GRC_shp/WDPA_WDOECM_Dec2022_Public_GRC_shp_0/WDPA_WDOECM_Dec2022_Public_GRC_shp-polygons.shp")

wdpa_gr_1 <- sf::st_read("~/Documents/spatial_data/WDPA_WDOECM_Dec2022_Public_GRC_shp/WDPA_WDOECM_Dec2022_Public_GRC_shp_1/WDPA_WDOECM_Dec2022_Public_GRC_shp-polygons.shp")

wdpa_gr_2 <- sf::st_read("~/Documents/spatial_data/WDPA_WDOECM_Dec2022_Public_GRC_shp/WDPA_WDOECM_Dec2022_Public_GRC_shp_2/WDPA_WDOECM_Dec2022_Public_GRC_shp-polygons.shp")

wdpa_gr <- rbind(wdpa_gr_0, wdpa_gr_1, wdpa_gr_2)

wdpa_crete <- st_intersection(wdpa_gr, crete_shp)

st_write(wdpa_crete,"../data/wdpa_crete/wdpa_crete.shp",
         layer_options = "ENCODING=UTF-8", delete_layer = TRUE)

#grid_iucn <- sf::st_read(dsn="~/Documents/AOOGrid_10x10kmshp/AOOGrid_10x10km.shp") %>% 
#    st_transform(., crs="WGS84")
##

arthropods <- readxl::read_excel("../data/arthropoda_crete_nhmc_for_analysis.xlsx")
# remove Opiliones because the dataset has some errors that are under review.
arthropods_kriti <- readxl::read_excel("../data/arthropoda_crete_nhmc_for_analysis.xlsx") %>% 
    filter(Order!="Opiliones") %>% 
    dplyr::select(-Ergasia)
# Data transformation for ConR package

locations_inland_df <- arthropods_occurrences %>%
    dplyr::rename(ddlon=logD, ddlat=latD, tax=subspeciesname) %>% 
    dplyr::select(-Order) %>%
    relocate(ddlat,ddlon, tax)

crete_spatial <- as(st_geometry(crete_polygon),"Spatial")  

eoo_results_list <- EOO.computing(locations_inland_df, 
                                  country_map=crete_spatial, 
                                  exclude.area=T,
                                  export_shp=T,
                                  write_shp=T)

eoo_results <- EOO.computing(locations_inland_df,export_shp=T, write_shp=T)

eoo_results <- read_delim("EOO.results.csv", delim=",", col_names=T)

# Raster spatial data

dem <- raster("~/Documents/spatial_data/EAA-DEM-EUD_CP-DEMS_5500015000-AA/EUD_CP-DEMS_5500015000-AA.tif")

crete_epsg <- st_transform(crete_shp, crs="EPSG:3035")
dem_crete <- crop(dem, crete_epsg)
wgs84 = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
dem_crete <- projectRaster(dem_crete, crs = wgs84, method = "ngb")
writeRaster(dem_crete, filename="../data/dem_crete/dem_crete.tif")
# ConR package returned an error with the exclution of the land
# so a new function named "eoo_calculation" is created
g2 <- ggplot() +
    geom_sf(crete_polygon, mapping=aes()) +
    geom_sf(locations_inland, mapping=aes(),color="blue", size=0.1, alpha=0.2) +
    geom_sf(crete_grid10m, mapping=aes(),color="red", alpha=0.2, size=0.1) +
    geom_sf(grid_10k_species_s, mapping=aes(fill=species), alpha=0.8, size=0.1) +
    scale_fill_gradient(low = "white", high = "black")+
    coord_sf(crs="WGS84") +
    theme_bw()

ggsave("tst.png", 
       plot =g2, 
       device = "png",
       width = 30,
       height = 30,
       units = "cm",
       dpi = 300 ,
       path = "../plots/")

# facet hot spot 
map_greece_plot_grid_hot_spot_facet <- ggplot()+
  geom_polygon(data = grid_10km_species_data,
               aes(x=long, y=lat,group = group, fill=hot_spot),
               lwd=0.12, color="orange")+
  geom_polygon(data = hellenic_borders_df,
               aes(x=long, y=lat,group = group),
               lwd=0.12,color="black")+
  geom_polygon(data = grid_10km_df,
               aes(x=long, y=lat,group = group),
               lwd=0.12, fill=NA, color="orange")+
  geom_point(data = locations,aes(x=logD, y=latD),size = 0.5)+
  ggtitle("hot spot")+
  scale_fill_gradient(low = "yellow", high = "red", na.value = NA)+ 
  labs(x="Longitude",y="Latitude")+
#  scale_fill_gradientn(colours = c("gray100",terrain.colors(10)),na.value =NA ,name="Number of species")+ #c("gray100","gray50","gray40","gray35","gray30","gray20","gray10","gray0")
  scale_x_continuous(breaks = seq(23,27,0.5),limits = c(23,27))+
  scale_y_continuous(breaks = seq(34,36,0.5),limits = c(34,37))+
  coord_map(xlim = c(23,27), ylim = c(34,37))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = c(0.87, 0.73),
        legend.text = element_text(size=9),
        legend.title = element_text(size=10))
  facet_wrap(~Order)
    
ggsave("map_greece_plot_grid_hot_spot.png", plot = map_greece_plot_grid_hot_spot, device = "png",width = 30,height = 30,units = "cm",dpi = 300 ,path = "../plots/")

map_greece_plot_grid_endemic <- ggplot()+
  geom_polygon(data = grid_10km_species_data,
               aes(x=long, y=lat,group = group, fill=endemic_species),
               lwd=0.12, color="orange")+
  geom_polygon(data = hellenic_borders_df,
               aes(x=long, y=lat,group = group),
               lwd=0.12,color="black")+
  geom_polygon(data = grid_10km_df,
               aes(x=long, y=lat,group = group),
               lwd=0.12, fill=NA, color="orange")+
  geom_point(data = locations,aes(x=logD, y=latD,color=Order),size = 0.2)+
  scale_fill_gradient(low = "yellow", high = "red", na.value = NA)+ 
  ggtitle("endemic")+
  labs(x="Longitude",y="Latitude")+
#  scale_fill_gradientn(colours = c("gray100",terrain.colors(10)),na.value =NA ,name="Number of species")+ #c("gray100","gray50","gray40","gray35","gray30","gray20","gray10","gray0")
  scale_x_continuous(breaks = seq(23,27,0.5),limits = c(23,27))+
  scale_y_continuous(breaks = seq(34,36,0.5),limits = c(34,37))+
  coord_map(xlim = c(23,27), ylim = c(34,37))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),legend.position = c(0.87, 0.73),legend.text = element_text(size=9),legend.title = element_text(size=10))
#geom_text(data = sisquoc, aes(label = paste("  ", as.character(name), sep="")), angle = 60, hjust = 0, color = "yellow")
    
ggsave("map_greece_plot_grid_endemic.png", plot = map_greece_plot_grid_endemic, device = "png",width = 30,height = 30,units = "cm",dpi = 300 ,path = "../plots/")

