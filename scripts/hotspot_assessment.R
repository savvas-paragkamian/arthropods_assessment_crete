#!/usr/bin/Rscript

## Script name: hotspot_assessment.R
##
## Purpose of script:
## Propose endemic and threatend species hotspots of the 
## Cretan endemic anthropods.
##
## How to run:
## Rscript hotspot_assessment.R
##
## Execution time: 
##
## Author: Savvas Paragkamian
##
## Date Created: 2024-01-14

library(tidyverse)
library(sf)
library(terra)
library(quadtree)
library(units)
source("functions.R")

# Load data
endemic_species <- read_delim("../results/endemic_species_assessment.tsv", delim="\t") 

locations_grid <- st_read("../results/locations_grid/locations_grid.shp")

# Adaptive resolution with quadtree

grid_10km <- terra::vect("../data/Greece_shapefile/gr_10km.shp") 
grid_10km_wgs <- terra::project(grid_10km,crete_shp)

### keep grid that overlaps with Crete
crete_bbox_polygon <- st_as_sf(st_as_sfc(st_bbox(crete_shp)))
crete_grid10 <- terra::crop(grid_10km_wgs, crete_bbox_polygon)

g_base + 
    geom_sf(sf::st_as_sf(crete_grid10), mapping=aes())

## Endemic hotspots

endemic_hotspots <- locations_grid %>%
    group_by(CELLCODE) %>%
    summarise(n_species=n()) %>%
    filter(n_species >= quantile(n_species, 0.90))

st_write(endemic_hotspots,
         "../results/endemic_hotspots/endemic_hotspots.shp", 
         append=F,
         delete_layer=T,
         delete_dsn = TRUE) 

endemic_hotspots_order <- locations_grid %>% 
    group_by(CELLCODE, Order) %>%
    summarise(n_species=n(), .groups="drop") %>%
    group_by(Order) %>%
    mutate(quant90= quantile(n_species, 0.90)) %>%
    filter(n_species >= quant90)

## threatspots

threatspots <- locations_grid %>%
    dplyr::select(CELLCODE,subspeciesname) %>%
    left_join(endemic_species, by=c("subspeciesname"="subspeciesname")) %>%
    st_drop_geometry() %>% # geometry must be dropped to run pivot wider
    dplyr::group_by(CELLCODE, paca) %>%
    dplyr::summarise(n_species = dplyr::n(), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = paca, values_from=n_species, values_fill=0)

threatspots_order <- locations_grid %>%
    dplyr::select(CELLCODE,subspeciesname) %>%
    left_join(endemic_species, by=c("subspeciesname"="subspeciesname")) %>%
    st_drop_geometry() %>% # geometry must be dropped to run pivot wider
    dplyr::group_by(CELLCODE, paca, Order) %>%
    dplyr::summarise(n_species = dplyr::n(), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = paca, values_from=n_species, values_fill=0)
# override the dataframe threatspots with a new one to include the 
# geometry of the grid.

threatspots <- locations_grid %>% 
    dplyr::select(CELLCODE) %>%
    distinct() %>%
    inner_join(threatspots, by=c("CELLCODE"="CELLCODE")) %>%
    mutate(LT.PT_pro=LT/PT) %>%
    dplyr::select(LT,PT, CELLCODE) %>%
    mutate(paca_threat=LT+PT) %>%
    filter(paca_threat>0)

threatspots_lt <- threatspots %>% 
    filter(paca_threat>= quantile(paca_threat,0.90))

st_write(threatspots,
         "../results/threatspots/threatspots.shp",
         append=F,
         delete_layer=T,
         delete_dsn = TRUE) 

## Per order
threatspots_order <- locations_grid %>% 
    dplyr::select(CELLCODE) %>%
    distinct() %>%
    inner_join(threatspots_order, by=c("CELLCODE"="CELLCODE")) %>%
    mutate(LT.PT_pro=LT/PT) %>%
    dplyr::select(LT,PT, CELLCODE, Order) %>%
    mutate(paca_threat=LT+PT) %>%
    filter(paca_threat>0)

threatspots_order_lt <- threatspots_order %>% 
    group_by(Order) %>%
    mutate(quant90 = quantile(paca_threat,0.90)) %>%
    filter(paca_threat>=quant90)

## Crete with endemic hotspots
##

g_e <- g_base +
    geom_sf(endemic_hotspots, mapping=aes(fill=n_species), alpha=0.3, size=0.1, na.rm = FALSE)+
    ggtitle("Endemic hotspots") +
    scale_fill_gradient(low = "yellow", high = "red", na.value = NA)

ggsave("../plots/crete-hotspots.png", plot=g_e, device="png")

g_e_order <- g_base +
    geom_sf(endemic_hotspots_order, mapping=aes(fill=n_species), alpha=0.3, size=0.1, na.rm = FALSE) +
    ggtitle("Endemic hotspots")+
    scale_fill_gradient(low = "yellow", high = "red", na.value = NA)+
    facet_wrap(vars(Order), ncol=4, scales = "fixed")

ggsave("../plots/crete-hotspots_order.png", 
       plot=g_e_order, 
       height = 20, 
       width = 50, 
       units="cm",
       device="png")

## Crete with threatspots
##
g_t <- g_base +
    geom_sf(threatspots_lt, mapping=aes(fill=paca_threat), alpha=0.3, size=0.1, na.rm = TRUE) +
    ggtitle("threatspots")+
    scale_fill_gradient(low = "yellow", high = "red", na.value = "transparent")+
    theme_bw()

ggsave("../plots/crete-threatspots.png", plot=g_t, device="png")


g_t_order <- g_base +
    geom_sf(threatspots_order_lt, mapping=aes(fill=paca_threat), alpha=0.3, size=0.1, na.rm = TRUE) +
    ggtitle("threatspots")+
    scale_fill_gradient(low = "yellow", high = "red", na.value = "transparent")+
    facet_wrap(vars(Order), ncol=4, scales = "fixed")

ggsave("../plots/crete-threatspots_order.png", 
       plot=g_t_order, 
       height = 20, 
       width = 50, 
       units="cm",
       device="png")

## Overlap

### Hotspots and threatspots
intersection_spots <- endemic_hotspots %>%
    st_drop_geometry() %>%
    inner_join(.,threatspots_lt, by=c("CELLCODE"="CELLCODE")) %>%
    st_as_sf()

g_e_t <- g_base +
    geom_sf(intersection_spots, mapping=aes(fill="red"), alpha=0.3, size=0.1, na.rm = TRUE) +
    ggtitle("Endemic hotspots and threatspots")+
#    scale_fill_gradient(low = "yellow", high = "red", na.value = "transparent")+
    theme_bw()

ggsave("../plots/crete-hotspots-threatspots.png", plot=g_e_t, device="png")

### With natura

endemic_hotspots_natura <- st_intersection(endemic_hotspots, natura_crete_land_sci)
sum(units::set_units(st_area(endemic_hotspots),km^2))
sum(units::set_units(st_area(endemic_hotspots_natura),km^2))

threatspots_natura <- st_intersection(threatspots_lt, natura_crete_land_sci)
sum(units::set_units(st_area(threatspots_lt),km^2))
sum(units::set_units(st_area(threatspots_natura),km^2))


### The threatened species that have AOO < 10% overlap with Natura2000.

species_10_natura <- endemic_species %>%
    filter(aoo_natura_percent<0.1 & threatened==T)

species_10_natura_l <- locations_grid %>%
    filter(subspeciesname %in% species_10_natura$subspeciesname) %>%
    group_by(CELLCODE) %>%
    summarise(n_species=n()) %>%
    filter(n_species>2)

species_10_natura_l_o <- locations_grid %>%
    filter(subspeciesname %in% species_10_natura$subspeciesname) %>%
    group_by(CELLCODE, Order) %>%
    summarise(n_species=n(), .groups="drop")

g_n_e <- g_base +
    geom_sf(species_10_natura_l, mapping=aes(fill=n_species), alpha=0.3, size=0.1, na.rm = FALSE)+
    ggtitle("Hotspots of AOO<10% overlap with Natura2000") +
    scale_fill_gradient(low = "yellow", high = "red", na.value = NA)

ggsave("../plots/hotspots_10_overlap_natura.png", plot=g_n_e, device="png")

g_e_order <- g_base +
    geom_sf(species_10_natura_l_o, mapping=aes(fill=n_species), alpha=0.3, size=0.1, na.rm = FALSE)+
    ggtitle("Hotspots of AOO<10% overlap with Natura2000") +
    scale_fill_gradient(low = "yellow", high = "red", na.value = NA)+
    facet_wrap(vars(Order), ncol=3, scales = "fixed")

ggsave("../plots/hotspots_10_overlap_natura_order.png", 
       plot=g_e_order, 
       height = 20, 
       width = 50, 
       units="cm",
       device="png")


