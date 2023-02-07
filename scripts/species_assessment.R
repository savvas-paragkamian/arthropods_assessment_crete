#!/usr/bin/Rscript

## Script name: species_assessment.R
##
## Purpose of script: IUCN criterion B and PACA 
## assessments of the cretan arthropoda occurrences as compiled from 
## Giannis Bolanakis and Apostolos Trichas from NHMC
##
## How to run:
## Rscript species_assessment.R
##
## Execution time: 22 minutes 
## for 7690 occurrences of 343 species
##
## Author: Savvas Paragkamian
##
## Date Created: 2022-12-22

library(tidyverse)
library(readxl)
library(sf)
library(units)
library(ConR)

source("functions.R")

# Load Data
## Endemic species occurrences in Crete
arthropods_kriti_endemic <- readxl::read_excel("../data/Data-ENDEMICS.xlsx") %>% 
    filter(Order!="Opiliones") %>% 
    dplyr::select(-Ergasia)

species_wo_coords <- arthropods_kriti_endemic[which(is.na(arthropods_kriti_endemic$latD)),]

## Sanity check on coordinates
print(paste0("there are ", nrow(species_wo_coords)," rows without coordinates"))

arthropods_occurrences <- arthropods_kriti_endemic %>% 
    dplyr::select(subspeciesname,latD, logD, Order, families) %>% 
    filter(latD<38, logD<30) %>%
    na.omit()

out_of_crete <- arthropods_kriti_endemic[which(!(arthropods_kriti_endemic$subspeciesname %in% arthropods_occurrences$subspeciesname)),]

locations_shp <- st_as_sf(arthropods_occurrences,
                          coords=c("logD", "latD"),
                          remove=TRUE,
                          crs="WGS84")

st_write(locations_shp, "../data/arthropods_occurrences/arthropods_occurrences.shp", layer_options = "ENCODING=UTF-8" , append=FALSE, delete_layer=T)

st_write(locations_shp, "../data/arthropods_occurrences/arthropods_occurrences.csv",
         layer_options = "ENCODING=UTF-8", 
         append=FALSE, 
         delete_layer=T, 
         delete_dsn = TRUE)
## Crete Spatial data

crete_shp <- sf::st_read("../data/crete/crete.shp")

### Natura2000

natura_crete <- sf::st_read("../data/natura2000/natura2000_crete.shp")

natura_crete_land <- st_intersection(natura_crete, crete_shp)

# split the SPA SCI

natura_crete_land_sci <- natura_crete_land %>% filter(SITETYPE=="B")

### World Database of Protected areas

wdpa_crete <- sf::st_read("../data/wdpa_crete/wdpa_crete.shp")

wdpa_crete_wildlife <- wdpa_crete %>% filter(DESIG_ENG=="Wildlife Refugee")

## here we calculate the area of each polygon of the region of Crete
## and subsequently we keep only the largest polygon e.g. Crete island only.

crete_polygon <- st_cast(crete_shp, "POLYGON") %>% 
    mutate(area_m2=as.numeric(st_area(.))) %>%
    filter(area_m2==max(area_m2))

# There is a difference in the points of the following objects
# some occurrences fall in the sea and others are on the peripheral
# islands of Crete.

locations_inland <- st_join(locations_shp, crete_shp, left=F)
locations_out <- st_difference(locations_shp, crete_shp)

print(paste0("these occurrences are in the sea", nrow(locations_out)))
# return the coordinates to a dataframe format

locations_inland_df <- locations_inland %>%
    mutate(ddlon = unlist(map(locations_inland$geometry,1)),
           ddlat = unlist(map(locations_inland$geometry,2))) %>%
    st_drop_geometry() %>% 
    dplyr::rename(tax=subspeciesname) %>%
    dplyr::select(ddlat, ddlon, tax)

## EEA reference grid
## https://www.eea.europa.eu/data-and-maps/data/eea-reference-grids-2

grid_10km <- sf::st_read(dsn="../data/Greece_shapefile/gr_10km.shp") %>%
    st_transform(., crs="WGS84")

## make a new grid using the st_make_grid
gr <- st_make_grid(crete_polygon,cellsize = 0.1)   
### keep grid that overlaps with Crete
crete_grid10 <- st_join(grid_10km, crete_shp, left=F)
#crete_iucn_grid10m <- st_join(grid_iucn, crete_polygon, left=F)

### Sites to locations. Each 10km grid is a location of a species
locations_grid <- st_join(crete_grid10, locations_inland, left=TRUE) %>%
    dplyr::select(CELLCODE, subspeciesname, Order, families) %>% 
    distinct() %>%
    na.omit() %>%
    group_by(subspeciesname) %>%
    mutate(n_locations=n()) %>% 
    ungroup()

st_write(locations_grid,
         "../results/locations_grid/locations_grid.shp", 
         append=FALSE, 
         delete_layer=T, 
         delete_dsn = TRUE)

## Here is Crete with all the sampling points
##
g_base <- g_base()

g <- g_base +
    geom_sf(locations_inland, mapping=aes(),color="blue", size=0.1, alpha=0.2) +
    geom_sf(locations_out, mapping=aes(),color="red", size=0.1, alpha=0.2) +
#    geom_sf(crete_grid10m, mapping=aes(),color="red", alpha=0.2, size=0.1) +
#    geom_sf(gr, mapping=aes(),color="orange", alpha=0.2, size=0.1) +
    coord_sf(crs="WGS84") +
    theme_bw()

ggsave("../plots/crete-occurrences.png", plot=g, device="png")

# Processing

## EOO
# eoo_calculation is a custom function that takes 3 inputs.
# the location shapefile, the polygon and a TRUE/FALSE value to 
# export or not to plot
eoo_results <- eoo_calculation(locations_inland, crete_shp,crete_shp,FALSE, "EOO")

write_delim(eoo_results, "../results/eoo_resuls.tsv", delim ="\t")

#eoo_results_df <- read_delim("../results/eoo_resuls.tsv", delim="\t")

### Natura overlap with eoo of species
eoo_natura <- eoo_calculation(locations_inland, crete_shp, natura_crete_land_sci, FALSE, "natura")
eoo_natura <- eoo_natura %>%
    dplyr::select(-c(n_sites,eoo)) %>%
    rename("eoo_natura"="eoo_overlap")
write_delim(eoo_natura, "../results/eoo_natura.tsv", delim="\t")

### Wildlife refugees overlap with EOO

eoo_wildlife <- eoo_calculation(locations_inland, crete_shp, wdpa_crete_wildlife, FALSE, "wildlife")

eoo_wildlife <- eoo_wildlife %>%
    dplyr::select(-c(n_sites, eoo)) %>%
    rename("eoo_wildlife"="eoo_overlap")
write_delim(eoo_wildlife, "../results/eoo_wildlife.tsv", delim="\t")
#eoo_wildlife_df <- read_delim("../results/eoo_wildlife.tsv", delim="\t")
## AOO
## see for bootstrap
AOO_endemic <- AOO.computing(locations_inland_df,Cell_size_AOO = 2, nbe.rep.rast.AOO=1000, export_shp = T)
# in AOO.computing when export_shp is T a list is returned.

if (is.list(AOO_endemic)) {
    AOO_endemic_df <- data.frame(AOO= AOO_endemic[[1]]) %>% rownames_to_column(var="subspeciesname")
} else {
    AOO_endemic_df <- tibble(subspeciesname=names(AOO_endemic),AOO=AOO_endemic, row.names=NULL)
}

# This function returns a tibble with the calculations as well a figure for each species.
#
aoo_overlap_natura_sci <- aoo_overlap(AOO_endemic,crete_shp, natura_crete_land_sci, T)

aoo_overlap_natura_sci <- aoo_overlap_natura_sci %>%
    rename("aoo_natura"="aoo_overlap" )

write_delim(aoo_overlap_natura_sci, "../results/aoo_endemic.tsv", delim="\t")

#aoo_overlap_natura_sci <- read_delim( "../results/aoo_endemic.tsv", delim="\t")
### AOO overlap with wildlife refugees
aoo_overlap_wildlife <- aoo_overlap(AOO_endemic,crete_shp, wdpa_crete_wildlife, T)

aoo_overlap_wildlife <- aoo_overlap_wildlife %>% 
    dplyr::select(-aoo) %>% 
    rename("aoo_wildlife"="aoo_overlap")

write_delim(aoo_overlap_wildlife, "../results/aoo_overlap_natura_sci.tsv", delim="\t")

## Preliminary Automated Conservation Assessments (PACA)

endemic_species <- locations_grid %>% 
    st_drop_geometry() %>%
    distinct(subspeciesname,Order,families,n_locations) %>% 
    ungroup() %>%
    left_join(eoo_results, by=c("subspeciesname"="subspeciesname")) %>%
    left_join(eoo_natura, by=c("subspeciesname"="subspeciesname")) %>%
    left_join(eoo_wildlife, by=c("subspeciesname"="subspeciesname")) %>%
    left_join(aoo_overlap_natura_sci, by=c("subspeciesname"="species")) %>%
    left_join(aoo_overlap_wildlife, by=c("subspeciesname"="species")) %>%
    mutate(potentially_VU=if_else(n_locations<=10 & (eoo<20000 | aoo<2000),TRUE,FALSE)) %>%
    mutate(potentially_EN=if_else(n_locations<=5 & (eoo<5000 | aoo<500),TRUE,FALSE)) %>%
    mutate(potentially_CR=if_else(n_locations==1 & (if_else(is.na(eoo),0,eoo)<100 | aoo<10),TRUE,FALSE)) %>%
    ungroup() %>%
    mutate(paca=if_else(potentially_CR==TRUE,"LT",
                        if_else(potentially_EN==TRUE,"LT",
                                if_else(potentially_VU==TRUE,"PT","LNT")))) %>%
    mutate(iucn=if_else(potentially_CR==TRUE,"CR",
                        if_else(potentially_EN==TRUE,"EN",
                                if_else(potentially_VU==TRUE,"VU","NT/LC")))) %>%
    mutate(threatened= if_else(paca %in% c("PT", "LT"),TRUE, FALSE)) %>%
    mutate(aoo_natura_percent=round(aoo_natura/aoo, digits=5))

write_delim(endemic_species, "../results/endemic_species_assessment.tsv", delim="\t") 

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


