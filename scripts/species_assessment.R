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
arthropods_kriti_endemic <- readxl::read_excel("../data/endemic-arthropods-crete.xlsx") %>% 
    filter(Order!="Opiliones") %>% 
    dplyr::select(-Ergasia)

species_wo_coords <- arthropods_kriti_endemic[which(is.na(arthropods_kriti_endemic$latD)),]

## Sanity check on coordinates
print(paste0("there are ", nrow(species_wo_coords)," rows without coordinates"))

arthropods_occurrences <- arthropods_kriti_endemic %>% 
    dplyr::select(subspeciesname,latD, logD, Order, families) %>% 
    filter(latD<38, logD<30) %>%
    na.omit()

arthropods_occurrences |> dplyr::select(-c(families,Order)) |> write_delim("../results/arthropods_occurrences.tsv",delim="\t")

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
## Source of the occurrences, NHMC collection and Bibliography

locations_source <- readxl::read_excel("../data/Data-ENDEMICS.xlsx", 
                                     col_types = c("text",
                                                   "text",
                                                   "text",
                                                   "text",
                                                   "text",
                                                   "numeric",
                                                   "numeric",
                                                   "date",
                                                   "text")) %>% 
    filter(Order!="Opiliones") %>%
    mutate(source=if_else(is.na(Ergasia),"Bibliography","NHMC")) %>% 
    distinct(source, latD, logD) %>%
    filter(latD<38, logD<30) %>%
    na.omit()

write_delim(locations_source, "../data/locations_source.tsv", delim="\t", col_names=T)

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
                                if_else(potentially_VU==TRUE,"PT","PNT")))) %>%
    mutate(iucn=if_else(potentially_CR==TRUE,"CR",
                        if_else(potentially_EN==TRUE,"EN",
                                if_else(potentially_VU==TRUE,"VU","NT/LC")))) %>%
    mutate(threatened= if_else(paca %in% c("PT", "LT"),TRUE, FALSE)) %>%
    mutate(aoo_natura_percent=round(aoo_natura/aoo, digits=5))

write_delim(endemic_species, "../results/endemic_species_assessment.tsv", delim="\t") 

