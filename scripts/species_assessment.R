#!/usr/bin/Rscript

library(tidyverse)
library(readxl)
library(rredlist)
library(taxize)
library(sf)
library(units)
#library(rgdal)
library(ConR)
library(vegan)

source("species_assessment_functions.R")

# Load Data
## Endemic species occurrences in Crete
arthropods_kriti_endemic <- readxl::read_excel("../data/Data-ENDEMICS.xlsx") %>% 
    filter(Order!="Opiliones") %>% 
    dplyr::select(-Ergasia)

species_wo_coords <- arthropods_kriti_endemic[which(is.na(arthropods_kriti_endemic$latD)),]

## Sanity check on coordinates
print(paste0("there are ", nrow(species_wo_coords)," rows without coordinates"))

arthropods_occurrences <- arthropods_kriti_endemic %>% 
    dplyr::select(subspeciesname,latD, logD, Order) %>% 
    filter(latD<38, logD<30) %>%
    na.omit()

out_of_crete <- arthropods_kriti_endemic[which(!(arthropods_kriti_endemic$subspeciesname %in% arthropods_occurrences$subspeciesname)),]

locations_shp <- st_as_sf(arthropods_occurrences,
                          coords=c("logD", "latD"),
                          remove=TRUE,
                          crs="WGS84")

## Crete Spatial data

# not run
#periphereies_shp <- sf::st_read("~/Documents/spatial_data/periphereies/periphereies.shp")

#### crete is the 12th region in this shapefile
#### https://geodata.gov.gr/en/dataset/28121643-d977-48eb-a8ca-a6fac6b4af6d/resource/7c80a2c1-93b7-4814-9fc4-245e775acaa6/download/periphereies.zip
#
#crete_shp <- periphereies_shp[12,] %>% st_transform(., "WGS84")
#
#st_write(crete_shp,"../data/crete/crete.shp",
#         layer_options = "ENCODING=UTF-8", delete_layer = TRUE)

crete_shp <- sf::st_read("../data/crete/crete.shp")

### Natura2000

# natura2000 shapefile downloaded from here https://www.eea.europa.eu/data-and-maps/data/natura-13 
#natura2000 <- sf::st_read("~/Downloads/Natura2000_end2021_Shapefile/Natura2000_end2021_epsg3035.shp") 
#                   %>% st_transform(., "WGS84")
#
#natura_crete <- st_crop(natura2000,
#                        y=st_bbox(c(xmin=23,ymin=34,xmax=27,ymax=36),
#                                  crs="WGS84"))
#
#st_write(natura_crete,"../data/natura2000/natura2000_crete.shp",
#         layer_options = "ENCODING=UTF-8", delete_layer = TRUE)

natura_crete <- sf::st_read("../data/natura2000/natura2000_crete.shp")

natura_crete_land <- st_intersection(natura_crete, crete_shp)

# split the SPA SCI

natura_crete_land_sci <- natura_crete_land %>% filter(SITETYPE=="B")

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

print(paste0("these occurrences are in the sea", locations_out))
# return the coordinates to a dataframe format

locations_inland_df <- locations_inland %>%
    mutate(ddlon = unlist(map(locations_inland$geometry,1)),
           ddlat = unlist(map(locations_inland$geometry,2))) %>%
    st_drop_geometry() %>% 
    dplyr::rename(tax=subspeciesname) %>%
    dplyr::select(ddlat, ddlon, tax)

## EEA reference grid
## https://www.eea.europa.eu/data-and-maps/data/eea-reference-grids-2
##

grid_10km <- sf::st_read(dsn="../data/Greece_shapefile/gr_10km.shp") %>%
    st_transform(., crs="WGS84")

#grid_iucn <- sf::st_read(dsn="~/Downloads/AOOGrid_10x10kmshp/AOOGrid_10x10km.shp") %>% 
#    st_transform(., crs="WGS84")

## make a new grid using the st_make_grid
gr <- st_make_grid(crete_polygon,cellsize = 0.1)   
### keep grid that overlaps with Crete
crete_grid10 <- st_join(grid_10km, crete_shp, left=F)
#crete_iucn_grid10m <- st_join(grid_iucn, crete_polygon, left=F)

### Sites to locations. Each 10km grid is a location of a species
locations_grid <- st_join(crete_grid10, locations_inland, left=TRUE) %>%
    dplyr::select(CELLCODE, subspeciesname, Order) %>% 
    distinct() %>%
    na.omit() %>%
    group_by(subspeciesname) %>%
    mutate(n_locations=n()) %>% 
    ungroup()

## Here is Crete with all the sampling points
##
g <- ggplot() +
    geom_sf(crete_shp, mapping=aes()) +
    geom_sf(locations_inland, mapping=aes(),color="blue", size=0.1, alpha=0.2) +
    geom_sf(locations_out, mapping=aes(),color="red", size=0.1, alpha=0.2) +
#    geom_sf(crete_grid10m, mapping=aes(),color="red", alpha=0.2, size=0.1) +
    geom_sf(natura_crete_land, mapping=aes(fill=SITETYPE, color=SITETYPE), alpha=0.2, size=0.1) +
#    geom_sf(gr, mapping=aes(),color="orange", alpha=0.2, size=0.1) +
    coord_sf(crs="WGS84") +
    theme_bw()

ggsave("../plots/crete-occurrences.png", plot=g, device="png")

# Processing

## EOO
# eoo_calculation is a custom function that takes 3 inputs.
# the location shapefile, the polygon and a TRUE/FALSE value to 
# export or not to plot
eoo_results <- eoo_calculation(locations_inland, crete_shp,"nothing",FALSE, "EOO")

eoo_results_df <- convert_ll_df(eoo_results) %>% as_tibble() 
colnames(eoo_results_df) <- c("subspeciesname", "n_sites", "eoo", "land_eoo")

eoo_results_df$n_sites <- as.numeric(eoo_results_df$n_sites)
eoo_results_df$eoo <- as.numeric(eoo_results_df$eoo)
eoo_results_df$land_eoo <- as.numeric(eoo_results_df$land_eoo)

write.table(eoo_results_df, file="../results/eoo_resuls.tsv", sep="\t")

### Natura overlap with eoo of species
eoo_natura <- eoo_calculation(locations_inland, crete_shp, natura_crete_land_sci, FALSE, "natura")

eoo_natura_df <- convert_ll_df(eoo_natura) %>% as_tibble()

colnames(eoo_natura_df) <- c("subspeciesname", "n_sites", "eoo", "natura_eoo")

eoo_natura_df$n_sites <- as.numeric(eoo_natura_df$n_sites)
eoo_natura_df$eoo <- as.numeric(eoo_natura_df$eoo)
eoo_natura_df$natura_eoo <- as.numeric(eoo_natura_df$natura_eoo)
write.table(eoo_natura_df, file="../results/eoo_natura.tsv", sep="\t")

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

write_delim(aoo_overlap_natura_sci, "../results/aoo_endemic.tsv", delim="\t")
## Preliminary Automated Conservation Assessments (PACA)

endemic_species <- locations_grid %>% 
    st_drop_geometry() %>%
    distinct(subspeciesname,Order,n_locations) %>% 
    ungroup() %>%
    left_join(eoo_natura_df, by=c("subspeciesname"="subspeciesname")) %>%
    left_join(aoo_overlap_natura_sci, by=c("subspeciesname"="species")) %>%
    mutate(Potentially_VU=if_else(n_locations<=10 & (eoo<20000 | aoo_area<2000),TRUE,FALSE)) %>%
    mutate(Potentially_EN=if_else(n_locations<=5 & (eoo<5000 | aoo_area<500),TRUE,FALSE)) %>%
    mutate(Potentially_CR=if_else(n_locations==1 & (eoo<100 | aoo_area<10),TRUE,FALSE)) %>%
    ungroup() %>%
    mutate(potential_status=if_else(
                                    Potentially_CR=="TRUE","Potentially_CR", 
                                    if_else(
                                            Potentially_EN=="TRUE","Potentially_EN",
                                            if_else(
                                                    Potentially_VU=="TRUE","Potentially_VU","FALSE"))))

write_delim(endemic_species, "../results/endemic_species_paca.tsv", delim="\t") 

endemic_species_threatened <- endemic_species %>% filter(potential_status!="FALSE") 

write_delim(endemic_species_threatened, "../results/endemic_species_paca_threatened.tsv", delim="\t") 

endemic_species_cr <- endemic_species %>% filter(potential_status=="Potentially_CR") 

write_delim(endemic_species_cr, "../results/endemic_species_cr.tsv", delim="\t") 

## Endemic hotspots and Threadspots








