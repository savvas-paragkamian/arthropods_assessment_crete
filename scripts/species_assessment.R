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
#st_write(locations_shp, "../data/arthropods_occurrences/arthropods_occurrences.shp", layer_options = "ENCODING=UTF-8" )
## Crete Spatial data

crete_shp <- sf::st_read("../data/crete/crete.shp")

### Natura2000

natura_crete <- sf::st_read("../data/natura2000/natura2000_crete.shp")

natura_crete_land <- st_intersection(natura_crete, crete_shp)

# split the SPA SCI

natura_crete_land_sci <- natura_crete_land %>% filter(SITETYPE=="B")

### World Database of Protected areas

wdpa_crete <- sf::st_read("../data/wdpa_crete/wdpa_crete.shp")

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

grid_10km <- sf::st_read(dsn="../data/Greece_shapefile/gr_10km.shp") %>%
    st_transform(., crs="WGS84")

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

eoo_results_df <- convert_ll_df(eoo_results) %>% as_tibble() 
colnames(eoo_results_df) <- c("subspeciesname", "n_sites", "eoo", "land_eoo")

eoo_results_df$n_sites <- as.numeric(eoo_results_df$n_sites)
eoo_results_df$eoo <- as.numeric(eoo_results_df$eoo)
eoo_results_df$land_eoo <- as.numeric(eoo_results_df$land_eoo)

write_delim(eoo_results_df, "../results/eoo_resuls.tsv", delim ="\t")

#eoo_results_df <- read_delim("../results/eoo_resuls.tsv", delim="\t")

### Natura overlap with eoo of species
eoo_natura <- eoo_calculation(locations_inland, crete_shp, natura_crete_land_sci, FALSE, "natura")

eoo_natura_df <- convert_ll_df(eoo_natura) %>% as_tibble()

colnames(eoo_natura_df) <- c("subspeciesname", "n_sites", "eoo", "natura_eoo")

eoo_natura_df$n_sites <- as.numeric(eoo_natura_df$n_sites)
eoo_natura_df$eoo <- as.numeric(eoo_natura_df$eoo)
eoo_natura_df$natura_eoo <- as.numeric(eoo_natura_df$natura_eoo)
write_delim(eoo_natura_df, "../results/eoo_natura.tsv", delim="\t")

#eoo_natura_df <- read_delim("../results/eoo_natura.tsv", delim="\t")
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
aoo_overlap_natura_sci <- read_delim( "../results/aoo_endemic.tsv", delim="\t")
## Preliminary Automated Conservation Assessments (PACA)

endemic_species <- locations_grid %>% 
    st_drop_geometry() %>%
    distinct(subspeciesname,Order,n_locations) %>% 
    ungroup() %>%
    left_join(eoo_natura_df, by=c("subspeciesname"="subspeciesname")) %>%
    left_join(aoo_overlap_natura_sci, by=c("subspeciesname"="species")) %>%
    mutate(potentially_VU=if_else(n_locations<=10 & (eoo<20000 | aoo_area<2000),TRUE,FALSE)) %>%
    mutate(potentially_EN=if_else(n_locations<=5 & (eoo<5000 | aoo_area<500),TRUE,FALSE)) %>%
    mutate(potentially_CR=if_else(n_locations==1 & (if_else(is.na(eoo),0,eoo)<100 | aoo_area<10),TRUE,FALSE)) %>%
    ungroup() %>%
    mutate(paca=if_else(potentially_CR==TRUE,"LT",
                        if_else(potentially_EN==TRUE,"LT",
                                if_else(potentially_VU==TRUE,"PT","FALSE")))) %>%
    mutate(threatened= if_else(paca=="FALSE",FALSE, TRUE))

write_delim(endemic_species, "../results/endemic_species_paca.tsv", delim="\t") 

## Endemic hotspots

endemic_hotspots <- locations_grid %>%
    group_by(CELLCODE) %>%
    summarise(n_species=n()) %>%
    filter(n_species >= quantile(n_species, 0.90))

st_write(endemic_hotspots, "../results/endemic_hotspots/endemic_hotspots.shp") 

endemic_hotspots_order <- locations_grid %>% 
    group_by(CELLCODE, Order) %>%
    summarise(n_species=n(), .groups="drop") %>%
    group_by(Order) %>%
    mutate(quant90= quantile(n_species, 0.90)) %>%
    filter(n_species >= quant90)

## Threadspots

threadspots <- locations_grid %>%
    dplyr::select(CELLCODE,subspeciesname) %>%
    left_join(endemic_species, by=c("subspeciesname"="subspeciesname")) %>%
    st_drop_geometry() %>% # geometry must be dropped to run pivot wider
    dplyr::group_by(CELLCODE, paca) %>%
    dplyr::summarise(n_species = dplyr::n(), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = paca, values_from=n_species, values_fill=0)

threadspots_order <- locations_grid %>%
    dplyr::select(CELLCODE,subspeciesname) %>%
    left_join(endemic_species, by=c("subspeciesname"="subspeciesname")) %>%
    st_drop_geometry() %>% # geometry must be dropped to run pivot wider
    dplyr::group_by(CELLCODE, paca, Order) %>%
    dplyr::summarise(n_species = dplyr::n(), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = paca, values_from=n_species, values_fill=0)
# override the dataframe threadspots with a new one to include the 
# geometry of the grid.

threadspots <- locations_grid %>% 
    dplyr::select(CELLCODE) %>%
    distinct() %>%
    inner_join(threadspots, by=c("CELLCODE"="CELLCODE")) %>%
    mutate(LT.PT_pro=LT/PT) %>%
    dplyr::select(LT,PT, CELLCODE) %>%
    mutate(paca_threat=LT+PT) %>%
    filter(paca_threat>0)

threadspots_lt <- threadspots %>% 
    filter(paca_threat>= quantile(paca_threat,0.90))

st_write(threadspots, "../results/threadspots/threadspots.shp") 

## Per order
threadspots_order <- locations_grid %>% 
    dplyr::select(CELLCODE) %>%
    distinct() %>%
    inner_join(threadspots_order, by=c("CELLCODE"="CELLCODE")) %>%
    mutate(LT.PT_pro=LT/PT) %>%
    dplyr::select(LT,PT, CELLCODE, Order) %>%
    mutate(paca_threat=LT+PT) %>%
    filter(paca_threat>0)

threadspots_order_lt <- threadspots_order %>% 
    group_by(Order) %>%
    mutate(quant90 = quantile(paca_threat,0.90)) %>%
    filter(paca_threat>=quant90)

## Crete with endemic hotspots
##

g_e <- g_base +
    geom_sf(endemic_hotspots, mapping=aes(fill=n_species), alpha=0.3, size=0.1, na.rm = FALSE) +
    ggtitle("Endemic hotspots")+
    scale_fill_gradient(low = "yellow", high = "red", na.value = NA)

ggsave("../plots/crete-hotspots.png", plot=g_e, device="png")

g_e_order <- g_base +
    geom_sf(endemic_hotspots_order, mapping=aes(fill=n_species), alpha=0.3, size=0.1, na.rm = FALSE) +
    ggtitle("Endemic hotspots")+
    scale_fill_gradient(low = "yellow", high = "red", na.value = NA)+
    facet_grid(rows=vars(Order))

ggsave("../plots/crete-hotspots_order.png", plot=g_e_order, device="png")

## Crete with threadspots
##
g_t <- g_base +
    geom_sf(threadspots_lt, mapping=aes(fill=paca_threat), alpha=0.3, size=0.1, na.rm = TRUE) +
    ggtitle("Threadspots")+
    scale_fill_gradient(low = "yellow", high = "red", na.value = "transparent")+
    theme_bw()

ggsave("../plots/crete-threadspots.png", plot=g_t, device="png")


g_t_order <- g_base +
    geom_sf(threadspots_order_lt, mapping=aes(fill=paca_threat), alpha=0.3, size=0.1, na.rm = TRUE) +
    ggtitle("Threadspots")+
    scale_fill_gradient(low = "yellow", high = "red", na.value = "transparent")+
    facet_grid(rows=vars(Order))

ggsave("../plots/crete-threadspots_order.png", plot=g_t_order, device="png")

## Overlap

### Hotspots and threadspots
intersection_spots <- endemic_hotspots %>%
    st_drop_geometry() %>%
    inner_join(.,threadspots_lt, by=c("CELLCODE"="CELLCODE")) %>%
    st_as_sf()

g_e_t <- g_base +
    geom_sf(intersection_spots, mapping=aes(fill="red"), alpha=0.3, size=0.1, na.rm = TRUE) +
    ggtitle("Endemic hotspots and Threadspots")+
#    scale_fill_gradient(low = "yellow", high = "red", na.value = "transparent")+
    theme_bw()

ggsave("../plots/crete-hotspots-threadspots.png", plot=g_e_t, device="png")

### With natura

endemic_hotspots_natura <- st_intersection(endemic_hotspots, natura_crete_land_sci)
sum(units::set_units(st_area(endemic_hotspots),km^2))
sum(units::set_units(st_area(endemic_hotspots_natura),km^2))

threadspots_natura <- st_intersection(threadspots_lt, natura_crete_land_sci)
sum(units::set_units(st_area(threadspots_lt),km^2))
sum(units::set_units(st_area(threadspots_natura),km^2))


