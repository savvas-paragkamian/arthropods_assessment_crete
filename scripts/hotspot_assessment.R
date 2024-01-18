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
library(ggnewscale)
library(sf)
library(terra)
library(quadtree)
library(WEGE) 
library(units)
source("functions.R")

# Load data
## Base
g_base <- g_base()
endemic_species <- read_delim("../results/endemic_species_assessment.tsv", delim="\t") 
crete_shp <- sf::st_read("../data/crete/crete.shp")

## Locations
locations_shp <- sf::st_read("../data/arthropods_occurrences/arthropods_occurrences.shp")
locations_shp <- locations_shp |>
    mutate(long = unlist(map(locations_shp$geometry,1)),
           lat = unlist(map(locations_shp$geometry,2)))
locations_grid <- st_read("../results/locations_grid/locations_grid.shp")
locations_spatial <- st_read("../results/locations_spatial/locations_spatial.shp")

## grid 1km
grid_1km_o <- st_read("../data/Greece_shapefile/gr_1km.shp")
crete_eea <- crete_shp |> st_transform(crs(grid_1km_o))
crete_1km <- st_read("../data/Greece_shapefile/gr_1km.shp") |>
    st_transform(., crs="WGS84") |>
    st_join(crete_shp, left=F)
#crete_bbox_polygon <- st_as_sf(st_as_sfc(st_bbox(crete_eea)))
#crete_grid1 <- st_crop(grid_1km, crete_bbox_polygon)

## grid 10km
grid_10km <- sf::st_read(dsn="../data/Greece_shapefile/gr_10km.shp") 
grid_10km_w <- grid_10km |> st_transform(., crs="WGS84")
crete_grid10 <- st_join(grid_10km_w, crete_shp, left=F)

### locations modifications
locations_inland <- st_join(locations_shp, crete_shp, left=F)
locations_inland <- locations_inland |> distinct(sbspcsn,long, lat, geometry)
locations_eea <- locations_inland |> st_transform(crs(grid_1km_o)) 
### 1km over shp

locations_1_grid <- st_join(crete_1km, locations_inland, left=F) |>
    dplyr::distinct(CELLCODE, sbspcsn, geometry) |>
    group_by(sbspcsn) |>
    summarise(geometry=st_union(geometry)) |>
    mutate(area=st_area(geometry))

### 1km over raster
crete_manual_grid <- rast(ncol=311,
          nrow=152,
          xmin=5523000,
          xmax=5834000,
          ymin=1400000,
          ymax=1552000,
          crs=crs(grid_1km_o))

values(crete_manual_grid) <- 0
crete_manual_grid_w <- terra::project(crete_manual_grid, "WGS84")
values(crete_manual_grid_w) <- 0


locations_eea <- locations_eea |> distinct(sbspcsn,long, lat, geometry) 
grid_occ_count <- rasterize(locations_eea, crete_manual_grid, fun='count')
#locations_inland$grid_occ <- extract(crete_manual_grid_w, locations_inland,cellnumbers=F, ID=F)

raster_df <- terra::as.data.frame(grid_occ_count, xy=TRUE, cells=TRUE)


# Adaptive resolution with quadtree
qt <- quadtree(grid_occ_count,
               min_cell_length=4000,
               max_cell_length=10000,
               split_threshold=10,
               split_if_any_na=T,
               split_method = "range")

qt_rast <- as_raster(qt)
qt_rast_df <- terra::as.data.frame(qt_rast, xy=TRUE, cells=TRUE) 

### plots
png(file="../plots/quadtree_crete.png",
    width = 40,
    height = 20,
    res=300,
    units = "cm",
    bg="white")
plot(qt,
     border_lwd = .01,
     xlim = c(5523000, 5834000),
     ylim = c(1420000,1548000),
     main = "expand raster dimensions")
plot(crete_eea, add=T, color="")
#plot(grid_occ_count, add=T)
plot(crete_grid10,border_lwd = .02, add=T, color="")
dev.off()


crete_sampling_1_grid <- ggplot() +
    geom_sf(crete_eea, mapping=aes()) +
    geom_tile(qt_rast_df, mapping=aes(x=x, y=y, fill=lyr.1)) + 
    scale_fill_gradientn(colours = c("black", "grey80")) +
    labs(fill = "Quadtrees")+
    new_scale_fill() +
    geom_tile(raster_df,mapping=aes(x=x, y=y, fill=count)) +
    scale_fill_gradientn(colours = c("purple", "orange")) +
    labs(fill = "1km grid",y = "Group") +
    geom_sf(locations_eea, mapping=aes(),color="red", size=0.01)+
#    geom_sf(locations_1_grid, mapping=aes(fill=n_occurrences), alpha=0.3,size=0.02)+
    ggtitle("Adaptive resolution")+
    theme_bw()

ggsave("../plots/crete_sampling_1_grid.png",
       plot=crete_sampling_1_grid,
       width=30,
       height=30,
       dpi=300,
       unit="cm",
       device="png")

# WEGE
## first create the input sf object of the species

# WEGE with EOO
#occurrences <- locations_inland |>
#    rename(subspeciesname=sbspcsn)

#species <- occurrences %>% st_drop_geometry() %>%
#    dplyr::select(subspeciesname) %>% 
#    distinct() %>%
#    pull()
#
### range as EOO, extend of occurrence
#range_list <- list()
#for(s in seq_along(species)){
#
#    n_occurrences <- occurrences %>%
#        filter(subspeciesname==species[s])
#
#    if (nrow(n_occurrences) > 2) {
#
#        species_convex <- st_convex_hull(st_union(n_occurrences))
#        range_list[[s]] <- st_sf(geom=species_convex, ID=species[s])
#    }
#}
#
#input <- do.call(rbind, range_list)

occurrences <- locations_1_grid |>
    rename(subspeciesname=sbspcsn)|>
    mutate(ID=subspeciesname)

input <- occurrences |>
    mutate(binomial=ID) |>
    left_join(endemic_species, by=c("ID"="subspeciesname")) |>
    mutate(category=gsub("NT/","",iucn))

## for each cellgrid calculate the WEGE

grid_for_wege <- crete_grid10

cellcodes <- grid_for_wege[[1]]
wege_l <- list()
for (c in seq_along(cellcodes)){

    target_area <- crete_grid10 |> filter(CELLCODE==cellcodes[c])
    temp <- get_wege(target_area,input,species = 'binomial',category = 'category')
    if (!is.null(temp)){
        wege_l[[c]] <- st_sf(st_as_sfc(target_area), wege=temp, CELLCODE=c)
    }
}

wege_results <- do.call(rbind, wege_l)

wege_results_f <- wege_results |> filter(wege>3)

wege_map <- g_base +
    geom_sf(wege_results_f, mapping=aes(fill=wege)) + 
    scale_fill_gradient(low = "yellow", high = "red", na.value = NA) +
    theme_bw()

ggsave("../plots/crete_wege_hotspots_map.png", plot=wege_map, device="png")


# Endemic hotspots

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


