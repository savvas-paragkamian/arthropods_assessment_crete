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
locations_grid <- st_read("../results/locations_grid/locations_grid.shp") |>
    rename("CELLCODE"="CELLCOD") |>
    rename("subspeciesname"="sbspcsn") |>
    rename("families"="familis")
locations_spatial <- st_read("../results/locations_spatial/locations_spatial.shp")

## Natura2000

natura_crete <- sf::st_read("../data/natura2000/natura2000_crete.shp")
wdpa_crete <- sf::st_read("../data/wdpa_crete/wdpa_crete.shp")

wdpa_crete_wildlife <- wdpa_crete |> filter(DESIG_ENG=="Wildlife Refugee")
natura_crete_land <- st_intersection(natura_crete, crete_shp)

# split the SPA SCI

natura_crete_land_sci <- natura_crete_land |> filter(SITETYPE=="B")

########################## Grids 1, 4, 8, 10 sq.km  ########################
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

locations_10_grid_samples <- st_join(crete_grid10, locations_inland, left=F) |>
    distinct(geometry,CELLCODE, long, lat) |>
    group_by(geometry,CELLCODE) |>
    summarise(n_samples=n(),.groups="keep")

locations_10_grid <- st_join(crete_grid10, locations_inland, left=F) |>
    dplyr::distinct(CELLCODE, sbspcsn, geometry) |>
    group_by(CELLCODE, geometry) |>
    summarise(n_species=n(), .groups="keep") |>
    mutate(area=st_area(geometry)) |>
    left_join(st_drop_geometry(locations_10_grid_samples))

### pearson correlation of sampling effort and species
cor(locations_10_grid$n_species, locations_10_grid$n_samples)
### locations modifications
locations_inland <- st_join(locations_shp, crete_shp, left=F)
locations_inland <- locations_inland |> distinct(sbspcsn,long, lat, geometry)
sampling_intensity <- locations_inland |>
    group_by(geometry) |>
    summarise(n_occurrences=n())
locations_eea <- locations_inland |> st_transform(crs(grid_1km_o)) 
### 1km over shp

locations_1_grid_samples <- st_join(crete_1km, locations_inland, left=F) |>
    distinct(geometry,CELLCODE, long, lat) |>
    group_by(geometry,CELLCODE) |>
    summarise(n_samples=n(),.groups="keep")

locations_1_grid <- st_join(crete_1km, locations_inland, left=F) |>
    dplyr::distinct(CELLCODE, sbspcsn, geometry) |>
    group_by(CELLCODE, geometry) |>
    summarise(n_species=n(), .groups="keep") |>
    mutate(area=st_area(geometry)) |>
    left_join(st_drop_geometry(locations_1_grid_samples))

### pearson correlation of sampling intensity and species in 1 sq. km
cor(locations_1_grid$n_species, locations_1_grid$n_samples) 

### occurrences of each species, species distribution
occurrences_1_grid <- st_join(crete_1km, locations_inland, left=F) |>
    dplyr::distinct(CELLCODE, sbspcsn, geometry) |>
    group_by(sbspcsn) |>
    summarise(geometry=st_union(geometry)) |>
    mutate(area=st_area(geometry))

### 8km over raster
crete_manual_grid8 <- rast(ncol=36,
          nrow=19,
          xmin=5522000,
          xmax=5834000,
          ymin=1400000,
          ymax=1552000,
          crs=crs(grid_1km_o))

values(crete_manual_grid8) <- 0

locations_grid_8_base <- st_as_sf(as.polygons(crete_manual_grid8,aggregate=F,values=T)) |> 
    st_transform(., crs="WGS84") 

locations_grid_8_base$CELLCODE <- seq(1:nrow(locations_grid_8_base))
locations_grid_8_base <- locations_grid_8_base |>
    st_join(locations_inland, left=F)

locations_8_grid_samples <- locations_grid_8_base |>
    distinct(geometry, long, lat, CELLCODE) |>
    group_by(geometry, CELLCODE) |>
    summarise(n_samples=n(),.groups="keep")

locations_8_grid <- locations_grid_8_base |>
    dplyr::distinct(sbspcsn, geometry, CELLCODE) |>
    group_by(geometry, CELLCODE) |>
    summarise(n_species=n(), .groups="keep") |>
    mutate(area=st_area(geometry)) |>
    left_join(st_drop_geometry(locations_8_grid_samples))

### pearson correlation of sampling effort and species
cor(locations_8_grid$n_species, locations_8_grid$n_samples)

### 4km over raster
crete_manual_grid4 <- rast(ncol=72,
          nrow=38,
          xmin=5522000,
          xmax=5834000,
          ymin=1400000,
          ymax=1552000,
          crs=crs(grid_1km_o))

values(crete_manual_grid4) <- 0

locations_grid_4_base <- st_as_sf(as.polygons(crete_manual_grid4,aggregate=F,values=T)) |> 
    st_transform(., crs="WGS84") 

locations_grid_4_base$CELLCODE <- seq(1:nrow(locations_grid_4_base))
locations_grid_4_base <- locations_grid_4_base |>
    st_join(locations_inland, left=F)

locations_4_grid_samples <- locations_grid_4_base |>
    distinct(geometry, long, lat, CELLCODE) |>
    group_by(geometry, CELLCODE) |>
    summarise(n_samples=n(),.groups="keep")

locations_4_grid <- locations_grid_4_base |>
    dplyr::distinct(sbspcsn, geometry, CELLCODE) |>
    group_by(geometry, CELLCODE) |>
    summarise(n_species=n(), .groups="keep") |>
    mutate(area=st_area(geometry)) |>
    left_join(st_drop_geometry(locations_4_grid_samples))

### pearson correlation of sampling effort and species
cor(locations_4_grid$n_species, locations_4_grid$n_samples)

### 1km over raster
crete_manual_grid <- rast(ncol=311,
          nrow=152,
          xmin=5523000,
          xmax=5832000,
          ymin=1400000,
          ymax=1552000,
          crs=crs(grid_1km_o))

values(crete_manual_grid) <- 0

locations_eea <- locations_eea |> distinct(sbspcsn,long, lat, geometry) 
locations_eea_sampling <- locations_eea |> distinct(geometry) 
grid_occ_count <- rasterize(locations_eea, crete_manual_grid, fun='count')
grid_sampling_count <- rasterize(locations_eea_sampling, crete_manual_grid, fun='count')
grid_occ_count_y <- subst(grid_occ_count, NA, 0)
grid_sampling_count_y <- subst(grid_sampling_count, NA, 0)

### WGS 84 of 1km grid EEA
crete_manual_grid_w <- terra::project(crete_manual_grid, "WGS84")
values(crete_manual_grid_w) <- 0

grid_occ_count_w <- rasterize(locations_inland, crete_manual_grid_w, fun='count')
grid_occ_count_w_df <- terra::as.data.frame(grid_occ_count_w, xy=TRUE, cells=TRUE)  


######################### Adaptive resolution with quadtree ######################
qt_sampling <- quadtree(grid_sampling_count_y,
               min_cell_length=2000,
               max_cell_length=10000,
               split_threshold=4,
               split_if_any_na=F,
               split_method = "range")
qt <- qt_sampling
#qt_species <- quadtree(grid_occ_count_y,
#               min_cell_length=2000,
#               max_cell_length=10000,
#               split_threshold=20,
#               split_if_any_na=F,
#               split_method = "range")

#qt_rast <- as_raster(qt)
#qt_rast_w <- project(qt_rast, crs(locations_inland))
#qt_rast_df <- terra::as.data.frame(qt_rast_w, xy=TRUE, cells=TRUE) 

# transform the QT object directly to dataframe to avoid errors from 
# the multiple transformation from terra to sf
qt_df <- as_data_frame(qt)

# function to create polygon from extend
f1 <- function(x) {st_polygon(list(rbind(c(x[['xmin']],x[['ymin']]),
                                          c(x[['xmax']],x[['ymin']]),
                                          c(x[['xmax']],x[['ymax']]),
                                          c(x[['xmin']],x[['ymax']]),
                                          c(x[['xmin']],x[['ymin']]))))
}

qt_df$geom <- apply(X=qt_df,1 , FUN=f1)

# create the sf polygon object from the dataframe
qt_sf <- qt_df |>
    mutate(geometry=st_sfc(geom)) |>
    st_as_sf() |>
    dplyr::select(-c(geom)) |>
    st_set_crs(st_crs(locations_eea))


qt_sf_sp <- qt_sf |>
    st_transform(crs(locations_inland)) |> 
    st_join(locations_inland, left=F) |>
    distinct(sbspcsn,geometry) |>
    group_by(geometry) |>
    mutate(n_species=n()) |>
    left_join(endemic_species, by=c("sbspcsn"="subspeciesname"))

### quadstrees plots
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
plot(vect(grid_10km),
     lwd = .05,
     add=T)
plot(crete_eea, add=T, color="")
#plot(locations_eea_sampling, add=T, color="#d55e00")
dev.off()

### quads and biodiversity
qt_sf_all <- qt_sf_sp |>
    group_by(geometry) |>
    summarise(n_species=n()) |>
    mutate(Order="All")

qt_sf_order <- qt_sf_sp |> 
    group_by(geometry, Order) |> 
    summarise(n_species=n(), .groups="keep")

qt_sf_n_order <- qt_sf_sp |> 
    distinct(geometry, Order) |> 
    group_by(geometry) |> 
    summarise(n_orders=n())

crete_quads_endemics_map <- ggplot() +
    geom_sf(crete_shp, mapping=aes()) +
    geom_sf(locations_grid, mapping=aes()) +
    geom_sf(qt_sf_all, mapping=aes(fill=n_species),
            alpha=0.8,
            colour="transparent",
            na.rm = F,
            show.legend=F) +
    scale_fill_gradient(low="#f0e442",
                        high="#d55e00",
                        guide = "colourbar")+
    coord_sf(crs="wgs84") +
    guides(fill = guide_colourbar(ticks = F,
                                  label = T,
                                  title="# endemics",
                                  title.vjust = 0.8,
                                  order = 1))+
    theme_bw()+
    theme(axis.title=element_blank(),
          axis.text=element_text(colour="black"),
          legend.title = element_text(size=8),
          legend.position = "bottom",
          legend.box.background = element_blank())

ggsave("../plots/crete_quads_endemics_map.png", 
       plot=crete_quads_endemics_map, 
       height = 10, 
       width = 20,
       dpi = 600, 
       units="cm",
       device="png")

crete_quads_endemics_orders_map <- ggplot() +
    geom_sf(crete_shp, mapping=aes()) +
    geom_sf(qt_sf_order, mapping=aes(fill=n_species),
            alpha=0.8,
            colour="transparent",
            na.rm = FALSE,
            show.legend=T) +
    scale_fill_gradient(low="#F0E442",
                        high="#D55E00",
                        guide = "colourbar")+
    coord_sf(crs="WGS84") +
    guides(fill = guide_colourbar(ticks = FALSE,
                                  label = TRUE,
                                  title="# endemics",
                                  title.vjust = 0.8,
                                  order = 1))+
    theme_bw()+
    theme(axis.title=element_blank(),
          axis.text=element_text(colour="black"),
          legend.title = element_text(size=8),
          legend.position = "bottom",
          legend.box.background = element_blank()) +
    facet_grid(Order ~ .)

ggsave("../plots/crete_quads_endemics_orders_map.png",
       plot=crete_quads_endemics_orders_map, 
       height = 60, 
       width = 20,
       dpi = 600, 
       units="cm",
       device="png")

crete_quads_orders_map <- ggplot() +
    geom_sf(crete_shp, mapping=aes()) +
    geom_sf(qt_sf_n_order, mapping=aes(fill=n_orders),
            alpha=0.8,
            colour="transparent",
            na.rm = F,
            show.legend=T) +
    scale_fill_gradient(low="#f0e442",
                        high="#d55e00",
                        guide = "colourbar")+
    coord_sf(crs="wgs84") +
    guides(fill = guide_colourbar(ticks = F,
                                  label = T,
                                  title="# orders",
                                  title.vjust = 0.8,
                                  order = 1))+
    theme_bw()+
    theme(axis.title=element_blank(),
          axis.text=element_text(colour="black"),
          legend.title = element_text(size=8),
          legend.position = "bottom",
          legend.box.background = element_blank())

ggsave("../plots/crete_quads_orders_map.png", 
       plot=crete_quads_orders_map, 
       height = 10, 
       width = 20,
       dpi = 600, 
       units="cm",
       device="png")


## Locations grid and quads

crete_grid10_q <- crete_grid10 |> 
    distinct(CELLCODE,geometry) |>
    mutate(grid_area=units::set_units(st_area(geometry),km^2)) |>
    st_intersection(qt_sf_all) |>
    mutate(ovelap_area=units::set_units(st_area(geometry),km^2)) |>
    group_by(CELLCODE) |> 
    summarise(ovelap_area=sum(ovelap_area)) |>
    filter(ovelap_area>units::set_units(50,km^2))

############################ Endemic hotspots ###############################

endemic_hotspots_quads <- qt_sf_all |>
    filter(n_species >= quantile(n_species, 0.95)) 

endemic_hotspots_1km <- locations_1_grid |>
    filter(n_species >= quantile(n_species, 0.95)) 

endemic_hotspots_4km <- locations_grid_4 |>
    filter(n_species >= quantile(n_species, 0.95)) 

endemic_hotspots_8km <- locations_8_grid |>
    filter(n_species >= quantile(n_species, 0.92)) 

endemic_hotspots <- locations_grid |>
    group_by(CELLCODE) |>
    summarise(n_species=n()) |>
    filter(n_species >= quantile(n_species, 0.90)) |>
    filter(CELLCODE %in% as.character(crete_grid10_q$CELLCODE))

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

crete_quads_sampling_grid <- ggplot() +
    geom_sf(crete_shp, mapping=aes()) +
    geom_sf(locations_1_grid, mapping=aes(fill=n_species), linewidth=0.01, alpha=0.3,size=0.02)+
    scale_fill_gradientn(colours = c("grey70", "purple")) +
    labs(fill = "1 sq. km grid\noccurrences") +
    new_scale_fill() +
    geom_sf(endemic_hotspots_quads, mapping=aes(fill=n_species), alpha=0.8, size=0.1, na.rm = FALSE)+
    scale_fill_gradientn(colours = c("grey80", "#D55E00")) +
    labs(fill = "Quadtrees hotspots")+
    new_scale_fill() +
    geom_sf(locations_inland, mapping=aes(),color="red", size=0.01)+
    new_scale_fill() +
    geom_sf(endemic_hotspots, mapping=aes(fill=n_species), alpha=0.3, size=0.1, na.rm = FALSE)+
    scale_fill_gradientn(colours = c("yellow", "red")) +
    labs(fill = "10 sq. km hotspots") +
    new_scale_fill() +
    geom_sf(endemic_hotspots_1km, mapping=aes(fill=n_species), alpha=0.3, size=0.1, na.rm = FALSE)+
    scale_fill_gradientn(colours = c("brown", "black")) +
    labs(fill = "1 sq. km hotspots") +
    new_scale_fill() +
    geom_sf(endemic_hotspots_4km, mapping=aes(fill=n_species), alpha=0.3, size=0.1, na.rm = FALSE)+
    scale_fill_gradientn(colours = c("cornsilk1", "cyan4")) +
    labs(fill = "4 sq. km hotspots") +
    new_scale_fill() +
    geom_sf(endemic_hotspots_8km, mapping=aes(fill=n_species), alpha=0.3, size=0.1, na.rm = FALSE)+
    scale_fill_gradientn(colours = c("deepskyblue", "dodgerblue4")) +
    labs(fill = "8 sq. km hotspots") +
    ggtitle("Multiple grid sizes of hotspots")+
    theme_bw()+
    theme(axis.title=element_blank(),
          axis.text=element_text(colour="black"),
          legend.title = element_text(size=8),
          legend.position = "bottom",
          legend.box.background = element_blank())


ggsave("../plots/crete_multiple_grids_hotspots.png",
       plot=crete_quads_sampling_grid,
       height = 20, 
       width = 40,
       dpi = 600, 
       unit="cm",
       device="png")
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

endemic_hotspots_natura <- st_intersection(endemic_hotspots, natura_crete_land_sci) |>
    mutate(natura_area=units::set_units(st_area(geometry),km^2)) |>
    group_by(CELLCODE) |>
    summarise(natura_area=sum(natura_area)) |>
    mutate(natura="yes")

locations_grid_yes <- endemic_hotspots |>
    distinct(CELLCODE,geometry) |>
    left_join(st_drop_geometry(endemic_hotspots_natura), by=c("CELLCODE"="CELLCODE")) |>
    mutate(classification="hotspot") 


locations_grid_natura <- locations_grid |>
    distinct(CELLCODE,geometry) |>
    filter(!c(CELLCODE %in% endemic_hotspots$CELLCODE)) |> 
    mutate(classification="no hotspot") |>
    st_intersection(natura_crete_land_sci) |>
    mutate(natura_area=units::set_units(st_area(geometry),km^2)) |>
    group_by(CELLCODE) |>
    summarise(natura_area=sum(natura_area)) |>
    mutate(natura="yes")


locations_grid_no <- locations_grid |>
    distinct(CELLCODE,geometry) |>
    left_join(st_drop_geometry(locations_grid_natura), by=c("CELLCODE"="CELLCODE")) |>
    mutate(classification="no hotspot") 

locations_all_classification <- rbind(locations_grid_no, locations_grid_yes) |>
    replace_na(list(natura="no", natura_area=units::set_units(0,km^2)))


chisq.test(locations_all_classification$natura,locations_all_classification$classification)
table(locations_all_classification$natura,locations_all_classification$classification)

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

occurrences <- occurrences_1_grid |>
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



