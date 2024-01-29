#!/usr/bin/Rscript

## Script name: spatial_analysis.R
##
## Purpose of script: Spatial analysis using data from the Copernicus system
## for elevation, habitats and corine land cover and land use change (HILDA+).
## Main output is an enriched locations file which is exported for further analysis.
##
## How to run:
## Rscript spacial_analysis.R
##
## Execution time: 4 minutes
##
## Author: Savvas Paragkamian
##
## Date Created: 2022-12-22

library(tidyverse)
library(sf)
library(terra)
library(units)
library(ggpubr)
source("functions.R")

# Load data
g_base <- g_base()
supplementary_material_1 <- readxl::read_excel("../data/Supplementary-material-1.xlsx", sheet="arthropods_occurrences")

arthropods_occurrences <- st_as_sf(supplementary_material_1,
                                   coords=c("decimalLongitude","decimalLatitude"),
                                   remove=F,
                                   crs="WGS84")
## Occurrences
locations_shp <- arthropods_occurrences |> 
    dplyr::select(-bibliographicCitation) |>
    distinct()

crete_shp <- sf::st_read("../data/crete/crete.shp")

endemic_species <- read_delim("../results/endemic_species_assessment.tsv", delim="\t")
endemic_hotspots <- st_read("../results/endemic_hotspots/endemic_hotspots.shp")
threatspots <- st_read("../results/threatspots/threatspots.shp")
threatspots_lt <- threatspots |> 
    filter(pc_thrt>= quantile(pc_thrt,0.90))

wege_results <- st_read("../results/wege_results/wege_results.shp")
natura_crete <- sf::st_read("../data/natura2000/natura2000_crete.shp")
wdpa_crete <- sf::st_read("../data/wdpa_crete/wdpa_crete.shp")

wdpa_crete_wildlife <- wdpa_crete |> filter(DESIG_ENG=="Wildlife Refugee")
natura_crete_land <- st_intersection(natura_crete, crete_shp)

# split the SPA SCI

natura_crete_land_sci <- natura_crete_land |> filter(SITETYPE=="B")

# Spatial data

locations_shp <- st_join(locations_shp, natura_crete_land_sci, left=T)
## CORINE Land Cover nomenclature
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
                threatspots = threatspots_lt,
                wege_kba = wege_results)

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

dem_crete_locations <- terra::extract(dem_crete, locations_shp, cellnumbers=F, ID=F)
locations_shp$elevation <- dem_crete_locations$dem_crete

dem_crete_df <- terra::as.data.frame(dem_crete, xy=TRUE, cells=TRUE)

g_dem <- g_base +
    geom_raster(dem_crete_df, mapping=aes(x=x, y=y, fill=dem_crete))+
    geom_sf(locations_shp, mapping=aes(),color="blue", size=0.1, alpha=0.2)

ggsave("../plots/crete_dem.png", plot=g_dem, device="png")


g_ele <- g_base + 
    geom_sf(locations_shp, mapping=aes(color=elevation), size=0.1, alpha=0.2) +
    scale_color_gradientn(colours = terrain.colors(5)) 

ggsave("../plots/crete_occurrences_dem.png", plot=g_ele, device="png")


## Crete land use change
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


#### Hilda analysis
hilda_path <- "../data/hildap_GLOB-v1.0_lulc-states_crete/"
hilda_id_names <- read_delim(paste0(hilda_path, "hilda_transitions_names.tsv", sep=""), delim="\t")
hilda_all <- list.files(hilda_path)
hilda_files <- hilda_all[grepl("*.tif", hilda_all)] 

#create_dir("../plots/hilda_crete")

for (i in 1:length(hilda_files)){
    print(i)
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
    
#    ggsave(paste0("../plots/hilda_crete/crete_",raster_name,"_map.png",sep=""),
#           plot=g_hilda_map,
#           height = 10, 
#           width = 20,
#           dpi = 300, 
#           units="cm",
#           device="png")
    
    hilda_sum <- zonal(cellSize(raster_file), raster_file, "sum") |> 
        mutate(area_m2=units::set_units(area,m^2)) |>
        mutate(area=units::set_units(area/10^6, km^2)) 
    
    hilda_sum <- hilda_sum |>
        mutate(hilda_id=as.character(hilda_sum[,1])) |>
        filter(hilda_sum[,1]>0) |>
        left_join(hilda_cat) 
    
    hilda_sum_g <- ggplot()+
        geom_col(hilda_sum,
                 mapping= aes(y=as.numeric(area),
                              x="",
                              fill = hilda_name),
                 position = position_stack(),
                 width = 0.2) +
        scale_fill_manual(values=hilda_cat_v) +
        scale_x_discrete(expand = expansion(add=c(0,0)))+
        scale_y_continuous(breaks=seq(0,9000,1000),
                           limits=c(0,8900),
                           expand = c(0,0))+
        ylab("Area sq. km") +
        xlab("") +
        theme_bw()+
        theme(legend.position='none',
              axis.ticks.x=element_blank(),
              panel.border = element_blank(),
              panel.grid.major = element_blank(), #remove major gridlines
              panel.grid.minor = element_blank()) #remove minor gridlines
    
#    ggsave(paste0("../plots/hilda_crete/crete_",raster_name,"_bar.png",sep=""),
#           plot=hilda_sum_g,
#           height = 10, 
#           width = 10,
#           dpi = 300, 
#           units="cm",
#           device="png")

    fig_hilda <- ggarrange(g_hilda_map,hilda_sum_g,
              labels = c("A", "B"),
              ncol = 2,
              nrow = 1,
              widths = c(0.85,0.15),
              font.label=list(color="black",size=15),
              common.legend = TRUE,
              legend="bottom") + bgcolor("white")
    
    ggsave(paste0("../plots/hilda_crete/crete_",raster_name,".png",sep=""), 
           plot=fig_hilda, 
           height = 10, 
           width = 25,
           dpi = 300, 
           units="cm",
           device="png")
}

## HILDA difference of land use, 1999-2019
##
hilda_1998 <- rast("../data/hildap_GLOB-v1.0_lulc-states_crete/crete_hilda_plus_1998_states_GLOB-v1-0_wgs84-nn.tif")

hilda_1998_locations <- terra::extract(hilda_1998, locations_shp, cellnumbers=F, ID=F)
colnames(hilda_1998_locations) <- c("hilda_1998")

hilda_2018 <- rast("../data/hildap_GLOB-v1.0_lulc-states_crete/crete_hilda_plus_2018_states_GLOB-v1-0_wgs84-nn.tif")
hilda_2018_locations <- terra::extract(hilda_2018, locations_shp, cellnumbers=F, ID=F)
colnames(hilda_2018_locations) <- c("hilda_2018")
## what is the transitions?
## the raster objects contain numeric values. The smart thing about this dataset
## is that I can use the numbers that are to show the category and create new
## numbers of the difference to symbolise the transitions.
## The oldest raster is the origin so transform it to the closest decade,
## 11 = 10, 22 = 20 etc. This can be accomplished with
## Modulus operation 44 - 44 %% 10
## Transform the latest raster to the units. so 44 = 4,
## Modulus operation 44 %% 10

hilda_1998_o <- app(hilda_1998, fun=function(i) i-i %% 10)
hilda_2018_o <- app(hilda_2018, fun=function(i) i %% 10)

hilda_1998_2018 <- hilda_1998_o + hilda_2018_o

hilda_1998_2018_sf <- st_as_sf(as.polygons(hilda_1998_2018,aggregate=F,values=T)) 

st_write(hilda_1998_2018_sf,
         "../results/hilda_1998_2018/hilda_1998_2018.shp", 
         append=F,
         delete_layer=T,
         delete_dsn = TRUE) 

hilda_transition <- terra::extract(hilda_1998_2018, locations_shp, cellnumbers=F, ID=F)
colnames(hilda_transition) <- c("hilda_transition")

hilda_all <- cbind(hilda_transition,hilda_2018_locations,hilda_1998_locations)

locations_shp <- cbind(locations_shp, hilda_all)

### Hilda summary
hilda_1998_2018_df <- terra::as.data.frame(hilda_1998_2018, xy=TRUE, cells=TRUE) |>
    filter(lyr.1>0) |>
    mutate(hilda_transition=as.character(lyr.1))
#    left_join(hilda_cat, by=c("hilda_id"="hilda_id"))

hilda_sum <- zonal(cellSize(hilda_1998_2018), hilda_1998_2018, "sum") |> 
    mutate(crete_m2=units::set_units(area,m^2)) |>
    mutate(crete=units::set_units(area/10^6, km^2)) 

hilda_sum <- hilda_sum |>
    mutate(hilda_id_transition=hilda_sum[,1]) |>
    filter(hilda_sum[,1]>0) |> 
    left_join(hilda_id_names, by=c("hilda_id_transition"="hilda_id"))

write_delim(hilda_sum, "../results/hilda_1998_2018.tsv", delim="\t")

## Natura2000 diffence of land use, 1998-2018
natura_crete_land_hilda <- terra::mask(hilda_1998_2018, natura_crete_land)

natura_crete_land_hilda_sum <- zonal(cellSize(natura_crete_land_hilda), natura_crete_land_hilda, "sum") |> 
    mutate(natura2000_m2=units::set_units(area,m^2)) |>
    mutate(natura2000=units::set_units(area/10^6, km^2)) |>
    rename(hilda_id_transition=lyr.1) |>
    filter(hilda_id_transition>0) |>
    left_join(hilda_id_names, by=c("hilda_id_transition"="hilda_id")) |>
    dplyr::select(-c(hilda_id_transition, area))

write_delim(natura_crete_land_hilda_sum, "../results/natura_crete_land_hilda_sum.tsv", delim="\t")

## Endemicity hotspots and hilda

endemic_hotspots_hilda <- terra::mask(hilda_1998_2018, endemic_hotspots)

endemic_hotspots_hilda_sum <- zonal(cellSize(endemic_hotspots_hilda), endemic_hotspots_hilda, "sum") |> 
    mutate(hotspots_m2=units::set_units(area,m^2)) |>
    mutate(hotspots=units::set_units(area/10^6, km^2)) |>
    rename(hilda_id_transition=lyr.1) |>
    filter(hilda_id_transition>0) |>
    left_join(hilda_id_names, by=c("hilda_id_transition"="hilda_id"))|>
    dplyr::select(-c(hilda_id_transition, area))

write_delim(endemic_hotspots_hilda_sum, "../results/endemic_hotspots_hilda_sum.tsv", delim="\t")

## WEGE KBAs and hilda

wege_kba_hilda <- terra::mask(hilda_1998_2018, wege_results)

wege_kba_hilda_sum <- zonal(cellSize(wege_kba_hilda), wege_kba_hilda, "sum") |> 
    mutate(wege_kba_m2=units::set_units(area,m^2)) |>
    mutate(wege_kba=units::set_units(area/10^6, km^2)) |>
    rename(hilda_id_transition=lyr.1) |>
    filter(hilda_id_transition>0) |>
    left_join(hilda_id_names, by=c("hilda_id_transition"="hilda_id")) |>
    dplyr::select(-c(hilda_id_transition, area))

write_delim(wege_kba_hilda_sum, "../results/wege_kba_hilda_sum.tsv", delim="\t")

## Threatspots and hilda
threatspots_hilda <- terra::mask(hilda_1998_2018, threatspots_lt)

threatspots_hilda_sum <- zonal(cellSize(threatspots_hilda), threatspots_hilda, "sum") |> 
    mutate(threatspots_m2=units::set_units(area,m^2)) |>
    mutate(threatspots=units::set_units(area/10^6, km^2)) |>
    rename(hilda_id_transition=lyr.1) |>
    filter(hilda_id_transition>0) |>
    left_join(hilda_id_names, by=c("hilda_id_transition"="hilda_id")) |>
    dplyr::select(-c(hilda_id_transition, area))

write_delim(threatspots_hilda_sum, "../results/threatspots_hilda_sum.tsv", delim="\t")

### all together
###
hilda_trans <- hilda_sum |>
    left_join(natura_crete_land_hilda_sum, by=c("hilda_name"="hilda_name")) |>
    left_join(endemic_hotspots_hilda_sum, by=c("hilda_name"="hilda_name")) |>
    left_join(wege_kba_hilda_sum, by=c("hilda_name"="hilda_name")) |>
    left_join(threatspots_hilda_sum, by=c("hilda_name"="hilda_name")) |>
    mutate_if(is.numeric, round, 0)

write_delim(hilda_trans, "../results/hilda_summary_transitions.tsv", delim="\t")

# Export of locations shapefile with all the spatial metadata
st_write(locations_shp,
         "../results/locations_spatial/locations_spatial.shp",
         layer_options = "ENCODING=UTF-8", 
         delete_layer=T, 
         delete_dsn = TRUE,
         append=TRUE)

write_delim(st_drop_geometry(locations_shp),
         "../results/locations_spatial/locations_spatial.tsv",
         delim="\t")

##################### create CELLCOD metadata #######################

locations_grid <- st_read("../results/locations_grid/locations_grid.shp") |>
    rename("CELLCODE"="CELLCOD") |>
    rename("scientificName"="scntfcN")

threatspots_df <- st_drop_geometry(wege_results) |> mutate(threatspot="threatspot") 
endemic_hotspots_df <- st_drop_geometry(endemic_hotspots) |>
    mutate(hotspot="hotspot") |>
    dplyr::select(-n_species)

locations_grid_d <- locations_grid |> distinct(CELLCODE,geometry)
eea_10_crete_metadata <- st_intersection(locations_grid_d, clc_crete_shp)

eea_10_crete_metadata_d <- eea_10_crete_metadata |>
    mutate(area = st_area(geometry)) |>
    group_by(LABEL2, CELLCODE) |>
    summarise(area_sum=sum(area), .groups="keep") |>
    st_drop_geometry() |> units::drop_units() |>
    pivot_wider(id_cols="CELLCODE", names_from="LABEL2", values_from="area_sum", values_fill=0)

grid_stats <- locations_grid |>
    st_drop_geometry() |>
    group_by(CELLCODE) |>
    summarise(n_species=n()) |>
    left_join(endemic_hotspots_df, by=c("CELLCODE"="CELLCODE")) |>
    left_join(threatspots_df, by=c("CELLCODE"="CELLCOD")) |>
    left_join(eea_10_crete_metadata_d)

grid_stats_t <- grid_stats |>
    mutate_at(c('hotspot','threatspot'), ~replace_na(.,"no")) |>
    mutate_at(c('LT','PT','pc_thrt','wege'), ~replace_na(.,0))

grid_stats_long <- grid_stats_t |>
    pivot_longer(!c(CELLCODE,hotspot,threatspot),
                 names_to="variables",
                 values_to="values")
