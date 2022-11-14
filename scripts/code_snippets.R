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

# ConR package returned an error with the exclution of the land
# so a new function named "eoo_calculation" is created
