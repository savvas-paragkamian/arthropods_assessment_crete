#!/usr/bin/Rscript

## Script name: get_species_info.R
##
## Purpose of script: use public databases to retrieve information of 
## species regarding their taxonomy, their global distribution and their 
## IUCN status
##
## Author: Savvas Paragkamian
##
## Date Created: 2022-12-22

library(sf)
library(tidyverse)
library(readxl)
library(rredlist)
library(taxize)
library(units)
library(vegan)


locations_shp <- sf::st_read("../data/arthropods_occurrences/arthropods_occurrences.shp")
locations_spatial <- sf::st_read("../results/locations_spatial/locations_spatial.shp")
locations_grid <- sf::st_read("../results/locations_grid/locations_grid.shp") 
crete_shp <- sf::st_read("../data/crete/crete.shp")
endemic_species <- read_delim("../results/endemic_species_assessment.tsv", delim="\t")



# IUCN dataset from website
iucn_arthropods <- read_delim("~/Downloads/redlist_species_data_e63da44c-39d3-46bc-819d-1a3203f508bd/points_data.csv",
                              delim=",") %>% 
    st_as_sf(coords=c("longitude", "latitude"),
             remove=FALSE,
             crs="WGS84")
iucn_arthropods_species <- iucn_arthropods %>% mutate(subspeciesname=if_else(subspecies=="<NULL>" | is.na(subspecies),
                                                         sci_name, 
                                                         paste(sci_name,subspecies))) %>%
                            distinct(sci_name,subspeciesname)


iucn_arthropods_species_u <- unique(iucn_arthropods_species$sci_name)

gbif_iucn <- get_gbifid(iucn_arthropods_species_u,ask=F)

iucn_arthropods_species_gbif <- data.frame(sci_name=iucn_arthropods_species_u, gbifid=gbif_iucn)

classification_iucn <- classification(iucn_arthropods_species_gbif$gbifid.ids, db = 'gbif')

classification_iucn_d <- do.call(rbind, classification_iucn) %>%
    rownames_to_column(var="gbif") %>% 
    mutate(gbif = gsub("\\.(.*)","", gbif)) %>%
    dplyr::select(-id) %>%
    group_by(gbif,name,rank) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::select(-n) %>%
    pivot_wider(names_from=rank, values_from=name) %>%
    dplyr::select(-`NA`) %>%
    filter(!is.na(gbif))

write_delim(classification_iucn_d, "../results/classification_iucn_gbif.tsv", delim="\t")

# Resolve names
## gnr_datasources() %>% filter(title=="GBIF Backbone Taxonomy") id=11
gnr_species <- gnr_resolve(endemic_species$subspeciesname)
gnr_species_gbif <- gnr_resolve(endemic_species$subspeciesname, data_source_ids=11)

write_delim(gnr_species, "../results/gnr_species_names.tsv", delim="\t")
# Get GBIF ids

res <- get_gbifid(endemic_species$subspeciesname,ask=F)
# Total: 343
# Found: 306

endemic_species$gbif <- as.numeric(res)

endemic_species_no_gbif <- endemic_species[which(is.na(endemic_species$gbif)),]

# Classification

classification_s <- classification(endemic_species$gbif, db = 'gbif')

classification_s_d <- do.call(rbind, classification_s) %>%
    rownames_to_column(var="gbif") %>% 
    mutate(gbif = gsub("\\.(.*)","", gbif)) %>%
    pivot_wider()

classification_s_d_w <- classification_s_d %>% dplyr::select(-id) %>%
    group_by(gbif, rank, name) %>% 
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::select(-n) %>%
    pivot_wider(names_from=rank, values_from=name) %>%
    mutate(gbif=as.numeric(gbif)) %>%
    dplyr::select(-`NA`) %>%
    filter(!is.na(gbif))

endemic_species_tax <- endemic_species %>%
    left_join(classification_s_d_w, by=c("gbif"="gbif"))
write_delim(endemic_species_tax, "../results/endemic_species_taxonomy.tsv", delim="\t")

endemic_species_tax <- read_delim("../results/endemic_species_taxonomy.tsv", delim="\t")

