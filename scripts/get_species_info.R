#!/usr/bin/Rscript

## Script name: code_snippets.R
##
## Purpose of script: here are some parts of code to crop spatial 
## objects around the area of Crete and other test code.
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
endemic_species <- read_delim("../results/endemic_species_paca.tsv", delim="\t")

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
