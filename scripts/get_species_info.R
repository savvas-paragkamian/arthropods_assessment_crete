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

greek_redlist_a <- read_delim("../data/greece-redlist/assessments.csv", delim=",")
greek_redlist_t <- read_delim("../data/greece-redlist/taxonomy.csv", delim=",") %>%
    dplyr::select(-scientificName)

greek_redlist <- greek_redlist_a %>%
    left_join(greek_redlist_t, by=c("internalTaxonId"="internalTaxonId")) %>%
    mutate(threatened=if_else(redlistCategory %in% c("Critically Endangered","Endangered","Vulnerable"), TRUE, FALSE)) %>%
    group_by(orderName, threatened) %>%
    summarise(n=n()) %>%
    mutate(source="greek_redlist")

endemic_greek_redlist_a <- read_delim("../data/endemic-greek-redlist/assessments.csv", delim=",")
endemic_greek_redlist_t <- read_delim("../data/endemic-greek-redlist/taxonomy.csv", delim=",") %>%
    dplyr::select(-scientificName)

endemic_greek_redlist <- endemic_greek_redlist_a %>%
    left_join(endemic_greek_redlist_t, by=c("internalTaxonId"="internalTaxonId")) %>%
    mutate(threatened=if_else(redlistCategory %in% c("Critically Endangered","Endangered","Vulnerable"), TRUE, FALSE)) %>%
    group_by(orderName, threatened) %>%
    summarise(n=n()) %>%
    mutate(source="endemic_greek_redlist")


endemic_crete_redlist_s <- endemic_species %>% 
    left_join(greek_redlist_a, by=c("subspeciesname"="scientificName")) %>%
    dplyr::select(subspeciesname,paca,iucn,redlistCategory) %>%
    na.omit()

write_delim(endemic_crete_redlist_s, "../results/endemic_crete_redlist.tsv", delim="\t")

endemic_crete_redlist <- greek_redlist_a %>% 
    filter(scientificName %in% endemic_species$subspeciesname) %>%
    left_join(greek_redlist_t, by=c("internalTaxonId"="internalTaxonId")) %>%
    mutate(threatened=if_else(redlistCategory %in% c("Critically Endangered","Endangered","Vulnerable"), TRUE, FALSE)) %>%
    group_by(orderName, threatened) %>%
    summarise(n=n()) %>%
    mutate(source="endemic_crete_redlist")


endemic_europe_redlist_a <- read_delim("../data/endemic-europe-redlist/assessments.csv", delim=",")
endemic_europe_redlist_t <- read_delim("../data/endemic-europe-redlist/taxonomy.csv", delim=",") %>%
    dplyr::select(-scientificName)

endemic_europe_redlist <- endemic_europe_redlist_a %>%
    left_join(endemic_europe_redlist_t, by=c("internalTaxonId"="internalTaxonId")) %>%
    mutate(threatened=if_else(redlistCategory %in% c("Critically Endangered","Endangered","Vulnerable"), TRUE, FALSE)) %>%
    group_by(orderName, threatened) %>%
    summarise(n=n()) %>%
    mutate(source="endemic_europe_redlist")

world_redlist_a <- read_delim("../data/world-redlist/assessments.csv", delim=",")
world_redlist_t <- read_delim("../data/world-redlist/taxonomy.csv", delim=",") %>%
    dplyr::select(-scientificName)

world_redlist <- world_redlist_a %>%
    left_join(world_redlist_t, by=c("internalTaxonId"="internalTaxonId")) %>%
    mutate(threatened=if_else(redlistCategory %in% c("Critically Endangered","Endangered","Vulnerable"), TRUE, FALSE)) %>%
    group_by(orderName, threatened) %>%
    summarise(n=n()) %>%
    mutate(source="world_redlist")

redlist_threatened <- rbind(endemic_crete_redlist,endemic_greek_redlist, endemic_europe_redlist, world_redlist) %>%
    rename(Order=orderName) %>% 
    mutate(Order=paste(substr(Order, 1, 1),substr(tolower(Order), 2, nchar(Order)),sep=""))

write_delim(redlist_threatened, "../results/redlist_threatened.tsv", delim="\t")

# clean the subspecies column because it has null value and dates.
iucn_arthropods_species <- iucn_arthropods %>% mutate(subspeciesname=if_else(grepl("^[[:digit:]]", subspecies) | subspecies %in% c("<NULL>","<Null>") | is.na(subspecies),sci_name, paste(sci_name,subspecies))) %>%
                            distinct(sci_name,subspeciesname)

write_delim(iucn_arthropods_species, "../results/iucn_arthropods_species", delim="\t")
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

