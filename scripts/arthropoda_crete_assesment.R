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

#arthropods <- readxl::read_excel("../data/arthropoda_crete_nhmc_for_analysis.xlsx")
## remove Opiliones because the dataset has some errors that are under review.
#arthropods_kriti <- readxl::read_excel("../data/arthropoda_crete_nhmc_for_analysis.xlsx") %>% 
#    filter(Order!="Opiliones") %>% 
#    dplyr::select(-Ergasia)

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

locations_shp <- st_as_sf(arthropods_occurrences,
                          coords=c("logD", "latD"),
                          remove=TRUE,
                          crs="WGS84")

## Crete Spatial data

periphereies_shp <- sf::st_read("~/Documents/spatial_data/periphereies/periphereies.shp")

#### crete is the 12th region in this shapefile
#### https://geodata.gov.gr/en/dataset/28121643-d977-48eb-a8ca-a6fac6b4af6d/resource/7c80a2c1-93b7-4814-9fc4-245e775acaa6/download/periphereies.zip

crete_shp <- periphereies_shp[12,] %>% st_transform(., "WGS84")

## here we calculate the area of each polygon of the region of Crete
## and subsequently we keep only the largest polygon e.g. Crete island only.

crete_polygon <- st_cast(crete_shp, "POLYGON") %>% 
    mutate(area_m2=as.numeric(st_area(.))) %>%
    filter(area_m2==max(area_m2))

locations_inland <- st_join(locations_shp, crete_polygon, left=F)

## EEA reference grid
## https://www.eea.europa.eu/data-and-maps/data/eea-reference-grids-2
##

grid_10km <- sf::st_read(dsn="../data/Greece_shapefile/gr_10km.shp") %>%
    st_transform(., crs="WGS84")

### keep grid that overlaps with Crete
crete_grid10m <- st_join(grid_10km, crete_polygon, left=F)

## Here is Crete with all the sampling points
##
g <- ggplot() +
    geom_sf(crete_polygon, mapping=aes()) +
    geom_sf(locations_inland, mapping=aes(),color="blue", size=0.1, alpha=0.2) +
    geom_sf(crete_grid10m, mapping=aes(),color="red", alpha=0.2, size=0.1) +
    coord_sf(crs="WGS84") +
    theme_bw()

ggsave("../plots/crete-occurrences.png", plot=g, device="png")

# Processing

## EOO
## Data transformation for ConR package

locations_inland_df <- arthropods_occurrences %>%
    dplyr::rename(ddlon=logD, ddlat=latD, tax=subspeciesname) %>% 
    dplyr::select(-Order) %>%
    relocate(ddlat,ddlon, tax)

crete_spatial <- as(crete_polygon, "Spatial")

### calculations

eoo_results_list <- EOO.computing(locations_inland_df, country_map=crete_spatial, export_shp=T,write_shp=T)

eoo_results <- read_delim("EOO.results.csv", delim=",", col_names=T)

## AOO
## see for bootstrap
AOO_endemic <- AOO.computing(locations_inland_df,Cell_size_AOO =2 )

AOO_endemic_df <- tibble(subspeciesname=names(AOO_endemic),AOO=AOO_endemic, row.names=NULL)

## Preliminary Automated Conservation Assessments (PACA)

endemic_species <- arthropods_occurrences %>% 
    group_by(subspeciesname,Order) %>% 
    summarize(n_locations=n()) %>% 
    ungroup() %>%
    left_join(AOO_endemic_df, by=c("subspeciesname"="subspeciesname")) %>%
    left_join(eoo_results, by=c("subspeciesname"="...1")) %>% 
    mutate(Potentially_VU=if_else(n_locations<10 & (EOO<20000 | AOO<2000),TRUE,FALSE)) %>%
    mutate(Potentially_EN=if_else(n_locations<5 & (EOO<5000 | AOO<500),TRUE,FALSE)) %>%
    mutate(Potentially_CR=if_else(n_locations==1 & (EOO<500 | AOO<10),TRUE,FALSE)) %>%
    ungroup() %>%
    mutate(potential_status=if_else(Potentially_CR=="TRUE","Potentially_CR", if_else(Potentially_EN=="TRUE","Potentially_EN",if_else(Potentially_VU=="TRUE","Potentially_VU","FALSE"))))

write_delim(endemic_species, "../data/endemic_species_iucn.tsv", delim="\t") 

endemic_species_threatened <- endemic_species %>% filter(potential_status!="FALSE") 

write_delim(endemic_species_threatened, "../data/endemic_species_threatened.tsv", delim="\t") 

endemic_species_cr <- endemic_species %>% filter(potential_status=="Potentially_CR") 

write_delim(endemic_species_cr, "../data/endemic_species_cr.tsv", delim="\t") 

## Raster analysis

over_grid_10k_all <- sp::over( x = locations_shp_all , y = grid_10km.wgs , fn = NULL)

species_over_grid_all <- over_grid_10k_all %>%
    mutate(id=seq(1:nrow(.))) %>% 
    left_join(locations_all, by=c("id"="id")) %>% 
    distinct(CELLCODE,subspeciesname) %>%
    group_by(CELLCODE) %>%
    summarize(total_species=n())


grid_10km_species_data_all <- grid_10km.wgs_data %>%
    left_join(species_over_grid_all, by=c("CELLCODE"="CELLCODE"))

### Endemics grid overlap
over_grid_10k <- sp::over( x = locations_shp , y = grid_10km.wgs , fn = NULL)

species_over_grid <- over_grid_10k %>%
    mutate(id=seq(1:nrow(.))) %>% 
    left_join(locations, by=c("id"="id")) %>% 
    distinct(CELLCODE,subspeciesname) %>%
    group_by(CELLCODE) %>%
    summarize(endemic_species=n())

grid_10km_species_data <- grid_10km.wgs_data %>%
    left_join(species_over_grid, by=c("CELLCODE"="CELLCODE")) %>%
    left_join(species_over_grid_all, by=c("CELLCODE"="CELLCODE")) %>% 
    mutate(hot_spot=if_else(total_species==1,0,endemic_species/total_species))

map_greece_plot_grid_endemic <- ggplot()+
  geom_polygon(data = grid_10km_species_data,
               aes(x=long, y=lat,group = group, fill=endemic_species),
               lwd=0.12, color="orange")+
  geom_polygon(data = hellenic_borders_df,
               aes(x=long, y=lat,group = group),
               lwd=0.12,color="black")+
  geom_polygon(data = grid_10km_df,
               aes(x=long, y=lat,group = group),
               lwd=0.12, fill=NA, color="orange")+
  geom_point(data = locations,aes(x=logD, y=latD,color=Order),size = 0.2)+
  scale_fill_gradient(low = "yellow", high = "red", na.value = NA)+ 
  ggtitle("endemic")+
  labs(x="Longitude",y="Latitude")+
#  scale_fill_gradientn(colours = c("gray100",terrain.colors(10)),na.value =NA ,name="Number of species")+ #c("gray100","gray50","gray40","gray35","gray30","gray20","gray10","gray0")
  scale_x_continuous(breaks = seq(23,27,0.5),limits = c(23,27))+
  scale_y_continuous(breaks = seq(34,36,0.5),limits = c(34,37))+
  coord_map(xlim = c(23,27), ylim = c(34,37))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),legend.position = c(0.87, 0.73),legend.text = element_text(size=9),legend.title = element_text(size=10))
#geom_text(data = sisquoc, aes(label = paste("  ", as.character(name), sep="")), angle = 60, hjust = 0, color = "yellow")
    
ggsave("map_greece_plot_grid_endemic.png", plot = map_greece_plot_grid_endemic, device = "png",width = 30,height = 30,units = "cm",dpi = 300 ,path = "../plots/")

# total
map_greece_plot_grid_total <- ggplot()+
  geom_polygon(data = grid_10km_species_data,
               aes(x=long, y=lat,group = group, fill=total_species),
               lwd=0.12, color="orange")+
  geom_polygon(data = hellenic_borders_df,
               aes(x=long, y=lat,group = group),
               lwd=0.12,color="black")+
  geom_polygon(data = grid_10km_df,
               aes(x=long, y=lat,group = group),
               lwd=0.12, fill=NA, color="orange")+
  geom_point(data = locations,aes(x=logD, y=latD,color=Order),size = 0.2)+
  scale_fill_gradient(low = "yellow", high = "red", na.value = NA)+ 
  ggtitle("total")+
  labs(x="Longitude",y="Latitude")+
#  scale_fill_gradientn(colours = c("gray100",terrain.colors(10)),na.value =NA ,name="Number of species")+ #c("gray100","gray50","gray40","gray35","gray30","gray20","gray10","gray0")
  scale_x_continuous(breaks = seq(23,27,0.5),limits = c(23,27))+
  scale_y_continuous(breaks = seq(34,36,0.5),limits = c(34,37))+
  coord_map(xlim = c(23,27), ylim = c(34,37))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),legend.position = c(0.87, 0.73),legend.text = element_text(size=9),legend.title = element_text(size=10))
#geom_text(data = sisquoc, aes(label = paste("  ", as.character(name), sep="")), angle = 60, hjust = 0, color = "yellow")
    
ggsave("map_greece_plot_grid_total.png", plot = map_greece_plot_grid_total, device = "png",width = 30,height = 30,units = "cm",dpi = 300 ,path = "../plots/")
### hot spots 

map_greece_plot_grid_hot_spot <- ggplot()+
  geom_polygon(data = grid_10km_species_data,
               aes(x=long, y=lat,group = group, fill=hot_spot),
               lwd=0.12, color="orange")+
  geom_polygon(data = hellenic_borders_df,
               aes(x=long, y=lat,group = group),
               lwd=0.12,color="black")+
  geom_polygon(data = grid_10km_df,
               aes(x=long, y=lat,group = group),
               lwd=0.12, fill=NA, color="orange")+
  geom_point(data = locations,aes(x=logD, y=latD,color=Order),size = 0.2)+
  ggtitle("hot spot")+
  scale_fill_gradient(low = "yellow", high = "red", na.value = NA)+ 
  labs(x="Longitude",y="Latitude")+
#  scale_fill_gradientn(colours = c("gray100",terrain.colors(10)),na.value =NA ,name="Number of species")+ #c("gray100","gray50","gray40","gray35","gray30","gray20","gray10","gray0")
  scale_x_continuous(breaks = seq(23,27,0.5),limits = c(23,27))+
  scale_y_continuous(breaks = seq(34,36,0.5),limits = c(34,37))+
  coord_map(xlim = c(23,27), ylim = c(34,37))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),legend.position = c(0.87, 0.73),legend.text = element_text(size=9),legend.title = element_text(size=10))
#geom_text(data = sisquoc, aes(label = paste("  ", as.character(name), sep="")), angle = 60, hjust = 0, color = "yellow")
    
ggsave("map_greece_plot_grid_hot_spot.png", plot = map_greece_plot_grid_hot_spot, device = "png",width = 30,height = 30,units = "cm",dpi = 300 ,path = "../plots/")

# facet hot spot 
map_greece_plot_grid_hot_spot_facet <- ggplot()+
  geom_polygon(data = grid_10km_species_data,
               aes(x=long, y=lat,group = group, fill=hot_spot),
               lwd=0.12, color="orange")+
  geom_polygon(data = hellenic_borders_df,
               aes(x=long, y=lat,group = group),
               lwd=0.12,color="black")+
  geom_polygon(data = grid_10km_df,
               aes(x=long, y=lat,group = group),
               lwd=0.12, fill=NA, color="orange")+
  geom_point(data = locations,aes(x=logD, y=latD),size = 0.5)+
  ggtitle("hot spot")+
  scale_fill_gradient(low = "yellow", high = "red", na.value = NA)+ 
  labs(x="Longitude",y="Latitude")+
#  scale_fill_gradientn(colours = c("gray100",terrain.colors(10)),na.value =NA ,name="Number of species")+ #c("gray100","gray50","gray40","gray35","gray30","gray20","gray10","gray0")
  scale_x_continuous(breaks = seq(23,27,0.5),limits = c(23,27))+
  scale_y_continuous(breaks = seq(34,36,0.5),limits = c(34,37))+
  coord_map(xlim = c(23,27), ylim = c(34,37))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = c(0.87, 0.73),
        legend.text = element_text(size=9),
        legend.title = element_text(size=10))
  facet_wrap(~Order)
    
ggsave("map_greece_plot_grid_hot_spot.png", plot = map_greece_plot_grid_hot_spot, device = "png",width = 30,height = 30,units = "cm",dpi = 300 ,path = "../plots/")



### only the top species
grid_10km_species_data_top <- grid_10km_species_data %>% filter(total_species>50)

map_greece_plot_grid_10_top <- ggplot()+
  geom_polygon(data = grid_10km_species_data_top,
               aes(x=long, y=lat,group = group, fill=total_species),
               lwd=0.12, color="orange")+
  geom_polygon(data = hellenic_borders_df,
               aes(x=long, y=lat,group = group),
               lwd=0.12,color="black")+
  geom_polygon(data = grid_10km_df,
               aes(x=long, y=lat,group = group),
               lwd=0.12, fill=NA, color="orange")+
  geom_point(data = locations,aes(x=logD, y=latD,color=Order),size = 0.2)+
  scale_fill_gradient(low = "yellow", high = "red", na.value = NA)+ 
  labs(x="Longitude",y="Latitude")+
#  scale_fill_gradientn(colours = c("gray100",terrain.colors(10)),na.value =NA ,name="Number of species")+ #c("gray100","gray50","gray40","gray35","gray30","gray20","gray10","gray0")
  scale_x_continuous(breaks = seq(23,27,0.5),limits = c(23,27))+
  scale_y_continuous(breaks = seq(34,36,0.5),limits = c(34,37))+
  coord_map(xlim = c(23,27), ylim = c(34,37))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),legend.position = c(0.87, 0.73),legend.text = element_text(size=9),legend.title = element_text(size=10))
#geom_text(data = sisquoc, aes(label = paste("  ", as.character(name), sep="")), angle = 60, hjust = 0, color = "yellow")
    
ggsave("map_greece_plot_grid_10_top.png", plot =map_greece_plot_grid_10_top, device = "png",width = 30,height = 30,units = "cm",dpi = 300 ,path = "../plots/")



### Threatened maps

locations_threatened <- endemic_species_threatened %>% 
    dplyr::select(!Order) %>%
    left_join(arthropods_kriti_endemic,by=c("subspeciesname"="subspeciesname")) %>%
    filter(latD<38, logD<30) %>%
    mutate(id=seq(1:nrow(.)))

write_delim(locations_threatened,"../data/locations_threatened.tsv",delim="\t")

locations_threatened_shp <- locations_threatened %>% 
    dplyr::select(subspeciesname,latD, logD) %>% 
    relocate(latD,logD,subspeciesname)

coordinates(locations_threatened_shp)<-~ logD+latD
proj4string(locations_threatened_shp) <- CRS("+proj=longlat +datum=WGS84")

over_grid_10k_thr <- sp::over( x = locations_threatened_shp, y = grid_10km.wgs , fn = NULL)

species_over_grid_thr <- over_grid_10k_thr %>%
    mutate(id=seq(1:nrow(.))) %>% 
    left_join(locations_threatened, by=c("id"="id")) %>% 
    distinct(CELLCODE,subspeciesname) %>%
    group_by(CELLCODE) %>%
    summarize(total_species=n())


grid_10km_species_data_thr <- grid_10km.wgs_data %>%
    left_join(species_over_grid_thr, by=c("CELLCODE"="CELLCODE"))


map_greece_plot_grid_10_thr <- ggplot()+
  geom_polygon(data = grid_10km_species_data_thr,
               aes(x=long, y=lat,group = group, fill=total_species),
               lwd=0.12, color="orange")+
  geom_polygon(data = hellenic_borders_df,
               aes(x=long, y=lat,group = group),
               lwd=0.12,color="black")+
  geom_polygon(data = grid_10km_df,
               aes(x=long, y=lat,group = group),
               lwd=0.12, fill=NA, color="orange")+
  geom_point(data = locations_threatened,aes(x=logD, y=latD,color=Order),size = 0.2)+
  scale_fill_gradient(low = "yellow", high = "red", na.value = NA)+ 
  ggtitle("Potentially Threatened")+
  labs(x="Longitude",y="Latitude")+
#  scale_fill_gradientn(colours = c("gray100",terrain.colors(10)),na.value =NA ,name="Number of species")+ #c("gray100","gray50","gray40","gray35","gray30","gray20","gray10","gray0")
  scale_x_continuous(breaks = seq(23,27,0.5),limits = c(23,27))+
  scale_y_continuous(breaks = seq(34,36,0.5),limits = c(34,37))+
  coord_map(xlim = c(23,27), ylim = c(34,37))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),legend.position = c(0.87, 0.73),legend.text = element_text(size=9),legend.title = element_text(size=10))
#geom_text(data = sisquoc, aes(label = paste("  ", as.character(name), sep="")), angle = 60, hjust = 0, color = "yellow")
    
ggsave("map_greece_plot_grid_10_thr.png", 
       plot =map_greece_plot_grid_10_thr, 
       device = "png",
       width = 30,
       height = 30,
       units = "cm",
       dpi = 300 ,
       path = "../plots/")

### only the top species
grid_10km_species_data_trh_5 <- grid_10km_species_data_thr %>% filter(total_species>15)

map_greece_plot_grid_10_thr_5 <- ggplot()+
  geom_polygon(data = grid_10km_species_data_trh_5,
               aes(x=long, y=lat,group = group, fill=total_species),
               lwd=0.12, color="orange")+
  geom_polygon(data = hellenic_borders_df,
               aes(x=long, y=lat,group = group),
               lwd=0.12,color="black")+
  geom_polygon(data = grid_10km_df,
               aes(x=long, y=lat,group = group),
               lwd=0.12, fill=NA, color="orange")+
  geom_point(data = locations_threatened,aes(x=logD, y=latD,color=Order),size = 0.2)+
  scale_fill_gradient(low = "yellow", high = "red", na.value = NA)+ 
  labs(x="Longitude",y="Latitude")+
  ggtitle("Potentially Threatened > 15 " )+
#  scale_fill_gradientn(colours = c("gray100",terrain.colors(10)),na.value =NA ,name="Number of species")+ #c("gray100","gray50","gray40","gray35","gray30","gray20","gray10","gray0")
  scale_x_continuous(breaks = seq(23,27,0.5),limits = c(23,27))+
  scale_y_continuous(breaks = seq(34,36,0.5),limits = c(34,37))+
  coord_map(xlim = c(23,27), ylim = c(34,37))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),legend.position = c(0.87, 0.73),legend.text = element_text(size=9),legend.title = element_text(size=10))
#geom_text(data = sisquoc, aes(label = paste("  ", as.character(name), sep="")), angle = 60, hjust = 0, color = "yellow")
    
ggsave("map_greece_plot_grid_10_thr_15.png", 
       plot =map_greece_plot_grid_10_thr_5, 
       device = "png",
       width = 30,
       height = 30,
       units = "cm",
       dpi = 300 ,
       path = "../plots/")



### Critically endangered maps

locations_cr <- endemic_species_cr %>% 
    dplyr::select(!Order) %>%
    left_join(arthropods_kriti_endemic,by=c("subspeciesname"="subspeciesname")) %>%
    dplyr::select(subspeciesname,latD, logD, Order) %>% 
    filter(latD<38, logD<30) %>%
    mutate("id"=seq(1:nrow(.)))

locations_cr_shp <- locations_cr %>% 
    relocate(latD,logD,subspeciesname)

coordinates(locations_cr_shp)<-~ logD+latD
proj4string(locations_cr_shp) <- CRS("+proj=longlat +datum=WGS84")

over_grid_10k_cr <- sp::over( x = locations_cr_shp, y = grid_10km.wgs , fn = NULL)

species_over_grid_cr <- over_grid_10k_cr %>%
    mutate(id=seq(1:nrow(.))) %>% 
    left_join(locations_cr, by=c("id"="id")) %>% 
    distinct(CELLCODE,subspeciesname) %>%
    group_by(CELLCODE) %>%
    summarize(total_species=n())


grid_10km_species_data_cr <- grid_10km.wgs_data %>%
    left_join(species_over_grid_cr, by=c("CELLCODE"="CELLCODE"))


map_greece_plot_grid_10_cr <- ggplot()+
  geom_polygon(data = grid_10km_species_data_cr,
               aes(x=long, y=lat,group = group, fill=total_species),
               lwd=0.12, color="orange")+
  geom_polygon(data = hellenic_borders_df,
               aes(x=long, y=lat,group = group),
               lwd=0.12,color="black")+
  geom_polygon(data = grid_10km_df,
               aes(x=long, y=lat,group = group),
               lwd=0.12, fill=NA, color="orange")+
  geom_point(data = locations_cr,aes(x=logD, y=latD,color=Order),size = 0.2)+
  ggtitle("Potentially CR")+
  scale_fill_gradient(low = "yellow", high = "red", na.value = NA)+ 
  labs(x="Longitude",y="Latitude")+
#  scale_fill_gradientn(colours = c("gray100",terrain.colors(10)),na.value =NA ,name="Number of species")+ #c("gray100","gray50","gray40","gray35","gray30","gray20","gray10","gray0")
  scale_x_continuous(breaks = seq(23,27,0.5),limits = c(23,27))+
  scale_y_continuous(breaks = seq(34,36,0.5),limits = c(34,37))+
  coord_map(xlim = c(23,27), ylim = c(34,37))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),legend.position = c(0.87, 0.73),legend.text = element_text(size=9),legend.title = element_text(size=10))
#geom_text(data = sisquoc, aes(label = paste("  ", as.character(name), sep="")), angle = 60, hjust = 0, color = "yellow")
    
ggsave("map_greece_plot_grid_10_cr.png", 
       plot =map_greece_plot_grid_10_cr, 
       device = "png",
       width = 30,
       height = 30,
       units = "cm",
       dpi = 300 ,
       path = "../plots/")



grid_10km_species_data_cr_5 <- grid_10km_species_data_cr %>% filter(total_species>4)

map_greece_plot_grid_10_cr_5 <- ggplot()+
  geom_polygon(data = grid_10km_species_data_cr_5,
               aes(x=long, y=lat,group = group, fill=total_species),
               lwd=0.12, color="orange")+
  geom_polygon(data = hellenic_borders_df,
               aes(x=long, y=lat,group = group),
               lwd=0.12,color="black")+
  geom_polygon(data = grid_10km_df,
               aes(x=long, y=lat,group = group),
               lwd=0.12, fill=NA, color="orange")+
  geom_point(data = locations_cr,aes(x=logD, y=latD,color=Order),size = 0.2)+
  scale_fill_gradient(low = "yellow", high = "red", na.value = NA)+ 
  labs(x="Longitude",y="Latitude")+
  ggtitle("Potentially CR > 4 " )+
#  scale_fill_gradientn(colours = c("gray100",terrain.colors(10)),na.value =NA ,name="Number of species")+ #c("gray100","gray50","gray40","gray35","gray30","gray20","gray10","gray0")
  scale_x_continuous(breaks = seq(23,27,0.5),limits = c(23,27))+
  scale_y_continuous(breaks = seq(34,36,0.5),limits = c(34,37))+
  coord_map(xlim = c(23,27), ylim = c(34,37))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),legend.position = c(0.87, 0.73),legend.text = element_text(size=9),legend.title = element_text(size=10))
#geom_text(data = sisquoc, aes(label = paste("  ", as.character(name), sep="")), angle = 60, hjust = 0, color = "yellow")
    
ggsave("map_greece_plot_grid_10_cr_5.png", 
       plot =map_greece_plot_grid_10_cr_5, 
       device = "png",
       width = 30,
       height = 30,
       units = "cm",
       dpi = 300 ,
       path = "../plots/")

