#!/usr/bin/Rscript

library(tidyverse)
library(readxl)
library(rredlist)
library(taxize)
library(rgdal)
library(ConR)
library(vegan)

###
###

#arthropods <- readxl::read_excel("../data/arthropoda_crete_nhmc_for_analysis.xlsx")

#arthropods_kriti <- readxl::read_excel("../data/arthropoda_crete_nhmc_for_analysis.xlsx") %>% filter(grepl("^Kriti",distribution)) %>% dplyr::select(-Ergasia)

arthropods_kriti_sheet <- readxl::read_excel("../data/arthropoda_crete_nhmc_for_analysis.xlsx", sheet="Endemics") %>% dplyr::select(-Ergasia)

# Spatial data
locations <- arthropods_kriti_sheet %>% dplyr::select(subspeciesname,latD, logD, Order) %>% 
    filter(latD<38, logD<30) %>%
    na.omit() %>%
    mutate(ID=seq(1:nrow(.)))

locations_shp <- locations %>% dplyr::select(!Order) %>% relocate(latD,logD,subspeciesname)
coordinates(locations_shp)<-~ logD+latD
proj4string(locations_shp) <- CRS("+proj=longlat +datum=WGS84")


hellenic_borders_shp <- rgdal::readOGR(dsn="/Users/savvas/Documents/spatial_data/hellenic_borders",layer="hellenic_borders",verbose=TRUE)

proj4string(hellenic_borders_shp) <- CRS("+proj=longlat +datum=WGS84") 

#CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")  # this is WGS84

hellenic_borders_df <- broom::tidy(hellenic_borders_shp)
bbox_hellenic_borders <- hellenic_borders_shp@bbox
bbox_hellenic_borders_lat <- bbox_hellenic_borders


map_greece_plot_lines <- ggplot()+
  geom_polygon(data = hellenic_borders_df,aes(x=long, y=lat,group = group),lwd=0.12,color="black")+
  geom_point(data = locations,aes(x=logD, y=latD,color=Order),size = 0.2)+
  labs(x="Longitude",y="Latitude")+
#  scale_fill_manual(values = c("chartreuse3","purple","cyan3","chocolate2"),labels = c("Natura2000 v30 SCI", "Natura2000 v30 SPA", "Natura2000 v30 SCISPA","Wildlife Refuge"),name="Protected areas")+
#  scale_color_manual(name="", values = c("Caves"="red"))+
  scale_x_continuous(breaks = seq(23,27,0.5),limits = c(23,27))+
  scale_y_continuous(breaks = seq(34,36,0.5),limits = c(34,37))+
  coord_map(xlim = c(23,27), ylim = c(34,37))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),legend.position = c(0.87, 0.73),legend.text = element_text(size=9),legend.title = element_text(size=10))
#geom_text(data = sisquoc, aes(label = paste("  ", as.character(name), sep="")), angle = 60, hjust = 0, color = "yellow")
    
ggsave("map_greece_plot_lines.png", plot = map_greece_plot_lines, device = "png",width = 30,height = 30,units = "cm",dpi = 300 ,path = "../plots/")



### eoo
# Run once on 19/11/2021
#EOO <- EOO.computing(locations_shp, country_map=hellenic_borders_shp, export_shp=T,write_shp=T)

eoo_results <- read_delim("../data/EOO.results.csv", delim=",", col_names=T)

Acanthopetalum_minotauri <- locations %>% filter(subspeciesname=="Carabus banoni")

Acanthopetalum_minotauri_plot <- ggplot()+
  geom_polygon(data = hellenic_borders_df,aes(x=long, y=lat,group = group),lwd=0.12,color="black")+
  geom_point(data = Acanthopetalum_minotauri ,aes(x=logD, y=latD,color=Order),size = 0.2)+
  labs(x="Longitude",y="Latitude")+
#  scale_fill_manual(values = c("chartreuse3","purple","cyan3","chocolate2"),labels = c("Natura2000 v30 SCI", "Natura2000 v30 SPA", "Natura2000 v30 SCISPA","Wildlife Refuge"),name="Protected areas")+
#  scale_color_manual(name="", values = c("Caves"="red"))+
  scale_x_continuous(breaks = seq(23,27,0.5),limits = c(23,27))+
  scale_y_continuous(breaks = seq(34,36,0.5),limits = c(34,37))+
  coord_map(xlim = c(23,27), ylim = c(34,37))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),legend.position = c(0.87, 0.73),legend.text = element_text(size=9),legend.title = element_text(size=10))
#geom_text(data = sisquoc, aes(label = paste("  ", as.character(name), sep="")), angle = 60, hjust = 0, color = "yellow")
    
ggsave("Carabus_banoni.png", plot = Acanthopetalum_minotauri_plot, device = "png",width = 30,height = 30,units = "cm",dpi = 300 ,path = "../plots/")


### Carabus banoni_EOO_poly
###

carabus_banoni_shp <- rgdal::readOGR(dsn="../shapesIUCN/",layer="Carabus banoni_EOO_poly",verbose=TRUE)

proj4string(carabus_banoni_shp) <- CRS("+proj=longlat +datum=WGS84") 

carabus_banoni <- broom::tidy(carabus_banoni_shp)
bbox_carabus_banoni <- carabus_banoni_shp@bbox
bbox_carabus_banoni_lat <- bbox_carabus_banoni
carabus_banoni_points <- locations %>% filter(subspeciesname=="Carabus banoni")

carabus_banoni_plot <- ggplot()+
  geom_polygon(data = hellenic_borders_df,aes(x=long, y=lat,group = group),lwd=0.12,color="black")+
  geom_polygon(data = carabus_banoni,aes(x=long, y=lat,group = group),alpha=0.5,lwd=0.12,color="orange")+
  geom_point(data = carabus_banoni_points,aes(x=logD, y=latD,color=Order),size = 0.2)+
  labs(x="Longitude",y="Latitude")+
#  scale_fill_manual(values = c("chartreuse3","purple","cyan3","chocolate2"),labels = c("Natura2000 v30 SCI", "Natura2000 v30 SPA", "Natura2000 v30 SCISPA","Wildlife Refuge"),name="Protected areas")+
#  scale_color_manual(name="", values = c("Caves"="red"))+
  scale_x_continuous(breaks = seq(23,27,0.5),limits = c(23,27))+
  scale_y_continuous(breaks = seq(34,36,0.5),limits = c(34,37))+
  coord_map(xlim = c(23,27), ylim = c(34,37))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),legend.position = c(0.87, 0.73),legend.text = element_text(size=9),legend.title = element_text(size=10))
#geom_text(data = sisquoc, aes(label = paste("  ", as.character(name), sep="")), angle = 60, hjust = 0, color = "yellow")
    
ggsave("Carabus_banoni.png", plot = carabus_banoni_plot, device = "png",width = 30,height = 30,units = "cm",dpi = 300 ,path = "../plots/")



## Grid 10km

grid_10km <- rgdal::readOGR(dsn="../data/Greece_shapefile/",layer="gr_10km", verbose=T, p4s="EPSG:3035")

#proj4string(grid_10km) <- CRS("EPSG:6258") # this is the EEA reference grid crs.
# see here https://www.eea.europa.eu/data-and-maps/data/eea-reference-grids-2 
grid_10km.wgs <- spTransform(grid_10km, CRS=CRS("+proj=longlat +datum=WGS84"))




over_grid_10k <- sp::over( x = locations_shp , y = grid_10km.wgs , fn = NULL)

species_over_grid <- over_grid_10k %>%
    mutate(ID=seq(1:nrow(.))) %>% 
    left_join(locations, by=c("ID"="ID")) %>% 
    distinct(CELLCODE,subspeciesname) %>%
    group_by(CELLCODE) %>%
    summarize(n_species=n())

grid_10km_df <- broom::tidy(grid_10km.wgs)

grid_10km.wgs_data <- grid_10km.wgs@data %>% 
        mutate(id=as.character(seq(0,nrow(.)-1))) %>%
        left_join(grid_10km_df, by=c("id"="id"))

grid_10km_species_data <- grid_10km.wgs_data %>%
    left_join(species_over_grid, by=c("CELLCODE"="CELLCODE"))

map_greece_plot_grid_10 <- ggplot()+
  geom_polygon(data = grid_10km_species_data,
               aes(x=long, y=lat,group = group, fill=n_species),
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
    
ggsave("map_greece_plot_grid_10.png", plot = map_greece_plot_grid_10, device = "png",width = 30,height = 30,units = "cm",dpi = 300 ,path = "../plots/")

grid_10km_species_data_top <- grid_10km_species_data %>% filter(n_species>50)

map_greece_plot_grid_10_top <- ggplot()+
  geom_polygon(data = grid_10km_species_data_top,
               aes(x=long, y=lat,group = group, fill=n_species),
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


## AOO
##
## see for bootstrap
AOO <- AOO.computing(locations_shp,Cell_size_AOO =2 )

