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

arthropods <- readxl::read_excel("../data/arthropoda_crete_nhmc_for_analysis.xlsx")



# Spatial data
locations <-  arthropods %>% dplyr::select(subspeciesname,latD, logD, Order) %>% 
    filter(latD<36, logD<30) %>%
    na.omit()

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

EOO <- EOO.computing(locations_shp, country_map=hellenic_borders_shp, export_shp=T)


eoo_results <- read_delim("EOO.results.csv", delim=",", col_names=T)

Acanthopetalum_minotauri <- locations %>% filter(subspeciesname=="Acanthopetalum minotauri")

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
    
ggsave("Acanthopetalum_minotauri.png", plot = Acanthopetalum_minotauri_plot, device = "png",width = 30,height = 30,units = "cm",dpi = 300 ,path = "../plots/")

