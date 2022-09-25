#!/usr/bin/Rscript

library(tidyverse)
library(readxl)
library(rredlist)
library(taxize)
library(rgdal)
library(ConR)
library(vegan)

#arthropods <- readxl::read_excel("../data/arthropoda_crete_nhmc_for_analysis.xlsx")
## remove Opiliones because the dataset has some errors that are under review.
arthropods_kriti <- readxl::read_excel("../data/arthropoda_crete_nhmc_for_analysis.xlsx") %>% 
    filter(Order!="Opiliones") %>% 
    dplyr::select(-Ergasia)

# Spatial data
locations_all <- arthropods_kriti %>% 
    dplyr::select(subspeciesname,latD, logD, Order) %>% 
    filter(latD<38, logD<30) %>%
    na.omit() %>%
    mutate(id=seq(1:nrow(.)))

locations_shp_all <- locations_all %>% 
    dplyr::select(!Order) %>% 
    relocate(latD,logD,subspeciesname)

coordinates(locations_shp_all)<-~ logD+latD
proj4string(locations_shp_all) <- CRS("+proj=longlat +datum=WGS84")

## endemic
arthropods_kriti_endemic <- readxl::read_excel("../data/arthropoda_crete_nhmc_for_analysis.xlsx", sheet="Endemics") %>% 
    filter(Order!="Opiliones") %>% 
    dplyr::select(-Ergasia)

# Spatial data
locations <- arthropods_kriti_endemic %>% 
    dplyr::select(subspeciesname,latD, logD, Order) %>% 
    filter(latD<38, logD<30) %>%
    na.omit() %>%
    mutate(id=seq(1:nrow(.)))

locations_shp <- locations %>% dplyr::select(!Order) %>% relocate(latD,logD,subspeciesname)
coordinates(locations_shp)<-~ logD+latD
proj4string(locations_shp) <- CRS("+proj=longlat +datum=WGS84")


hellenic_borders_shp <- rgdal::readOGR(dsn="~//Documents/spatial_data/hellenic_borders",layer="hellenic_borders",verbose=TRUE)

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

grid_10km_df <- broom::tidy(grid_10km.wgs)

grid_10km.wgs_data <- grid_10km.wgs@data %>% 
        mutate(id=as.character(seq(0,nrow(.)-1))) %>%
        left_join(grid_10km_df, by=c("id"="id"))
### all species grid overlap

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



## AOO
##
## see for bootstrap
AOO_endemic <- AOO.computing(locations_shp,Cell_size_AOO =2 )

AOO_endemic_df <- data.frame(subspeciesname=names(AOO_endemic),AOO=AOO_endemic, row.names=NULL)

endemic_species <- locations %>% 
    group_by(subspeciesname,Order) %>% 
    summarize(n_locations=n()) %>% 
    left_join(AOO_endemic_df, by=c("subspeciesname"="subspeciesname")) %>%
    left_join(eoo_results, by=c("subspeciesname"="X1")) %>% 
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

