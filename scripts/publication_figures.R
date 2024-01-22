#!/usr/bin/Rscript

## Script name: publications_figures.R
##
## Purpose of script: create figures for the publication of
## The conservation status of Cretan Endemic Arthropods under 
## Natura 2000 network
## How to run:
## Rscript publication_figures.R
##
## Execution time: 5 minutes 
##
## Author: Savvas Paragkamian
##
## Date Created: 2023-02-22

# load packages and functions
library(tidyverse)
library(ggrepel)
library(ggnewscale)
library(scales)
library(ggpubr)
library(sf)
library(jpeg)
library(raster)
source("functions.R")

g_base <- g_base()

# load data
locations_shp <- sf::st_read("../data/arthropods_occurrences/arthropods_occurrences.shp")
locations_source <- read_delim("../data/locations_source.tsv", delim="\t", col_names=T) %>%
        st_as_sf(coords=c("logD","latD"),
             remove=F,
             crs="WGS84")

locations_spatial <- sf::st_read("../results/locations_spatial/locations_spatial.shp")
locations_grid <- sf::st_read("../results/locations_grid/locations_grid.shp") 
crete_shp <- sf::st_read("../data/crete/crete.shp")
crete_peaks <- read_delim("../data/crete_mountain_peaks.csv", delim=";", col_names=T) %>% 
    st_as_sf(coords=c("X", "Y"),
             remove=F,
             crs="WGS84")
endemic_species <- read_delim("../results/endemic_species_assessment.tsv", delim="\t")
clc_crete_shp <- st_read("../data/clc_crete_shp/clc_crete_shp.shp")
natura_crete <- sf::st_read("../data/natura2000/natura2000_crete.shp")
wdpa_crete <- sf::st_read("../data/wdpa_crete/wdpa_crete.shp")
natura_crete_land <- st_intersection(natura_crete, crete_shp)

# raster DEM hangling
dem_crete <- raster("../data/dem_crete/dem_crete.tif")
dem_crete_pixel <- as(dem_crete, "SpatialPixelsDataFrame")
dem_crete_df <- as.data.frame(dem_crete_pixel) %>% filter(dem_crete>0)


# split the SPA SCI

natura_crete_land_sci <- natura_crete_land %>% filter(SITETYPE=="B")
## Hotspots and threatspots
endemic_hotspots <- st_read("../results/endemic_hotspots/endemic_hotspots.shp")
threatspots <- st_read("../results/threatspots/threatspots.shp")
threatspots_lt <- threatspots %>% 
    filter(pc_thrt>= quantile(pc_thrt,0.90))

locations_inland <- st_join(locations_shp, crete_shp, left=F)
# IUCN
redlist_orders <- read_delim("../results/redlist_threatened.tsv", delim="\t")

# Colorblind palette
palette.colors(palette = "Okabe-Ito") 
# Crete figures

fig1ab <- readJPEG("../figures/Fig1ab.jpg")

g_fig1ab <- ggplot() +
    background_image(fig1ab)+
    theme(plot.margin = margin(t=0.5, l=0.5, r=0.5, b=0, unit = "cm"))
## Fig1c

crete_base <- ggplot() +
    geom_sf(crete_shp, mapping=aes()) +
    geom_raster(dem_crete_df, mapping=aes(x=x, y=y, fill=dem_crete))+
    scale_fill_gradientn(guide = "colourbar",colours = c("snow3","#F0E442","#D55E00","#CC79A7"),
                        breaks = c(100, 800, 1500, 2400),
                        labels = c(100, 800, 1500, 2400))+
    geom_sf(natura_crete_land_sci,
            mapping=aes(colour="Natura2000 HSD"),
            linewidth=0.4,
            alpha=1,
            fill=NA,
            show.legend=T) +
    scale_colour_manual(values = c("Natura2000 HSD" = "#56B4E9"),
                        guide = guide_legend(override.aes = list(linetype="solid",shape = NA)),
                        name="")+
    new_scale_color()+
    geom_point(locations_source,
            mapping=aes(x=logD, y=latD, color=source, shape=source),
            size=1.8,
            alpha=0.8,
            show.legend=T) +
    geom_label(data = crete_peaks, 
               mapping=aes(x = X, y = Y, label = name),
               size = 1.5,
               nudge_x = 0.05,
               nudge_y=0.05,
               label.padding = unit(0.1, "lines"))+ 
    scale_colour_manual(values = c("NHMC" = "#009E73", 
                                   "Bibliography" = "#999999"),
                        guide = guide_legend(override.aes = list(size= c(3,3),linetype = c("blank", "blank"))),
                        name = "Sampling") +
    scale_shape_manual(values = c("NHMC" = 17,
                                   "Bibliography" = 4),
                        name = "Sampling") +
    guides(fill = guide_colourbar(ticks = FALSE,
                                  label = TRUE,
                                  title="Elevation",
                                  title.vjust = 0.8),
           colour = guide_legend())+
    coord_sf(crs="WGS84") +
    theme_bw()+
    theme(axis.title=element_blank(),
          axis.text=element_text(colour="black"),
          legend.title = element_text(size=8),
          legend.position = "bottom",
          legend.box.background = element_blank())


ggsave("../figures/Fig1c.tiff", 
       plot=crete_base, 
       height = 10, 
       width = 20,
       dpi = 600, 
       units="cm",
       device="tiff")

ggsave("../figures/Fig1c.png", 
       plot=crete_base, 
       height = 10, 
       width = 20,
       dpi = 600, 
       units="cm",
       device="png")

## Combine
fig1 <- ggarrange(g_fig1ab,crete_base,
#          align = "hv",
          nrow = 2,
          legend="right") + bgcolor("white")

ggsave("../figures/Fig1.tiff", 
       plot=fig1, 
       height = 20, 
       width = 25,
       dpi = 600, 
       units="cm",
       device="tiff")

ggsave("../figures/Fig1.png", 
       plot=fig1, 
       height = 20, 
       width = 25,
       dpi = 600, 
       units="cm",
       device="png")

ggsave("../figures/Fig1.pdf", 
       plot=fig1, 
       height = 20, 
       width = 25,
       dpi = 600, 
       units="cm",
       device="pdf")

ggsave("../figures/Fig1-small.png", 
       plot=fig1, 
       height = 20, 
       width = 25,
       dpi = 300, 
       units="cm",
       device="png")

# Figure 2
## fig2a
crete_hotspot <- ggplot() +
    geom_sf(crete_shp, mapping=aes()) +
    geom_sf(natura_crete_land_sci,
            mapping=aes(fill="Natura2000 HSD"),
            alpha=1,
            colour="transparent",
            show.legend=T) +
    scale_fill_manual(values = c("Natura2000 HSD" = "#56B4E9"),
                      guide = guide_legend(title=""))+
    new_scale_fill() +
    geom_sf(endemic_hotspots, mapping=aes(fill=n_species),
            alpha=0.6,
            colour="transparent",
            na.rm = FALSE,
            show.legend=T) +
    scale_fill_gradient(low="#F0E442",
                        high="#D55E00",
                        guide = "colourbar")+
    geom_sf(crete_peaks,
            mapping=aes(),
            colour="#D55E00",
            size=1,
            alpha=1,
            show.legend=F) +
    geom_label(data = crete_peaks, 
               mapping=aes(x = X, y = Y, label = name),
               size = 1.5,
               nudge_x = 0.05,
               nudge_y=0.05, label.padding = unit(0.1, "lines"))+ 
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
          legend.box.background = element_blank())

ggsave("../figures/Fig2a.tiff", 
       plot=crete_hotspot, 
       height = 10, 
       width = 20,
       dpi = 600, 
       units="cm",
       device="tiff")

ggsave("../figures/Fig2a.png", 
       plot=crete_hotspot, 
       height = 10, 
       width = 20,
       dpi = 600, 
       units="cm",
       device="png")

## fig2b
crete_threat <- ggplot() +
    geom_sf(crete_shp, mapping=aes()) +
    geom_sf(natura_crete_land_sci,
            mapping=aes(fill="Natura2000 HSD"),
            alpha=1,
            colour="transparent",
            show.legend=T) +
    scale_fill_manual(values = c("Natura2000 HSD" = "#56B4E9"), 
                      guide = guide_legend(title=""))+
    new_scale_fill() +
    geom_sf(threatspots_lt, mapping=aes(fill=pc_thrt),
            alpha=0.6,
            colour="transparent",
            size=0.1,
            na.rm = FALSE,
            show.legend=T) +
    scale_fill_gradient(low="#E69F00",
                        high="#CC79A7",
                        breaks = c(25,30,35,40,45,50),
                        labels = c(25,30,35,40,45,50),
                        guide = "colourbar")+
    geom_sf(crete_peaks,
            mapping=aes(),
            colour="#D55E00",
            size=1,
            alpha=1,
            show.legend=F) +
    geom_label(data = crete_peaks, 
               mapping=aes(x = X, y = Y, label = name),
               size = 1.5,
               nudge_x = 0.05,
               nudge_y=0.05, label.padding = unit(0.1, "lines"))+ 
    coord_sf(crs="WGS84") +
    guides(fill = guide_colourbar(ticks = FALSE,
                                  label = TRUE,
                                  title="# threatened",
                                  title.vjust = 0.8,
                                  order = 1))+
    theme_bw()+
    theme(axis.title=element_blank(),
          axis.text=element_text(colour="black"),
          legend.title = element_text(size=8),
          legend.position = "bottom",
          legend.box.background = element_blank())

ggsave("../figures/Fig2b.tiff", 
       plot=crete_threat, 
       height = 10, 
       width = 20,
       dpi = 600, 
       units="cm",
       device="tiff")

ggsave("../figures/Fig2b.png", 
       plot=crete_threat, 
       height = 10, 
       width = 20,
       dpi = 600, 
       units="cm",
       device="png")

## fig2c
species_10_natura <- endemic_species %>%
    mutate(aoo_natura_percent=round(aoo_natura/aoo, digits=4)) %>%
    mutate(aoo_natura_relative=round(1-abs(aoo_natura-aoo)/aoo, digits=4)) %>%
    filter(aoo_natura_relative<0.1 & threatened==T)

species_10_natura_l <- locations_grid %>%
    filter(sbspcsn %in% species_10_natura$subspeciesname) %>%
    group_by(CELLCOD) %>%
    summarise(n_species=n()) %>%
    filter(n_species>2, CELLCOD!="10kmE570N150")
species_10_natura_l_o <- locations_grid %>%
    filter(sbspcsn %in% species_10_natura$subspeciesname) %>%
    group_by(CELLCOD, Order) %>%
    summarise(n_species=n(), .groups="drop")

table(species_10_natura$Order)


crete_aoo <- ggplot() +
    geom_sf(crete_shp, mapping=aes()) +
    geom_sf(natura_crete_land_sci,
            mapping=aes(fill="Natura2000 HSD"),
            colour="transparent",
            show.legend=T) +
    scale_fill_manual(values = c("Natura2000 HSD" = "#56B4E9"),
                      guide = guide_legend(title = ""))+
    new_scale_fill() +
    geom_sf(species_10_natura_l, mapping=aes(fill=n_species),
            alpha=0.6,
            colour="transparent",
            size=0.1,
            na.rm = FALSE,
            show.legend=T) +
    scale_fill_gradient(low="#999999",
                        high="#E69F00",
                       # breaks = c(3,6,9,12),
                       # labels = c(3,6,9,12),
                        guide = guide_colourbar(override.aes = list(alpha=1),
                                                ticks = FALSE,
                                                label = TRUE,
                                                title="# taxa \n(AOO < 10% in N2K)",
                                                title.vjust = 0.8,
                                                order = 1))+
    geom_sf(crete_peaks,
            mapping=aes(),
            colour="#D55E00",
            size=1,
            alpha=1,
            show.legend=F) +
    geom_label(data = crete_peaks,
               mapping=aes(x = X, y = Y, label = name),
               size = 1.5,
               nudge_x = 0.05,
               nudge_y=0.05, label.padding = unit(0.1, "lines"))+
    coord_sf(crs="WGS84") +
    theme_bw()+
    theme(axis.title=element_blank(),
          axis.text=element_text(colour="black"),
          legend.title = element_text(size=8),
          legend.position = "bottom",
          legend.box.background = element_blank())

ggsave("../figures/Fig2c.tiff",
       plot=crete_aoo,
       height = 10,
       width = 20,
       dpi = 600,
       units="cm",
       device="tiff")

ggsave("../figures/Fig2c.png",
       plot=crete_aoo,
       height = 10,
       width = 20,
       dpi = 600,
       units="cm",
       device="png")

## fig2d

crete_corine <- ggplot() +
    geom_sf(crete_shp, mapping=aes()) +
    geom_sf(clc_crete_shp,
            mapping=aes(fill=LABEL1),
            alpha=1,
            colour="transparent",
            show.legend=T) +
    geom_sf(natura_crete_land_sci,
            mapping=aes(color="Natura2000 HSD"),
            linewidth=0.6,
            fill=NA,
            alpha=1,
            show.legend=T) +
    geom_sf(crete_peaks,
            mapping=aes(),
            color = "#D55E00",
            size=1,
            alpha=1,
            show.legend=F) +
    geom_label(data = crete_peaks, 
               mapping=aes(x = X, y = Y, label = name),
               size = 1.5,
               nudge_x = 0.05,
               nudge_y=0.05, label.padding = unit(0.1, "lines"))+ 
    scale_fill_manual(values = c("Artificial surfaces"="#000000",
                                 "Agricultural areas"="#E69F00",
                                 "Forest and semi natural areas" = "#009E73",
                                 "Water bodies" = "#0072B2",
                                 "Natura2000 HSD"=NA),
                      guide = "legend") +
    scale_colour_manual(values = c("Natura2000 HSD" = "#56B4E9"),
                        guide = "legend") +
    guides(fill = guide_legend(override.aes = list(color = "transparent", alpha=1) ),
           colour = guide_legend(override.aes = list(alpha=1, fill="transparent") ) )+
    coord_sf(crs="WGS84") +
    theme_bw()+
    theme(axis.title=element_blank(),
          axis.text=element_text(colour="black"),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.box.background = element_blank(),
          legend.key.size = unit(6, "mm"), 
          legend.text=element_text(size=8))

ggsave("../figures/Fig2d.tiff", 
       plot=crete_corine, 
       height = 10, 
       width = 20,
       dpi = 600, 
       units="cm",
       device="tiff")

ggsave("../figures/Fig2d.png", 
       plot=crete_corine, 
       height = 10, 
       width = 20,
       dpi = 600, 
       units="cm",
       device="png")


fig2 <- ggarrange(crete_hotspot,crete_threat,crete_aoo,crete_corine,
          labels = c("A", "B", "C", "D"),
          align = "hv",
          widths = c(1,1,1,0.6),
          ncol = 1,
          nrow = 4,
          font.label=list(color="black",size=22),
          legend="right") + bgcolor("white")

ggsave("../figures/Fig2.tiff", 
       plot=fig2, 
       height = 40, 
       width = 30,
       dpi = 600, 
       units="cm",
       device="tiff")

ggsave("../figures/Fig2.png", 
       plot=fig2, 
       height = 40, 
       width = 30,
       dpi = 600, 
       units="cm",
       device="png")

ggsave("../figures/Fig2.pdf", 
       plot=fig2, 
       height = 40, 
       width = 30,
       dpi = 600, 
       units="cm",
       device="pdf")

ggsave("../figures/Fig2-small.png", 
       plot=fig2, 
       height = 40, 
       width = 30,
       dpi = 300, 
       units="cm",
       device="png")

#figure 3

#fig3a
aoo_dist <- endemic_species %>%
    pivot_longer(cols=c(aoo,eoo,n_locations)) %>%
    dplyr::select(subspeciesname, Order, name, value) %>%
    filter(value>0) %>%
    mutate(Order=gsub("Lepidoptera", "Lepidoptera\n(Geometrid moths)", Order))


fig3a <- ggplot() +
    geom_boxplot(aoo_dist,
                 mapping=aes(x=Order, y=value,color=name),
                 outlier.size = 0) +
    geom_point(aoo_dist, 
               mapping=aes(x=Order, y=value, color=name, shape=name),
               position=position_jitterdodge(0.3)) + 
    geom_vline(xintercept = seq(0.5, length(aoo_dist$Order), by = 1), 
               color="gray", 
               linewidth=.5,
               alpha=.5) + # # set vertical lines between x groups
    scale_y_continuous(trans='log10', name = "Value",
                     breaks=trans_breaks('log10', function(x) 10^x),
                     labels=trans_format('log10', math_format(10^.x))) + 
    scale_shape_manual(values = c(2,1,3),
                       labels=c("AOO in sq. km","EOO in sq. km", "# locations"),
                       name="Quantity")+
    scale_color_manual(values=c("gray15", "gray45", "gray65"),
                       labels=c("AOO in sq. km","EOO in sq. km", "# locations"),
                       name="Quantity")+
    scale_fill_manual(values=c("gray15", "gray45", "gray65"),
                       labels=c("AOO in sq. km","EOO in sq. km", "# locations"),
                       name="Quantity")+
    theme_bw()+
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(face="bold",angle = 90, hjust = 0),
          axis.text = element_text(size=13), 
          axis.title.x=element_blank(),
          axis.title.y=element_text(face="bold", size=13),
          legend.position = c(0.88, 0.1))

ggsave("../figures/fig3a.png", 
       plot=fig3a, 
       device="png", 
       height = 20, 
       width = 23, 
       units="cm")

#fig3b
# Overlap of hotspots
endemic_hotspots_o <- locations_grid %>% 
    filter(CELLCOD %in% endemic_hotspots$CELLCODE) %>%
    distinct(CELLCOD, Order)

heatmaps_hotspots <- heatmaps(endemic_hotspots_o)

ggsave("../figures/fig3aa.png",
       plot = heatmaps_hotspots[[2]],
       width = 25,
       height = 25,
       units='cm', 
       device = "png",
       dpi = 300)

# Overlap of threatspots
threatspots_o <- locations_grid %>% 
    filter(CELLCOD %in% threatspots_lt$CELLCOD) %>%
    distinct(CELLCOD, Order)

heatmaps_threatspots <- heatmaps(threatspots_o) 

heatmap_sort <- heatmaps_threatspots[[1]] %>%
    mutate(from=gsub("Lepidoptera", "Lepidoptera\n(Geometrid moths)", from)) %>%
    mutate(to=gsub("Lepidoptera", "Lepidoptera\n(Geometrid moths)",to))


order_cell_long <- heatmap_sort %>%
    mutate(count=if_else(from==to,0,count))

order_cell_long_t <- order_cell_long %>%
    filter(count!=0)

diagonal <- heatmap_sort %>%
    filter(from==to)


fig3b <- ggplot()+
      geom_tile(data=order_cell_long,
                aes(x=from, y=to,fill=count),
                color="white",
                alpha=1,
                show.legend = T)+
      geom_point(data=diagonal,
                 aes(x=from, y=to),
                 colour="lightyellow4",
                 size=1,
                 show.legend = F)+
      geom_text(data=order_cell_long_t,
                aes(x=from, y=to, label=count),
                size=4) +
      scale_fill_gradient(low="gray87",
                          high="#0072B2",
                          limits=c(1, max(order_cell_long$count)),
                          na.value="white",
                          guide = guide_legend(override.aes = list(alpha=1),
                                                ticks = FALSE,
                                                label = TRUE,
                                                title="# threat-spots",
                                                title.vjust = 0.8,
                                                order = 1))+
      scale_y_discrete(limits = rev)+
      xlab("") +
      ylab("")+
      theme_bw()+
      theme(
            panel.border=element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor=element_blank(),
            axis.text.x = element_text(face="bold",angle = 90, hjust = 0),
            axis.text.y = element_text(face="bold"),
            axis.text = element_text(size=13), 
            axis.title.x=element_blank(),
            axis.title.y=element_text(face="bold", size=13),
            legend.position = c(.90, .83),
            legend.title=element_text(size=9))

ggsave("../figures/fig3b.png",
       plot = fig3b,
       width = 20,
       height = 20,
       units='cm', 
       device = "png",
       dpi = 300)


#fig3 <- ggarrange(fig3a, fig3b,
#          labels = c("A", "B"),
#          align = "hv",
#          widths = c(0.8,1),
#          ncol = 2,
#          nrow = 1,
#          font.label=list(color="black",size=22)) + bgcolor("white")

fig3 <- ggarrange(fig3a, fig3b,
          labels = c("A", "B"),
          align = "hv",
          widths = c(1,1),
          ncol = 1,
          nrow = 2,
          font.label=list(color="black",size=22)) + bgcolor("white")

ggsave("../figures/Fig3.tiff", 
       plot=fig3, 
       height = 40, 
       width = 24,
       dpi = 600, 
       units="cm",
       device="tiff")

ggsave("../figures/Fig3.png", 
       plot=fig3, 
       height = 40, 
       width = 24,
       dpi = 600, 
       units="cm",
       device="png")

ggsave("../figures/Fig3.pdf", 
       plot=fig3, 
       height = 40, 
       width = 24,
       dpi = 600, 
       units="cm",
       device="pdf")

ggsave("../figures/Fig3-small.png", 
       plot=fig3, 
       height = 40, 
       width = 24,
       dpi = 300, 
       units="cm",
       device="png")

## Supplementary Figure 1
orders <- unique(endemic_species$Order)

redlist_threatened <- redlist_orders %>%
    group_by(Order, source) %>%
    mutate(total=sum(n), proportion = round(n/total,digits=2)) %>%
    ungroup() %>%
    filter(Order %in% orders)

redlist_threatened$source <- factor(redlist_threatened$source,
                                    levels=c("endemic_crete_redlist",
                                             "endemic_greek_redlist",
                                             "endemic_europe_redlist",
                                             "world_redlist"))

redlist_threatened$label <- factor(redlist_threatened$source,
                                    labels=c("Cretan endemics IUCN",
                                             "Greek endemic IUCN",
                                             "Europe endemic IUCN",
                                             "World Red IUCN"))

figS1 <- ggplot() +
    geom_col(redlist_threatened,
             mapping=aes(x=Order, y=n, fill=threatened),
             position = position_stack()) +
    geom_text_repel(data=redlist_threatened,
              aes(x=Order,y=n, group=threatened,
                  label = paste(n," (",proportion,")", sep="")),
              size=3,
              position = position_stack(vjust = .6),
              direction="y",
              min.segment.length = 0) +
    scale_fill_manual(values=c("FALSE"="#56B4E9","TRUE"="#D55E00"),
                      labels=c("Not Threatened","Threatened")) +
    theme_bw()+
    ylab("# of taxa") + 
    facet_grid(label ~ ., scales="free")+
    scale_y_continuous(breaks = scales::pretty_breaks(6), limits = c(0, NA))+
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(face="bold",size = 13,angle = 90, hjust = 0),
          axis.text = element_text(size=13), 
          axis.title.y=element_text(face="bold",size = 13),
          axis.title.x=element_blank(),
          legend.position = c(0.90, 0.9),
          strip.text = element_text(size = 12, face="bold"),
          legend.title=element_blank())

ggsave("../figures/figS1.png", 
       plot=figS1, 
       device="png", 
       height = 23, 
       width = 23, 
       units="cm")

## Supplementary Figure 2

order_aoo <- endemic_species %>%
    mutate(aoo_natura_relative=round(1-abs(aoo_natura-aoo)/aoo, digits=4)) %>%
    group_by(Order) %>%
    mutate(average=mean(aoo_natura_relative), std=sd(aoo_natura_relative)) %>%
    mutate(Order=gsub("Lepidoptera", "Lepidoptera\n(Geometrid moths)", Order))

figS2 <- ggplot() +
    geom_boxplot(order_aoo,
                 mapping=aes(x=Order, y=aoo_natura_relative),
                 outlier.size = 0) +
    geom_jitter(order_aoo, 
               mapping=aes(x=Order, y=aoo_natura_relative)) + 
    geom_vline(xintercept = seq(0.5, length(order_aoo$Order), by = 1), 
               color="gray", 
               linewidth=.5, 
               alpha=.5) + # # set vertical lines between x groups
    labs(y="Proportion of AOO overlap with N2K")+
    theme_bw()+
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 0,face="bold"),
          axis.text = element_text(size=13), 
          axis.title.x=element_blank(),
          axis.title.y=element_text(face="bold"),
          legend.position = c(0.85, 0.1))

ggsave("../figures/figS2.png", 
       plot=figS2, 
       device="png", 
       height = 20, 
       width = 23, 
       units="cm")

## Supplementary Figure box plot


##################### create CELLCOD metadata #######################
threatspots_df <- st_drop_geometry(threatspots) %>% mutate(threatspot="threatspot") 
endemic_hotspots_df <- st_drop_geometry(endemic_hotspots) %>%
    mutate(hotspot="hotspot") %>%
    dplyr::select(-n_species)

locations_grid_d <- locations_grid %>% distinct(CELLCOD,geometry)
eea_10_crete_metadata <- st_intersection(locations_grid_d, clc_crete_shp)

eea_10_crete_metadata_d <- eea_10_crete_metadata %>%
    mutate(area = st_area(geometry)) %>%
    group_by(LABEL2, CELLCOD) %>%
    summarise(area_sum=sum(area), .groups="keep") %>%
    st_drop_geometry() %>% units::drop_units() %>%
    pivot_wider(id_cols="CELLCOD", names_from="LABEL2", values_from="area_sum", values_fill=0)

grid_stats <- locations_grid %>%
    st_drop_geometry() %>%
    group_by(CELLCOD) %>%
    summarise(n_species=n()) %>%
    left_join(endemic_hotspots_df, by=c("CELLCOD"="CELLCODE")) %>%
    left_join(threatspots_df) %>%
    left_join(eea_10_crete_metadata_d)

grid_stats_t <-grid_stats %>%
    mutate_at(c('hotspot','threatspot'), ~replace_na(.,"no")) %>%
    mutate_at(c('LT','PT','pc_thrt'), ~replace_na(.,0))

grid_label2_boxplot <- ggplot()+
    geom_boxplot(grid_stats_long_label2,
                 mapping=aes(x= variables, y=values),
                 outlier.size = 0) +
#    geom_point(locations_grid_d, 
#               mapping=aes(x=group, y=n_species),
#               position=position_jitterdodge(0.3)) + 
    scale_y_continuous(name = "Area in m2",
                       labels = scales::comma) + 
#    scale_shape_manual(values = c(2,1,3),
#                       labels=c("AOO in sq. km","EOO in sq. km", "# locations"),
#                       name="Quantity")+
#    scale_color_manual(values=c("gray15", "gray45", "gray65"),
#                       labels=c("AOO in sq. km","EOO in sq. km", "# locations"),
#                       name="Quantity")+
#    scale_fill_manual(values=c("gray15", "gray45", "gray65"),
#                       labels=c("AOO in sq. km","EOO in sq. km", "# locations"),
#                       name="Quantity")+
    theme_bw()+
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(face="bold",angle = 90, hjust = 0),
          axis.text = element_text(size=10), 
          axis.title.x=element_blank(),
          axis.title.y=element_text(face="bold", size=13),
          legend.position = c(0.88, 0.1)) +
    facet_wrap(~ hotspot)

ggsave("../figures/fig_grid_labal2_boxplot.png", 
       plot=grid_label2_boxplot, 
       device="png", 
       height = 20, 
       width = 23, 
       units="cm")


