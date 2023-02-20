#!/usr/bin/Rscript

# load packages and functions
library(tidyverse)
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
            mapping=aes(colour="Natura2000 SAC"),
            linewidth=0.4,
            alpha=1,
            fill=NA,
            show.legend=T) +
    scale_colour_manual(values = c("Natura2000 SAC" = "#56B4E9"),
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
            mapping=aes(fill="Natura2000 SAC"),
            alpha=1,
            colour="transparent",
            show.legend=T) +
    scale_fill_manual(values = c("Natura2000 SAC" = "#56B4E9"),
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
            mapping=aes(fill="Natura2000 SAC"),
            alpha=1,
            colour="transparent",
            show.legend=T) +
    scale_fill_manual(values = c("Natura2000 SAC" = "#56B4E9"), 
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
    filter(aoo_natura_percent<0.15 & threatened==T)

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
            mapping=aes(fill="Natura2000 SAC"),
            colour="transparent",
            show.legend=T) +
    scale_fill_manual(values = c("Natura2000 SAC" = "#56B4E9"),
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
                        breaks = c(3,6,9,12),
                        labels = c(3,6,9,12),
                        guide = guide_colourbar(override.aes = list(alpha=1),
                                                ticks = FALSE,
                                                label = TRUE,
                                                title="# taxa \n(AOO < 15% in N2K)",
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
            mapping=aes(color="Natura2000 SAC"),
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
                                 "Natura2000 SAC"=NA),
                      guide = "legend") +
    scale_colour_manual(values = c("Natura2000 SAC" = "#56B4E9"),
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

heatmap_sort <- heatmaps_threatspots[[1]]

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
                                                title="# threatspots",
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
            axis.text = element_text(size=13), 
            axis.text.x = element_text(angle = 90, hjust = 0),
            legend.position = c(.90, .83),
            legend.title=element_text(size=9))

ggsave("../figures/fig3b.png",
       plot = fig3b,
       width = 20,
       height = 20,
       units='cm', 
       device = "png",
       dpi = 300)

#fig3b
aoo_dist <- endemic_species %>%
    pivot_longer(cols=c(aoo,eoo,n_locations)) %>%
    dplyr::select(subspeciesname, Order, name, value) %>%
    filter(value>0)

fig3a <-ggplot() +
    geom_point(aoo_dist, 
               mapping=aes(x=Order, y=value, color=name, shape=name),
               position=position_jitterdodge(0.3)) +
#    stat_summary(fun=mean, geom="pointrange", color="red")+
    scale_y_continuous(trans='log10',
                     breaks=trans_breaks('log10', function(x) 10^x),
                     labels=trans_format('log10', math_format(10^.x))) + 
    scale_shape_manual(values = c(2,1,3),
                       labels=c("AOO","EOO", "# locations"),
                       name="Quantity")+
    scale_color_manual(values=c("gray48", "gray60", "gray40"),
                       labels=c("AOO","EOO", "# locations"),
                       name="Quantity")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust = 0),
          axis.text = element_text(size=13), 
          axis.title.x=element_blank(),
          legend.position = c(0.85, 0.1))

ggsave("../figures/fig3a.png", 
       plot=fig3a, 
       device="png", 
       height = 20, 
       width = 23, 
       units="cm")

fig3 <- ggarrange(fig3a, fig3b,
          labels = c("A", "B"),
          align = "hv",
          widths = c(0.8,1),
          ncol = 2,
          nrow = 1,
          font.label=list(color="black",size=22)) + bgcolor("white")

ggsave("../figures/Fig3.tiff", 
       plot=fig3, 
       height = 20, 
       width = 43,
       dpi = 600, 
       units="cm",
       device="tiff")

ggsave("../figures/Fig3.png", 
       plot=fig3, 
       height = 20, 
       width = 43,
       dpi = 600, 
       units="cm",
       device="png")

ggsave("../figures/Fig3.pdf", 
       plot=fig3, 
       height = 20, 
       width = 43,
       dpi = 600, 
       units="cm",
       device="pdf")

ggsave("../figures/Fig3-small.png", 
       plot=fig3, 
       height = 20, 
       width = 43,
       dpi = 300, 
       units="cm",
       device="png")
