#!/usr/bin/Rscript

aoo_overlap <- function(aoo_shp, baseline_map, overlap_area, plots){

    aoo <- aoo_shp[[2]]
    species_aoo <- list()
#    overlap_area_c <- st_combine(overlap_area)
    if (plots){
    }

    for (a in seq_along(aoo)){
        
        species <- names(aoo[a])
        aoo_sf <- st_as_sf(aoo[[a]])
        aoo_area <- sum(st_area(aoo_sf))
        print(aoo_area)
        # Calculate the overlap with protected area

        aoo_overlap <- st_intersection(aoo_sf,overlap_area)

        aoo_overlap_area <- sum(st_area(aoo_overlap))
        print(aoo_overlap_area)
        species_row <- cbind(species, aoo_area, aoo_overlap_area)

        species_aoo[[a]]<- species_row

        if (plots) {
            
            if (is(baseline_map, "sf")==TRUE){
        
                g <- g_species(aoo_sf, overlap_area, baseline_map)

                if (as.numeric(aoo_overlap_area) > 0){
                    g <- g + 
                        geom_sf(aoo_overlap, mapping=aes(fill="orange"),color="orange", alpha=0.2, size=0.2)
                }
                ggsave(paste0("../plots/aoo/","aoo_", as.character(species),".png"), plot=g, device="png")
            }
        }
    }

    species_aoo_df <- as_tibble(do.call(rbind, species_aoo)) 
    species_aoo_df$aoo_area <- as.numeric(species_aoo_df$aoo_area)
    species_aoo_df$aoo_overlap_area <- as.numeric(species_aoo_df$aoo_overlap_area)

    return(species_aoo_df)
}

eoo_calculation <- function(occurrences,baseline_map, overlap_area, plots, prefix) {
    
# eoo_calculation is a custom function that takes 3 inputs.
# the location shapefile, the polygon and a TRUE/FALSE value to 
# export or not to plot
    species <- occurrences %>% st_drop_geometry() %>%
        dplyr::select(subspeciesname) %>% 
        distinct() %>%
        pull()

    calculations <- list(list())

    for(s in seq_along(species)){

        n_occurrences <- occurrences %>%
            filter(subspeciesname==species[s])

        calculations[[s]] <- eoo_single(n_occurrences,baseline_map, overlap_area, plots, prefix)
        
    }

    #calculations_df <- as.data.frame(do.call(rbind, calculations)) 
    return(calculations)
}

#eoo_protected <- function(taxon_occurrences,overlap_area, plots) {
#
#    calculations <- eoo_calculation(taxon_occurrences,overlap_area, plots)
#
#}

eoo_single <- function(taxon_occurrences,baseline_map, overlap_area, plots, prefix) {
    
    prefix <- as.character(prefix)
    rows <- nrow(taxon_occurrences)
    species <- taxon_occurrences %>% st_drop_geometry() %>%
        dplyr::select(subspeciesname) %>% 
        distinct() %>%
        pull()

    calculations <- list()
    print(as.character(species))
    calculations[[1]] <- as.character(species)
    calculations[[2]] <- rows

    if (length(species) > 1){
        print(paste0(" has more than one taxon."))

    } else if (rows < 3){
        print(paste0(species, " has ", rows, 
                         " occurrences. Moving to the next taxon"))
        calculations[[3]] <- NA
        calculations[[4]] <- NA

    } else {
        # union the points of each species
        # and then calculate the convex hull of the points of 
        # each species
        species_convex <- st_convex_hull(st_union(taxon_occurrences))
        print(st_area(species_convex))

        # here we check whether the EOO of a species forms a polygon.
        # There are some rare cases that it forms a line. When that
        # happens the overlap cannot be calculated.

        if (st_geometry_type(species_convex)!="POLYGON"){
            calculations[[3]] <- 0
            calculations[[4]] <- 0

        } else {
            
            calculations[[3]] <- st_area(species_convex)

            if (is(overlap_area, "sf")==TRUE){

                species_convex_overlap <- st_intersection(species_convex,overlap_area)
                calculations[[4]] <- st_area(species_convex_overlap)

            } else {
                calculations[[4]] <- NA
            }
        }
        
        if (plots==TRUE){
            g <- g_species(taxon_occurrences,species_convex, baseline_map)
            
            if (is(baseline_map, "sf")==TRUE){
        
                g <- g + 
                    geom_sf(overlap_area, mapping=aes(),color="green", alpha=0.2, size=0.2)
        
            }

            ggsave(paste0("../plots/polygons/",prefix, "_", as.character(species),".png"), plot=g, device="png")

        }
    }

    return(calculations)
}

g_species <- function(occurrences,convex, baseline_map){
    
#    convex <- st_convex_hull(st_union(occurrences))
    name <- unique(as.character(occurrences$subspeciesname))
    name_ <- gsub(" ", "_", name)

    g <- ggplot() +
        geom_sf(baseline_map, mapping=aes()) +
        geom_sf(occurrences, mapping=aes(),color="blue", size=0.1, alpha=0.2) +
        geom_sf(convex, mapping=aes(),color="red", alpha=0.2, size=0.1) +
        coord_sf(crs="WGS84") +
        ggtitle(name)+
        theme_bw()+
        theme(plot.title = element_text(face = "italic"))
    
        return(g)

}

convert_ll_df <- function(list_of_l_of_l){
    
    list_of_l <- lapply(list_of_l_of_l, unlist)
    df <- data.frame(do.call(rbind, list_of_l))

    return(df)

}
