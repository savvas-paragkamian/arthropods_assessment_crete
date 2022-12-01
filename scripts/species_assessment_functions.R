#!/usr/bin/Rscript

eoo_calculation <- function(occurrences,land_area, plots) {
    
# eoo_calculation is a custom function that takes 3 inputs.
# the location shapefile, the polygon and a TRUE/FALSE value to 
# export or not to plot
    species <- unique(st_drop_geometry(occurrences)[,1])
    calculations <- matrix(ncol=4,nrow=length(species))

    for(s in seq_along(species)){

        calculations[s,1] <- species[s]
        n_occurrences <- occurrences %>%
            filter(subspeciesname==species[s])
        
        rows <- nrow(n_occurrences)
        calculations[s,2] <- rows
    
        if (rows < 3){
            # the EOO is not suitable below 3 occurrences of a taxon

            print(paste0(species[s], " has ", rows, 
                         " occurrences. Moving to the next taxon"))
            calculations[s,3] <- NA
            calculations[s,4] <- NA
            next
        } else { 
            print(paste0(species[s], " has ", rows, " occurrences. Proceeding with calculations"))
            # union the points of each species
            # and then calculate the convex hull of the points of 
            # each species
            species_convex <- st_convex_hull(st_union(n_occurrences))
            species_convex_land <- st_intersection(species_convex,land_area)
            calculations[s,3] <- st_area(species_convex)
            calculations[s,4] <- st_area(species_convex_land)

            if (plots==TRUE){
                g_species(n_occurrences,species_convex, land_area)

            }
        }
    }

    calculations <- as_tibble(calculations)
    colnames(calculations) <- c("subspeciesname", "occurrences","area_convex", "area_convex_overlap" )
    calculations$occurrences <- as.numeric(calculations$occurrences)
    calculations$area_convex <- as.numeric(calculations$area_convex)
    calculations$area_land <- as.numeric(calculations$area_land)
    return(calculations)
}

eoo_single <- function(taxon_occurrences,overlap_area, plots) {
    
    rows <- nrow(taxon_occurrences)
    species <- unique(st_drop_geometry(ss)[,1])
    calculations <- list()
    calculations[2] <- rows

    if (length(species) > 1){
        print(paste0(" has more than one taxon."))

    } else if (rows < 3){
        print(paste0(species, " has ", rows, 
                         " occurrences. Moving to the next taxon"))
        calculations[3] <- NA
        calculations[4] <- NA

    } else {
        # union the points of each species
        # and then calculate the convex hull of the points of 
        # each species
        
        calculations[1] <- species

        species_convex <- st_convex_hull(st_union(taxon_occurrences))
        species_convex_overlap <- st_intersection(species_convex,overlap_area)
        print(st_area(species_convex))
        calculations[3] <- st_area(species_convex)
        calculations[4] <- st_area(species_convex_overlap)
        
        if (plots==TRUE){
            g_species(taxon_occurrences,species_convex, overlap_area)
        }
        return(calculations)
    }

}

g_species <- function(occurrences,convex, land_area){
    
#    convex <- st_convex_hull(st_union(occurrences))
    name <- unique(as.character(occurrences$subspeciesname))
    name_ <- gsub(" ", "_", name)

    g <- ggplot() +
        geom_sf(land_area, mapping=aes()) +
        geom_sf(occurrences, mapping=aes(),color="blue", size=0.1, alpha=0.2) +
        geom_sf(convex, mapping=aes(),color="red", alpha=0.2, size=0.1) +
        coord_sf(crs="WGS84") +
        ggtitle(name)+
        theme_bw()+
        theme(plot.title = element_text(face = "italic"))
    
    ggsave(paste0("../plots/polygons/",name_,"_convex.png"), plot=g, device="png")

}
