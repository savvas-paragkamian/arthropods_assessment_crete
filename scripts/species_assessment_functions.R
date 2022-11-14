#!/usr/bin/Rscript

eoo_calculation <- function(occurrences,land_area) {
    
    species <- unique(occurrences$subspeciesname)

    for(s in seq_along(species)){

        n_occurrences <- occurrences %>%
            filter(subspeciesname==species[s])
        
        rows <- nrow(n_occurrences)
    
        if (rows <=3 ){

            print(paste0(species[s], " has ", rows, 
                         " occurrences. Moving to the next taxon"))
            next
        } else { 
            print(paste0(species[s], " has ", rows, " occurrences. Proceeding with 
                         calculations"))
            # union the points of each species
            # and then calculate the convex hull of the points of 
            # each species
            species_convex <- st_convex_hull(st_union(n_occurrences))
            print(class(species_convex))
            print(st_area(species_convex))
            species_convex_land <- st_intersection(species_convex,land_area)
            print(st_area(species_convex_land))
        }
    }
}
