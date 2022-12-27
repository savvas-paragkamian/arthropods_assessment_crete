# Arhropoda species assessment in Crete island

## Structure of the repo

```
|-- data
|-- results
|-- plots
|-- scripts

```
First the script `species_assessment.R` is executed. Then the 
`spatial_analysis.R` script. These 2 scripts generate all the 
results that are needed for the `species_assessment_statistict.Rmd`
R markdown file that generates a report. The `functions.R` is a 
collection of functions required in all scripts.

There is an additional script, `code_snippets.R`, which contains 
code and data manipulation of large files, e.g. Natura2000 sites, 
elevation maps of Europe. This script is not meant to be executed
but to guide to the original versions of the data as downloaded from 
the sources.

## Occurrence Data

Occurrences of endemic arthropod species in the island of Crete were 
compiled by Natural History Museaum of Crete.


## Assessment

IUCN

Preliminary Automated Conservation Assessments (PACA)

## Hotspots

Using the [European Environment Agency reference grid](https://www.eea.europa.eu/data-and-maps/data/eea-reference-grids-2) 
of 10km X 10km we defined a location as a cell grid that a species has occurrered. 
The endemic hotspots and the threatened hotspots are the cell grids 
the 10% highest endemic species and threatened species, respectively.

## Protected areas

Natura2000 SCI (habitats directive)

Wildlife refugees.


## Spatial data

Copernicus 

Elevation and slope

Corine 

